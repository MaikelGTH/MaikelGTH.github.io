"""
ULTIMATE SHAFT SOLVER v1.5
Cálculo estático y a fatiga de ejes, bulones y tornillos (Shigley / Norton).

Cambios respecto a v1.0:
- Aluminio 6082-T6 con propiedades reales (EN 755-2) y Se' = 0.4·Sut a 5e8 ciclos (Norton).
- kb continuo (kb = 1 para d <= 7.62 mm) y diámetro efectivo d_e = 0.37·d en piezas que no giran.
- Axial con su propio Kf y con kc = 0.85 dentro del Von Mises alternante (Shigley, ec. 7-7).
- Media de compresión tratada correctamente (no penaliza la fatiga).
- FoS estático sin Kf (material dúctil), con la carga de pico nominal.
- Basquin con f(Sut) de Shigley en acero; σm >= Sut se declara fallo estático.
- Optimización con brentq en un intervalo acotado y tres modos geométricos:
  espesor fijo, relación d_int/D fija y mínima masa con restricciones (v1.2).
- Base de materiales ampliada: aceros, tornillería 8.8/10.9/12.9, inox 304,
  aluminios 6061/6082/7020/7075 y titanio Ti-6Al-4V (v1.3).
- Gráficas: curva S-N, criterios estáticos (σ-τ) y diagrama de Haigh con tu punto (v1.4).
- Flexión media con signo negativo: en flexión siempre hay una fibra traccionada, así que
  se usa |M_m| (antes un M_m < 0 se trataba como compresión y daba un FoS a fatiga optimista) (v1.5).
"""
import math
from scipy.optimize import brentq, minimize_scalar

# ==============================================================================
#  1. TEXTOS (ES / EN)
# ==============================================================================
T = {
    'es': {
        'ask_mat': ">> Elige material (número): ",
        'mat_list': "\nMATERIALES DISPONIBLES:",
        'mat_ninf_note': "  N_inf = ciclos a partir de los que se considera vida infinita. (*) Aluminio: referencia, no hay límite real.",
        'pl_sn': "Curva S-N", 'pl_n': "Ciclos N", 'pl_s': "Tensión alternante eq. (MPa)",
        'pl_you': "Tu punto", 'pl_static': "Criterios estáticos (pico nominal, sin Kf)",
        'pl_fat': "Diagrama de Haigh (fatiga, con Kf)", 'pl_safe': "Zona segura",
        'pl_sigma': "σ normal (MPa)", 'pl_tau': "τ cortante (MPa)",
        'pl_sm': "σ' media (MPa)", 'pl_sa': "σ' alternante (MPa)", 'pl_load': "Recta de carga",
        'pl_saved': ">> Gráficas guardadas en: {}",
        'pl_nolib': ">> AVISO: no se pueden dibujar las gráficas. Instala las librerías con:  pip install matplotlib numpy",
        'note_al':  "    > Nota: el aluminio no tiene límite de fatiga; Se' es la resistencia a {:.0e} ciclos.",
        'note_ti':  "    > Nota: el titanio es muy sensible a las entallas (q ≈ 1). Usa Kf ≈ Kt.",
        'note_bolt':"    > Nota: en tornillos la rosca manda. Para un cálculo fino usa Se de rosca (Shigley tabla 8-17) o VDI 2230.",
        'note_ss':  "    > Nota: Se' del inoxidable es una aproximación (0.45·Sut).",
        'mat_src':  "    > Fuente datos: {}",
        'ask_len': ">> Longitud / brazo de palanca (mm): ",
        'ask_D':   ">> Diámetro exterior (mm): ",
        'ask_t':   ">> Espesor de pared (mm) [0 o D/2 = MACIZO]: ",
        'ask_rot': ">> ¿La pieza GIRA bajo la carga de flexión? (s/n): ",
        'ask_Kf':  "\n>> Kf flexión (ej. 1.7): ",
        'ask_Kfa': ">> Kf axial (Intro = igual que flexión): ",
        'ask_Kfs': ">> Kfs torsión (ej. 1.5): ",
        'err_num': ">> Error: introduce un número válido.",
        'err_opt': "Error: selecciona 1 o 2.",
        'err_neg': ">> Error: el valor debe ser positivo.",

        'data_for': "\n--- DATOS DE: {} ---",
        'val_max':  "1. Valor MÁXIMO (pico, en valor absoluto) en {}: ",
        'skip':     "   >> {} es 0. Omitiendo...",
        'cycle':    "2. ¿Cómo es el ciclo?",
        'type_p':   "   'p' PULSANTE (0 a máx) | 'r' ALTERNANTE PURO (±máx) | 'v' VARIABLE (máx y media): ",
        'auto_p':   "   [AUTO] Media = {:.2f} | Alternante = {:.2f}",
        'val_mean': "3. Valor MEDIO en {} (negativo = compresión): ",
        'err_type': "   >> Error: escribe 'p', 'r' o 'v'.",
        'err_mean': "   >> Error: |media| no puede superar el máximo.",

        'flex_conf': "   CONFIGURACIÓN DE FLEXIÓN",
        'flex_opt1': " 1. Modo ESTRIBO/VOLADIZO (fuerza en punta -> M = F·L)",
        'flex_opt2': " 2. Modo GENERAL/EJE (momento flector directo)",
        'choose':    ">> Elige opción (1 o 2): ",
        'lbl_f_tip': "FUERZA EN PUNTA", 'lbl_m_dir': "MOMENTO FLECTOR",
        'lbl_axil':  "FUERZA AXIL",     'lbl_tor':   "MOMENTO TORSOR",
        'src_tip':   "Fuerza en punta", 'src_dir':   "Momento directo",

        'rep_title': "INFORME TÉCNICO - {}",
        'sec_geo':   "[1] GEOMETRÍA Y CARGAS",
        'geo_D':     "    > Diámetro ext:    {:.3f} mm",
        'geo_solid': "    > Tipo:            MACIZO (d.int = 0)",
        'geo_hollow':"    > Espesor:         {:.2f} mm (d.int = {:.2f} mm)",
        'geo_L':     "    > Longitud:        {} mm",
        'geo_W':     "    > Peso estimado:   {:.1f} g",
        'geo_rot':   "    > Giro:            {}",
        'geo_K':     "    > Concentradores:  Kf flex={} | Kf axial={} | Kfs tors={}",
        'load_det':  "\n    > CARGAS (origen flexión: {}):",
        'col_type': "TIPO", 'col_peak': "PICO", 'col_mean': "MEDIA", 'col_alt': "ALTERNANTE",
        'row_flex': "M flex (N·m)", 'row_ax': "F axial (N)", 'row_tor': "Torsor (N·m)",
        'yes': "Sí (rotativo)", 'no': "No (estático, d_e = 0.37·d)",

        'sec_mat':     "[2] MATERIAL Y LÍMITE DE FATIGA",
        'mat_name':    "    > Material: {}",
        'mat_sy':      "    > Sy:  {} MPa",
        'mat_sut':     "    > Sut: {} MPa",
        'mat_se_base': "    > Se' (base): {:.1f} MPa @ {:.0e} ciclos",
        'mat_marin':   "    > Factores de Marin:",
        'mat_se_real': "    > Se (corregido): {:.1f} MPa",

        'sec_stress':  "[3] TENSIONES (VON MISES)",
        'vm_peak':     "    > VM pico nominal (sin Kf, estática): {:.2f} MPa",
        'vm_m':        "    > VM medio (con Kf, fatiga):          {:.2f} MPa",
        'vm_a':        "    > VM alternante (con Kf, fatiga):     {:.2f} MPa",
        'comp_note':   "    > Nota: media de compresión -> no penaliza la fatiga.",

        'sec_fos':   "[4] FACTORES DE SEGURIDAD",
        'fos_stat':  "    [ESTÁTICA] Von Mises vs Sy (sin Kf):  {:.3f}",
        'fos_lang':  "    [FLUENCIA 1er ciclo] Langer (con Kf): {:.3f}",
        'fos_good':  "    [FATIGA]   Goodman modificado:        {:.3f}",
        'fos_asme':  "    [FATIGA]   ASME elíptico:             {:.3f}",

        'sec_life':  "[5] VIDA A FATIGA (Basquin + Goodman)",
        'life_rev':  "    > Tensión eq. totalmente alternante: {:.2f} MPa",
        'life_N':    "    > Ciclos (N): {}",
        'life_inf_st': "Infinita (σ_rev <= Se)",
        'life_inf_al': "> {:.0e} (el aluminio NO tiene límite de fatiga real)",
        'life_low':  "< 1.000 ciclos (zona de bajo número de ciclos)",
        'life_stat': "FALLO ESTÁTICO (σm >= Sut)",

        'opt_ask':  "¿Optimizar el diámetro? (s/n): ",
        'opt_mode': "¿Criterio ESTÁTICO (e) o FATIGA (f)?: ",
        'opt_targ': "FoS objetivo: ",
        'opt_run':  "\nOptimizando...",
        'opt_geo':  "\n¿Cómo varía la geometría?\n 1. Espesor de pared FIJO (tubo comercial; un macizo sigue macizo)\n 2. Relación d.int/D FIJA (escalado geométrico)\n 3. MÍNIMA MASA (D y espesor libres, con restricciones)\n>> Opción (1/2/3): ",
        'opt_tmin': ">> Espesor mínimo fabricable (mm) [Intro = 1.0]: ",
        'opt_dmax': ">> D máximo por espacio disponible (mm) [Intro = 1000]: ",
        'opt_dt':   ">> Relación D/t máxima (Recomendable < 20). Por encima de unos 30, habría que comprobar el pandeo local, y el código no lo calcula. [Intro = 20]: ",
        'geo_dt':   "    > Relación D/t:    {:.2f} (Recomendable < 20)",
        'dt_warn':  "    > AVISO: D/t > 30. Habría que comprobar el pandeo local, y el código no lo calcula.",
        'opt_fail': "\n>> No hay solución con esas restricciones para FoS = {}. Revisa cargas, objetivo o límites.",
        'opt_ok':   "\n>> Solución: D = {:.3f} mm | espesor = {} | masa = {:.1f} g | FoS = {:.3f}",
        'solid_w':  "macizo",
    },
    'en': {
        'ask_mat': ">> Choose material (number): ",
        'mat_list': "\nAVAILABLE MATERIALS:",
        'mat_ninf_note': "  N_inf = cycles beyond which life is taken as infinite. (*) Aluminium: reference only, no true limit.",
        'pl_sn': "S-N curve", 'pl_n': "Cycles N", 'pl_s': "Equivalent alternating stress (MPa)",
        'pl_you': "Your point", 'pl_static': "Static criteria (nominal peak, no Kf)",
        'pl_fat': "Haigh diagram (fatigue, with Kf)", 'pl_safe': "Safe zone",
        'pl_sigma': "Normal σ (MPa)", 'pl_tau': "Shear τ (MPa)",
        'pl_sm': "Mean σ' (MPa)", 'pl_sa': "Alternating σ' (MPa)", 'pl_load': "Load line",
        'pl_saved': ">> Plots saved to: {}",
        'pl_nolib': ">> WARNING: plots cannot be drawn. Install the libraries with:  pip install matplotlib numpy",
        'note_al':  "    > Note: aluminium has no endurance limit; Se' is the strength at {:.0e} cycles.",
        'note_ti':  "    > Note: titanium is highly notch-sensitive (q ≈ 1). Use Kf ≈ Kt.",
        'note_bolt':"    > Note: in bolts the thread governs. For accurate work use thread Se (Shigley table 8-17) or VDI 2230.",
        'note_ss':  "    > Note: stainless Se' is an approximation (0.45·Sut).",
        'mat_src':  "    > Data source: {}",
        'ask_len': ">> Length / lever arm (mm): ",
        'ask_D':   ">> Outer diameter (mm): ",
        'ask_t':   ">> Wall thickness (mm) [0 or D/2 = SOLID]: ",
        'ask_rot': ">> Does the part ROTATE under the bending load? (y/n): ",
        'ask_Kf':  "\n>> Kf bending (e.g. 1.7): ",
        'ask_Kfa': ">> Kf axial (Enter = same as bending): ",
        'ask_Kfs': ">> Kfs torsion (e.g. 1.5): ",
        'err_num': ">> Error: please enter a valid number.",
        'err_opt': "Error: select 1 or 2.",
        'err_neg': ">> Error: value must be positive.",

        'data_for': "\n--- DATA FOR: {} ---",
        'val_max':  "1. MAXIMUM value (peak, absolute) in {}: ",
        'skip':     "   >> {} is 0. Skipping...",
        'cycle':    "2. Load cycle type?",
        'type_p':   "   'p' PULSATING (0 to max) | 'r' FULLY REVERSED (±max) | 'f' FLUCTUATING (max & mean): ",
        'auto_p':   "   [AUTO] Mean = {:.2f} | Alternating = {:.2f}",
        'val_mean': "3. MEAN value in {} (negative = compression): ",
        'err_type': "   >> Error: type 'p', 'r' or 'f'.",
        'err_mean': "   >> Error: |mean| cannot exceed the maximum.",

        'flex_conf': "   BENDING CONFIGURATION",
        'flex_opt1': " 1. CANTILEVER/PEG mode (tip force -> M = F·L)",
        'flex_opt2': " 2. GENERAL/SHAFT mode (direct bending moment)",
        'choose':    ">> Choose option (1 or 2): ",
        'lbl_f_tip': "TIP FORCE",       'lbl_m_dir': "BENDING MOMENT",
        'lbl_axil':  "AXIAL FORCE",     'lbl_tor':   "TORQUE",
        'src_tip':   "Tip force",       'src_dir':   "Direct moment",

        'rep_title': "TECHNICAL REPORT - {}",
        'sec_geo':   "[1] GEOMETRY & LOADS",
        'geo_D':     "    > Outer dia:       {:.3f} mm",
        'geo_solid': "    > Type:            SOLID (inner dia = 0)",
        'geo_hollow':"    > Wall thickness:  {:.2f} mm (inner dia = {:.2f} mm)",
        'geo_L':     "    > Length:          {} mm",
        'geo_W':     "    > Est. weight:     {:.1f} g",
        'geo_rot':   "    > Rotation:        {}",
        'geo_K':     "    > Concentrations:  Kf bend={} | Kf axial={} | Kfs tors={}",
        'load_det':  "\n    > LOADS (bending source: {}):",
        'col_type': "TYPE", 'col_peak': "PEAK", 'col_mean': "MEAN", 'col_alt': "ALTERNATING",
        'row_flex': "Bend. M (N·m)", 'row_ax': "Axial F (N)", 'row_tor': "Torque (N·m)",
        'yes': "Yes (rotating)", 'no': "No (static, d_e = 0.37·d)",

        'sec_mat':     "[2] MATERIAL & ENDURANCE LIMIT",
        'mat_name':    "    > Material: {}",
        'mat_sy':      "    > Sy:  {} MPa",
        'mat_sut':     "    > Sut: {} MPa",
        'mat_se_base': "    > Se' (base): {:.1f} MPa @ {:.0e} cycles",
        'mat_marin':   "    > Marin factors:",
        'mat_se_real': "    > Se (corrected): {:.1f} MPa",

        'sec_stress':  "[3] STRESSES (VON MISES)",
        'vm_peak':     "    > Peak nominal VM (no Kf, static): {:.2f} MPa",
        'vm_m':        "    > Mean VM (with Kf, fatigue):      {:.2f} MPa",
        'vm_a':        "    > Alternating VM (with Kf, fatigue): {:.2f} MPa",
        'comp_note':   "    > Note: compressive mean -> does not penalise fatigue.",

        'sec_fos':   "[4] SAFETY FACTORS",
        'fos_stat':  "    [STATIC]  Von Mises vs Sy (no Kf):   {:.3f}",
        'fos_lang':  "    [1st-cycle YIELD] Langer (with Kf):  {:.3f}",
        'fos_good':  "    [FATIGUE] Modified Goodman:          {:.3f}",
        'fos_asme':  "    [FATIGUE] ASME elliptic:             {:.3f}",

        'sec_life':  "[5] FATIGUE LIFE (Basquin + Goodman)",
        'life_rev':  "    > Equivalent fully reversed stress: {:.2f} MPa",
        'life_N':    "    > Cycles (N): {}",
        'life_inf_st': "Infinite (σ_rev <= Se)",
        'life_inf_al': "> {:.0e} (aluminium has NO true endurance limit)",
        'life_low':  "< 1,000 cycles (low-cycle region)",
        'life_stat': "STATIC FAILURE (σm >= Sut)",

        'opt_ask':  "Optimise the diameter? (y/n): ",
        'opt_mode': "STATIC (s) or FATIGUE (f) criterion?: ",
        'opt_targ': "Target FoS: ",
        'opt_run':  "\nOptimising...",
        'opt_geo':  "\nHow does the geometry vary?\n 1. FIXED wall thickness (commercial tube; a solid stays solid)\n 2. FIXED inner/outer diameter ratio (geometric scaling)\n 3. MINIMUM MASS (D and thickness free, with constraints)\n>> Option (1/2/3): ",
        'opt_tmin': ">> Minimum manufacturable thickness (mm) [Enter = 1.0]: ",
        'opt_dmax': ">> Maximum D for available space (mm) [Enter = 1000]: ",
        'opt_dt':   ">> Maximum D/t ratio (Recommended < 20). Above about 30, local buckling should be checked, and the code does not calculate it. [Enter = 20]: ",
        'geo_dt':   "    > D/t ratio:       {:.2f} (Recommended < 20)",
        'dt_warn':  "    > WARNING: D/t > 30. Local buckling should be checked, and the code does not calculate it.",
        'opt_fail': "\n>> No solution with these constraints for FoS = {}. Check loads, target or limits.",
        'opt_ok':   "\n>> Solution: D = {:.3f} mm | thickness = {} | mass = {:.1f} g | FoS = {:.3f}",
        'solid_w':  "solid",
    },
}

# ==============================================================================
#  2. MATERIALES
# ==============================================================================
# Valores de diseño (mín. de norma o típicos, ver "src"). Unidades: MPa, kg/mm3.
# Se_prime: resistencia a fatiga del material base (probeta pulida, flexión rotativa)
#           a N_e ciclos. Si es None, se estima con la regla de su familia.
# ka: 'machined' -> 4.51·Sut^-0.265 (Shigley, mecanizado) | número fijo.
def _m(es, en, fam, Sy, Sut, E, rho, Se_prime, N_e, ka, src):
    return {"name_es": es, "name_en": en, "family": fam, "S_y": Sy, "S_ut": Sut, "E": E,
            "rho": rho, "Se_prime": Se_prime, "N_e": N_e, "ka": ka, "src": src}

MATERIALS = {
    # --- Aceros ---
    'c45':      _m("Acero C45 normalizado (+N)", "C45 steel, normalised (+N)", "steel",
                   340, 620, 210000, 7.85e-6, None, 1e6, "machined", "EN 10083-2, d<=16 mm, mín."),
    '42crmo4':  _m("Acero 42CrMo4 bonificado (+QT)", "42CrMo4 steel, Q&T (+QT)", "steel",
                   750, 1000, 210000, 7.85e-6, None, 1e6, "machined", "EN 10083-3, 16<d<=40 mm, mín."),
    'bolt_88':  _m("Tornillo clase 8.8", "Bolt, property class 8.8", "bolt",
                   640, 800, 210000, 7.85e-6, None, 1e6, "machined", "ISO 898-1, nominal"),
    'bolt_109': _m("Tornillo clase 10.9", "Bolt, property class 10.9", "bolt",
                   940, 1040, 210000, 7.85e-6, None, 1e6, "machined", "ISO 898-1, nominal"),
    'bolt_129': _m("Tornillo clase 12.9", "Bolt, property class 12.9", "bolt",
                   1100, 1220, 210000, 7.85e-6, None, 1e6, "machined", "ISO 898-1, nominal"),
    'aisi304':  _m("Inoxidable AISI 304 recocido", "Stainless AISI 304, annealed", "stainless",
                   205, 515, 193000, 8.00e-6, 0.45 * 515, 1e7, "machined", "ASTM A276, mín.; Se' aprox. 0.45·Sut"),
    # --- Aluminios (sin límite de fatiga: Se' a 5e8 ciclos) ---
    'al6061':   _m("Aluminio 6061-T6", "Aluminium 6061-T6", "aluminum",
                   276, 310, 68900, 2.70e-6, 96.5, 5e8, 0.85, "ASM, típicos; Sf' R.R. Moore @5e8"),
    'al6082':   _m("Aluminio 6082-T6", "Aluminium 6082-T6", "aluminum",
                   250, 290, 70000, 2.70e-6, 95.0, 5e8, 0.85, "EN 755-2, mín.; Sf' típico @5e8"),
    'al7020':   _m("Aluminio 7020-T6", "Aluminium 7020-T6", "aluminum",
                   290, 350, 70000, 2.78e-6, None, 5e8, 0.85, "EN 755-2, mín.; Sf' regla de Norton"),
    'al7075':   _m("Aluminio 7075-T6", "Aluminium 7075-T6", "aluminum",
                   503, 572, 71700, 2.81e-6, 159.0, 5e8, 0.85, "ASM, típicos; Sf' R.R. Moore @5e8"),
    # --- Titanio ---
    'ti64':     _m("Titanio Ti-6Al-4V (Grado 5) recocido", "Titanium Ti-6Al-4V (Grade 5), annealed", "titanium",
                   880, 950, 113800, 4.43e-6, 510.0, 1e7, "machined", "ASM, típicos; Se' probeta lisa R=-1 @1e7"),
}

FATIGUE_OPTS = {"kd_temperatura": 1.0, "ke_fiabilidad": 0.814, "k_extra": 1.0}  # ke 0.814 = 99 %
KC_AXIAL = 0.85
S1000_RATIO_AL = 0.9


# ==============================================================================
#  3. MOTOR DE CÁLCULO (sin input/print: función pura)
# ==============================================================================
def geometria(d_ext, espesor):
    """Devuelve d_int, área, módulo resistente W y momento polar J (mm)."""
    if espesor <= 0 or espesor >= d_ext / 2:
        d_int = 0.0
    else:
        d_int = d_ext - 2 * espesor
    area = math.pi / 4 * (d_ext**2 - d_int**2)
    W = math.pi * (d_ext**4 - d_int**4) / (32 * d_ext)
    J = math.pi * (d_ext**4 - d_int**4) / 32
    return d_int, area, W, J


def factor_kb(d):
    """Shigley ec. 6-20 en mm, continuo y <= 1."""
    if d <= 7.62:
        return 1.0
    if d <= 51:
        return 1.24 * d ** -0.107
    return 1.51 * d ** -0.157


def factor_f(S_ut, family):
    """Fracción de Sut a 1000 ciclos (Shigley fig. 6-18)."""
    if family not in ("steel", "bolt"):
        return S1000_RATIO_AL   # aprox. 0.9·Sut (Norton) para no férreos e inoxidable
    sut_kpsi = S_ut / 6.894757
    if sut_kpsi < 70:
        return 0.9
    sut_kpsi = min(sut_kpsi, 200)
    return 1.06 - 2.8e-3 * sut_kpsi + 6.9e-6 * sut_kpsi**2


def limite_fatiga(mat, d_ext, gira, solo_axial):
    """Se' y factores de Marin."""
    if mat["Se_prime"] is not None:
        S_prime = mat["Se_prime"]
    elif mat["family"] in ("steel", "bolt", "stainless"):
        S_prime = 0.5 * mat["S_ut"] if mat["S_ut"] <= 1400 else 700.0
    elif mat["family"] == "aluminum":
        S_prime = 0.4 * mat["S_ut"] if mat["S_ut"] < 330 else 130.0   # Norton, @5e8
    else:
        S_prime = 0.5 * mat["S_ut"]

    ka = 4.51 * mat["S_ut"] ** -0.265 if mat["ka"] == "machined" else mat["ka"]
    ka = min(ka, 1.0)

    if solo_axial:
        kb = 1.0
    else:
        d_e = d_ext if gira else 0.37 * d_ext
        kb = factor_kb(d_e)

    kc = KC_AXIAL if solo_axial else 1.0   # combinadas: el 0.85 axial va dentro del VM
    kd, ke, kx = FATIGUE_OPTS["kd_temperatura"], FATIGUE_OPTS["ke_fiabilidad"], FATIGUE_OPTS["k_extra"]
    Se = S_prime * ka * kb * kc * kd * ke * kx
    return S_prime, (ka, kb, kc, kd, ke, kx), Se


def calcular(d_ext, espesor, cargas, kf_f, kf_a, kf_t, mat_key, gira=True):
    mat = MATERIALS[mat_key]
    d_int, area, W, J = geometria(d_ext, espesor)
    r = d_ext / 2

    M_p, M_m, M_a = cargas['flex']
    F_p, F_m, F_a = cargas['axil']
    T_p, T_m, T_a = cargas['tor']

    # Tensiones nominales. Flexión: la fibra crítica es la traccionada -> |M_m|.
    # Axial: el signo de la media se conserva (negativo = compresión).
    sb_m, sb_a = abs(M_m) / W, M_a / W
    sx_ax_m, sx_ax_a = F_m / area, F_a / area
    t_m, t_a = T_m * r / J, T_a * r / J

    solo_axial = (M_a == 0 and M_m == 0 and T_a == 0 and T_m == 0 and (F_a or F_m))

    # --- Fatiga (con Kf). Axial alternante dividido entre 0.85 (Shigley 7-7) ---
    kc_ax_in_vm = 1.0 if solo_axial else KC_AXIAL
    sx_a = kf_f * abs(sb_a) + kf_a * abs(sx_ax_a) / kc_ax_in_vm
    sx_m = kf_f * sb_m + kf_a * sx_ax_m
    vm_a = math.sqrt(sx_a**2 + 3 * (kf_t * t_a)**2)

    compresion = sx_m < 0
    # Media de compresión: no penaliza la fatiga; solo cuenta la torsión media
    vm_m = math.sqrt((0 if compresion else sx_m)**2 + 3 * (kf_t * t_m)**2)

    # --- Estática (sin Kf, material dúctil): pico nominal ---
    sx_peak = abs(sb_m) + abs(sb_a) + abs(sx_ax_m) + abs(sx_ax_a)
    t_peak = abs(t_m) + abs(t_a)
    vm_peak = math.sqrt(sx_peak**2 + 3 * t_peak**2)

    S_prime, factors, Se = limite_fatiga(mat, d_ext, gira, solo_axial)
    Sy, Sut = mat["S_y"], mat["S_ut"]

    inf = float('inf')
    fos_est = Sy / vm_peak if vm_peak > 0 else inf
    vm_m_full = math.sqrt(sx_m**2 + 3 * (kf_t * t_m)**2)
    fos_lang = Sy / (vm_a + vm_m_full) if (vm_a + vm_m_full) > 0 else inf

    term_good = vm_a / Se + vm_m / Sut
    fos_good = 1 / term_good if term_good > 0 else inf
    term_asme = (vm_a / Se)**2 + (vm_m / Sy)**2
    fos_asme = 1 / math.sqrt(term_asme) if term_asme > 0 else inf

    # --- Vida (Basquin con Goodman) ---
    N_e = mat["N_e"]
    life = {'type': None, 'N': None, 'Sf_check': None}
    S_1000 = factor_f(Sut, mat["family"]) * Sut
    if vm_m >= Sut:
        sig_rev = float('nan')
        life['type'] = 'static'
    else:
        sig_rev = vm_a / (1 - vm_m / Sut)
        if sig_rev <= Se:
            life['type'] = 'inf'
        elif sig_rev >= S_1000:
            life['type'] = 'low'
        else:
            b = math.log10(S_1000 / Se) / math.log10(1e3 / N_e)
            a = S_1000 / 1e3**b
            N = (sig_rev / a) ** (1 / b)
            life.update(type='finite', N=N, Sf_check=a * N**b)

    return {
        'd_int': d_int, 'vm_peak': vm_peak, 'vm_m': vm_m, 'vm_a': vm_a, 'compresion': compresion,
        'Se': Se, 'S_prime': S_prime, 'N_e': N_e, 'factors': factors,
        'fos': {'static': fos_est, 'langer': fos_lang, 'goodman': fos_good, 'asme': fos_asme},
        'sig_rev': sig_rev, 'life': life,
        # datos para las gráficas
        'sx_peak': sx_peak, 't_peak': t_peak, 'S_1000': S_1000, 'S_y': Sy, 'S_ut': Sut,
        'family': mat["family"], 'vm_m_full': vm_m_full,
    }


def fos_criterio(D, t, cargas, kf_f, kf_a, kf_t, mat_key, gira, modo):
    """FoS que manda según el criterio: estático, o el más restrictivo de fatiga."""
    r = calcular(D, t, cargas, kf_f, kf_a, kf_t, mat_key, gira)['fos']
    return r['static'] if modo == 'stat' else min(r['goodman'], r['asme'])


def _resolver_D(esp_de_D, D_min, D_max, args, target):
    """D mínimo que cumple FoS >= target con una ley t = esp_de_D(D). FoS es monótono en D."""
    g = lambda D: fos_criterio(D, esp_de_D(D), *args) - target
    try:
        if g(D_min) >= 0:
            return D_min
        if g(D_max) < 0:
            return None
        return brentq(g, D_min, D_max, xtol=1e-4)
    except ValueError:
        return None


def optimizar(geo_modo, espesor_ini, D_ini, cargas, kf_f, kf_a, kf_t, mat_key, gira, target, modo,
              D_max=1000.0, t_min=1.0, DT_max=20.0):
    """
    geo_modo:
      'fixed_t' -> espesor de pared fijo (tubo comercial). Si es macizo, sigue macizo.
      'ratio'   -> relación d_int/D fija (escalado geométrico).
      'mass'    -> D y t libres: mínima masa con t >= t_min, D <= D_max y D/t <= DT_max.
    Devuelve (D, t) o None si no hay solución con las restricciones.
    """
    args = (cargas, kf_f, kf_a, kf_t, mat_key, gira, modo)
    macizo = espesor_ini <= 0 or espesor_ini >= D_ini / 2

    if geo_modo == 'fixed_t' or (geo_modo == 'ratio' and macizo):
        if macizo:
            D = _resolver_D(lambda D: D / 2, 0.5, D_max, args, target)
            return None if D is None else (D, D / 2)
        D = _resolver_D(lambda D: espesor_ini, 2 * espesor_ini + 0.5, D_max, args, target)
        return None if D is None else (D, espesor_ini)

    if geo_modo == 'ratio':
        k = (D_ini - 2 * espesor_ini) / D_ini        # d_int / D
        esp = lambda D: D * (1 - k) / 2
        D = _resolver_D(esp, 0.5, D_max, args, target)
        return None if D is None else (D, esp(D))

    # --- 'mass': barrido en t, para cada t el D mínimo; nos quedamos con la menor área ---
    def area(D, t):
        return geometria(D, t)[1]

    candidatos = []
    # Macizo como candidato
    Ds = _resolver_D(lambda D: D / 2, 0.5, D_max, args, target)
    if Ds is not None:
        candidatos.append((area(Ds, Ds / 2), Ds, Ds / 2))
    # Tubos
    t_hi = D_max / 2
    n = 120
    for i in range(n + 1):
        t = t_min * (t_hi / t_min) ** (i / n)        # barrido logarítmico
        D_lo = max(2 * t + 0.5, 0.5)
        if D_lo > D_max:
            break
        D_hi = min(D_max, DT_max * t)
        if D_lo > D_hi:
            continue
        D = _resolver_D(lambda D, t=t: t, D_lo, D_hi, args, target)
        if D is not None and D / t <= DT_max + 1e-9:
            candidatos.append((area(D, t), D, t))
    if not candidatos:
        return None
    A_best, D, t = min(candidatos)

    # Refinado local en t alrededor del mejor punto del barrido
    if 2 * t < D:
        paso = (t_hi / t_min) ** (1 / n)
        def area_de_t(tt):
            D_lo, D_hi = 2 * tt + 0.5, min(D_max, DT_max * tt)
            if D_lo > D_hi:
                return float('inf')
            Dt = _resolver_D(lambda D, tt=tt: tt, D_lo, D_hi, args, target)
            return float('inf') if Dt is None else area(Dt, tt)
        r = minimize_scalar(area_de_t, bounds=(max(t_min, t / paso), t * paso), method='bounded',
                            options={'xatol': 1e-4})
        if r.success and r.fun < A_best:
            tt = r.x
            D = _resolver_D(lambda D: tt, 2 * tt + 0.5, min(D_max, DT_max * tt), args, target)
            t = tt
    return D, t


# ==============================================================================
#  4. ENTRADA DE DATOS
# ==============================================================================
def get_float(msg, lang, positivo=False, defecto=None):
    while True:
        raw = input(msg).replace(',', '.').strip()
        if raw == '' and defecto is not None:
            return defecto
        try:
            v = float(raw)
            if positivo and v <= 0:
                print(T[lang]['err_neg']); continue
            return v
        except ValueError:
            print(T[lang]['err_num'])


def elegir_material(lang):
    keys = list(MATERIALS)
    print(T[lang]['mat_list'])
    for i, k in enumerate(keys, 1):
        m = MATERIALS[k]
        nombre = m["name_es"] if lang == 'es' else m["name_en"]
        n_txt = f"{m['N_e']:.0e}" + ("*" if m['family'] == 'aluminum' else " ")
        print(f"  {i:>2}. {nombre:<40} Sy={m['S_y']:>5.0f}  Sut={m['S_ut']:>5.0f} MPa  N_inf={n_txt}")
    print(T[lang]['mat_ninf_note'])
    while True:
        raw = input(T[lang]['ask_mat']).strip()
        if raw.isdigit() and 1 <= int(raw) <= len(keys):
            return keys[int(raw) - 1]
        print(T[lang]['err_num'])


def solicitar_carga(nombre, unidad, lang):
    print(T[lang]['data_for'].format(nombre))
    val_max = abs(get_float(T[lang]['val_max'].format(unidad), lang))
    if val_max == 0:
        print(T[lang]['skip'].format(nombre))
        return 0.0, 0.0, 0.0
    print(T[lang]['cycle'])
    while True:
        tipo = input(T[lang]['type_p']).lower().strip()
        if tipo == 'p':
            m, a = val_max / 2, val_max / 2
        elif tipo == 'r':
            m, a = 0.0, val_max
        elif tipo in ('v', 'f'):
            m = get_float(T[lang]['val_mean'].format(unidad), lang)
            if abs(m) > val_max:
                print(T[lang]['err_mean']); continue
            a = val_max - abs(m)
        else:
            print(T[lang]['err_type']); continue
        print(T[lang]['auto_p'].format(m, a))
        return val_max, m, a


def gestionar_flexion(L, lang):
    print(f"\n{'=' * 60}\n{T[lang]['flex_conf']}\n{'=' * 60}")
    print(T[lang]['flex_opt1']); print(T[lang]['flex_opt2'])
    while True:
        op = input(T[lang]['choose']).strip()
        if op == '1':
            p, m, a = solicitar_carga(T[lang]['lbl_f_tip'], "N", lang)
            return (p * L, m * L, a * L), T[lang]['src_tip']
        if op == '2':
            p, m, a = solicitar_carga(T[lang]['lbl_m_dir'], "N·mm", lang)
            return (p, m, a), T[lang]['src_dir']
        print(T[lang]['err_opt'])


# ==============================================================================
#  5. INFORME
# ==============================================================================
def imprimir_informe(D, t, L, cargas, src, kf_f, kf_a, kf_t, mat_key, gira, lang):
    tx = T[lang]
    res = calcular(D, t, cargas, kf_f, kf_a, kf_t, mat_key, gira)
    mat = MATERIALS[mat_key]
    d_int = res['d_int']
    peso = math.pi / 4 * (D**2 - d_int**2) * L * mat["rho"] * 1000
    nombre = mat["name_es"] if lang == 'es' else mat["name_en"]

    print(f"\n{'#' * 80}\n       {tx['rep_title'].format(nombre)}\n{'#' * 80}")
    print(f"\n{tx['sec_geo']}\n{'-' * 30}")
    print(tx['geo_D'].format(D))
    if d_int == 0:
        print(tx['geo_solid'])
    else:
        print(tx['geo_hollow'].format(t, d_int))
        print(tx['geo_dt'].format(D / t))
        if D / t > 30:
            print(tx['dt_warn'])
    print(tx['geo_L'].format(L)); print(tx['geo_W'].format(peso))
    print(tx['geo_rot'].format(tx['yes'] if gira else tx['no']))
    print(tx['geo_K'].format(kf_f, kf_a, kf_t))

    print(tx['load_det'].format(src))
    fmt = "      {:<15} | {:>12} | {:>12} | {:>12}"
    print(fmt.format(tx['col_type'], tx['col_peak'], tx['col_mean'], tx['col_alt']))
    print("      " + "-" * 60)
    for lbl, tup, k in ((tx['row_flex'], cargas['flex'], 1e-3), (tx['row_ax'], cargas['axil'], 1.0),
                        (tx['row_tor'], cargas['tor'], 1e-3)):
        print(fmt.format(lbl, *(f"{x * k:.2f}" for x in tup)))

    print(f"\n{tx['sec_mat']}\n{'-' * 30}")
    print(tx['mat_name'].format(nombre)); print(tx['mat_sy'].format(mat["S_y"])); print(tx['mat_sut'].format(mat["S_ut"]))
    print(tx['mat_se_base'].format(res['S_prime'], res['N_e']))
    ka, kb, kc, kd, ke, kx = res['factors']
    print(tx['mat_marin'])
    print(f"      ka={ka:.3f} | kb={kb:.3f} | kc={kc} | kd={kd} | ke={ke} | kx={kx}")
    print(tx['mat_se_real'].format(res['Se']))
    print(tx['mat_src'].format(mat["src"]))
    fam = mat["family"]
    if fam == "aluminum":
        print(tx['note_al'].format(res['N_e']))
    elif fam == "titanium":
        print(tx['note_ti'])
    elif fam == "bolt":
        print(tx['note_bolt'])
    elif fam == "stainless":
        print(tx['note_ss'])

    print(f"\n{tx['sec_stress']}\n{'-' * 30}")
    print(tx['vm_peak'].format(res['vm_peak'])); print(tx['vm_m'].format(res['vm_m'])); print(tx['vm_a'].format(res['vm_a']))
    if res['compresion']:
        print(tx['comp_note'])

    fo = res['fos']
    print(f"\n{tx['sec_fos']}\n{'-' * 30}")
    print(tx['fos_stat'].format(fo['static'])); print(tx['fos_lang'].format(fo['langer']))
    print(tx['fos_good'].format(fo['goodman'])); print(tx['fos_asme'].format(fo['asme']))

    print(f"\n{tx['sec_life']}\n{'-' * 30}")
    life = res['life']
    if life['type'] == 'static':
        print(tx['life_N'].format(tx['life_stat']))
    else:
        print(tx['life_rev'].format(res['sig_rev']))
        if life['type'] == 'inf':
            txt = tx['life_inf_al'].format(res['N_e']) if mat["family"] == "aluminum" else tx['life_inf_st']
            print(tx['life_N'].format(txt))
        elif life['type'] == 'low':
            print(tx['life_N'].format(tx['life_low']))
        else:
            print(tx['life_N'].format(f"{int(life['N']):,d}"))
    print(f"{'#' * 80}\n")
    return res



# ==============================================================================
#  5b. GRÁFICAS (S-N, criterios estáticos y diagrama de Haigh)
# ==============================================================================
def graficar(res, D, t, nombre_mat, lang, titulo_extra="", mostrar=True, archivo=None):
    tx = T[lang]
    try:
        import matplotlib
        import matplotlib.pyplot as plt
        import numpy as np
    except ImportError:
        print(tx['pl_nolib'])
        return None
    Sy, Sut, Se, S1k, Ne = res['S_y'], res['S_ut'], res['Se'], res['S_1000'], res['N_e']
    fo = res['fos']
    C_SAFE, C_LINE2, C_LINE3, C_PT = "#2e7d5b", "#1f4fd1", "#8a5cc2", "#c8102e"

    fig, axs = plt.subplots(1, 3, figsize=(16, 5.2))
    geo = f"D = {D:.2f} mm, " + ("macizo" if res['d_int'] == 0 else f"t = {t:.2f} mm")
    fig.suptitle(f"{nombre_mat} | {geo} {titulo_extra}", fontsize=12, fontweight='bold')

    # ---------- 1) S-N ----------
    ax = axs[0]
    N_end = Ne * 100
    b = math.log10(S1k / Se) / math.log10(1e3 / Ne)
    a = S1k / 1e3**b
    Nn = np.logspace(3, math.log10(Ne), 200)
    ax.loglog(Nn, a * Nn**b, color=C_SAFE, lw=2, label="Basquin")
    ax.loglog([Ne, N_end], [Se, Se], color=C_SAFE, lw=2, ls='--' if res['family'] == 'aluminum' else '-')
    ax.fill_between(np.r_[Nn, Ne, N_end], 1, np.r_[a * Nn**b, Se, Se], color=C_SAFE, alpha=0.08)
    ax.axvline(Ne, color='grey', lw=0.8, ls=':')
    ax.text(Ne * 1.15, Se * 1.08, f"N_inf = {Ne:.0e}", fontsize=8, color='grey')
    life = res['life']
    sr = res['sig_rev']
    if life['type'] == 'finite':
        ax.loglog(life['N'], sr, 'o', color=C_PT, ms=9, zorder=5, label=f"{tx['pl_you']}: N = {life['N']:.2e}")
    elif life['type'] == 'inf':
        ax.loglog(Ne * 10, sr, 'o', color=C_PT, ms=9, zorder=5, label=f"{tx['pl_you']}: σ_rev = {sr:.0f} ≤ Se")
    elif life['type'] == 'low':
        ax.loglog(1e3, sr, 'o', color=C_PT, ms=9, zorder=5, label=f"{tx['pl_you']}: N < 1e3")
    ax.set_xlim(1e3, N_end)
    ax.set_ylim(max(1, min(Se, sr if sr == sr else Se) * 0.3), max(S1k, sr if sr == sr else S1k) * 1.4)
    ax.set_title(tx['pl_sn'] + f" (Se = {Se:.0f} MPa)")
    ax.set_xlabel(tx['pl_n']); ax.set_ylabel(tx['pl_s'])
    ax.grid(True, which='both', alpha=0.25); ax.legend(fontsize=8, loc='lower left')

    # ---------- 2) Estáticos: plano σ-τ ----------
    ax = axs[1]
    sg = np.linspace(-Sy, Sy, 400)
    t_vm = np.sqrt(np.clip(Sy**2 - sg**2, 0, None) / 3)
    t_tr = np.sqrt(np.clip(Sy**2 - sg**2, 0, None)) / 2
    ax.fill_between(sg, -t_vm, t_vm, color=C_SAFE, alpha=0.12, label=f"{tx['pl_safe']} (Von Mises)")
    ax.plot(sg, t_vm, color=C_SAFE, lw=2); ax.plot(sg, -t_vm, color=C_SAFE, lw=2)
    ax.plot(sg, t_tr, color=C_LINE2, lw=1.5, ls='--', label="Tresca"); ax.plot(sg, -t_tr, color=C_LINE2, lw=1.5, ls='--')
    sp, tp = res['sx_peak'], res['t_peak']
    ax.plot(sp, tp, 'o', color=C_PT, ms=9, zorder=5, label=f"{tx['pl_you']}: FoS = {fo['static']:.2f}")
    ax.plot([0, sp], [0, tp], color=C_PT, lw=1, ls=':')
    lim = max(Sy, sp * 1.1)
    ax.set_xlim(-lim, lim); ax.set_ylim(-lim / 1.6, lim / 1.6)
    ax.axhline(0, color='grey', lw=0.6); ax.axvline(0, color='grey', lw=0.6)
    ax.set_title(tx['pl_static'] + f"\nSy = {Sy:.0f} MPa")
    ax.set_xlabel(tx['pl_sigma']); ax.set_ylabel(tx['pl_tau'])
    ax.grid(True, alpha=0.25); ax.legend(fontsize=8, loc='lower right')

    # ---------- 3) Fatiga: diagrama de Haigh ----------
    ax = axs[2]
    sm = np.linspace(0, Sut, 400)
    good = Se * (1 - sm / Sut)
    asme = np.where(sm <= Sy, Se * np.sqrt(np.clip(1 - (sm / Sy)**2, 0, None)), np.nan)
    lang_y = np.clip(Sy - sm, 0, None)
    safe = np.minimum(good, lang_y)
    ax.fill_between(sm, 0, safe, color=C_SAFE, alpha=0.12, label=f"{tx['pl_safe']} (Goodman + Langer)")
    ax.plot(sm, good, color=C_SAFE, lw=2, label=f"Goodman (FoS = {fo['goodman']:.2f})")
    ax.plot(sm, asme, color=C_LINE2, lw=1.8, ls='--', label=f"ASME (FoS = {fo['asme']:.2f})")
    ax.plot(sm, lang_y, color=C_LINE3, lw=1.2, ls='-.', label=f"Langer (FoS = {fo['langer']:.2f})")
    vm, va = res['vm_m'], res['vm_a']
    ax.plot(vm, va, 'o', color=C_PT, ms=9, zorder=5, label=tx['pl_you'])
    if vm > 0 or va > 0:
        k = fo['goodman'] if fo['goodman'] != float('inf') else 1
        ax.plot([0, vm * k], [0, va * k], color=C_PT, lw=1, ls=':', label=tx['pl_load'])
    ax.set_xlim(0, max(Sut, vm * 1.1) * 1.02); ax.set_ylim(0, max(Se, Sy * 0.5, va * 1.2) * 1.1)
    ax.set_title(tx['pl_fat'] + f"\nSe = {Se:.0f} | Sy = {Sy:.0f} | Sut = {Sut:.0f} MPa")
    ax.set_xlabel(tx['pl_sm']); ax.set_ylabel(tx['pl_sa'])
    ax.grid(True, alpha=0.25); ax.legend(fontsize=8, loc='upper right')

    fig.tight_layout()
    import os
    ruta = os.path.abspath(archivo or "informe.png")
    fig.savefig(ruta, dpi=130)
    print(tx['pl_saved'].format(ruta))
    if mostrar:
        interactivo = matplotlib.get_backend().lower() not in ('agg', 'pdf', 'svg', 'ps', 'cairo', 'template')
        if interactivo:
            try:
                plt.show(block=False)
                plt.pause(0.5)
                return fig
            except Exception:
                pass
        _abrir_archivo(ruta)   # sin ventana disponible: abre el PNG con el visor del sistema
    return fig


def _abrir_archivo(ruta):
    import os, sys, subprocess
    try:
        if sys.platform.startswith('win'):
            os.startfile(ruta)
        elif sys.platform == 'darwin':
            subprocess.Popen(['open', ruta])
        else:
            subprocess.Popen(['xdg-open', ruta], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    except Exception:
        pass

# ==============================================================================
#  6. PROGRAMA PRINCIPAL
# ==============================================================================
def main():
    print("\n" + "-" * 61 + "\n   ULTIMATE SHAFT SOLVER v1.5 (STEEL · ALUMINIUM · TITANIUM)\n" + "-" * 61)
    while True:
        l_in = input("Select Language / Seleccione Idioma (EN / ES): ").lower().strip()
        if l_in in ('en', 'english'): lang = 'en'; break
        if l_in in ('es', 'español', 'espanol', 'spanish'): lang = 'es'; break
    tx = T[lang]

    mat_key = elegir_material(lang)

    L = get_float(tx['ask_len'], lang, positivo=True)
    D = get_float(tx['ask_D'], lang, positivo=True)
    t = get_float(tx['ask_t'], lang)
    gira = input(tx['ask_rot']).lower().strip() in ('s', 'si', 'sí', 'y', 'yes')

    flex, src = gestionar_flexion(L, lang)
    cargas = {'flex': flex,
              'axil': solicitar_carga(tx['lbl_axil'], "N", lang),
              'tor': solicitar_carga(tx['lbl_tor'], "N·mm", lang)}

    kf_f = get_float(tx['ask_Kf'], lang, positivo=True)
    kf_a = get_float(tx['ask_Kfa'], lang, positivo=True, defecto=kf_f) if any(cargas['axil']) else kf_f
    kf_t = get_float(tx['ask_Kfs'], lang, positivo=True)

    res0 = imprimir_informe(D, t, L, cargas, src, kf_f, kf_a, kf_t, mat_key, gira, lang)
    nombre_mat = MATERIALS[mat_key]["name_es" if lang == 'es' else "name_en"]
    graficar(res0, D, t, nombre_mat, lang, "", archivo="informe_inicial.png")

    if input(tx['opt_ask']).lower().strip() in ('s', 'si', 'sí', 'y', 'yes'):
        while True:
            g_in = input(tx['opt_geo']).strip()
            if g_in in ('1', '2', '3'):
                geo_modo = {'1': 'fixed_t', '2': 'ratio', '3': 'mass'}[g_in]; break
        t_min, D_max, DT_max = 1.0, 1000.0, 20.0
        if geo_modo == 'mass':
            t_min = get_float(tx['opt_tmin'], lang, positivo=True, defecto=1.0)
            D_max = get_float(tx['opt_dmax'], lang, positivo=True, defecto=1000.0)
            DT_max = get_float(tx['opt_dt'], lang, positivo=True, defecto=20.0)
            if DT_max > 30:
                print(tx['dt_warn'].strip())
        modo = 'stat' if input(tx['opt_mode']).lower().strip() in ('e', 's') else 'fat'
        target = get_float(tx['opt_targ'], lang, positivo=True)
        print(tx['opt_run'])
        sol = optimizar(geo_modo, t, D, cargas, kf_f, kf_a, kf_t, mat_key, gira, target, modo,
                        D_max=D_max, t_min=t_min, DT_max=DT_max)
        if sol is None:
            print(tx['opt_fail'].format(target))
        else:
            D_opt, t_opt = sol
            res1 = imprimir_informe(D_opt, t_opt, L, cargas, src, kf_f, kf_a, kf_t, mat_key, gira, lang)
            graficar(res1, D_opt, t_opt, nombre_mat, lang, "(optimizado)" if lang == 'es' else "(optimised)",
                     archivo="informe_optimizado.png")
            d_int, area, _, _ = geometria(D_opt, t_opt)
            masa = area * L * MATERIALS[mat_key]["rho"] * 1000
            val = fos_criterio(D_opt, t_opt, cargas, kf_f, kf_a, kf_t, mat_key, gira, modo)
            esp_txt = tx['solid_w'] if d_int == 0 else f"{t_opt:.2f} mm"
            print(tx['opt_ok'].format(D_opt, esp_txt, masa, val))
    try:
        import matplotlib.pyplot as plt
        if plt.get_fignums() and plt.get_backend().lower() not in ('agg', 'pdf', 'svg', 'ps', 'cairo', 'template'):
            print(">> Cierra las ventanas de las gráficas para terminar." if lang == 'es' else ">> Close the plot windows to finish.")
            plt.show()   # mantiene abiertas las ventanas hasta que las cierres
    except ImportError:
        pass

if __name__ == "__main__":
    main()
