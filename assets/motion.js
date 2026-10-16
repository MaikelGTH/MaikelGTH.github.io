/* Movimiento al desplazarse. Se carga en <head> (sin defer) para poner la
   clase .js-motion antes del primer pintado y evitar parpadeos.
   - Elementos con [data-reveal] aparecen al entrar en pantalla (escalonados
     dentro de .projects, .grid, .tl y .two).
   - La línea de tiempo (.tl) se dibuja según avanza el scroll.
   - La cabecera gana fondo translúcido y un filete de progreso de lectura.
   Con prefers-reduced-motion no se anima nada y todo se ve desde el inicio. */
(function () {
  var root = document.documentElement;
  var reduce = window.matchMedia('(prefers-reduced-motion: reduce)').matches;
  var motion = !reduce && 'IntersectionObserver' in window;
  if (motion) root.classList.add('js-motion');

  // Llegada con transición entre páginas: lo visible se muestra ya, sin
  // esperar al revelado, para que el título pueda transformarse.
  window.addEventListener('pagereveal', function (e) {
    if (!motion || !e.viewTransition) return;
    root.classList.add('vt-instant');
    [].forEach.call(document.querySelectorAll('[data-reveal]'), function (el) {
      var r = el.getBoundingClientRect();
      if (r.top < window.innerHeight && r.bottom > 0) el.classList.add('in');
    });
    e.viewTransition.finished.finally(function () { root.classList.remove('vt-instant'); });
  });

  document.addEventListener('DOMContentLoaded', function () {
    var header = document.querySelector('header.site');
    var bar = null;
    if (header) {
      bar = document.createElement('div');
      bar.className = 'progress';
      bar.setAttribute('aria-hidden', 'true');
      header.appendChild(bar);
    }
    var tls = motion ? [].slice.call(document.querySelectorAll('.tl')) : [];

    if (motion) {
      // Retardo escalonado para hermanos dentro del mismo contenedor
      [].forEach.call(document.querySelectorAll('.projects, .grid, .tl, .two, .links'), function (box) {
        var i = 0;
        [].forEach.call(box.children, function (el) {
          if (el.hasAttribute('data-reveal')) el.style.setProperty('--d', Math.min(i++, 4) * 90 + 'ms');
        });
      });
      var io = new IntersectionObserver(function (entries) {
        entries.forEach(function (e) {
          if (e.isIntersecting) { e.target.classList.add('in'); io.unobserve(e.target); }
        });
      }, { rootMargin: '0px 0px -8% 0px', threshold: 0.12 });
      var items = document.querySelectorAll('[data-reveal]');
      [].forEach.call(items, function (el) { io.observe(el); });
      // Red de seguridad: si algo falla, que nada quede oculto
      setTimeout(function () {
        [].forEach.call(items, function (el) {
          if (el.getBoundingClientRect().top < window.innerHeight) el.classList.add('in');
        });
      }, 2500);
    }

    var ticking = false;
    function update() {
      ticking = false;
      var y = window.scrollY, max = root.scrollHeight - window.innerHeight;
      if (header) header.classList.toggle('scrolled', y > 8);
      if (bar) bar.style.setProperty('--progress', max > 0 ? Math.min(y / max, 1).toFixed(4) : 0);
      var mark = window.innerHeight * 0.62;
      tls.forEach(function (tl) {
        var r = tl.getBoundingClientRect();
        var p = Math.max(0, Math.min(1, (mark - r.top) / r.height));
        tl.style.setProperty('--tl', p.toFixed(4));
        [].forEach.call(tl.children, function (li) {
          li.classList.toggle('lit', li.getBoundingClientRect().top + 12 < mark);
        });
      });
    }
    function onScroll() { if (!ticking) { ticking = true; requestAnimationFrame(update); } }
    window.addEventListener('scroll', onScroll, { passive: true });
    window.addEventListener('resize', onScroll);
    update();
  });
})();
