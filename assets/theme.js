/* Tema claro / oscuro. Se carga en <head> para evitar el parpadeo.
   Por defecto sigue al sistema; el botón de la cabecera lo fija y se recuerda. */
(function () {
  var KEY = 'tema';
  var mq = window.matchMedia('(prefers-color-scheme: dark)');
  function guardado() { try { return localStorage.getItem(KEY); } catch (e) { return null; } }
  function actual() { return guardado() || (mq.matches ? 'dark' : 'light'); }
  function aplicar(t) {
    document.documentElement.setAttribute('data-theme', t);
    var b = document.querySelector('.theme-btn');
    if (b) {
      b.textContent = t === 'dark' ? '☀' : '☾';
      b.setAttribute('aria-label', t === 'dark' ? 'Cambiar a tema claro' : 'Cambiar a tema oscuro');
      b.title = b.getAttribute('aria-label');
    }
    window.dispatchEvent(new Event('themechange'));
  }
  document.documentElement.setAttribute('data-theme', actual());
  mq.addEventListener('change', function () { if (!guardado()) aplicar(actual()); });
  document.addEventListener('DOMContentLoaded', function () {
    var bar = document.querySelector('header.site .wrap');
    if (!bar) return;
    var b = document.createElement('button');
    b.type = 'button'; b.className = 'theme-btn';
    b.addEventListener('click', function () {
      var t = actual() === 'dark' ? 'light' : 'dark';
      try { localStorage.setItem(KEY, t); } catch (e) {}
      aplicar(t);
    });
    bar.appendChild(b);
    aplicar(actual());
  });
})();
