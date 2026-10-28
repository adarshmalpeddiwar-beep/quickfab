// QuickFab mobile service worker: offline app shell + Android share target.
const VERSION = '/*VERSION*/';
const CACHE = 'quickfab-mobile-' + VERSION;
const SHELL = [
  './', './index.html', './manifest.webmanifest', './occt-import-js.js', './occt-import-js.wasm',
  './vendor/three.module.min.js', './vendor/OrbitControls.js', './vendor/RoomEnvironment.js',
  './slider-block-pa6.step', './as1-oc-214.stp', './sg1-c5-214.stp',
  './icons/icon-192.png', './icons/icon-512.png', './icons/apple-touch-icon.png'
];

self.addEventListener('install', e => {
  e.waitUntil(caches.open(CACHE).then(c => c.addAll(SHELL)).then(() => self.skipWaiting()));
});

self.addEventListener('activate', e => {
  e.waitUntil(caches.keys().then(keys => Promise.all(keys.filter(k => k.startsWith('quickfab-mobile-') && k !== CACHE).map(k => caches.delete(k)))).then(() => self.clients.claim()));
});

self.addEventListener('fetch', e => {
  const req = e.request, url = new URL(req.url);
  // Android share sheet: stash the shared STEP file, then open the app on it
  if (req.method === 'POST' && url.pathname.endsWith('/share-target')) {
    e.respondWith((async () => {
      try {
        const form = await req.formData(); const f = form.get('file');
        if (f) { const c = await caches.open('quickfab-share'); await c.put('shared-file', new Response(f, { headers: { 'x-name': encodeURIComponent(f.name || 'shared.step') } })); }
      } catch (err) {}
      return Response.redirect('./?shared=1', 303);
    })());
    return;
  }
  if (req.method !== 'GET' || url.origin !== location.origin) return;
  // pages: network first so updates arrive, cache when offline
  if (req.mode === 'navigate') {
    e.respondWith(fetch(req).then(r => { const copy = r.clone(); caches.open(CACHE).then(c => c.put('./index.html', copy)); return r; }).catch(() => caches.match('./index.html')));
    return;
  }
  // everything else: cache first
  e.respondWith(caches.match(req, { ignoreSearch: true }).then(hit => hit || fetch(req).then(r => { if (r.ok) { const copy = r.clone(); caches.open(CACHE).then(c => c.put(req, copy)); } return r; })));
});
