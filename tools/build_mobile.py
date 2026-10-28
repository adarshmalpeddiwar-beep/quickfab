#!/usr/bin/env python3
"""Build the QuickFab phone app (PWA) into mobile/.

The phone app is a separate entity from web/ with its own UI (src/mobile.html), but it shares
the planning engine (src/engine.js) and the stock-to-part simulator and case-study code, which are
copied out of src/app.html at build time so both apps always run the same machining logic.
web/ is not touched by this script.

Usage: python3 tools/build_mobile.py [path-to-three.js-r160-checkout]
"""
import hashlib, json, pathlib, shutil, sys

root = pathlib.Path(__file__).resolve().parent.parent
src, out, web = root / 'src', root / 'mobile', root / 'web'
three = pathlib.Path(sys.argv[1]) if len(sys.argv) > 1 else None

app = (src / 'app.html').read_text()
engine = (src / 'engine.js').read_text()
page = (src / 'mobile.html').read_text()

# shared slices of the web app
sim = app[app.index('const METAL'):app.index('\nlet sim;')]
case = app[app.index('function renderCase()'):app.index('function renderGeo()')]
page = page.replace('<script>/*ENGINE*/</script>', '<script>\n' + engine + '\n</script>')
page = page.replace('/*SIM*/', sim).replace('/*CASE*/', case)
for marker in ('/*ENGINE*/', '/*SIM*/', '/*CASE*/'):
    assert marker not in page, marker

out.mkdir(exist_ok=True)
(out / 'vendor').mkdir(exist_ok=True)
(out / 'icons').mkdir(exist_ok=True)
(out / 'index.html').write_text(page)

# runtime files shared with the web build (STEP reader + samples)
for f in ['occt-import-js.js', 'occt-import-js.wasm', 'license.occt.txt', 'license.occt-import-js.txt',
          'slider-block-pa6.step', 'as1-oc-214.stp', 'sg1-c5-214.stp']:
    shutil.copyfile(web / f, out / f)

# three.js r160, vendored so the installed app works offline
if three:
    shutil.copyfile(three / 'build/three.module.min.js', out / 'vendor/three.module.min.js')
    shutil.copyfile(three / 'examples/jsm/controls/OrbitControls.js', out / 'vendor/OrbitControls.js')
    shutil.copyfile(three / 'examples/jsm/environments/RoomEnvironment.js', out / 'vendor/RoomEnvironment.js')
    shutil.copyfile(three / 'LICENSE', out / 'vendor/LICENSE.three.txt')
for f in ['three.module.min.js', 'OrbitControls.js', 'RoomEnvironment.js']:
    assert (out / 'vendor' / f).exists(), 'missing vendor/' + f + ' (pass the three.js r160 checkout path)'

manifest = {
    "name": "QuickFab — STEP to machining quote",
    "short_name": "QuickFab",
    "description": "Upload a STEP file, get a process plan, an Indian-market machining quote and a stock-to-part simulation.",
    "id": "./", "start_url": "./", "scope": "./",
    "display": "standalone", "orientation": "portrait",
    "background_color": "#f2f2f7", "theme_color": "#f2f2f7",
    "categories": ["business", "productivity", "utilities"],
    "icons": [
        {"src": "icons/icon-192.png", "sizes": "192x192", "type": "image/png"},
        {"src": "icons/icon-512.png", "sizes": "512x512", "type": "image/png"},
        {"src": "icons/icon-maskable-512.png", "sizes": "512x512", "type": "image/png", "purpose": "maskable"}
    ],
    "share_target": {"action": "./share-target", "method": "POST", "enctype": "multipart/form-data",
                     "params": {"files": [{"name": "file", "accept": [".step", ".stp", "application/step", "model/step", "application/octet-stream"]}]}},
    "file_handlers": [{"action": "./", "accept": {"model/step": [".step", ".stp"]}}]
}
(out / 'manifest.webmanifest').write_text(json.dumps(manifest, indent=2))

# service worker, versioned by the content it caches
h = hashlib.sha1()
for p in sorted(out.rglob('*')):
    if p.is_file() and p.name != 'sw.js':
        h.update(p.read_bytes())
version = h.hexdigest()[:10]
(out / 'sw.js').write_text((src / 'mobile-sw.js').read_text().replace('/*VERSION*/', version))
print('wrote mobile/ —', len(page), 'bytes index.html, cache version', version)
