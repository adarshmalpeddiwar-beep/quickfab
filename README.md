# QuickFab — Stock-to-Part estimator

Upload a STEP file of a machined part and get, in the browser:

- **Geometry analysis**: faces classified as flat / cylinder / cone / free-form; holes (through, blind, counterbored, tap-drill sizes), internal corner radii, slots; rotational vs prismatic test.
- **Process planning**: best orientation, number of setups, material removed per setup, unreachable regions, deep-hole handling.
- **Route comparison** at the chosen quantity: block → VMC 3-axis, 5-axis, laser/waterjet profile → VMC, wire EDM, bar → CNC lathe (+VMC), turn-mill.
- **Stock sizing** to standard Indian plate thicknesses and bar diameters.
- **Costing** at Indian job-shop rates: material ₹/kg, machine-hour rates, setup + CAM amortised over the batch, deburr/inspection, surface finish, overhead, margin, GST.
- **Stock-to-part simulation** (three.js): vise setups and flips, collision-checked roughing and finishing toolpaths on a voxel model, drilling cycles, lathe turning/boring/part-off, laser and wire-EDM profiles.
- **Editable market rates** (kept in the browser's localStorage), with JSON export/import.

Case study: `web/slider-block-pa6.step`, a cast-nylon linear slider block, loads by default.

There are two apps with the same features:

| | Folder | Live (GitHub Pages) |
|---|---|---|
| **Web** (desktop layout) | `web/` | `https://adarshmalpeddiwar-beep.github.io/quickfab/web/` |
| **Phone app** (installable PWA, iOS-style UI, works offline) | `mobile/` | `https://adarshmalpeddiwar-beep.github.io/quickfab/mobile/` |

The site root (`index.html`) sends phones to `mobile/` and computers to `web/`.

### Install the phone app

- **iPhone / iPad (Safari):** open the `mobile/` link → Share → **Add to Home Screen**.
- **Android (Chrome):** open the link → tap **Install** in the app (or menu ⋮ → **Install app**). Once installed, QuickFab also appears in the Android share sheet for `.step` files.

After the first launch the app, the STEP reader and the sample parts are cached, so it opens and quotes without a network connection. Rate edits on the phone are stored on that phone, separately from the web version.

## Run

It is a static site. Serve the `web/` folder with any static server (the STEP reader loads a `.wasm` file, so `file://` will not work):

```bash
cd web
python3 -m http.server 8080
# open http://localhost:8080
```

three.js r160 is loaded from jsDelivr, so the web page needs internet access for the 3D view. The phone app bundles three.js in `mobile/vendor/`.

Serve `mobile/` the same way to try the phone app; the service worker only runs on `localhost` or HTTPS.

## Layout

```
src/engine.js       geometry analysis, process planning, costing, rate tables (no dependencies; browser + Node)
src/app.html        UI, rate editor, case-study card, three.js simulation (engine is inlined at build)
tools/build.py      inlines src/engine.js into src/app.html -> web/index.html
tools/engine-test.js headless run of the planner on the sample parts
web/                deployable static site (index.html, occt-import-js reader + wasm, sample STEP files)
src/mobile.html     phone app UI (tab bar: Part, Simulate, Quote, Rates)
src/mobile-sw.js    phone app service worker (offline cache, Android share target)
tools/build_mobile.py builds mobile/ from src/mobile.html + the engine + the simulator and case-study code in src/app.html
mobile/             deployable PWA (index.html, manifest, sw.js, icons, vendored three.js, STEP reader, samples)
```

The phone app has its own interface but runs the same engine and simulator as the web app: `build_mobile.py` copies them in at build time, so a fix in `src/engine.js` or the simulator in `src/app.html` reaches both after rebuilding both.

After editing anything in `src/`:

```bash
python3 tools/build.py
python3 tools/build_mobile.py /path/to/three.js-r160   # path needed only to refresh mobile/vendor
node tools/engine-test.js
```

## Rates

Seeded October 2026 from Indian supplier and job-work listings, for example: aluminium 6061 bar ₹250–375/kg, SS 304 ₹150–325/kg, EN8 ₹55–130/kg, cast nylon ~₹360/kg; VMC 3-axis ₹300–600/h, 5-axis ₹1,500–4,500/h, CNC turning ₹250–500/h; laser ₹100–500/m; wire EDM ₹0.12–0.30/mm². Surface-finish rates are estimates. Edit them in `DEFAULT_RATES` in `src/engine.js` or in the page.

## Limits

- Process recognition is rule-based geometry, not a trained model.
- Threads are not stored in STEP; holes at tap-drill diameters are assumed tapped (toggle in the UI).
- Sheet-metal bending is not modelled.
- The animation approximates machining on a voxel grid; it is not CAM toolpath or G-code output.

## Third-party

- [occt-import-js](https://github.com/kovacsv/occt-import-js) (Open CASCADE, LGPL-2.1) — STEP reading in WebAssembly. Licences in `web/license.*.txt`.
- [three.js](https://threejs.org) (MIT), loaded from CDN on the web; bundled in `mobile/vendor/` for the phone app.
- `web/as1-oc-214.stp`, `web/sg1-c5-214.stp`: CAx-IF test models, shipped with occt-import-js.
