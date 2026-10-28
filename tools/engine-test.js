// Headless check of the planner: node tools/engine-test.js
// Reads the sample STEP files with occt-import-js and prints features, stock and routes with costs.
const fs = require('fs'), path = require('path');
const web = path.join(__dirname, '..', 'web');
const occtF = require(path.join(web, 'occt-import-js.js'));
const M = require(path.join(__dirname, '..', 'src', 'engine.js'));
const P = { linearUnit: 'millimeter', linearDeflectionType: 'bounding_box_ratio', linearDeflection: 0.0015, angularDeflection: 0.25 };
(async () => {
  const occt = await occtF({ locateFile: f => path.join(web, f) });
  const as1 = occt.ReadStepFile(fs.readFileSync(path.join(web, 'as1-oc-214.stp')), P);
  const sg1 = occt.ReadStepFile(fs.readFileSync(path.join(web, 'sg1-c5-214.stp')), P);
  const sld = occt.ReadStepFile(fs.readFileSync(path.join(web, 'slider-block-pa6.step')), P);
  const parts = [['Slider block', sld.meshes[0], 'pa6'], ['Plate', as1.meshes[10], 'al6061'], ['L-bracket', as1.meshes[9], 'al6061'], ['Bolt', as1.meshes[3], 'ms'], ['Flanged sleeve', sg1.meshes[0], 'al6061']];
  for (const [name, mesh, mat] of parts) {
    const { S, Pl } = M.runAll([mesh], mat, 10);
    const f = S.features;
    console.log(`\n${name}: ${S.A.bbox.size.map(v => v.toFixed(1)).join(' x ')} mm, ${(S.A.volume / 1000).toFixed(1)} cm3, ${f.faces} faces, ${S.turned ? 'rotational' : 'prismatic'}, ${f.holes.length} holes`);
    for (const r of Pl.routes) console.log(`  ${r.rank === 0 ? '*' : ' '} ${r.name.padEnd(38)} Rs ${r.cost.perPart.toFixed(0).padStart(6)} / part  ${r.cost.cycleMin.toFixed(1)} min  ${r.setups.length} setups  stock: ${r.stock.desc}`);
  }
})();
