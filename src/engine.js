/* Stock-to-Part engine: geometry analysis, process planning, costing, voxel simulation data.
   Pure JS, no dependencies. Works in browser (window.MFG) and Node (module.exports). */
(function (root) {
'use strict';

// ---------- small vector helpers ----------
const sub = (a, b) => [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
const add = (a, b) => [a[0] + b[0], a[1] + b[1], a[2] + b[2]];
const mul = (a, s) => [a[0] * s, a[1] * s, a[2] * s];
const dot = (a, b) => a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
const cross = (a, b) => [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]];
const len = a => Math.hypot(a[0], a[1], a[2]);
const norm = a => { const l = len(a) || 1; return [a[0] / l, a[1] / l, a[2] / l]; };
const AXES = [[1, 0, 0], [0, 1, 0], [0, 0, 1]];
const AXN = ['X', 'Y', 'Z'];

function anyPerp(u) {
  const t = Math.abs(u[0]) < 0.9 ? [1, 0, 0] : [0, 1, 0];
  return norm(cross(u, t));
}

// Jacobi eigen-decomposition of symmetric 3x3; returns {vals asc, vecs (columns as arrays)}
function eigSym3(m) {
  const a = [m[0].slice(), m[1].slice(), m[2].slice()];
  const v = [[1, 0, 0], [0, 1, 0], [0, 0, 1]];
  for (let sweep = 0; sweep < 50; sweep++) {
    let off = Math.abs(a[0][1]) + Math.abs(a[0][2]) + Math.abs(a[1][2]);
    if (off < 1e-14) break;
    for (let p = 0; p < 2; p++) for (let q = p + 1; q < 3; q++) {
      if (Math.abs(a[p][q]) < 1e-16) continue;
      const th = (a[q][q] - a[p][p]) / (2 * a[p][q]);
      const t = Math.sign(th || 1) / (Math.abs(th) + Math.sqrt(th * th + 1));
      const c = 1 / Math.sqrt(t * t + 1), s = t * c;
      for (let k = 0; k < 3; k++) { const akp = a[k][p], akq = a[k][q]; a[k][p] = c * akp - s * akq; a[k][q] = s * akp + c * akq; }
      for (let k = 0; k < 3; k++) { const apk = a[p][k], aqk = a[q][k]; a[p][k] = c * apk - s * aqk; a[q][k] = s * apk + c * aqk; }
      for (let k = 0; k < 3; k++) { const vkp = v[k][p], vkq = v[k][q]; v[k][p] = c * vkp - s * vkq; v[k][q] = s * vkp + c * vkq; }
    }
  }
  const idx = [0, 1, 2].sort((i, j) => a[i][i] - a[j][j]);
  return { vals: idx.map(i => a[i][i]), vecs: idx.map(i => [v[0][i], v[1][i], v[2][i]]) };
}

function solve3(A, b) {
  const M = A.map((r, i) => [...r, b[i]]);
  for (let c = 0; c < 3; c++) {
    let p = c; for (let r = c + 1; r < 3; r++) if (Math.abs(M[r][c]) > Math.abs(M[p][c])) p = r;
    [M[c], M[p]] = [M[p], M[c]];
    if (Math.abs(M[c][c]) < 1e-18) return null;
    for (let r = 0; r < 3; r++) if (r !== c) { const f = M[r][c] / M[c][c]; for (let k = c; k < 4; k++) M[r][k] -= f * M[c][k]; }
  }
  return [M[0][3] / M[0][0], M[1][3] / M[1][1], M[2][3] / M[2][2]];
}

// ---------- mesh preparation ----------
// accepts occt-import-js mesh objects (one or several, merged into one body)
function prepareMeshes(meshes) {
  const P = [], VN = [], faces = [];
  for (const m of meshes) {
    const pos = m.attributes.position.array, idx = m.index.array;
    const nor = m.attributes.normal ? m.attributes.normal.array : null;
    const ntri = idx.length / 3;
    const bf = (m.brep_faces && m.brep_faces.length) ? m.brep_faces : [{ first: 0, last: ntri - 1 }];
    const base = P.length / 9;
    for (let t = 0; t < ntri; t++) for (let k = 0; k < 3; k++) { const vi = idx[t * 3 + k] * 3; P.push(pos[vi], pos[vi + 1], pos[vi + 2]); if (nor) VN.push(nor[vi], nor[vi + 1], nor[vi + 2]); else VN.push(0, 0, 0); }
    for (const f of bf) faces.push({ first: base + f.first, last: base + f.last });
  }
  const T = buildTriData(new Float64Array(P), faces); T.vn = new Float32Array(VN); return T;
}

function buildTriData(pos, faces) {
  const n = pos.length / 9;
  const nrm = new Float64Array(n * 3), area = new Float64Array(n), cen = new Float64Array(n * 3), faceOf = new Int32Array(n).fill(-1);
  for (let t = 0; t < n; t++) {
    const o = t * 9;
    const a = [pos[o], pos[o + 1], pos[o + 2]], b = [pos[o + 3], pos[o + 4], pos[o + 5]], c = [pos[o + 6], pos[o + 7], pos[o + 8]];
    const cr = cross(sub(b, a), sub(c, a)); const l = len(cr);
    area[t] = l / 2;
    if (l > 0) { nrm[t * 3] = cr[0] / l; nrm[t * 3 + 1] = cr[1] / l; nrm[t * 3 + 2] = cr[2] / l; }
    cen[t * 3] = (a[0] + b[0] + c[0]) / 3; cen[t * 3 + 1] = (a[1] + b[1] + c[1]) / 3; cen[t * 3 + 2] = (a[2] + b[2] + c[2]) / 3;
  }
  faces.forEach((f, i) => { for (let t = f.first; t <= f.last && t < n; t++) faceOf[t] = i; });
  // unassigned triangles → own face
  let extra = null;
  for (let t = 0; t < n; t++) if (faceOf[t] < 0) { if (!extra) { extra = { first: t, last: t, loose: true }; faces.push(extra); } faceOf[t] = faces.length - 1; }
  return { pos, nrm, area, cen, faceOf, faces, n };
}

// ---------- face classification ----------
function faceTris(T, fi) { const f = T.faces[fi]; const r = []; for (let t = f.first; t <= f.last && t < T.n; t++) if (T.faceOf[t] === fi) r.push(t); return r; }
function triN(T, t) { return [T.nrm[t * 3], T.nrm[t * 3 + 1], T.nrm[t * 3 + 2]]; }
function triC(T, t) { return [T.cen[t * 3], T.cen[t * 3 + 1], T.cen[t * 3 + 2]]; }
function triV(T, t, k) { const o = t * 9 + k * 3; return [T.pos[o], T.pos[o + 1], T.pos[o + 2]]; }

function classifyFace(T, fi, diag) {
  const tris = faceTris(T, fi);
  let A = 0; const mn = [0, 0, 0];
  const C = [[0, 0, 0], [0, 0, 0], [0, 0, 0]];
  const verts = [];
  for (const t of tris) {
    const a = T.area[t]; if (a <= 0) continue; A += a;
    const n = triN(T, t);
    for (let i = 0; i < 3; i++) { mn[i] += a * n[i]; for (let j = 0; j < 3; j++) C[i][j] += a * n[i] * n[j]; }
    for (let k = 0; k < 3; k++) verts.push(triV(T, t, k));
  }
  const face = { id: fi, tris, area: A, type: 'free' };
  if (A <= 0) { face.type = 'degenerate'; return face; }
  const mnl = len(mn) / A;
  const nMean = norm(mn);
  // plane?
  let planar = mnl > 0.9995;
  if (planar) { for (const t of tris) if (T.area[t] > 0 && dot(triN(T, t), nMean) < 0.9990) { planar = false; break; } }
  if (planar) {
    face.type = 'plane'; face.normal = nMean;
    let d = 0; for (const t of tris) d += T.area[t] * dot(triC(T, t), nMean); face.offset = d / A;
    return face;
  }
  for (let i = 0; i < 3; i++) for (let j = 0; j < 3; j++) C[i][j] /= A;
  const E = eigSym3(C);
  // cylinder: normals perpendicular to axis
  if (E.vals[0] < 2e-3) {
    const u = norm(E.vecs[0]);
    let ok = true; for (const t of tris) if (T.area[t] > 0 && Math.abs(dot(triN(T, t), u)) > 0.06) { ok = false; break; }
    if (ok) {
      const fit = fitCircle(verts, u);
      if (fit && fit.r > 0 && fit.rms / fit.r < 0.02 && fit.r < diag * 5) {
        Object.assign(face, cylInfo(T, tris, verts, u, fit));
        face.type = 'cyl';
        return face;
      }
    }
  }
  // cone: normals at constant angle to axis (axis ~ mean normal)
  if (mnl > 0.15) {
    const u = nMean; let s = 0, s2 = 0, w = 0;
    for (const t of tris) { const c = dot(triN(T, t), u); s += T.area[t] * c; s2 += T.area[t] * c * c; w += T.area[t]; }
    const m = s / w, sd = Math.sqrt(Math.max(0, s2 / w - m * m));
    if (sd < 0.02 && m < 0.995) {
      // axis line: least squares point closest to all normal lines
      const pt = closestPointToLines(T, tris);
      if (pt) { face.type = 'cone'; face.axis = u; face.point = pt; face.cosAngle = m; face.vrad = verts.map(v => len(cross(sub(v, pt), u))); return face; }
    }
  }
  return face;
}

function closestPointToLines(T, tris) {
  const A = [[0, 0, 0], [0, 0, 0], [0, 0, 0]], b = [0, 0, 0];
  for (const t of tris) {
    const n = triN(T, t), c = triC(T, t), w = T.area[t];
    for (let i = 0; i < 3; i++) for (let j = 0; j < 3; j++) { const m = (i === j ? 1 : 0) - n[i] * n[j]; A[i][j] += w * m; b[i] += w * m * c[j]; }
  }
  return solve3(A, b);
}

function fitCircle(verts, u) {
  const e1 = anyPerp(u), e2 = cross(u, e1);
  let Sxx = 0, Sxy = 0, Syy = 0, Sx = 0, Sy = 0, Sxz = 0, Syz = 0, Sz = 0, N = verts.length;
  const pts = verts.map(v => [dot(v, e1), dot(v, e2)]);
  for (const [x, y] of pts) { const z = x * x + y * y; Sxx += x * x; Sxy += x * y; Syy += y * y; Sx += x; Sy += y; Sxz += x * z; Syz += y * z; Sz += z; }
  const sol = solve3([[Sxx, Sxy, Sx], [Sxy, Syy, Sy], [Sx, Sy, N]], [-Sxz, -Syz, -Sz]);
  if (!sol) return null;
  const cx = -sol[0] / 2, cy = -sol[1] / 2, r2 = cx * cx + cy * cy - sol[2];
  if (!(r2 > 0)) return null;
  const r = Math.sqrt(r2);
  let e = 0; for (const [x, y] of pts) { const d = Math.hypot(x - cx, y - cy) - r; e += d * d; }
  return { cx, cy, r, rms: Math.sqrt(e / N), e1, e2 };
}

function cylInfo(T, tris, verts, u, fit) {
  const { e1, e2, cx, cy, r } = fit;
  // axis point: (cx,cy) in plane, t=0
  const point = add(mul(e1, cx), mul(e2, cy));
  let s = 0;
  for (const t of tris) { const c = triC(T, t); const rv = sub(sub(c, point), mul(u, dot(sub(c, point), u))); s += T.area[t] * Math.sign(dot(triN(T, t), rv)); }
  const ang = verts.map(v => Math.atan2(dot(v, e2) - cy, dot(v, e1) - cx)).sort((a, b) => a - b);
  let gap = 0; for (let i = 1; i < ang.length; i++) gap = Math.max(gap, ang[i] - ang[i - 1]);
  gap = Math.max(gap, ang[0] + 2 * Math.PI - ang[ang.length - 1]);
  const ts = verts.map(v => dot(sub(v, point), u));
  return { axis: u, point, r, internal: s < 0, span: 2 * Math.PI - gap, t0: Math.min(...ts), t1: Math.max(...ts) };
}

// canonicalize axis direction (for grouping)
function canonAxis(u) {
  const ax = [Math.abs(u[0]), Math.abs(u[1]), Math.abs(u[2])];
  const k = ax[0] >= ax[1] && ax[0] >= ax[2] ? 0 : ax[1] >= ax[2] ? 1 : 2;
  return u[k] < 0 ? mul(u, -1) : u.slice();
}
function lineDist(p, u, q) { const d = sub(q, p); return len(sub(d, mul(u, dot(d, u)))); }

// ---------- part-level analysis ----------
function analyzePart(T) {
  // bbox, volume, area
  const mn = [Infinity, Infinity, Infinity], mx = [-Infinity, -Infinity, -Infinity];
  for (let i = 0; i < T.pos.length; i += 3) for (let a = 0; a < 3; a++) { const v = T.pos[i + a]; if (v < mn[a]) mn[a] = v; if (v > mx[a]) mx[a] = v; }
  const size = sub(mx, mn), diag = len(size), ctr = mul(add(mn, mx), 0.5);
  let vol = 0, area = 0;
  for (let t = 0; t < T.n; t++) {
    const o = t * 9;
    const a = [T.pos[o] - ctr[0], T.pos[o + 1] - ctr[1], T.pos[o + 2] - ctr[2]], b = [T.pos[o + 3] - ctr[0], T.pos[o + 4] - ctr[1], T.pos[o + 5] - ctr[2]], c = [T.pos[o + 6] - ctr[0], T.pos[o + 7] - ctr[1], T.pos[o + 8] - ctr[2]];
    vol += dot(a, cross(b, c)) / 6; area += T.area[t];
  }
  const flipped = vol < 0; vol = Math.abs(vol);
  if (flipped) { for (let i = 0; i < T.nrm.length; i++) T.nrm[i] = -T.nrm[i]; }

  const faces = T.faces.map((_, i) => classifyFace(T, i, diag));
  // merge coaxial cylinders
  const cyls = [];
  for (const f of faces) if (f.type === 'cyl') {
    const u = canonAxis(f.axis);
    const tShift = dot(u, f.axis) < 0;
    const t0 = tShift ? -f.t1 : f.t0, t1 = tShift ? -f.t0 : f.t1;
    let g = cyls.find(c => Math.abs(dot(c.axis, u)) > 0.9995 && lineDist(c.point, c.axis, f.point) < Math.max(0.02 * f.r, 0.01) && Math.abs(c.r - f.r) < Math.max(0.005 * f.r, 0.005) && c.internal === f.internal && t0 <= c.t1 + 0.01 && t1 >= c.t0 - 0.01);
    if (g) { g.span += f.span; g.t0 = Math.min(g.t0, t0); g.t1 = Math.max(g.t1, t1); g.faces.push(f.id); g.area += f.area; }
    else cyls.push({ axis: u, point: f.point, r: f.r, internal: f.internal, span: f.span, t0, t1, faces: [f.id], area: f.area });
  }
  // normalise cylinder point to t=0 projection
  for (const c of cyls) { c.span = Math.min(c.span, 2 * Math.PI); c.full = c.span > 5.4; }

  const planes = faces.filter(f => f.type === 'plane');
  const cones = faces.filter(f => f.type === 'cone');
  const free = faces.filter(f => f.type === 'free');
  const freeArea = free.reduce((s, f) => s + f.area, 0);

  // holes: full internal cylinders; group into stacks by axis line
  const holeCyl = cyls.filter(c => c.internal && c.full);
  const stacks = [];
  for (const c of holeCyl) {
    let s = stacks.find(s => Math.abs(dot(s.axis, c.axis)) > 0.9995 && lineDist(s.point, s.axis, c.point) < Math.max(0.02 * c.r, 0.01));
    if (!s) { s = { axis: c.axis, point: c.point, segs: [] }; stacks.push(s); }
    s.segs.push(c);
  }
  for (const s of stacks) {
    s.segs.sort((a, b) => a.t0 - b.t0);
    s.t0 = Math.min(...s.segs.map(c => c.t0)); s.t1 = Math.max(...s.segs.map(c => c.t1));
    const capAt = (tEnd) => {
      for (const f of planes) {
        if (Math.abs(Math.abs(dot(f.normal, s.axis)) - 1) > 1e-3) continue;
        // face lies at tEnd?
        const c0 = triC(T, f.tris[0]); if (Math.abs(dot(sub(c0, s.point), s.axis) - tEnd) > 0.02 + 1e-3 * diag) continue;
        let inside = true; const rEnd = Math.max(...s.segs.map(c => c.r)) * 1.02 + 0.01;
        for (const t of f.tris) for (let k = 0; k < 3; k++) if (lineDist(s.point, s.axis, triV(T, t, k)) > rEnd) { inside = false; break; }
        if (inside) return 'flat';
      }
      for (const f of cones) {
        if (Math.abs(dot(f.axis, s.axis)) < 0.999 || lineDist(s.point, s.axis, f.point) > 0.05 + 0.02 * s.segs[0].r) continue;
        const ts = []; for (const t of f.tris) for (let k = 0; k < 3; k++) ts.push(dot(sub(triV(T, t, k), s.point), s.axis));
        const lo = Math.min(...ts), hi = Math.max(...ts);
        const minR = Math.min(...f.vrad);
        if (minR < 0.05 * s.segs[0].r + 0.02 && ((tEnd <= s.t0 + 1e-6 && hi <= s.t0 + 0.05 && lo < s.t0 - 0.01) || (tEnd >= s.t1 - 1e-6 && lo >= s.t1 - 0.05 && hi > s.t1 + 0.01))) return 'drill-point';
      }
      return null;
    };
    s.cap0 = capAt(s.t0); s.cap1 = capAt(s.t1);
    s.through = !s.cap0 && !s.cap1;
    s.dMin = 2 * Math.min(...s.segs.map(c => c.r)); s.dMax = 2 * Math.max(...s.segs.map(c => c.r));
    s.depth = s.t1 - s.t0;
    s.cutDepth = s.segs.reduce((a, c) => a + (c.t1 - c.t0), 0);
    s.counterbored = new Set(s.segs.map(c => c.r.toFixed(2))).size > 1;
    s.interrupted = s.segs.length > 1 && !s.counterbored;
  }

  // internal corner radii (partial internal cylinders) → smallest tool
  const cornerR = cyls.filter(c => c.internal && !c.full).map(c => c.r);
  const slots = cyls.filter(c => c.internal && !c.full && c.span > 2.9 && c.span < 3.4);

  return { T, faces, cyls, planes, cones, free, freeArea, stacks, cornerR, slots, bbox: { mn, mx, size, diag, ctr }, volume: vol, area, flipped };
}

// ---------- rotational test ----------
function rotationalTest(A) {
  const { T, faces, cyls, bbox } = A;
  const cands = [];
  for (const c of cyls) if (c.span > 3) cands.push({ axis: c.axis, point: c.point });
  for (let k = 0; k < 3; k++) cands.push({ axis: AXES[k], point: bbox.ctr });
  let best = null;
  for (const cd of cands) {
    const u = cd.axis, p = cd.point;
    let side = 0, comp = 0, extCylLen = 0, perpArea = 0;
    for (const f of faces) {
      if (f.type === 'degenerate') continue;
      if (f.type === 'plane' && Math.abs(Math.abs(dot(f.normal, u)) - 1) < 2e-3) { perpArea += f.area; continue; }
      side += f.area;
      let ok = false;
      if (f.type === 'cyl') ok = Math.abs(dot(f.axis, u)) > 0.9995 && lineDist(p, u, f.point) < 0.01 * f.r + 0.01;
      else if (f.type === 'cone') ok = Math.abs(dot(f.axis, u)) > 0.999 && lineDist(p, u, f.point) < 0.02 * bbox.diag;
      else if (f.type === 'free') {
        let ca = 0; for (const t of f.tris) {
          const n = triN(T, t), c = triC(T, t), nx = cross(n, u), l = len(nx);
          const rc = lineDist(p, u, c);
          if (l < 1e-3 || Math.abs(dot(sub(c, p), nx)) / l < Math.max(0.08 * rc, 0.003 * bbox.diag)) ca += T.area[t];
        }
        ok = ca > 0.9 * f.area;
      }
      if (ok) comp += f.area;
    }
    for (const c of cyls) if (!c.internal && c.span > 5.4 && Math.abs(dot(c.axis, u)) > 0.9995 && lineDist(p, u, c.point) < 0.01 * c.r + 0.01) extCylLen += c.t1 - c.t0;
    const frac = side > 0 ? comp / side : 0;
    if (!best || frac + (extCylLen > 0 ? 0.01 : 0) > best.frac + (best.extCylLen > 0 ? 0.01 : 0)) best = { axis: u, point: p, frac, extCylLen, side, comp };
  }
  // length along axis
  let tMin = Infinity, tMax = -Infinity, rMax = 0;
  for (let i = 0; i < T.pos.length; i += 3) {
    const v = [T.pos[i], T.pos[i + 1], T.pos[i + 2]]; const t = dot(sub(v, best.point), best.axis);
    tMin = Math.min(tMin, t); tMax = Math.max(tMax, t); rMax = Math.max(rMax, lineDist(best.point, best.axis, v));
  }
  best.L = tMax - tMin; best.tMin = tMin; best.rMax = rMax;
  best.isRot = best.frac >= 0.55 && best.extCylLen >= 0.15 * best.L;
  return best;
}

// ---------- frames ----------
// frame: p' = R (p - o)  rows of R are new axes
function makeFrame(R, o) { return { R, o }; }
function toFrame(F, p) { const d = sub(p, F.o); return [dot(F.R[0], d), dot(F.R[1], d), dot(F.R[2], d)]; }
function dirToFrame(F, u) { return [dot(F.R[0], u), dot(F.R[1], u), dot(F.R[2], u)]; }
function frameTris(T, F) {
  const out = new Float64Array(T.pos.length);
  for (let i = 0; i < T.pos.length; i += 3) { const q = toFrame(F, [T.pos[i], T.pos[i + 1], T.pos[i + 2]]); out[i] = q[0]; out[i + 1] = q[1]; out[i + 2] = q[2]; }
  return out;
}

// ---------- standard stock ----------
const BAR_D = [6, 8, 10, 12, 14, 16, 18, 20, 22, 25, 28, 30, 32, 35, 38, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100, 110, 120, 125, 130, 140, 150, 160, 180, 200, 220, 250, 280, 300, 350, 400];
const PLATE_T = [1, 1.5, 2, 3, 4, 5, 6, 8, 10, 12, 16, 20, 25, 30, 32, 35, 40, 45, 50, 60, 65, 70, 80, 90, 100, 110, 120, 130, 150, 160, 180, 200, 250, 300];
const nextStd = (list, v) => list.find(x => x >= v - 1e-6) || Math.ceil(v / 10) * 10;
const roundUp = (v, s) => Math.ceil(v / s - 1e-9) * s;

// ---------- voxel model ----------
// Builds grid over stock box in frame coords. state: 0 outside stock, 1 removable, 2 part, 3 hole (drilled)
function voxelize(fpos, stockMin, stockMax, opts) {
  const size = sub(stockMax, stockMin);
  const maxCells = opts.maxCells || 110000;
  let v = Math.max(Math.max(...size) / (opts.res || 60), Math.cbrt(size[0] * size[1] * size[2] / maxCells));
  const nx = Math.max(1, Math.round(size[0] / v)), ny = Math.max(1, Math.round(size[1] / v)), nz = Math.max(1, Math.round(size[2] / v));
  const vx = size[0] / nx, vy = size[1] / ny, vz = size[2] / nz;
  const N = nx * ny * nz; const state = new Uint8Array(N);
  const id = (i, j, k) => (k * ny + j) * nx + i;
  // bin triangles in XY
  const ntri = fpos.length / 9; const bins = new Array(nx * ny);
  for (let t = 0; t < ntri; t++) {
    const o = t * 9;
    const x0 = Math.min(fpos[o], fpos[o + 3], fpos[o + 6]), x1 = Math.max(fpos[o], fpos[o + 3], fpos[o + 6]);
    const y0 = Math.min(fpos[o + 1], fpos[o + 4], fpos[o + 7]), y1 = Math.max(fpos[o + 1], fpos[o + 4], fpos[o + 7]);
    const i0 = Math.max(0, Math.floor((x0 - stockMin[0]) / vx - 0.5)), i1 = Math.min(nx - 1, Math.ceil((x1 - stockMin[0]) / vx - 0.5));
    const j0 = Math.max(0, Math.floor((y0 - stockMin[1]) / vy - 0.5)), j1 = Math.min(ny - 1, Math.ceil((y1 - stockMin[1]) / vy - 0.5));
    for (let j = j0; j <= j1; j++) for (let i = i0; i <= i1; i++) { const b = j * nx + i; (bins[b] || (bins[b] = [])).push(t); }
  }
  const inStock = opts.inStock || (() => true);
  const colIntervals = new Array(nx * ny);
  for (let j = 0; j < ny; j++) for (let i = 0; i < nx; i++) {
    const x = stockMin[0] + (i + 0.5) * vx + 1.234e-5 * vx, y = stockMin[1] + (j + 0.5) * vy + 2.718e-5 * vy;
    const hits = [];
    for (const t of (bins[j * nx + i] || [])) {
      const o = t * 9;
      const ax = fpos[o], ay = fpos[o + 1], bx = fpos[o + 3], by = fpos[o + 4], cx = fpos[o + 6], cy = fpos[o + 7];
      const d = (by - cy) * (ax - cx) + (cx - bx) * (ay - cy); if (Math.abs(d) < 1e-14) continue;
      const l1 = ((by - cy) * (x - cx) + (cx - bx) * (y - cy)) / d, l2 = ((cy - ay) * (x - cx) + (ax - cx) * (y - cy)) / d, l3 = 1 - l1 - l2;
      if (l1 < 0 || l2 < 0 || l3 < 0) continue;
      hits.push(l1 * fpos[o + 2] + l2 * fpos[o + 5] + l3 * fpos[o + 8]);
    }
    hits.sort((a, b) => a - b);
    const iv = []; for (let h = 0; h + 1 < hits.length; h += 2) iv.push([hits[h], hits[h + 1]]);
    colIntervals[j * nx + i] = iv;
    for (let k = 0; k < nz; k++) {
      const z = stockMin[2] + (k + 0.5) * vz;
      const p = [x, y, z];
      if (!inStock(p)) { state[id(i, j, k)] = 0; continue; }
      let inside = false; for (const [a, b] of iv) if (z >= a && z <= b) { inside = true; break; }
      state[id(i, j, k)] = inside ? 2 : 1;
    }
  }
  return { nx, ny, nz, vx, vy, vz, min: stockMin, max: stockMax, state, id, colIntervals,
    center: (i, j, k) => [stockMin[0] + (i + 0.5) * vx, stockMin[1] + (j + 0.5) * vy, stockMin[2] + (k + 0.5) * vz] };
}

// mark hole voxels; holes given in frame coords {point, axis, r, t0, t1, hid}
function markHoles(G, holes) {
  const owner = new Int16Array(G.state.length).fill(-1);
  for (const h of holes) {
    const ends = [add(h.point, mul(h.axis, h.t0)), add(h.point, mul(h.axis, h.t1))];
    const lo = [0, 1, 2].map(a => Math.min(ends[0][a], ends[1][a]) - h.r), hi = [0, 1, 2].map(a => Math.max(ends[0][a], ends[1][a]) + h.r);
    const r0 = [Math.floor((lo[0] - G.min[0]) / G.vx), Math.floor((lo[1] - G.min[1]) / G.vy), Math.floor((lo[2] - G.min[2]) / G.vz)];
    const r1 = [Math.ceil((hi[0] - G.min[0]) / G.vx), Math.ceil((hi[1] - G.min[1]) / G.vy), Math.ceil((hi[2] - G.min[2]) / G.vz)];
    for (let k = Math.max(0, r0[2]); k < Math.min(G.nz, r1[2]); k++) for (let j = Math.max(0, r0[1]); j < Math.min(G.ny, r1[1]); j++) for (let i = Math.max(0, r0[0]); i < Math.min(G.nx, r1[0]); i++) {
      const q = G.id(i, j, k); if (G.state[q] !== 1) continue;
      const c = G.center(i, j, k), d = sub(c, h.point), t = dot(d, h.axis);
      if (t < h.t0 - 0.5 * G.vz || t > h.t1 + 0.5 * G.vz) continue;
      if (len(sub(d, mul(h.axis, t))) <= h.r) { G.state[q] = 3; owner[q] = h.hid; }
    }
  }
  G.holeOwner = owner;
}

// accessibility: bit flags 1:+Z 2:-Z 4:+X 8:-X 16:+Y 32:-Y
const DIRS = [{ k: '+Z', a: 2, s: 1, bit: 1 }, { k: '-Z', a: 2, s: -1, bit: 2 }, { k: '+X', a: 0, s: 1, bit: 4 }, { k: '-X', a: 0, s: -1, bit: 8 }, { k: '+Y', a: 1, s: 1, bit: 16 }, { k: '-Y', a: 1, s: -1, bit: 32 }];
function accessibility(G) {
  const acc = new Uint8Array(G.state.length);
  const n = [G.nx, G.ny, G.nz];
  for (const D of DIRS) {
    const a = D.a, b = (a + 1) % 3, c = (a + 2) % 3;
    for (let ib = 0; ib < n[b]; ib++) for (let ic = 0; ic < n[c]; ic++) {
      let blocked = false;
      for (let s = 0; s < n[a]; s++) {
        const ia = D.s > 0 ? n[a] - 1 - s : s;
        const ijk = [0, 0, 0]; ijk[a] = ia; ijk[b] = ib; ijk[c] = ic;
        const q = G.id(ijk[0], ijk[1], ijk[2]);
        if (G.state[q] === 2) blocked = true;
        else if (!blocked && G.state[q] !== 0) acc[q] |= D.bit;
      }
    }
  }
  return acc;
}

// assign removable voxels to setups (greedy, in priority order); returns setups with voxel lists ordered for tool path
function assignSetups(G, acc, priority) {
  const owner = new Int8Array(G.state.length).fill(-1);
  const used = [];
  // first pass: count reachable per dir
  const remaining = [];
  for (let q = 0; q < G.state.length; q++) if (G.state[q] === 1) remaining.push(q);
  let left = remaining;
  const order = priority.slice();
  const setups = [];
  // always take directions in priority while they remove something; side dirs only if voxels unreachable from earlier dirs
  for (const dk of order) {
    const D = DIRS.find(d => d.k === dk);
    const mine = [], rest = [];
    for (const q of left) (acc[q] & D.bit ? mine : rest).push(q);
    if (mine.length) { setups.push({ dir: dk, vox: mine }); mine.forEach(q => owner[q] = setups.length - 1); }
    left = rest;
  }
  return { setups, unreachable: left, owner };
}

// order voxels for a setup: layer by layer from the tool side, serpentine within layer
function orderPath(G, vox, dir) {
  const D = DIRS.find(d => d.k === dir); const a = D.a, b = (a + 1) % 3, c = (a + 2) % 3;
  const n = [G.nx, G.ny, G.nz];
  const coords = vox.map(q => { const i = q % G.nx, j = Math.floor(q / G.nx) % G.ny, k = Math.floor(q / (G.nx * G.ny)); return [q, [i, j, k]]; });
  coords.sort((A, B) => {
    const la = D.s > 0 ? n[a] - 1 - A[1][a] : A[1][a], lb = D.s > 0 ? n[a] - 1 - B[1][a] : B[1][a];
    if (la !== lb) return la - lb;
    if (A[1][b] !== B[1][b]) return A[1][b] - B[1][b];
    const flip = A[1][b] % 2 ? -1 : 1;
    return flip * (A[1][c] - B[1][c]);
  });
  return coords.map(x => x[0]);
}

// ---------- default rates (India, Oct 2026) ----------
const DEFAULT_RATES = {
  updated: '2026-10-23',
  materials: [
    // price ₹/kg for bar and plate, density g/cc, machinability (time factor vs Al 6061), laser/waterjet factor, laser max thickness mm, finish options
    { id: 'al6061', name: 'Aluminium 6061-T6', group: 'Aluminium', bar: 330, plate: 340, density: 2.70, mi: 1.0, cutF: 1.4, laserMax: 15 },
    { id: 'al6082', name: 'Aluminium 6082-T6', group: 'Aluminium', bar: 320, plate: 330, density: 2.71, mi: 1.0, cutF: 1.4, laserMax: 15 },
    { id: 'al7075', name: 'Aluminium 7075-T651', group: 'Aluminium', bar: 650, plate: 620, density: 2.81, mi: 1.15, cutF: 1.4, laserMax: 12 },
    { id: 'ms', name: 'Mild steel IS 2062 / EN1A', group: 'Steel', bar: 72, plate: 68, density: 7.85, mi: 2.0, cutF: 1.0, laserMax: 25 },
    { id: 'en8', name: 'EN8 (C45) carbon steel', group: 'Steel', bar: 75, plate: 80, density: 7.85, mi: 2.3, cutF: 1.0, laserMax: 25 },
    { id: 'en24', name: 'EN24 alloy steel', group: 'Steel', bar: 100, plate: 110, density: 7.85, mi: 3.0, cutF: 1.05, laserMax: 20 },
    { id: 'ss304', name: 'Stainless steel 304', group: 'Stainless', bar: 230, plate: 220, density: 8.0, mi: 3.6, cutF: 1.6, laserMax: 20 },
    { id: 'ss316', name: 'Stainless steel 316', group: 'Stainless', bar: 330, plate: 320, density: 8.0, mi: 4.0, cutF: 1.6, laserMax: 20 },
    { id: 'brass', name: 'Brass IS 319 (free cutting)', group: 'Copper alloy', bar: 650, plate: 680, density: 8.5, mi: 0.8, cutF: 2.0, laserMax: 6 },
    { id: 'cu', name: 'Copper ETP', group: 'Copper alloy', bar: 1400, plate: 1400, density: 8.94, mi: 1.6, cutF: 2.2, laserMax: 6 },
    { id: 'ti5', name: 'Titanium Grade 5 (Ti-6Al-4V)', group: 'Titanium', bar: 2200, plate: 2400, density: 4.43, mi: 5.5, cutF: 2.5, laserMax: 0 },
    { id: 'pa6', name: 'Cast nylon PA6', group: 'Engineering plastic', bar: 360, plate: 360, density: 1.15, mi: 0.6, cutF: 1.0, laserMax: 0, plastic: true },
    { id: 'pom', name: 'Acetal / Delrin (POM-C)', group: 'Engineering plastic', bar: 320, plate: 340, density: 1.41, mi: 0.5, cutF: 1.0, laserMax: 0, plastic: true }
  ],
  machines: { // ₹ per hour
    vmc3: { name: 'VMC 3-axis', rate: 450 },
    vmc4: { name: 'VMC 4-axis', rate: 650 },
    vmc5: { name: '5-axis machining centre', rate: 1500 },
    lathe: { name: 'CNC lathe', rate: 350 },
    turnmill: { name: 'Turn-mill (live tooling)', rate: 750 },
    labour: { name: 'Deburr / inspection labour', rate: 250 },
    program: { name: 'CAM programming', rate: 600 }
  },
  cutting: {
    sawPerCut: 30, sawPerCm2: 0.5,
    laserBasePerM: 25, laserPerMmPerM: 10, laserPierce: 3,
    waterjetBasePerM: 60, waterjetPerMmPerM: 18,
    wireEdmPerMm2: 0.18, wireEdmStart: 40
  },
  process: { // base values for Aluminium 6061 (scaled by material machinability)
    millMRR: 30, millFinish: 45, turnMRR: 55, turnFinish: 90, drillFeed: 250, tapPerHole: 0.6, reamPerHole: 0.5,
    toolChange: 0.15, loadVmc: 3, loadLathe: 1.5,
    setupVmc: 45, setupLathe: 30, setupEdm: 30, setupLaser: 15, programPerSetup: 30,
    deburrBase: 4, deburrPerFeature: 0.4, freeformFactor: 3.0, partOffMin: 0.4
  },
  overheads: { overheadPct: 15, marginPct: 20, gstPct: 18 },
  finishes: [
    { id: 'none', name: 'As machined', perDm2: 0, perKg: 0, min: 0 },
    { id: 'anodise', name: 'Clear anodise (Al)', perDm2: 8, perKg: 0, min: 300 },
    { id: 'hardanod', name: 'Hard anodise (Al)', perDm2: 20, perKg: 0, min: 600 },
    { id: 'zinc', name: 'Zinc plating (steel)', perDm2: 0, perKg: 60, min: 250 },
    { id: 'blackox', name: 'Black oxide (steel)', perDm2: 0, perKg: 80, min: 250 },
    { id: 'passiv', name: 'Passivation (SS)', perDm2: 0, perKg: 100, min: 300 },
    { id: 'powder', name: 'Powder coat', perDm2: 10, perKg: 0, min: 400 }
  ]
};

const TOL = {
  general: { name: 'General (ISO 2768-m, ±0.1)', finish: 1.0, holes: 'drill', inspect: 1.0 },
  precision: { name: 'Precision (±0.02)', finish: 1.4, holes: 'ream', inspect: 1.5 },
  tight: { name: 'Tight (±0.01)', finish: 1.9, holes: 'ream', inspect: 2.2 }
};

const TAP_DRILL = [[2.5, 'M3'], [3.3, 'M4'], [4.2, 'M5'], [5.0, 'M6'], [6.8, 'M8'], [8.5, 'M10'], [10.2, 'M12'], [14.0, 'M16'], [17.5, 'M20']];
function likelyTap(d) { const m = TAP_DRILL.find(([dd]) => Math.abs(dd - d) < 0.12); return m ? m[1] : null; }

// ---------- study: everything geometric, independent of material & quantity ----------
function study(T, opts = {}) {
  const A = analyzePart(T);
  const rot = rotationalTest(A);
  const res = { A, rot };
  res.prismatic = studyPrismatic(A);
  if (rot.isRot) res.turned = studyTurned(A, rot);
  res.features = featureSummary(A);
  return res;
}

function featureSummary(A) {
  const holes = A.stacks.map((s, i) => ({ i, d: s.dMin, dMax: s.dMax, depth: s.depth, cutDepth: s.cutDepth, segs: s.segs.length, interrupted: s.interrupted, through: s.through, cbore: s.counterbored, tap: likelyTap(s.dMin), axis: s.axis }));
  return {
    faces: A.faces.length,
    planar: A.planes.length, cyl: A.faces.filter(f => f.type === 'cyl').length, cone: A.cones.length, free: A.free.length,
    holes, minCornerR: A.cornerR.length ? Math.min(...A.cornerR) : null, slots: A.slots.length,
    freeFrac: A.freeArea / A.area
  };
}

function studyPrismatic(A) {
  // try each axis as the machine Z and keep the orientation with the lowest machining-effort proxy
  let best = null;
  for (let k = 0; k < 3; k++) {
    const r = studyPrismaticK(A, k);
    if (!best || r.effort < best.effort) best = r;
  }
  return best;
}

function studyPrismaticK(A, bestK) {
  const { bbox, stacks, planes } = A;
  const perm = bestK === 2 ? [0, 1, 2] : bestK === 0 ? [1, 2, 0] : [2, 0, 1];
  const R = perm.map(i => AXES[i]);
  const F = makeFrame(R, bbox.mn);
  const fpos = frameTris(A.T, F);
  const sz = [bbox.size[perm[0]], bbox.size[perm[1]], bbox.size[perm[2]]];
  // block stock: +2 mm per side on X/Y (rounded up to 5), thickness next std plate ≥ T + 2
  const stock = [roundUp(sz[0] + 4, 5), roundUp(sz[1] + 4, 5), nextStd(PLATE_T, sz[2] + 2)];
  const off = [(stock[0] - sz[0]) / 2, (stock[1] - sz[1]) / 2, (stock[2] - sz[2]) / 2];
  const sMin = [-off[0], -off[1], -off[2]], sMax = [sz[0] + off[0], sz[1] + off[1], sz[2] + off[2]];
  // holes in frame
  const holes = [];
  stacks.forEach((s, hid) => {
    for (const c of s.segs) {
      const p = toFrame(F, c.point), u = dirToFrame(F, s.axis);
      holes.push({ hid, point: p, axis: u, r: c.r, t0: c.t0, t1: c.t1 });
    }
  });
  const G = voxelize(fpos, sMin, sMax, { res: 64 });
  markHoles(G, holes);
  const acc = accessibility(G);
  // hole directions in frame
  const holeInfo = stacks.map((s, hid) => {
    const u = dirToFrame(F, s.axis); const ak = [0, 1, 2].find(a => Math.abs(u[a]) > 0.999);
    const end0 = toFrame(F, add(s.segs[0].point, mul(s.axis, s.t0))), end1 = toFrame(F, add(s.segs[0].point, mul(s.axis, s.t1)));
    // opening: end without cap. from end → tool comes from outside along -u (at end1) or +u (at end0)
    const opens = [];
    if (!s.cap1) opens.push({ at: end1, toward: mul(u, -1) }); // tool moves along -u entering at end1; tool side = +u
    if (!s.cap0) opens.push({ at: end0, toward: u });
    let dirs = opens.map(o => { const side = mul(o.toward, -1); const k = ak === undefined ? null : (side[ak] > 0 ? '+' : '-') + AXN[ak]; return { ...o, dir: k }; });
    return { hid, u, ak, dirs, d: s.dMin, dMax: s.dMax, depth: s.depth, cutDepth: s.cutDepth, cbore: s.counterbored, through: s.through, segs: s.segs.map(c => ({ r: c.r, t0: c.t0, t1: c.t1 })), point: toFrame(F, s.segs[0].point) };
  });
  // setups: priority +Z, -Z, sides
  const pri = ['+Z', '-Z', '+X', '-X', '+Y', '-Y'];
  const asg = assignSetups(G, acc, pri);
  const setupDirs = new Set(asg.setups.map(s => s.dir));
  // holes choose setup
  const angledHoles = [];
  for (const h of holeInfo) {
    if (!h.dirs.length) { h.dirs = [{ at: null, dir: null }]; }
    let pick = h.dirs.find(d => d.dir && setupDirs.has(d.dir)) || h.dirs.find(d => d.dir) || null;
    if (!pick) { angledHoles.push(h); h.setup = 'angled'; continue; }
    h.setup = pick.dir; h.entry = pick;
    setupDirs.add(pick.dir);
  }
  const order = pri.filter(d => setupDirs.has(d));
  const vol = (n) => n * G.vx * G.vy * G.vz;
  const setups = order.map(d => {
    const s = asg.setups.find(x => x.dir === d);
    let reach = 1;
    if (s && s.vox.length) {
      const D = DIRS.find(x => x.k === d), a = D.a, nA = [G.nx, G.ny, G.nz][a], vA = [G.vx, G.vy, G.vz][a];
      let f = 0; for (const q of s.vox) { const ijk = [q % G.nx, Math.floor(q / G.nx) % G.ny, Math.floor(q / (G.nx * G.ny))]; const depth = (D.s > 0 ? nA - 1 - ijk[a] : ijk[a]) * vA; f += 1 + Math.min(5, Math.pow(depth / 40, 2)); }
      reach = f / s.vox.length;
    }
    return { dir: d, voxels: s ? orderPath(G, s.vox, d) : [], volume: s ? vol(s.vox.length) : 0, reach, holes: holeInfo.filter(h => h.setup === d).map(h => h.hid) };
  });
  // drop voxel-noise setups (faceting slivers) that carry no holes; fold their voxels into the last kept setup
  const noiseCut = Math.max(0.003 * A.volume, 4 * G.vx * G.vy * G.vz);
  const kept = setups.filter((x, i) => i < 2 || x.holes.length || x.volume > noiseCut);
  const dropped = setups.filter(x => !kept.includes(x));
  if (dropped.length && kept.length) { const tgt = kept[kept.length - 1]; for (const d of dropped) { tgt.voxels = tgt.voxels.concat(d.voxels); tgt.volume += d.volume; } }
  setups.length = 0; setups.push(...kept);
  // correct volumes to exact removal
  const stockVol = stock[0] * stock[1] * stock[2];
  const holeVol = holeInfo.reduce((s, h) => s + h.segs.reduce((a, c) => a + Math.PI * c.r * c.r * (c.t1 - c.t0), 0), 0);
  const exactRemove = Math.max(0, stockVol - A.volume - holeVol);
  const voxRemove = setups.reduce((s, x) => s + x.volume, 0) + vol(asg.unreachable.length);
  const scale = voxRemove > 0 ? exactRemove / voxRemove : 0;
  setups.forEach(s => s.volume *= scale);
  const unreachableVol = vol(asg.unreachable.length) * scale;

  // silhouette & 2.5D through-profile test (for laser / EDM routes)
  let silCells = 0, boundary = 0, fullBoundary = 0; const fullT = sz[2];
  const silAt = (i, j) => i >= 0 && j >= 0 && i < G.nx && j < G.ny && G.colIntervals[j * G.nx + i].length > 0;
  for (let j = 0; j < G.ny; j++) for (let i = 0; i < G.nx; i++) {
    if (!silAt(i, j)) continue; silCells++;
    const edge = !silAt(i + 1, j) || !silAt(i - 1, j) || !silAt(i, j + 1) || !silAt(i, j - 1);
    if (edge) { boundary++; const iv = G.colIntervals[j * G.nx + i]; const span = iv.reduce((s, [a, b]) => s + (b - a), 0); if (span > 0.97 * fullT) fullBoundary++; }
  }
  let edgeLen = 0;
  for (let j = 0; j < G.ny; j++) for (let i = 0; i < G.nx; i++) { if (!silAt(i, j)) continue; if (!silAt(i + 1, j)) edgeLen += G.vy; if (!silAt(i - 1, j)) edgeLen += G.vy; if (!silAt(i, j + 1)) edgeLen += G.vx; if (!silAt(i, j - 1)) edgeLen += G.vx; }
  const perimeter = edgeLen * Math.PI / 4;
  const silArea = silCells * G.vx * G.vy;
  let walls25 = true;
  for (const p of planes) { const n = dirToFrame(F, p.normal); if (!(Math.abs(Math.abs(n[2]) - 1) < 2e-3 || Math.abs(n[2]) < 2e-3)) { walls25 = false; break; } }
  for (const c of A.cyls) { const u = dirToFrame(F, c.axis); if (Math.abs(Math.abs(u[2]) - 1) > 1e-3) { walls25 = false; break; } }
  if (A.freeArea / A.area > 0.03) walls25 = false;
  const throughProfile = walls25 && boundary > 0 && fullBoundary / boundary > 0.9;
  const tilted = planes.filter(p => { const n = dirToFrame(F, p.normal); return n.every(x => Math.abs(x) < 0.999) && Math.abs(n[2]) > 2e-3; });
  const tiltedArea = tilted.reduce((s, p) => s + p.area, 0);

  // finishing area per setup: part area split by normal direction affinity
  const finishArea = A.area;

  const effort = setups.reduce((a, x) => a + x.volume / 1000 / 30 * x.reach, 0) + setups.length * 6 + angledHoles.length * 10 + unreachableVol / 1000 / 5;
  return { effort, upAxis: AXN[bestK], F, perm, R, size: sz, stock, sMin, sMax, G, holes, holeInfo, setups, unreachable: asg.unreachable, unreachableVol, angledHoles,
    removal: exactRemove, stockVol, perimeter, silArea, throughProfile, walls25, tiltedArea, finishArea };
}

function studyTurned(A, rot) {
  const { T } = A;
  const u = rot.axis;
  const e2 = anyPerp(u), e3 = cross(u, e2);
  const o = add(rot.point, mul(u, rot.tMin));
  const F = makeFrame([u, e2, e3], o);
  const L = rot.L;
  const NB = Math.max(60, Math.min(400, Math.round(L / 0.25)));
  const db = L / NB;
  const rOut = new Float64Array(NB), rIn = new Float64Array(NB).fill(Infinity);
  const fpos = frameTris(T, F);
  for (let t = 0; t < T.n; t++) {
    const o9 = t * 9;
    const P = [[fpos[o9], fpos[o9 + 1], fpos[o9 + 2]], [fpos[o9 + 3], fpos[o9 + 4], fpos[o9 + 5]], [fpos[o9 + 6], fpos[o9 + 7], fpos[o9 + 8]]];
    const n = dirToFrame(F, triN(T, t));
    const c = [(P[0][0] + P[1][0] + P[2][0]) / 3, (P[0][1] + P[1][1] + P[2][1]) / 3, (P[0][2] + P[1][2] + P[2][2]) / 3];
    const rc = Math.hypot(c[1], c[2]); const radialN = rc > 1e-9 ? (n[1] * c[1] + n[2] * c[2]) / rc : 0;
    const x0 = Math.min(P[0][0], P[1][0], P[2][0]), x1 = Math.max(P[0][0], P[1][0], P[2][0]);
    const b0 = Math.max(0, Math.floor(x0 / db)), b1 = Math.min(NB - 1, Math.floor(x1 / db));
    for (let b = b0; b <= b1; b++) {
      const xc = Math.min(Math.max((b + 0.5) * db, x0), x1);
      // intersect triangle with plane x = xc → segment endpoints
      const pts = [];
      for (let e = 0; e < 3; e++) {
        const p = P[e], q = P[(e + 1) % 3];
        if ((p[0] - xc) * (q[0] - xc) <= 0 && p[0] !== q[0]) { const s = (xc - p[0]) / (q[0] - p[0]); pts.push([p[1] + s * (q[1] - p[1]), p[2] + s * (q[2] - p[2])]); }
        else if (p[0] === xc) pts.push([p[1], p[2]]);
      }
      if (!pts.length) continue;
      for (const [y, z] of pts) {
        const r = Math.hypot(y, z);
        if (r > rOut[b]) rOut[b] = r;
        if (radialN < -0.3 && r < rIn[b]) rIn[b] = r;
      }
    }
  }
  for (let b = 0; b < NB; b++) if (!isFinite(rIn[b]) || rIn[b] >= rOut[b] - 1e-6) rIn[b] = 0;
  // fill gaps in rOut
  for (let b = 1; b < NB; b++) if (rOut[b] === 0) rOut[b] = rOut[b - 1];
  const rMax = Math.max(...rOut);
  const barD = nextStd(BAR_D, 2 * rMax + (2 * rMax < 50 ? 2 : 3));
  const barFed = barD <= 52;
  const cutLen = L + 2 + 2 + (barFed ? 3 : 0) + (barFed ? 0 : 15);
  const R = barD / 2;
  let envVol = 0, boreVol = 0, odRemove = 0;
  for (let b = 0; b < NB; b++) { envVol += Math.PI * (rOut[b] ** 2 - rIn[b] ** 2) * db; boreVol += Math.PI * rIn[b] ** 2 * db; odRemove += Math.PI * (R * R - rOut[b] ** 2) * db; }
  const faceRemove = Math.PI * R * R * 4; // two 2 mm facing passes
  const stockVol = Math.PI * R * R * cutLen;
  // secondary features: anything between envelope and the part, excluding coaxial holes
  const sMin = [0, -R, -R], sMax = [L, R, R];
  const inEnv = p => { const b = Math.min(NB - 1, Math.max(0, Math.floor(p[0] / db))); const r = Math.hypot(p[1], p[2]); return r <= rOut[b] && r >= rIn[b]; };
  const G = voxelize(fpos, sMin, sMax, { res: 56, inStock: inEnv });
  // off-axis holes
  const holes = [], holeInfo = [];
  A.stacks.forEach((s, hid) => {
    const uu = dirToFrame(F, s.axis);
    const coax = Math.abs(Math.abs(uu[0]) - 1) < 1e-3 && lineDist(toFrame(F, s.point), [1, 0, 0], [0, 0, 0]) < 0.02 * s.dMin + 0.02;
    if (coax) return;
    for (const c of s.segs) holes.push({ hid, point: toFrame(F, c.point), axis: uu, r: c.r, t0: c.t0, t1: c.t1 });
    const end0 = toFrame(F, add(s.segs[0].point, mul(s.axis, s.t0))), end1 = toFrame(F, add(s.segs[0].point, mul(s.axis, s.t1)));
    holeInfo.push({ hid, u: uu, d: s.dMin, dMax: s.dMax, depth: s.depth, cutDepth: s.cutDepth, cbore: s.counterbored, through: s.through, entry: !s.cap1 ? { at: end1, toward: mul(uu, -1) } : { at: end0, toward: uu }, segs: s.segs.map(c => ({ r: c.r, t0: c.t0, t1: c.t1 })), point: toFrame(F, s.segs[0].point) });
  });
  markHoles(G, holes);
  const acc = accessibility(G);
  const asg = assignSetups(G, acc, ['+Z', '-Z', '+Y', '-Y', '+X', '-X']);
  const vol = (n) => n * G.vx * G.vy * G.vz;
  const secSetups = asg.setups.map(s => ({ dir: s.dir, voxels: orderPath(G, s.vox, s.dir), volume: vol(s.vox.length) }));
  const secVolume = Math.max(0, envVol - A.volume - holes.reduce((s, h) => s + Math.PI * h.r * h.r * (h.t1 - h.t0), 0));
  const vsum = secSetups.reduce((s, x) => s + x.volume, 0);
  if (vsum > 0) secSetups.forEach(s => s.volume *= Math.max(secVolume, 0) / vsum);
  // drop tiny secondary setups (voxel noise from faceting) when less than 0.5% of part volume
  const secReal = secSetups.filter(s => s.volume > 0.005 * A.volume && s.voxels.length > 3);
  // OD profile length (for finishing)
  let odArea = 0, idArea = 0;
  for (let b = 0; b < NB; b++) { odArea += 2 * Math.PI * rOut[b] * db; if (rIn[b] > 0) idArea += 2 * Math.PI * rIn[b] * db; }
  // count steps (shoulders)
  let steps = 0; for (let b = 1; b < NB; b++) if (Math.abs(rOut[b] - rOut[b - 1]) > 0.2) steps++;
  const hasBore = boreVol > 1;
  const throughBore = hasBore && rIn[0] > 0 && rIn[NB - 1] > 0;
  return { F, L, NB, db, rOut: Array.from(rOut), rIn: Array.from(rIn), rMax, barD, barFed, cutLen, stockVol, envVol, odRemove, boreVol, faceRemove, G, holes, holeInfo,
    secSetups: secReal, odArea, idArea, steps, hasBore, throughBore, minBore: Math.min(...Array.from(rIn).filter(r => r > 0)) * 2 };
}

// ---------- planning & costing ----------
function plan(S, matId, q, rates, options = {}) {
  rates = rates || DEFAULT_RATES;
  const mat = rates.materials.find(m => m.id === matId) || rates.materials[0];
  const tol = Object.assign({}, TOL[options.tol || 'general'], { noTap: !!options.noTap });
  const fin = rates.finishes.find(f => f.id === (options.finish || 'none')) || rates.finishes[0];
  const routes = [];
  routes.push(routeBlockVMC(S, mat, rates, tol, 'vmc3'));
  const P = S.prismatic;
  if (P.setups.length >= 3 || P.unreachableVol > 0.02 * S.A.volume || P.angledHoles.length || P.tiltedArea > 0.05 * S.A.area || S.features.freeFrac > 0.1) routes.push(routeBlockVMC(S, mat, rates, tol, 'vmc5'));
  if (P.throughProfile) {
    const r = routeProfileCut(S, mat, rates, tol); if (r) routes.push(r);
    routes.push(routeWireEdm(S, mat, rates, tol));
  }
  if (S.turned) {
    routes.push(routeTurned(S, mat, rates, tol, false));
    if (S.turned.secSetups.length || S.turned.holeInfo.length) routes.push(routeTurned(S, mat, rates, tol, true));
  }
  for (const r of routes) costRoute(r, S, mat, q, rates, tol, fin);
  routes.sort((a, b) => a.cost.perPart - b.cost.perPart);
  routes.forEach((r, i) => r.rank = i);
  return { mat, tol, fin, q, routes, best: routes[0] };
}

const mins = (vol_mm3, rate_cm3) => vol_mm3 / 1000 / rate_cm3;

function holeOps(holes, mat, rates, tol, machine, setupLabel) {
  const pr = rates.process; const ops = [];
  for (const h of holes) {
    const d = h.d, depth = h.cutDepth || h.depth;
    const tap = likelyTap(d);
    const ld = depth / d;
    const big = d > 16 && ld <= 3;          // short large hole: helical interpolation
    const deepBig = d > 16 && ld > 3;       // long large hole: drill + bore/ream
    const feed = pr.drillFeed / mat.mi * (d < 4 ? 0.6 : 1);
    const peck = ld > 8 ? 2.2 : ld > 5 ? 1.7 : ld > 3 ? 1.35 : 1;
    let t, label;
    if (big) {
      const v = h.segs.reduce((s, c) => s + Math.PI * c.r * c.r * (c.t1 - c.t0), 0);
      t = mins(v, pr.millMRR / mat.mi * 0.5) + 0.3; label = `Interpolate bore Ø${d.toFixed(1)} × ${depth.toFixed(1)}`;
    } else if (deepBig) {
      t = (h.depth + 5) / (feed * 0.7) * peck + (h.depth + 5) / (feed * 1.5) + 2 * pr.toolChange;
      label = `Drill Ø${(Math.floor(d) - 0.5).toFixed(1)} + bore/ream Ø${d.toFixed(2)} × ${h.depth.toFixed(1)}`;
    } else {
      t = (h.depth + 3) / feed * peck + 0.15; label = h.segs.length > 1 && !h.cbore ? `Drill Ø${d.toFixed(1)} through ${h.segs.length} walls (${h.segs.map(c => (c.t1 - c.t0).toFixed(0)).join(' + ')} mm cut, ${h.depth.toFixed(0)} mm reach)` : `Drill Ø${d.toFixed(1)} × ${h.depth.toFixed(1)}`;
      if (h.cbore) t += 0.4;
    }
    if (tol.holes === 'ream' && !tap && !deepBig) t += pr.reamPerHole;
    label += (h.through ? ' thru' : '') + (ld > 5 ? ` (L/D ${ld.toFixed(1)}, peck)` : '');
    ops.push({ kind: big ? 'bore' : 'drill', machine, setup: setupLabel, hid: h.hid, d, depth: h.depth, min: t, label });
    if (tap && !tol.noTap) ops.push({ kind: 'tap', machine, setup: setupLabel, hid: h.hid, d, depth, min: pr.tapPerHole * Math.sqrt(mat.mi), label: `Tap ${tap} (tap-drill size detected)`, optional: true });
  }
  return ops;
}

function routeBlockVMC(S, mat, rates, tol, mc) {
  const P = S.prismatic, pr = rates.process;
  const five = mc === 'vmc5';
  const r = { id: five ? 'block-5ax' : 'block-vmc', name: five ? 'Block → 5-axis machining' : 'Block → VMC 3-axis', kind: 'prismatic', machineKey: mc, ops: [], notes: [] };
  r.stock = { form: 'plate', desc: `Plate cut to ${P.stock[0]} × ${P.stock[1]} × ${P.stock[2]} mm`, dims: P.stock, volume: P.stockVol };
  r.ops.push({ kind: 'saw', machine: 'saw', label: `Saw-cut blank ${P.stock[0]} × ${P.stock[1]} from ${P.stock[2]} mm plate`, cutArea: (P.stock[1] * P.stock[2] + P.stock[0] * P.stock[2]) / 100 });
  let setups = P.setups.map(s => ({ ...s }));
  if (five) {
    // 3+2 positioning: all non-bottom setups merge into op 1, bottom in op 2
    const top = setups.filter(s => s.dir !== '-Z'), bot = setups.filter(s => s.dir === '-Z');
    const merged = [];
    if (top.length) merged.push({ dir: '+Z', reach: Math.min(...top.map(s => s.reach || 1)), merged: top.map(s => s.dir), voxels: top.flatMap(s => s.voxels), volume: top.reduce((a, s) => a + s.volume, 0) + P.unreachableVol, holes: top.flatMap(s => s.holes).concat(P.angledHoles.map(h => h.hid)) });
    if (bot.length) merged.push(...bot);
    setups = merged;
  } else if (P.unreachableVol > 0.02 * S.A.volume) r.notes.push(`${(P.unreachableVol / 1000).toFixed(1)} cm³ is not reachable from any of the six axis directions (undercut). A 5-axis or EDM step is needed for that region.`);
  if (!five && P.angledHoles.length) r.notes.push(`${P.angledHoles.length} hole(s) are at a compound angle; a fixture or 4/5-axis is needed.`);
  r.setups = setups.map((s, i) => ({ n: i + 1, dir: s.dir, merged: s.merged }));
  const totalFinishArea = P.finishArea;
  setups.forEach((s, i) => {
    const lbl = `Setup ${i + 1} (${five && s.merged ? '5-axis' : 'tool from ' + s.dir})`;
    if (s.volume > Math.max(500, 0.002 * S.A.volume) || (i < 2 && s.volume > 0)) {
      r.ops.push({ kind: 'rough', machine: mc, setup: i, dir: s.dir, label: `Face + rough mill ${(s.volume / 1000).toFixed(1)} cm³`, min: mins(s.volume, pr.millMRR / mat.mi) * (s.reach || 1) + pr.toolChange, volume: s.volume, voxels: s.voxels });
      const wts = setups.map(x => (x.dir === '+Z' || x.dir === '-Z' || x.merged) ? 1 : (x.volume > 0.05 * P.removal ? 0.6 : 0.15));
      const share = totalFinishArea * wts[i] / wts.reduce((a, b) => a + b, 0);
      const freeExtra = S.features.freeFrac * share * (pr.freeformFactor - 1);
      r.ops.push({ kind: 'finish', machine: mc, setup: i, dir: s.dir, label: `Finish walls & floors${S.features.minCornerR ? `, Ø${Math.max(2, Math.min(12, 2 * S.features.minCornerR * 0.8)).toFixed(0)} end mill for corners` : ''}`, min: (share + freeExtra) / 100 / (pr.millFinish / Math.pow(mat.mi, 0.7)) * tol.finish + pr.toolChange });
    }
    const hs = P.holeInfo.filter(h => s.holes.includes(h.hid));
    r.ops.push(...holeOps(hs, mat, rates, tol, mc, i));
    r.ops.forEach(o => { if (o.setup === i && !o.setupLabel) o.setupLabel = lbl; });
  });
  const deep = P.holeInfo.filter(h => (h.cutDepth || h.depth) / h.d > 5);
  if (deep.length) r.notes.push(`${deep.length} hole(s) are deeper than 5×D (up to ${Math.max(...deep.map(h => h.depth / h.d)).toFixed(1)}×D). Use long-series drills with peck cycles, or drill halfway from each side.`);
  if (S.features.minCornerR && S.features.minCornerR < 1) r.notes.push(`Smallest internal corner radius is R${S.features.minCornerR.toFixed(2)}; needs a very small cutter or wire EDM for sharp corners.`);
  return r;
}

function routeProfileCut(S, mat, rates, tol) {
  const P = S.prismatic, pr = rates.process;
  const T = P.size[2];
  const std = nextStd(PLATE_T, T);
  const exact = Math.abs(std - T) < 0.05;
  const laser = mat.laserMax >= std;
  if (!laser && std > 100) return null;
  const r = { id: 'profile', name: (laser ? 'Laser' : 'Waterjet') + ' profile → VMC finish', kind: 'prismatic', machineKey: 'vmc3', ops: [], notes: [] };
  const bw = P.size[0] + 10, bh = P.size[1] + 10;
  r.stock = { form: 'plate', desc: `${std} mm plate, nested area ${bw.toFixed(0)} × ${bh.toFixed(0)} mm`, dims: [bw, bh, std], volume: bw * bh * std };
  // through-Z holes laser-cut if d >= max(T,3) and general tolerance
  const laserHoles = [], vmcHoles = [];
  for (const h of P.holeInfo) {
    const zAx = h.ak === 2;
    if (zAx && h.through && h.segs.length === 1 && h.d >= Math.max(std, 3) && tol.holes === 'drill' && laser) laserHoles.push(h); else vmcHoles.push(h);
  }
  const perim = P.perimeter + laserHoles.reduce((s, h) => s + Math.PI * h.d, 0);
  const c = rates.cutting;
  const perM = laser ? (c.laserBasePerM + c.laserPerMmPerM * std) * mat.cutF : (c.waterjetBasePerM + c.waterjetPerMmPerM * std) * Math.max(1, mat.cutF * 0.8);
  r.ops.push({ kind: 'laser', machine: laser ? 'laser' : 'waterjet', label: `${laser ? 'Fibre laser' : 'Abrasive waterjet'} cut outer profile${laserHoles.length ? ` + ${laserHoles.length} hole(s)` : ''}, ${(perim / 1000).toFixed(2)} m`, cost: perim / 1000 * perM + (1 + laserHoles.length) * c.laserPierce, perimeter: perim, holes: laserHoles.map(h => h.hid) });
  // VMC: remaining voxels inside silhouette (pockets/steps) and faces if plate thicker than part
  const G = P.G; let inner = 0;
  const silAt = (i, j) => G.colIntervals[j * G.nx + i].length > 0;
  for (const s of P.setups) for (const qv of s.voxels) { const i = qv % G.nx, j = Math.floor(qv / G.nx) % G.ny; if (silAt(i, j)) inner++; }
  const totalVox = P.setups.reduce((a, s) => a + s.voxels.length, 0) || 1;
  const vRem = P.removal * inner / totalVox * (exact ? 0.6 : 1);
  const setupsNeeded = new Set(vmcHoles.map(h => h.setup)); if (vRem > 0.01 * S.A.volume || !exact) { setupsNeeded.add('+Z'); if (!exact) setupsNeeded.add('-Z'); }
  const order = ['+Z', '-Z', '+X', '-X', '+Y', '-Y'].filter(d => setupsNeeded.has(d));
  r.setups = order.map((d, i) => ({ n: i + 1, dir: d }));
  order.forEach((d, i) => {
    const lbl = `Setup ${i + 1} (tool from ${d})`;
    const share = d === '+Z' ? 0.7 : d === '-Z' ? 0.3 : 0;
    if (share > 0 && vRem > 0) r.ops.push({ kind: 'rough', machine: 'vmc3', setup: i, setupLabel: lbl, dir: d, label: `${exact ? 'Mill pockets/steps' : 'Face to thickness + mill pockets'} ${(vRem * share / 1000).toFixed(1)} cm³`, min: mins(vRem * share, pr.millMRR / mat.mi) + pr.toolChange });
    if (share > 0) r.ops.push({ kind: 'finish', machine: 'vmc3', setup: i, setupLabel: lbl, dir: d, label: 'Finish floors & cut edges', min: P.finishArea * share * 0.5 / 100 / (pr.millFinish / Math.pow(mat.mi, 0.7)) * tol.finish + pr.toolChange });
    r.ops.push(...holeOps(vmcHoles.filter(h => h.setup === d), mat, rates, tol, 'vmc3', i).map(o => ({ ...o, setupLabel: lbl })));
  });
  if (!order.length) r.notes.push('Part comes off the profile cutter finished; only deburring remains.');
  r.notes.push(laser ? 'Laser edges carry a heat-affected zone and slight taper; precision edges are re-milled.' : 'Waterjet edges show slight taper on thick plate.');
  return r;
}

function routeWireEdm(S, mat, rates, tol) {
  const P = S.prismatic, pr = rates.process, c = rates.cutting;
  const r = { id: 'wire-edm', name: 'Block → Wire EDM profile', kind: 'prismatic', machineKey: 'vmc3', ops: [], notes: [] };
  const T = P.stock[2];
  r.stock = { form: 'plate', desc: `Plate cut to ${P.stock[0]} × ${P.stock[1]} × ${T} mm`, dims: P.stock, volume: P.stockVol };
  r.ops.push({ kind: 'saw', machine: 'saw', label: `Saw-cut blank ${P.stock[0]} × ${P.stock[1]} from ${T} mm plate`, cutArea: (P.stock[1] * T + P.stock[0] * T) / 100 });
  r.setups = [{ n: 1, dir: '+Z' }, { n: 2, dir: '-Z' }, { n: 3, dir: 'EDM' }];
  r.ops.push({ kind: 'finish', machine: 'vmc3', setup: 0, setupLabel: 'Setup 1 (tool from +Z)', dir: '+Z', label: 'Face top to thickness', min: P.stock[0] * P.stock[1] / 100 / (pr.millFinish / Math.pow(mat.mi, 0.7)) + mins(P.stock[0] * P.stock[1] * (T - P.size[2]) / 2, pr.millMRR / mat.mi) });
  r.ops.push({ kind: 'finish', machine: 'vmc3', setup: 1, setupLabel: 'Setup 2 (tool from -Z)', dir: '-Z', label: 'Face bottom to thickness', min: P.stock[0] * P.stock[1] / 100 / (pr.millFinish / Math.pow(mat.mi, 0.7)) + mins(P.stock[0] * P.stock[1] * (T - P.size[2]) / 2, pr.millMRR / mat.mi) });
  const startHoles = P.holeInfo.filter(h => h.ak === 2 && h.through && h.d > 3).length;
  r.ops.push({ kind: 'drill', machine: 'vmc3', setup: 0, setupLabel: 'Setup 1 (tool from +Z)', label: `Drill ${startHoles} wire start hole(s)`, min: startHoles * ((P.size[2] + 3) / (rates.process.drillFeed / mat.mi) + 0.15) });
  const cutLen = P.perimeter + P.holeInfo.filter(h => h.ak === 2 && h.through).reduce((s, h) => s + Math.PI * h.d, 0);
  const mm2 = cutLen * P.size[2];
  r.ops.push({ kind: 'edm', machine: 'wedm', setup: 2, setupLabel: 'Setup 3 (wire EDM)', label: `Wire-cut profile + thru features, ${(mm2 / 100).toFixed(0)} cm² cut area`, cost: mm2 * c.wireEdmPerMm2 * (mat.mi > 3 ? 1.15 : 1) + (1 + startHoles) * c.wireEdmStart, cutArea: mm2, holes: P.holeInfo.filter(h => h.ak === 2 && h.through).map(h => h.hid) });
  const blind = P.holeInfo.filter(h => !(h.ak === 2 && h.through));
  r.ops.push(...holeOps(blind, mat, rates, tol, 'vmc3', 0).map(o => ({ ...o, setupLabel: 'Setup 1 (tool from +Z)' })));
  r.notes.push('Wire EDM holds ±0.005 mm and sharp internal corners; slowest option on thick sections.');
  return r;
}

function routeTurned(S, mat, rates, tol, turnMill) {
  const Tn = S.turned, pr = rates.process;
  const mc = turnMill ? 'turnmill' : 'lathe';
  const r = { id: turnMill ? 'turn-mill' : 'lathe', name: turnMill ? 'Bar → Turn-mill (one machine)' : 'Bar → CNC lathe' + (Tn.secSetups.length || Tn.holeInfo.length ? ' + VMC' : ''), kind: 'turned', machineKey: mc, ops: [], notes: [] };
  r.stock = { form: 'bar', desc: `Ø${Tn.barD} round bar, ${Tn.cutLen.toFixed(0)} mm per part${Tn.barFed ? ' (bar-fed)' : ' incl. 15 mm chuck grip'}`, dims: [Tn.barD, Tn.cutLen], volume: Tn.stockVol };
  if (!Tn.barFed) r.ops.push({ kind: 'saw', machine: 'saw', label: `Saw-cut Ø${Tn.barD} billet to ${Tn.cutLen.toFixed(0)} mm`, cutArea: Math.PI * (Tn.barD / 2) ** 2 / 100 });
  r.setups = [{ n: 1, dir: 'Lathe op 1' }];
  const L1 = 'Setup 1 (lathe chuck)';
  const mrr = pr.turnMRR / mat.mi;
  r.ops.push({ kind: 'face', machine: mc, setup: 0, setupLabel: L1, label: 'Face front', min: mins(Tn.faceRemove / 2, mrr) + 0.1 });
  r.ops.push({ kind: 'turn-rough', machine: mc, setup: 0, setupLabel: L1, label: `Rough turn OD${Tn.steps ? `, ${Tn.steps} shoulder(s)` : ''}, ${(Tn.odRemove / 1000).toFixed(1)} cm³`, min: mins(Tn.odRemove, mrr) + pr.toolChange });
  r.ops.push({ kind: 'turn-finish', machine: mc, setup: 0, setupLabel: L1, label: 'Finish turn OD profile', min: Tn.odArea / 100 / (pr.turnFinish / Math.pow(mat.mi, 0.7)) * tol.finish + pr.toolChange });
  if (Tn.hasBore) {
    const dDrill = Math.max(3, Math.min(Tn.minBore * 0.9, 30));
    r.ops.push({ kind: 'turn-drill', machine: mc, setup: 0, setupLabel: L1, label: `Centre drill + drill Ø${dDrill.toFixed(1)}${Tn.throughBore ? ' thru' : ''}`, min: (Tn.L + 3) / (pr.drillFeed / mat.mi * 0.8) + pr.toolChange, dDrill });
    r.ops.push({ kind: 'turn-bore', machine: mc, setup: 0, setupLabel: L1, label: `Bore ID, ${(Tn.boreVol / 1000).toFixed(1)} cm³`, min: mins(Math.max(0, Tn.boreVol - Math.PI * (dDrill / 2) ** 2 * Tn.L), mrr * 0.6) + Tn.idArea / 100 / (pr.turnFinish * 0.7 / Math.pow(mat.mi, 0.7)) * tol.finish + pr.toolChange });
  }
  if (Tn.barFed) r.ops.push({ kind: 'part', machine: mc, setup: 0, setupLabel: L1, label: 'Part off', min: pr.partOffMin * Math.sqrt(mat.mi) * Tn.barD / 25 });
  else { r.setups.push({ n: 2, dir: 'Lathe op 2 (flipped)' }); r.ops.push({ kind: 'face', machine: mc, setup: 1, setupLabel: 'Setup 2 (flip, face back)', label: 'Flip, face back & remove grip', min: mins(Math.PI * (Tn.barD / 2) ** 2 * 17, mrr) + 0.2 }); }
  const secMc = turnMill ? 'turnmill' : 'vmc3';
  let sIdx = r.setups.length;
  const secDirs = Tn.secSetups.map(s => s.dir);
  const holeDirs = Tn.holeInfo.map(h => { const s = mul(h.entry.toward, -1); const k = [1, 2].sort((a, b) => Math.abs(s[b]) - Math.abs(s[a]))[0]; return (s[k] > 0 ? '+' : '-') + AXN[k]; });
  Tn.holeInfo.forEach((h, i) => h.setupDir = holeDirs[i]);
  const dirs = [...new Set([...secDirs, ...holeDirs])];
  dirs.forEach((d) => {
    const lbl = turnMill ? `Live tooling (C-axis, from ${d})` : `Setup ${sIdx + 1} (VMC, tool from ${d})`;
    if (!turnMill) r.setups.push({ n: sIdx + 1, dir: d });
    const s = Tn.secSetups.find(x => x.dir === d);
    if (s) {
      r.ops.push({ kind: 'rough', machine: secMc, setup: turnMill ? 0 : sIdx, setupLabel: lbl, dir: d, label: `Mill flats/slots ${(s.volume / 1000).toFixed(1)} cm³`, min: mins(s.volume, pr.millMRR / mat.mi * (turnMill ? 0.7 : 1)) + pr.toolChange * 2, voxels: s.voxels, sec: true });
    }
    const hs = Tn.holeInfo.filter(h => h.setupDir === d);
    r.ops.push(...holeOps(hs, mat, rates, tol, secMc, turnMill ? 0 : sIdx).map(o => ({ ...o, setupLabel: lbl, sec: true })));
    if (!turnMill) sIdx++;
  });
  if (turnMill && dirs.length) r.notes.push('Live tooling and C-axis do the cross features without re-clamping, so runout between turned and milled features is minimal.');
  if (!turnMill && dirs.length) r.notes.push('Cross features need a V-block or collet block fixture on the VMC.');
  return r;
}

function machineRate(rates, key) { return (rates.machines[key] || { rate: 0 }).rate; }

function costRoute(r, S, mat, q, rates, tol, fin) {
  const pr = rates.process, c = rates.cutting;
  // material
  const stockMass = r.stock.volume / 1e6 * mat.density; // kg (mm³→cm³ /1000 then g→kg /1000)
  const price = r.stock.form === 'bar' ? mat.bar : mat.plate;
  const material = stockMass * price * 1.05; // 5% kerf, offcut & handling
  let cutting = 0, machining = 0, setupBatch = 0;
  const lines = [];
  const setupKeys = new Set();
  for (const o of r.ops) {
    if (o.kind === 'saw') { o.cost = c.sawPerCut + c.sawPerCm2 * o.cutArea * Math.sqrt(mat.mi); cutting += o.cost; o.min = o.min || 2; }
    else if (o.kind === 'laser' || o.kind === 'edm') { cutting += o.cost; }
    else { const rate = machineRate(rates, o.machine); o.rate = rate; o.cost = o.min / 60 * rate; machining += o.cost; setupKeys.add(o.machine + ':' + (o.setup ?? 0)); }
  }
  // setups (batch) + load/unload (per part)
  let load = 0;
  for (const k of setupKeys) {
    const [mc] = k.split(':'); const rate = machineRate(rates, mc);
    const lathe = mc === 'lathe' || mc === 'turnmill';
    setupBatch += (lathe ? pr.setupLathe : pr.setupVmc) / 60 * rate + pr.programPerSetup / 60 * machineRate(rates, 'program');
    load += (lathe ? pr.loadLathe : pr.loadVmc) / 60 * rate;
  }
  if (r.ops.some(o => o.kind === 'edm')) setupBatch += pr.setupEdm / 60 * 500;
  if (r.ops.some(o => o.kind === 'laser')) setupBatch += pr.setupLaser / 60 * 600;
  const nFeat = S.A.stacks.length + (S.features.slots || 0) + (S.turned ? S.turned.steps : 0);
  const deburr = (pr.deburrBase + pr.deburrPerFeature * nFeat) * tol.inspect / 60 * machineRate(rates, 'labour');
  // finish/treatment
  let treat = 0;
  if (fin && fin.id !== 'none') {
    const dm2 = S.A.area / 1e4, kg = S.A.volume / 1e6 * mat.density;
    treat = Math.max(fin.perDm2 * dm2 + fin.perKg * kg, fin.min / Math.max(q, 1));
  }
  const setupPer = setupBatch / Math.max(1, q);
  const direct = material + cutting + machining + load + setupPer + deburr + treat;
  const overhead = (machining + load + setupPer + deburr) * rates.overheads.overheadPct / 100;
  const margin = (direct + overhead) * rates.overheads.marginPct / 100;
  const perPart = direct + overhead + margin;
  const gst = perPart * rates.overheads.gstPct / 100;
  const cycleMin = r.ops.reduce((s, o) => s + (o.kind === 'saw' || o.kind === 'laser' || o.kind === 'edm' ? 0 : (o.min || 0)), 0);
  r.cost = { stockMass, price, material, cutting, machining, load, setupBatch, setupPer, deburr, treat, overhead, margin, perPart, gst, perPartGst: perPart + gst, batch: perPart * q, batchGst: (perPart + gst) * q, cycleMin, partMass: S.A.volume / 1e6 * mat.density, utilisation: S.A.volume / r.stock.volume };
  return r;
}

// one-call convenience
function runAll(meshes, matId, q, rates, options) {
  const T = prepareMeshes(meshes);
  const S = study(T);
  const Pl = plan(S, matId, q, rates, options);
  return { T, S, Pl };
}

const api = { studyPrismaticK, prepareMeshes, buildTriData, study, plan, runAll, DEFAULT_RATES, TOL, BAR_D, PLATE_T, likelyTap, toFrame, dirToFrame, DIRS, vec: { sub, add, mul, dot, cross, len, norm } };
if (typeof module !== 'undefined' && module.exports) module.exports = api; else root.MFG = api;
})(typeof window !== 'undefined' ? window : globalThis);
