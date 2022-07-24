var k1;
var g3d, g2dprojection, multigraph;

var quassCoef = [];

function preload() {
  quassCoef = loadStrings("quass.txt");
}

function setup() {
  for(let i = 0; i < quassCoef.length; i++) {quassCoef[i] = parseFloat(quassCoef[i]);}
  
  k1 = new Knot(11);
  createCanvas(1200, 850);
  g3d = createGraphics(600, 600, WEBGL);
  g3d.setAttributes('antialias', true);
  g3d.setAttributes('alpha', true);
  g3d.setAttributes('depth', false); // don't set this to true!!!!!!!!!!
  g3d.setAttributes('stencil', true);
  g3d.setAttributes('premultipliedAlpha', true);
  g3d.setAttributes('preserveDrawingBuffer', true);
  g3d.setAttributes('perPixelLighting', true);
  g3d.camera(200.0, 0.0, 0.0, // eyeX, eyeY, eyeZ
    0.0, 0.0, 0.0, // centerX, centerY, centerZ
    0.0, 1.0, 0.0);
  // g3d.debugMode();
  
  g2dprojection = createGraphics(600, 600);
  multigraph = createGraphics(1200, 250);
}

var cameraRotateX = 0;
var cameraRotateY = 0;
var cameraSensitivity = 2*3.1415926535897932 / 600;

function draw() {
  let rate = 3;
  let step = 0.02;
  
  let camera_dX = mouseX - pmouseX;
  let camera_dY = mouseY - pmouseY;
  if (600 < mouseX && mouseX < 1200 && mouseIsPressed) {
    cameraRotateX += camera_dX * cameraSensitivity;
    cameraRotateY += camera_dY * cameraSensitivity;
    cameraRotateY = constrain(cameraRotateY, -HALF_PI, HALF_PI);
  }
  g3d.rotateZ(cameraRotateY);
  g3d.rotateY(cameraRotateX);
  g3d.background(200);
  // g3d.ortho();
  drawAxes(100);
  k1.display(40);
  for (let i = 0; i < rate; i++) {
    k1.flow(step);
  }
  
  g3d.rotateY(-cameraRotateX);
  g3d.rotateZ(-cameraRotateY);
  
  g2dprojection.background(255);
  projectKnot();
  
  multigraph.background(255);
  multigraph.stroke(255, 0, 0);
  graphEnergy(k1);
  
  image(g3d, 600, 0);
  image(g2dprojection, 0, 0);
  image(multigraph, 0, 600);
  
      console.log(k1.history.energy[k1.history.energy.length-1]);
}

function graphEnergy(k) {
  multigraph.strokeWeight(2);
  let h = k.history.energy;
  let scaleH = 1/5.0;
  if (h.length > 1) {
    for (let i = 1; i < h.length; i++) {
      let a = h[i-1];
      let b = h[i];
      multigraph.line(i-1, 250-a*scaleH, i, 250-b*scaleH);
    }
  }
  multigraph.stroke(240, 0, 0, 50);
  multigraph.line(0, 250-189.28*scaleH, 1200, 250-189.28*scaleH);
  multigraph.stroke(0, 240, 0, 50);
  multigraph.line(0, 250-33.141*scaleH, 1200, 250-33.141*scaleH);
}

function drawAxes(scale) {
  g3d.stroke(255, 0, 0);
  g3d.line(-scale, 0, 0, scale, 0, 0);
  g3d.stroke(0, 255, 0);
  g3d.line(0, -scale, 0, 0, scale, 0);
  g3d.stroke(0, 0, 255);
  g3d.line(0, 0, -scale, 0, 0, scale);
}

function drawArrowH(tail_, head_, scale) {
  let tail = p5.Vector.mult(tail_, scale);
  let head = p5.Vector.mult(head_, scale);
  g3d.line(tail.x, tail.y, tail.z, head.x, head.y, head.z);
  g3d.push();
  g3d.translate(head);
  let d = p5.Vector.sub(head, tail);
  let h2 = createVector(d.y, d.z);
  let h1 = createVector(d.x, h2.mag());
  g3d.rotateX(PI/2);
  g3d.rotateX(h2.heading() - PI/2);
  g3d.rotateZ(h1.heading() - PI/2);
  g3d.rotateY(millis() / 200);
  g3d.stroke(210, 255);
  g3d.fill(50, 255, 230);
  g3d.cone(2, 4, 5, 1, false);
  g3d.pop();
  g3d.line(tail.x, tail.y, tail.z, head.x, head.y, head.z);
}

function drawArrow(node, vector, scale) {
  let v = p5.Vector.add(node, vector);
  drawArrowH(node, v, scale);
}

function projectKnot() { // find the projection of a point onto a plane
  g2dprojection.stroke(0); // https://math.stackexchange.com/questions/100761
  g2dprojection.strokeWeight(2);
  let n = k1.system.nodes.length;
  let projNodes = [];
  for (let i = 0; i < n; i++) {
    let a = k1.system.nodes[i];
    a = orthoSync(a); // also has x coordinate
    let aX = a.z * -200 + 300;
    let aY = a.y * 200 + 300;
    let aZ = a.x;
    append(projNodes, createVector(aX, aY, aZ));
    // g2dprojection.text(i, aX, aY);
  }
  
  g2dprojection.stroke(0);
  let crossingInf = createArray(n, n); // crossing information
  for (let i = 0; i < n; i++) {
    let a = projNodes[i];
    let b = projNodes[(i+1)%n];
    for (let j = i + 2; j < (i == 0 ? n - 1 : n); j++) {
      let c = projNodes[j];
      let d = projNodes[(j+1)%n];
      let inter = intersect2d(a, b, c, d);
      if(inter.intersect) {
        let aVect = p5.Vector.lerp(projNodes[i], projNodes[(i+1)%n], inter.lerp);
        let bVect = p5.Vector.lerp(projNodes[j], projNodes[(j+1)%n], inter.lerp2);
        if (aVect.z > bVect.z) {
          crossingInf[i][j] = {crossing: true, orientation: false, cut: inter.lerp2};
          crossingInf[j][i] = {crossing: true, orientation: true, cut: inter.lerp};
        } else {
          crossingInf[i][j] = {crossing: true, orientation: true, cut: inter.lerp2};
          crossingInf[j][i] = {crossing: true, orientation: false, cut: inter.lerp};
        }
      } else {
      crossingInf[i][j] = {crossing: false, orientation: null, cut: NaN};
      crossingInf[j][i] = {crossing: false, orientation: null, cut: NaN};
      }
    }
  }
  
  // draw crossings
  for (let i = 0; i < n; i++) {
    let a = projNodes[i];
    let b = projNodes[(i+1)%n];
    let cuts = [];
    for (let j = 0; j < n; j++) {
      if (i != j && i != (j+1)%n && j != (i+1)%n && crossingInf[i][j].crossing && crossingInf[i][j].orientation) {
        append(cuts, crossingInf[i][j].cut);
      }
    }
    drawCrossing(a, b, cuts.sort());
  }
}

function drawCrossing(a, b, cuts) {
  let cutRadius = 6;
  let boundaries = [0];
  let segment = p5.Vector.sub(b, a).mag();
  for (let i = 0; i < cuts.length; i++) {
    let minBound = cuts[i] - cutRadius/segment;
    let maxBound = cuts[i] + cutRadius/segment;
    append(boundaries, minBound);
    append(boundaries, maxBound);
  }
  append(boundaries, 1);
  stroke(0);
  strokeWeight(2);
  for (let i = 0; i < boundaries.length; i+=2) { // the length of boundaries is twice that of cuts
    let x = p5.Vector.lerp(a, b, boundaries[i]);
    let y = p5.Vector.lerp(a, b, boundaries[i+1]);
    if (boundaries[i] < boundaries[i+1]) {
      g2dprojection.line(x.x, x.y, y.x, y.y);
    }
  }
}

function orthoSync(v_) { // sync projection angle with camera angle 
  /*
  g3d.rotateZ(cameraRotateY);
  g3d.rotateY(cameraRotateX);
  */
  let v;
  let cX = cameraRotateX;
  let cY = cameraRotateY;
  
  v = v_.copy();
  v.set(v.z*sin(cX) + v.x*cos(cX), v.y, v.z*cos(cX) - v.x*sin(cX));
  v.set(v.x*cos(cY) - v.y*sin(cY), v.x*sin(cY) + v.y*cos(cY), v.z);
  return v;
}

function intersection(a, b, c, d) { // one sided implementation, only checks for the h in (ab, cd)
  let e = p5.Vector.sub(b, a);
  let f = p5.Vector.sub(d, c);
  let p = createVector(-e.y, e.x);
  let dot = p5.Vector.dot(f, p);
  let u = p5.Vector.dot(p5.Vector.sub(a, c), p);
  let h = u / dot;
  return {intersect: (h > 0 && h < 1), lerp: h};
} // https://stackoverflow.com/a/563275

function intersect2d(a_, b_, c_, d_) { // symmetric implementation, checks for the g in (cd, ab) too
  a = createVector(a_.x, a_.y);
  b = createVector(b_.x, b_.y);
  c = createVector(c_.x, c_.y);
  d = createVector(d_.x, d_.y); // 2d-ify
  let x = intersection(a, b, c, d);
  let y = intersection(c, d, a, b);
  return {intersect: x.intersect && y.intersect, lerp: x.lerp, lerp2: y.lerp};
} // https://stackoverflow.com/a/1201356

function createArray(length) { // create n-dimensional array
  let arr = new Array(length || 0);
  let i = length;
  if (arguments.length > 1) {
    var args = Array.prototype.slice.call(arguments, 1);
    while(i--) {arr[length-1 - i] = createArray.apply(this, args);}
  }
  return arr;
} // https://stackoverflow.com/a/966938

function keyPressed() {
  if (keyCode === 32) {
    k1.drawForces = !k1.drawForces;
  }
  debugger;
}
