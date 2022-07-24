class Knot {
  constructor(numberOfNodes) {
    this.system = {nodes: [], velocities: [], forces: []};
    this.params = {ck: 0.1, cr: 1, damping: 0.4};
    this.history = {frames: 1200, energy: []};
    
    for (let i = 0; i < numberOfNodes; i++) {
      let t = TWO_PI / numberOfNodes * i;
      //append(this.system.nodes, createVector((2 + cos(2*t)) * cos(3*t), (2 + cos(2*t)) * sin(3*t), sin(4*t)));
      append(this.system.nodes, createVector((sin(t) + 2*sin(2*t))/4, (cos(t) - 2*cos(2*t))/4, -sin(3*t)/4));
      //append(this.system.nodes, createVector(0, sin(t), cos(t)));
      //append(this.system.nodes, p5.Vector.random3D().mult(1));
      this.system.nodes[i].add(createVector(random(-0.3, 0.3), random(-0.3, 0.3), random(-0.3, 0.3)));
      append(this.system.velocities, createVector(0, 0, 0));
      append(this.system.forces, createVector(0, 0, 0));
    }
    
    this.drawForces = false;
    this.showNodeLabels = true;
  }
  
  display(scale) {
    let n = this.system.nodes.length;
    for (let i = 0; i < n; i++) {
      let head = this.system.nodes[i];
      let tail = this.system.nodes[(i + 1)%n];
      g3d.stroke(120);
      g3d.line(scale * head.x, scale * head.y, scale * head.z,
        scale * tail.x, scale * tail.y, scale * tail.z);
      g3d.stroke(50, 210, 210);
      g3d.fill(50, 255, 230);
      if (this.drawForces) {drawArrow(head, p5.Vector.mult(this.system.forces[i], 1), scale);}
    }
  }
  
  flow(dt) { // rk4
    let active = true;
    this.measureSegLengths();
    this.recenter();
    if (active) {
      /* pseudocode for runge-kutta
      let fA = f(this.system);
      let fB = f(this.system + dt/2 * fA);
      let fC = f(this.system + dt/2 * fB);
      let fD = f(this.system + dt * fC);
      this.system = this.system + dt/6 * (fA + 2*fB + 2*fC + fD); 
      */// actual code
      let fA = totalDiff(this.system, this.params);
      let fB = totalDiff(addLinear(this.system, scaleLinear(fA, dt/2)), this.params);
      let fC = totalDiff(addLinear(this.system, scaleLinear(fB, dt/2)), this.params);
      let fD = totalDiff(addLinear(this.system, scaleLinear(fC, dt)), this.params);
      this.system = addLinear(this.system, scaleLinear(addLinear(addLinear(addLinear(fA, scaleLinear(fB, 2)), scaleLinear(fC, 2)), fD), dt/6));
    }
    
    append(this.history.energy, this.moebiusEnergy()); // update knot frames (as a movie)
    while(this.history.energy.length > this.history.frames) {this.history.energy.shift();}
  }
  
  recenter() {
    let minX = Number.MAX_VALUE;
    let minY = Number.MAX_VALUE;
    let minZ = Number.MAX_VALUE; // the minimum of 0 numbers is infty
    let maxX = -Number.MAX_VALUE;
    let maxY = -Number.MAX_VALUE;
    let maxZ = -Number.MAX_VALUE; // the maximum of 0 numbers is -infty
    let n = this.system.nodes.length;
    for (let i = 0; i < n; i++) {
      let thisNode = this.system.nodes[i];
      if (thisNode.x > maxX) {maxX = thisNode.x;}
      if (thisNode.y > maxY) {maxY = thisNode.y;}
      if (thisNode.z > maxZ) {maxZ = thisNode.z;}
      if (thisNode.x < minX) {minX = thisNode.x;}
      if (thisNode.y < minY) {minY = thisNode.y;}
      if (thisNode.z < minZ) {minZ = thisNode.z;}
    }
    let centerX = maxX + minX; centerX /= 2;
    let centerY = maxY + minY; centerY /= 2;
    let centerZ = maxZ + minZ; centerZ /= 2;
    for (let i = 0; i < n; i++) {
      this.system.nodes[i].sub(createVector(centerX, centerY, centerZ));
    }
  }
  
  measureSegLengths() {
    this.cumSegLengths = [];
    let ns = this.system.nodes;
    let x = 0;
    for (let i = 0; i < ns.length; i++) {
      x += p5.Vector.dist(ns[i], ns[(i+1)%(ns.length)]);
      append(this.cumSegLengths, x);
    }
  }
  
  moebiusEnergy() {
    let result = 0;
    let ns = this.system.nodes;
    let n = ns.length;
    for (let i = 0; i < n; i++) {
      for (let j = 0; j < n; j++) {
        if (i != j && !equimod(i, j+1, n) && !equimod(i+1, j, n)) { // i not next to j
          result += p5.Vector.dist(ns[i], ns[(i+1)%n]) *
                    p5.Vector.dist(ns[j], ns[(j+1)%n]) *
                    pow(segmentDistance(ns[i], ns[(i+1)%n], ns[j], ns[(j+1)%n]), -2);
        }
      }
    }
    return result;
  }
}

function totalDiff(system, params) { // derivatives of the system
  let n = system.nodes.length;
  let f = [];
  system.forces = totalTension(system.nodes, params.cr);
  let repulsion = totalRepulsion(system.nodes, params.ck);
  for (let i = 0; i < n; i++) {system.forces[i].add(repulsion[i]);}
  for (let i = 0; i < n; i++) {system.forces[i].add(p5.Vector.mult(system.velocities[i], -params.damping));}
  for (let i = 0; i < n; i++) {append(f, createVector(0, 0, 0));}
  let result = {nodes: system.velocities, velocities: system.forces, forces: f}; // forces = 0 because no 3rd order derivatives
  return result;
}

function totalTension(nodes, cr) {
  let n = nodes.length;
  let forces = [];
  for (let i = 0; i < n; i++) {
    append(forces, createVector(0, 0, 0));
  }
  for (let i = 0; i < n; i++) {
    for (let j = 0; j < n; j++) {
      if (equimod(i, j-1, n) || equimod(i, j+1, n)) { // i is next to j
        let left = nodes[i];
        let right = nodes[j];
        let mag = p5.Vector.dist(left, right) - 1;
        let dir = p5.Vector.sub(right, left).normalize();
        forces[i].add(p5.Vector.mult(dir, mag * cr));
        forces[j].add(p5.Vector.mult(dir, -mag * cr));
      }
    }
  }
  return forces;
}

function totalRepulsion(nodes, ck) {
  let n = nodes.length;
  let forces = [];
  for (let i = 0; i < n; i++) {
    append(forces, createVector(0, 0, 0));
  }
  for (let i = 0; i < n; i++) {
    for (let j = 0; j < n; j++) {
      if (i != j && !equimod(i, j+1, n) && !equimod(i+1, j, n)) { // i not next to j
      //if (i == 0 && j == 2) { // debug
        let a = nodes[i];
        let b = nodes[mod(i+1, n)];
        let c = nodes[j];
        let d = nodes[mod(j+1, n)];
        
        let re = (lerp) => p5.Vector.lerp(c, d, lerp);
        let res = (lerp) => lineSegmentDist(re(lerp), a, b);
        let rest = (lerp) => p5.Vector.sub(re(lerp), res(lerp).point);
        
        let l_ = (lerp) => p5.Vector.mult(p5.Vector.normalize(rest(lerp)), pow(res(lerp).mag, -2) * (1-lerp) * ck);
        let r_ = (lerp) => p5.Vector.mult(p5.Vector.normalize(rest(lerp)), pow(res(lerp).mag, -2) * lerp * ck);
                
        let leftTorque = quassVector(l_, 0, 1, 5); // these lambda expressions look so spaghetti because i am bad at js lambdas 
        let rightTorque = quassVector(r_, 0, 1, 5);
        /*
        if (i == 0 && j == 2) {
          drawArrow(_, leftTorque, 40);
          drawArrow(_, rightTorque, 40);
        }
        */
        forces[j].add(leftTorque);
        forces[mod(j+1, n)].add(rightTorque);
      }
    }
  }
  return forces;
}

function scaleLinear(system, scalar) {
  let result = {nodes: [], velocities: [], forces: []};
  for (let i = 0; i < system.nodes.length; i++) {
    append(result.nodes, p5.Vector.mult(system.nodes[i], scalar));
    append(result.velocities, p5.Vector.mult(system.velocities[i], scalar));
    append(result.forces, p5.Vector.mult(system.forces[i], scalar));
  }
  return result;
}
  
function addLinear(systemA, systemB) {
  let result = {nodes: [], velocities: [], forces: []};
  for (let i = 0; i < systemA.nodes.length; i++) {
    append(result.nodes, p5.Vector.add(systemA.nodes[i], systemB.nodes[i]));
    append(result.velocities, p5.Vector.add(systemA.velocities[i], systemB.velocities[i]));
    append(result.forces, p5.Vector.add(systemA.forces[i], systemB.forces[i]));
  }
  return result;
}

function lineSegmentDist(p, a, b) {
  let u = p5.Vector.sub(b, a);
  let v = p5.Vector.sub(p, a);
  let u_ = p5.Vector.normalize(u);
  let dot = p5.Vector.dot(u_, v);
  dot = constrain(dot, 0, u.mag());
  let closest = p5.Vector.mult(u_, dot);
  let result = p5.Vector.add(a, closest);
  let dist = p5.Vector.dist(result, p);
  return {point: result, mag: dist};
}

function segmentAngle(a, x, b) {
  return p5.Vector.sub(a, x).angleBetween(p5.Vector.sub(b, x));
}

function mod(x, y) {
  return (x % y < 0 ? x % y + y : x % y);
}

function equimod(x, y, m) {
  return mod(x, m) === mod(y, m);
}

function quassVector(f, a, b, order) { // gaussian quadrature for vector valued functions
  let m = createVector(0, 0, 0);
  b = (b - a) / 2;
  let w = quassCoef[order*order];
  m.add(f(a + b).mult(w));
  for (let i = order*order + 1; i < (order + 1)*(order + 1); i+=2) {
    let x_ = quassCoef[i];
    let w_ = quassCoef[i+1];
    m.add(f(a + b + b * x_).mult(w_));
    m.add(f(a + b - b * x_).mult(w_));
  }
  return m.mult(b);
}

function quass(f, a, b, order) { // gaussian quadrature
  let m = 0;
  b = (b - a) / 2;
  let w = quassCoef[order*order];
  m += f(a + b) * w;
  for (let i = order*order + 1; i < (order + 1)*(order + 1); i+=2) {
    let x_ = quassCoef[i];
    let w_ = quassCoef[i+1];
    m += f(a + b + b * x_) * w_;
    m += f(a + b - b * x_) * w_;
  }
  
  return m * b;
}

function quass2(f, a1, b1, a2, b2, order) { // gaussian quadrature, double integral
  let f_ = (x) => ((y) => f(x, y)); // fix x, integrate over y
  let I = (x) => quass(f_(x), a2(x), b2(x), order); // inner integral's bound might depend on x, that's why a2 and b2 are also lambdas
  let result = quass(I, a1, b1, order);
  return result;
}

function segmentDistance(a1, a2, b1, b2) {
  let u = p5.Vector.sub(a2, a1);
  let v = p5.Vector.sub(b2, b1);
  let w = p5.Vector.sub(a1, b1);
  let a = p5.Vector.dot(u, u);
  let b = p5.Vector.dot(u, v);
  let c = p5.Vector.dot(v, v);
  let d = p5.Vector.dot(u, w);
  let e = p5.Vector.dot(v, w);
  let D = a * c - b * b;
  let sc, sN, sD = D;
  let tc, tN, tD = D;
  let epsilon = 0.0001;
  if (D < epsilon)
  {
    sN = 0;
    sD = 1;
    tN = e;
    tD = c;
  } else {
    sN = (b * e - c * d);
    tN = (a * e - b * d);
    if (sN < 0) {
      sN = 0;
      tN = e;
      tD = c;
    } else if (sN > sD) {
      sN = sD;
      tN = e + b;
      tD = c;
    }
  }
  if (tN < 0) {
    tN = 0;
    if (-d < 0) {
      sN = 0;
    } else if (-d > a) {
      sN = sD;
    } else {
      sN = -d;
      sD = a;
    }
  } else if (tN > tD) {
    tN = tD;
    if ((-d + b) < 0) {
      sN = 0;
    } else if ((-d + b) > a) {
      sN = sD;
    } else {
      sN = (-d + b);
      sD = a;
    }
  }
  if (abs(sN) < epsilon) {
    sc = 0;
  } else {
    sc = sN / sD;
  }
  if (abs(tN) < epsilon) {
    tc = 0;
  } else {
    tc = tN / tD;
  }
  let dP = p5.Vector.add(w, p5.Vector.sub(p5.Vector.mult(u, sc), p5.Vector.mult(v, tc)));
  let distance = sqrt(p5.Vector.dot(dP, dP));
  return distance;
}
