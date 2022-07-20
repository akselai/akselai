class Knot {
  constructor(numberOfNodes) {
    /* params
    float ck = 1; // coefficient of repulsion between segments, so they don't unknot
    float cr = 2; // coefficient of restoring to unit length tendency, so they don't explode to infinity
    float mass = 1; // per unit length
    float damping = 0.99;
    */
    
    this.system = {nodes: [], velocities: [], forces: []};
    this.params = {ck: 0.1, cr: 1.0, damping: 0.4};
    for (let i = 0; i < numberOfNodes; i++) {
      let t = TWO_PI / numberOfNodes * i;
      //append(this.system.nodes, createVector((2 + cos(2*t)) * cos(3*t), (2 + cos(2*t)) * sin(3*t), sin(4*t)));
      append(this.system.nodes, createVector((sin(t) + 2*sin(2*t))/4, (cos(t) - 2*cos(2*t))/4, -sin(3*t)/4));
      //append(this.system.nodes, createVector(0, sin(t), cos(t)));
      //append(this.system.nodes, p5.Vector.random3D().mult(1));
      append(this.system.velocities, createVector(0, 0, 0));
      append(this.system.forces, createVector(0, 0, 0));
    }
    /*
    append(this.system.nodes, createVector(0, -1, -1));
    append(this.system.nodes, createVector(0, 1, -1));
    append(this.system.nodes, createVector(sin(millis() / 1000), cos(millis() / 1000), 1));
    append(this.system.nodes, createVector(-sin(millis() / 1000), -cos(millis() / 1000), 1));
    */
    
    this.drawForces = false;
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
      if (this.drawForces) {drawArrow(this.system.nodes[i], p5.Vector.mult(this.system.forces[i], 1), scale);}
    }
  }
  
  flow(dt) { // rk4
    let active = true;
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
      if (i != j && !equimod(i, j+1, n) && !equimod(i+1, j, n)) {
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
        
        let leftTorque = quass5(l_, 0, 1); // these lambda expressions look so spaghetti because i am bad at js lambdas 
        let rightTorque = quass5(r_, 0, 1);
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

function quass5(f, a, b) { // gaussian quadrature, 5 points
  let m = createVector(0, 0, 0);
  b = (b - a) / 2;
  m.add(f(a + b).mult(0.6222222222222222222222));
  m.add(f(a + b + b * 0.9258200997725514615666).mult(0.197979797979797979798));
  m.add(f(a + b + b * -0.9258200997725514615666).mult(0.197979797979797979798));
  m.add(f(a + b + b * 0.5773502691896257645092).mult(0.4909090909090909090909));
  m.add(f(a + b + b * -0.5773502691896257645092).mult(0.4909090909090909090909));
  return m.mult(b);
}

function quass11(f, a, b) { // gaussian quadrature, 11 points
  let m = createVector(0, 0, 0);
  b = (b - a) / 2;
  m.add(f(a + b).mult(0.2829874178574912132043));
  m.add(f(a + b + b * 0.9840853600948424644962).mult(0.0425820367510818328645));
  m.add(f(a + b + b * -0.9840853600948424644962).mult(0.0425820367510818328645));
  m.add(f(a + b + b * 0.9061798459386639927976).mult(0.1152333166224733940246));
  m.add(f(a + b + b * -0.9061798459386639927976).mult(0.1152333166224733940246));
  m.add(f(a + b + b * 0.7541667265708492204408).mult(0.1868007965564926574678));
  m.add(f(a + b + b * -0.7541667265708492204408).mult(0.1868007965564926574678));
  m.add(f(a + b + b * 0.5384693101056830910363).mult(0.2410403392286475866999));
  m.add(f(a + b + b * -0.5384693101056830910363).mult(0.2410403392286475866999));
  m.add(f(a + b + b * 0.2796304131617831934135).mult(0.2728498019125589223410));
  m.add(f(a + b + b * -0.2796304131617831934135).mult(0.2728498019125589223410));
  return m.mult(b);
}
