let logo = [];

function preload() {
  logo = loadStrings("logo.txt");
}

var p;
var colors;

function setup() {
  noCanvas();
  p = createP("");
  colors = ["<span style=\"color:cyan;\">", "<span style=\"color:magenta;\">", "<span style=\"color:yellow;\">"]
}

function draw() {
  let buffer = (millis()%7000<3500?slideDigit(sin(millis()/300)+1):noVariation());
  p.html(buffer);
  p.style('font-family', 'Courier New');
  p.style('font-size', '20px');
  p.position(0, 0);
}

let backDigit = 0;
function noVariation() {
  let buffer = "";
  for (let i = 0; i < logo.length; i++) {
    for (let j = 0; j < logo[i].length; j++) {
      if (logo[i].charAt(j) == ' ') {
        if(millis()%100 <= 10) {backDigit = round(random(-0.4999, 9.4999));}
        buffer += backDigit + '';
      } else {
        buffer += colors[round(random(-0.4999, 2.4999))] + round(random(-0.4999, 9.4999)) + "</span>";
      }
    }
    buffer += '<br>';
  }
  return "<span style=\"letter-spacing: 5px; background-color:#444444;\">" + buffer + "</span>";
}

function slideDigit(mag10) {
  let seed = (millis()/628)%65536;
  let rng64 = (n)=>(q=n>>8^(n&=255))/2^n<<7^n^57460^q%2*40564; // mario 64 rng (https://codegolf.stackexchange.com/a/250175)
  let buffer = "";
  for (let i = 0; i < logo.length; i++) {
    let bufferS = "";
    let parallel = "";
    for (let j = 0; j < logo[i].length; j++) {
      if (logo[i].charAt(j) == ' ') {
        parallel += (parallel.length > 78 ? '' : '0');
        bufferS += (parallel.length > 78 ? '' : '0');
      } else {
        let d = round(random(-0.4999, pow(10, mag10)-0.5001)) + '';
        parallel += d;
        bufferS += colors[rng64(seed)%3] + '' + d.slice(0, (parallel.length > 78 ? 78 - parallel.length : d.length)) + "</span>";
      }
      seed = rng64(seed);
    }
    buffer += bufferS
    buffer += '<br>';
  }
  return "<span style=\"letter-spacing: 5px; background-color:#444444;\">" + buffer + "</apn>";
}
