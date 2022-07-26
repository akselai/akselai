var logo;

function preload() {
  logo = [
    "             kk                                   llll                    ii  ",
    "             kk                                   llll                    ii  ",
    "             kk                                     ll                        ",
    "             kk                                     ll                        ",
    "  aaaaaa     kk    kk     ssssssss     eeeeee       ll       aaaaaa     iiii  ",
    "  aaaaaaa    kk   kkk    sssssssss    eeeeeeee      ll       aaaaaaa    iiii  ",
    "       aaa   kk  kkk    sss          eee    eee     ll            aaa     ii  ",
    "        aa   kk kkk     sss          ee      ee     ll             aa     ii  ",
    "  aaaaaaaa   kkkkk       sssssss     eeeeeeeeee     ll       aaaaaaaa     ii  ",
    " aaaaaaaaa   kkkkk        sssssss    eeeeeeeeee     ll      aaaaaaaaa     ii  ",
    "aaa     aa   kk kkk            sss   ee             ll     aaa     aa     ii  ",
    "aaa     aa   kk  kkk           sss   eee            ll     aaa     aa     ii  ",
    " aaaaaaaaa   kk   kkk   sssssssss     eeeeeee     llllll    aaaaaaaaa   iiiiii",
    "  aaaaaaaa   kk    kk   ssssssss       eeeeee     llllll     aaaaaaaa   iiiiii",
  ];
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
  // p.style('text-align', 'center');
  p.style('font-family', 'Courier New');
  p.style('font-size', '20px');
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
        /*
        let str = (parallel.length > 78 ? '' : '0'); 
        parallel += str;
        bufferS += str;
        */
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
