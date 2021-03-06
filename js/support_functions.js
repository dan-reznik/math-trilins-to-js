// support fn
const sqrt3 = Math.sqrt(3.0);
const phi = (Math.sqrt(5)+1.0)/2;
const cPi3 = Math.cos(Math.PI/3);
const sPi3 = Math.sin(Math.PI/3);
const cPi6 = Math.cos(Math.PI/6);
const sPi6 = Math.sin(Math.PI/6);
const Pi = Math.PI;

function oneOv(x) { return 1.0/x; }

function triAreaHeron(a, b, c) {
  let s = (a + b + c)/2.0;
  return sqrt(s*(s-a)*(s-b)*(s-c));
}

function lawOfCosines(a, b, c) { return (b*b+c*c-a*a)/(2.0*b*c); }

function getSin(theCos) { return sqrt(1.0 - theCos*theCos); }

function cosDoubleAngle(cosAngle) { return 2*cosAngle*cosAngle - 1.0; }
function sinDoubleAngle(sinAngle,cosAngle) { return 2*sinAngle*cosAngle; }

function cosHalfAngle(cosAngle) { return sqrt((1.0+cosAngle)/2); }
function sinHalfAngle(cosAngle) { return sqrt((1.0-cosAngle)/2); }

function getSinApmB(sa, sb, ca, cb) { return [sa*cb + sb*ca, sa*cb - sb*ca]; }

function getSinApmB1(sa, sb, ca, cb) { return sa*cb + sb*ca; }
function getSinApmB2(sa, sb, ca, cb) { return sa*cb - sb*ca; }

function cosThirdAngle(c) {
  let th = Math.acos(c);
  return Math.cos(th/3);
}

function sinThirdAngle(s) {
  let th = Math.asin(s);
  return Math.sin(th/3);
}

function sinCosTripleAngle(s, c, s2, c2) {
  let c3 = c2*c - s2*s;
  let s3 = s2*c + s*c2;
  return [s3, c3];
}

function cosTripleAngle(s, c, s2, c2) {
  return c2*c - s2*s;
}

function sinTripleAngle(s, c, s2, c2) {
  return s2*c + s*c2;
}

function getJ(R,SW) {
  return Math.sqrt(9*R*R-2*SW);
}

function getE(S,SW) {
  // return (S*SW*SW)/Math.sqrt(S*S+SW*SW);
  // moses: Sin[w]=S/Sqrt[S^2+SW^2]
  // since e = Sqrt[1-4 Sin[w]^2], then:
  // e = Sqrt[1-4/(1+(SW/S)^2)]
  let ratio=SW/S;
  return Math.sqrt(1-4/(1+ratio*ratio));
}

function getCotPrime(a,b,c) {
  return 1/Math.tan((2*a*Math.PI)/(a+b+c));
}