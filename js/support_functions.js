// support fn

function trilin_to_cartesian(
  [A, B, C], //(* vertices *)
  [a, b, c], // (* side lengths *)
  [x, y, z]) //(* trilinears *)
{
  //let denom = a * x + b * y + c * z;
  let v = [a * x, b * y, c * z];
  let denom = sum(v);
  return [dot(v, [A[0], B[0], C[0]]) / denom,
    dot(v, [A[1], B[1], C[1]]) / denom
  ];
}

function get_Xn_low(orbit, sides, fn_trilin) {
  return fn_trilin(orbit,sides);
}

function get_fn_trilin(n) {
  let fn_name = sprintf("trilin_X%d",n);
  return window[fn_name];
}

function get_Xn(orbit, sides, n) {
  return get_Xn_low(orbit, sides, get_fn_trilin(n));
}

const cPi3 = Math.cos(Math.PI/3);
const sPi3 = Math.sin(Math.PI/3);
const cPi6 = Math.cos(Math.PI/6);
const sPi6 = Math.sin(Math.PI/6);

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
