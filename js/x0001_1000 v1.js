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

function get_Xn_low_bary(orbit, fn_bary) {
  return fn_bary(orbit);
}

function get_fn_bary(n) {
  let fn_name = sprintf("bary_X%d",n);
  return window[fn_name];
}

function get_Xn_bary(orbit, n) {
  return get_Xn_low_bary(orbit, get_fn_bary(n));
}

// vs = vertices (2-vectors), bs = barycentrics (scalars)
function barys_to_cartesian(vs, bs) {
  const bs_sum = sum(bs);
  const vs_scaled = vs.map((v,i)=>vscale(v,bs[i]));
  const vs_scaled_norm = vscale(vsum3(...vs_scaled),1/bs_sum);
  return vs_scaled_norm;
}

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
  return Math.cot((2*a*Math.PI)/(a+b+c));
}

function barys_X1(orbit) {
   /* begin vars */

   /* end vars */
   let v1 = a;
   let v2 = b;
   let v3 = c;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X2(orbit) {
   /* begin vars */

   /* end vars */
   let v1 = 1;
   let v2 = 1;
   let v3 = 1;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X3(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   /* end vars */
   let v1 = a2*(a2-b2-c2);
   let v2 = b2*(-a2+b2-c2);
   let v3 = c2*(-a2-b2+c2);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X4(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   /* end vars */
   let v1 = (a2+b2-c2)*(a2-b2+c2);
   let v2 = (a2+b2-c2)*(-a2+b2+c2);
   let v3 = (a2-b2+c2)*(-a2+b2+c2);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X5(orbit) {
   /* begin vars */
   let a2=a*a;
   let c2=c*c;
   let b2=b*b;
   let a4=a2*a2;
   let c4=c2*c2;
   let b4=b2*b2;
   /* end vars */
   let v1 = -(a2*b2)+b4-a2*c2-2*b2*c2+c4;
   let v2 = a4-a2*b2-2*a2*c2-b2*c2+c4;
   let v3 = a4-2*a2*b2+b4-a2*c2-b2*c2;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X6(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   /* end vars */
   let v1 = a2;
   let v2 = b2;
   let v3 = c2;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X7(orbit) {
   /* begin vars */

   /* end vars */
   let v1 = (a+b-c)*(a-b+c);
   let v2 = (a+b-c)*(-a+b+c);
   let v3 = (a-b+c)*(-a+b+c);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X8(orbit) {
   /* begin vars */

   /* end vars */
   let v1 = -a+b+c;
   let v2 = a-b+c;
   let v3 = a+b-c;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X9(orbit) {
   /* begin vars */

   /* end vars */
   let v1 = a*(a-b-c);
   let v2 = b*(-a+b-c);
   let v3 = c*(-a-b+c);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X10(orbit) {
   /* begin vars */

   /* end vars */
   let v1 = b+c;
   let v2 = a+c;
   let v3 = a+b;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X11(orbit) {
   /* begin vars */

   /* end vars */
   let v1 = (b-c)*(b-c)*(-a+b+c);
   let v2 = (-a+c)*(-a+c)*(a-b+c);
   let v3 = (a-b)*(a-b)*(a+b-c);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X12(orbit) {
   /* begin vars */

   /* end vars */
   let v1 = (a+b-c)*(a-b+c)*(b+c)*(b+c);
   let v2 = (a+b-c)*(a+c)*(a+c)*(-a+b+c);
   let v3 = (a+b)*(a+b)*(a-b+c)*(-a+b+c);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X13(orbit) {
   /* begin vars */
   let area=triAreaHeron(a,b,c);
   let S=2*area;
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   /* end vars */
   let v1 = 3*(a2+b2-c2)*(a2-b2+c2)+4*S*(S+a2*sqrt3);
   let v2 = 3*(a2+b2-c2)*(-a2+b2+c2)+4*S*(S+b2*sqrt3);
   let v3 = 3*(a2-b2+c2)*(-a2+b2+c2)+4*S*(S+c2*sqrt3);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X14(orbit) {
   /* begin vars */
   let area=triAreaHeron(a,b,c);
   let S=2*area;
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   /* end vars */
   let v1 = 3*(a2+b2-c2)*(a2-b2+c2)-4*S*(-S+a2*sqrt3);
   let v2 = 3*(a2+b2-c2)*(-a2+b2+c2)-4*S*(-S+b2*sqrt3);
   let v3 = 3*(a2-b2+c2)*(-a2+b2+c2)-4*S*(-S+c2*sqrt3);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X15(orbit) {
   /* begin vars */
   let area=triAreaHeron(a,b,c);
   let c2=c*c;
   let b2=b*b;
   let S=2*area;
   let a2=a*a;
   /* end vars */
   let v1 = a2*(-2*S+(a2-b2-c2)*sqrt3);
   let v2 = b2*(-2*S+(-a2+b2-c2)*sqrt3);
   let v3 = c2*(-2*S+(-a2-b2+c2)*sqrt3);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X16(orbit) {
   /* begin vars */
   let area=triAreaHeron(a,b,c);
   let c2=c*c;
   let b2=b*b;
   let S=2*area;
   let a2=a*a;
   /* end vars */
   let v1 = a2*(2*S+(a2-b2-c2)*sqrt3);
   let v2 = b2*(2*S+(-a2+b2-c2)*sqrt3);
   let v3 = c2*(2*S+(-a2-b2+c2)*sqrt3);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X17(orbit) {
   /* begin vars */
   let area=triAreaHeron(a,b,c);
   let S=2*area;
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   /* end vars */
   let v1 = (a2+b2-c2+2*S*sqrt3)*(a2-b2+c2+2*S*sqrt3);
   let v2 = (a2+b2-c2+2*S*sqrt3)*(-a2+b2+c2+2*S*sqrt3);
   let v3 = (a2-b2+c2+2*S*sqrt3)*(-a2+b2+c2+2*S*sqrt3);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X18(orbit) {
   /* begin vars */
   let area=triAreaHeron(a,b,c);
   let S=2*area;
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   /* end vars */
   let v1 = (-a2+b2-c2+2*S*sqrt3)*(-a2-b2+c2+2*S*sqrt3);
   let v2 = (a2-b2-c2+2*S*sqrt3)*(-a2-b2+c2+2*S*sqrt3);
   let v3 = (a2-b2-c2+2*S*sqrt3)*(-a2+b2-c2+2*S*sqrt3);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X19(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   /* end vars */
   let v1 = a*(a2+b2-c2)*(a2-b2+c2);
   let v2 = b*(a2+b2-c2)*(-a2+b2+c2);
   let v3 = c*(a2-b2+c2)*(-a2+b2+c2);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X20(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let c4=c2*c2;
   let b4=b2*b2;
   let a4=a2*a2;
   /* end vars */
   let v1 = -3*a4+2*a2*b2+b4+2*a2*c2-2*b2*c2+c4;
   let v2 = a4+2*a2*b2-3*b4-2*a2*c2+2*b2*c2+c4;
   let v3 = a4-2*a2*b2+b4+2*a2*c2+2*b2*c2-3*c4;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X21(orbit) {
   /* begin vars */

   /* end vars */
   let v1 = a*(a+b)*(a-b-c)*(a+c);
   let v2 = b*(a+b)*(-a+b-c)*(b+c);
   let v3 = c*(a+c)*(-a-b+c)*(b+c);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X22(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let c4=c2*c2;
   let b4=b2*b2;
   let a4=a2*a2;
   /* end vars */
   let v1 = a2*(a4-b4-c4);
   let v2 = b2*(-a4+b4-c4);
   let v3 = c2*(-a4-b4+c4);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X23(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let c4=c2*c2;
   let b4=b2*b2;
   let a4=a2*a2;
   /* end vars */
   let v1 = a2*(a4-b4+b2*c2-c4);
   let v2 = b2*(-a4+b4+a2*c2-c4);
   let v3 = c2*(-a4+a2*b2-b4+c4);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X24(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let c4=c2*c2;
   let b4=b2*b2;
   let a4=a2*a2;
   /* end vars */
   let v1 = a2*(a2+b2-c2)*(a2-b2+c2)*(a4-2*a2*b2+b4-2*a2*c2+c4);
   let v2 = b2*(a2+b2-c2)*(-a2+b2+c2)*(a4-2*a2*b2+b4-2*b2*c2+c4);
   let v3 = c2*(a2-b2+c2)*(-a2+b2+c2)*(a4+b4-2*a2*c2-2*b2*c2+c4);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X25(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   /* end vars */
   let v1 = a2*(a2+b2-c2)*(a2-b2+c2);
   let v2 = b2*(a2+b2-c2)*(-a2+b2+c2);
   let v3 = c2*(a2-b2+c2)*(-a2+b2+c2);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X26(orbit) {
   /* begin vars */
   let a2=a*a;
   let c2=c*c;
   let c4=c2*c2;
   let c6=c2*c4;
   let b2=b*b;
   let b4=b2*b2;
   let b6=b2*b4;
   let a4=a2*a2;
   let a6=a2*a4;
   let c8=c2*c6;
   let b8=b2*b6;
   let a8=a2*a6;
   /* end vars */
   let v1 = a2*(a8-2*a6*b2+2*a2*b6-b8-2*a6*c2+2*b6*c2-2*b4*c4+2*a2*c6+2*b2*c6-c8);
   let v2 = b2*(-a8+2*a6*b2-2*a2*b6+b8+2*a6*c2-2*b6*c2-2*a4*c4+2*a2*c6+2*b2*c6-c8);
   let v3 = c2*(-a8+2*a6*b2-2*a4*b4+2*a2*b6-b8+2*a6*c2+2*b6*c2-2*a2*c6-2*b2*c6+c8);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X27(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   /* end vars */
   let v1 = (a+b)*(a+c)*(a2+b2-c2)*(a2-b2+c2);
   let v2 = (a+b)*(b+c)*(a2+b2-c2)*(-a2+b2+c2);
   let v3 = (a+c)*(b+c)*(a2-b2+c2)*(-a2+b2+c2);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X28(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   /* end vars */
   let v1 = a*(a+b)*(a+c)*(a2+b2-c2)*(a2-b2+c2);
   let v2 = b*(a+b)*(b+c)*(a2+b2-c2)*(-a2+b2+c2);
   let v3 = c*(a+c)*(b+c)*(a2-b2+c2)*(-a2+b2+c2);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X29(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   /* end vars */
   let v1 = (a+b)*(a-b-c)*(a+c)*(a2+b2-c2)*(a2-b2+c2);
   let v2 = (a+b)*(-a+b-c)*(b+c)*(a2+b2-c2)*(-a2+b2+c2);
   let v3 = (a+c)*(-a-b+c)*(b+c)*(a2-b2+c2)*(-a2+b2+c2);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X30(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let c4=c2*c2;
   let b4=b2*b2;
   let a4=a2*a2;
   /* end vars */
   let v1 = -2*a4+a2*b2+b4+a2*c2-2*b2*c2+c4;
   let v2 = a4+a2*b2-2*b4-2*a2*c2+b2*c2+c4;
   let v3 = a4-2*a2*b2+b4+a2*c2+b2*c2-2*c4;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X31(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let c3=c2*c;
   let b3=b2*b;
   let a3=a2*a;
   /* end vars */
   let v1 = a3;
   let v2 = b3;
   let v3 = c3;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X32(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let c4=c2*c2;
   let b4=b2*b2;
   let a4=a2*a2;
   /* end vars */
   let v1 = a4;
   let v2 = b4;
   let v3 = c4;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X33(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   /* end vars */
   let v1 = a*(a-b-c)*(a2+b2-c2)*(a2-b2+c2);
   let v2 = b*(-a+b-c)*(a2+b2-c2)*(-a2+b2+c2);
   let v3 = c*(-a-b+c)*(a2-b2+c2)*(-a2+b2+c2);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X34(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   /* end vars */
   let v1 = a*(a+b-c)*(a-b+c)*(a2+b2-c2)*(a2-b2+c2);
   let v2 = b*(a+b-c)*(-a+b+c)*(a2+b2-c2)*(-a2+b2+c2);
   let v3 = c*(a-b+c)*(-a+b+c)*(a2-b2+c2)*(-a2+b2+c2);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X35(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   /* end vars */
   let v1 = a2*(a2-b2-b*c-c2);
   let v2 = b2*(-a2+b2-a*c-c2);
   let v3 = c2*(-a2-a*b-b2+c2);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X36(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   /* end vars */
   let v1 = a2*(a2-b2+b*c-c2);
   let v2 = b2*(-a2+b2+a*c-c2);
   let v3 = c2*(-a2+a*b-b2+c2);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X37(orbit) {
   /* begin vars */

   /* end vars */
   let v1 = a*(b+c);
   let v2 = b*(a+c);
   let v3 = (a+b)*c;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X38(orbit) {
   /* begin vars */
   let a2=a*a;
   let c2=c*c;
   let b2=b*b;
   /* end vars */
   let v1 = a*(b2+c2);
   let v2 = b*(a2+c2);
   let v3 = (a2+b2)*c;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X39(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   /* end vars */
   let v1 = a2*(b2+c2);
   let v2 = b2*(a2+c2);
   let v3 = (a2+b2)*c2;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X40(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let c3=c2*c;
   let b3=b2*b;
   let a3=a2*a;
   /* end vars */
   let v1 = a*(a3+a2*b-a*b2-b3+a2*c-2*a*b*c+b2*c-a*c2+b*c2-c3);
   let v2 = b*(-a3-a2*b+a*b2+b3+a2*c-2*a*b*c+b2*c+a*c2-b*c2-c3);
   let v3 = c*(-a3+a2*b+a*b2-b3-a2*c-2*a*b*c-b2*c+a*c2+b*c2+c3);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X41(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let c3=c2*c;
   let b3=b2*b;
   let a3=a2*a;
   /* end vars */
   let v1 = a3*(a-b-c);
   let v2 = b3*(-a+b-c);
   let v3 = (-a-b+c)*c3;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X42(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   /* end vars */
   let v1 = a2*(b+c);
   let v2 = b2*(a+c);
   let v3 = (a+b)*c2;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X43(orbit) {
   /* begin vars */

   /* end vars */
   let v1 = a*(a*b+a*c-b*c);
   let v2 = b*(a*b-a*c+b*c);
   let v3 = c*(-(a*b)+a*c+b*c);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X44(orbit) {
   /* begin vars */

   /* end vars */
   let v1 = a*(2*a-b-c);
   let v2 = b*(-a+2*b-c);
   let v3 = c*(-a-b+2*c);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X45(orbit) {
   /* begin vars */

   /* end vars */
   let v1 = a*(a-2*b-2*c);
   let v2 = b*(-2*a+b-2*c);
   let v3 = c*(-2*a-2*b+c);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X46(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let c3=c2*c;
   let b3=b2*b;
   let a3=a2*a;
   /* end vars */
   let v1 = a*(a3+a2*b-a*b2-b3+a2*c+b2*c-a*c2+b*c2-c3);
   let v2 = b*(-a3-a2*b+a*b2+b3+a2*c+b2*c+a*c2-b*c2-c3);
   let v3 = c*(-a3+a2*b+a*b2-b3-a2*c-b2*c+a*c2+b*c2+c3);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X47(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let c3=c2*c;
   let b3=b2*b;
   let c4=c2*c2;
   let b4=b2*b2;
   let a4=a2*a2;
   let a3=a2*a;
   /* end vars */
   let v1 = a3*(a4-2*a2*b2+b4-2*a2*c2+c4);
   let v2 = b3*(a4-2*a2*b2+b4-2*b2*c2+c4);
   let v3 = c3*(a4+b4-2*a2*c2-2*b2*c2+c4);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X48(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let c3=c2*c;
   let b3=b2*b;
   let a3=a2*a;
   /* end vars */
   let v1 = a3*(a2-b2-c2);
   let v2 = b3*(-a2+b2-c2);
   let v3 = (-a2-b2+c2)*c3;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X49(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let c4=c2*c2;
   let b4=b2*b2;
   let a4=a2*a2;
   /* end vars */
   let v1 = a4*(a2-b2-c2)*(a4-2*a2*b2+b4-2*a2*c2-b2*c2+c4);
   let v2 = b4*(-a2+b2-c2)*(a4-2*a2*b2+b4-a2*c2-2*b2*c2+c4);
   let v3 = (-a2-b2+c2)*c4*(a4-a2*b2+b4-2*a2*c2-2*b2*c2+c4);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X50(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let c4=c2*c2;
   let b4=b2*b2;
   let a4=a2*a2;
   /* end vars */
   let v1 = a4*(a2-b2-b*c-c2)*(a2-b2+b*c-c2);
   let v2 = b4*(-a2+b2-a*c-c2)*(-a2+b2+a*c-c2);
   let v3 = (-a2-a*b-b2+c2)*(-a2+a*b-b2+c2)*c4;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X51(orbit) {
   /* begin vars */
   let a2=a*a;
   let c2=c*c;
   let b2=b*b;
   let a4=a2*a2;
   let c4=c2*c2;
   let b4=b2*b2;
   /* end vars */
   let v1 = a2*(a2*b2-b4+a2*c2+2*b2*c2-c4);
   let v2 = b2*(-a4+a2*b2+2*a2*c2+b2*c2-c4);
   let v3 = c2*(-a4+2*a2*b2-b4+a2*c2+b2*c2);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X52(orbit) {
   /* begin vars */
   let a2=a*a;
   let c2=c*c;
   let b2=b*b;
   let a4=a2*a2;
   let c4=c2*c2;
   let b4=b2*b2;
   /* end vars */
   let v1 = a2*(a2*b2-b4+a2*c2+2*b2*c2-c4)*(a4-2*a2*b2+b4-2*a2*c2+c4);
   let v2 = b2*(-a4+a2*b2+2*a2*c2+b2*c2-c4)*(a4-2*a2*b2+b4-2*b2*c2+c4);
   let v3 = c2*(-a4+2*a2*b2-b4+a2*c2+b2*c2)*(a4+b4-2*a2*c2-2*b2*c2+c4);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X53(orbit) {
   /* begin vars */
   let a2=a*a;
   let c2=c*c;
   let b2=b*b;
   let a4=a2*a2;
   let c4=c2*c2;
   let b4=b2*b2;
   /* end vars */
   let v1 = (a2+b2-c2)*(a2-b2+c2)*(a2*b2-b4+a2*c2+2*b2*c2-c4);
   let v2 = (a2+b2-c2)*(-a2+b2+c2)*(-a4+a2*b2+2*a2*c2+b2*c2-c4);
   let v3 = (a2-b2+c2)*(-a2+b2+c2)*(-a4+2*a2*b2-b4+a2*c2+b2*c2);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X54(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let c4=c2*c2;
   let b4=b2*b2;
   let a4=a2*a2;
   /* end vars */
   let v1 = a2*(a4-2*a2*b2+b4-a2*c2-b2*c2)*(a4-a2*b2-2*a2*c2-b2*c2+c4);
   let v2 = b2*(a4-2*a2*b2+b4-a2*c2-b2*c2)*(-(a2*b2)+b4-a2*c2-2*b2*c2+c4);
   let v3 = c2*(-(a2*b2)+b4-a2*c2-2*b2*c2+c4)*(a4-a2*b2-2*a2*c2-b2*c2+c4);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X55(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   /* end vars */
   let v1 = a2*(a-b-c);
   let v2 = b2*(-a+b-c);
   let v3 = (-a-b+c)*c2;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X56(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   /* end vars */
   let v1 = a2*(a+b-c)*(a-b+c);
   let v2 = b2*(a+b-c)*(-a+b+c);
   let v3 = (a-b+c)*(-a+b+c)*c2;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X57(orbit) {
   /* begin vars */

   /* end vars */
   let v1 = a*(a+b-c)*(a-b+c);
   let v2 = b*(a+b-c)*(-a+b+c);
   let v3 = c*(a-b+c)*(-a+b+c);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X58(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   /* end vars */
   let v1 = a2*(a+b)*(a+c);
   let v2 = (a+b)*b2*(b+c);
   let v3 = (a+c)*(b+c)*c2;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X59(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   /* end vars */
   let v1 = a2*(a-b)*(a-b)*(a-c)*(a-c)*(a+b-c)*(a-b+c);
   let v2 = (-a+b)*(-a+b)*b2*(b-c)*(b-c)*(a+b-c)*(-a+b+c);
   let v3 = (-a+c)*(-a+c)*(-b+c)*(-b+c)*(a-b+c)*(-a+b+c)*c2;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X60(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   /* end vars */
   let v1 = a2*(a+b)*(a+b)*(a-b-c)*(a+c)*(a+c);
   let v2 = (a+b)*(a+b)*b2*(-a+b-c)*(b+c)*(b+c);
   let v3 = (a+c)*(a+c)*(-a-b+c)*(b+c)*(b+c)*c2;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X61(orbit) {
   /* begin vars */
   let area=triAreaHeron(a,b,c);
   let S=2*area;
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   /* end vars */
   let v1 = a2*(a2-b2-c2-2*S*sqrt3);
   let v2 = b2*(-a2+b2-c2-2*S*sqrt3);
   let v3 = c2*(-a2-b2+c2-2*S*sqrt3);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X62(orbit) {
   /* begin vars */
   let area=triAreaHeron(a,b,c);
   let S=2*area;
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   /* end vars */
   let v1 = a2*(a2-b2-c2+2*S*sqrt3);
   let v2 = b2*(-a2+b2-c2+2*S*sqrt3);
   let v3 = c2*(-a2-b2+c2+2*S*sqrt3);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X63(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   /* end vars */
   let v1 = a*(a2-b2-c2);
   let v2 = b*(-a2+b2-c2);
   let v3 = c*(-a2-b2+c2);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X64(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let c4=c2*c2;
   let b4=b2*b2;
   let a4=a2*a2;
   /* end vars */
   let v1 = a2*(a4-2*a2*b2+b4+2*a2*c2+2*b2*c2-3*c4)*(a4+2*a2*b2-3*b4-2*a2*c2+2*b2*c2+c4);
   let v2 = b2*(a4-2*a2*b2+b4+2*a2*c2+2*b2*c2-3*c4)*(-3*a4+2*a2*b2+b4+2*a2*c2-2*b2*c2+c4);
   let v3 = c2*(-3*a4+2*a2*b2+b4+2*a2*c2-2*b2*c2+c4)*(a4+2*a2*b2-3*b4-2*a2*c2+2*b2*c2+c4);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X65(orbit) {
   /* begin vars */

   /* end vars */
   let v1 = a*(a+b-c)*(a-b+c)*(b+c);
   let v2 = b*(a+b-c)*(a+c)*(-a+b+c);
   let v3 = (a+b)*c*(a-b+c)*(-a+b+c);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X66(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let c4=c2*c2;
   let b4=b2*b2;
   let a4=a2*a2;
   /* end vars */
   let v1 = (a4+b4-c4)*(a4-b4+c4);
   let v2 = (a4+b4-c4)*(-a4+b4+c4);
   let v3 = (a4-b4+c4)*(-a4+b4+c4);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X67(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let c4=c2*c2;
   let b4=b2*b2;
   let a4=a2*a2;
   /* end vars */
   let v1 = (a4-a2*b2+b4-c4)*(a4-b4-a2*c2+c4);
   let v2 = (a4-a2*b2+b4-c4)*(-a4+b4-b2*c2+c4);
   let v3 = (a4-b4-a2*c2+c4)*(-a4+b4-b2*c2+c4);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X68(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let c4=c2*c2;
   let b4=b2*b2;
   let a4=a2*a2;
   /* end vars */
   let v1 = (a2-b2-c2)*(a4-2*a2*b2+b4-2*b2*c2+c4)*(a4+b4-2*a2*c2-2*b2*c2+c4);
   let v2 = (-a2+b2-c2)*(a4-2*a2*b2+b4-2*a2*c2+c4)*(a4+b4-2*a2*c2-2*b2*c2+c4);
   let v3 = (-a2-b2+c2)*(a4-2*a2*b2+b4-2*a2*c2+c4)*(a4-2*a2*b2+b4-2*b2*c2+c4);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X69(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   /* end vars */
   let v1 = -a2+b2+c2;
   let v2 = a2-b2+c2;
   let v3 = a2+b2-c2;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X70(orbit) {
   /* begin vars */
   let c2=c*c;
   let c4=c2*c2;
   let c6=c2*c4;
   let b2=b*b;
   let b4=b2*b2;
   let b6=b2*b4;
   let a2=a*a;
   let a4=a2*a2;
   let a6=a2*a4;
   let c8=c2*c6;
   let b8=b2*b6;
   let a8=a2*a6;
   /* end vars */
   let v1 = (a8-2*a6*b2+2*a4*b4-2*a2*b6+b8-2*a6*c2-2*b6*c2+2*a2*c6+2*b2*c6-c8)*(a8-2*a6*b2+2*a2*b6-b8-2*a6*c2+2*b6*c2+2*a4*c4-2*a2*c6-2*b2*c6+c8);
   let v2 = (a8-2*a6*b2+2*a4*b4-2*a2*b6+b8-2*a6*c2-2*b6*c2+2*a2*c6+2*b2*c6-c8)*(-a8+2*a6*b2-2*a2*b6+b8+2*a6*c2-2*b6*c2+2*b4*c4-2*a2*c6-2*b2*c6+c8);
   let v3 = (a8-2*a6*b2+2*a2*b6-b8-2*a6*c2+2*b6*c2+2*a4*c4-2*a2*c6-2*b2*c6+c8)*(-a8+2*a6*b2-2*a2*b6+b8+2*a6*c2-2*b6*c2+2*b4*c4-2*a2*c6-2*b2*c6+c8);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X71(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   /* end vars */
   let v1 = a2*(b+c)*(a2-b2-c2);
   let v2 = b2*(a+c)*(-a2+b2-c2);
   let v3 = (a+b)*c2*(-a2-b2+c2);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X72(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   /* end vars */
   let v1 = a*(b+c)*(a2-b2-c2);
   let v2 = b*(a+c)*(-a2+b2-c2);
   let v3 = (a+b)*c*(-a2-b2+c2);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X73(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   /* end vars */
   let v1 = a2*(a+b-c)*(a-b+c)*(b+c)*(a2-b2-c2);
   let v2 = b2*(a+b-c)*(a+c)*(-a+b+c)*(-a2+b2-c2);
   let v3 = (a+b)*(a-b+c)*(-a+b+c)*c2*(-a2-b2+c2);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X74(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let c4=c2*c2;
   let b4=b2*b2;
   let a4=a2*a2;
   /* end vars */
   let v1 = a2*(a4-2*a2*b2+b4+a2*c2+b2*c2-2*c4)*(a4+a2*b2-2*b4-2*a2*c2+b2*c2+c4);
   let v2 = b2*(a4-2*a2*b2+b4+a2*c2+b2*c2-2*c4)*(-2*a4+a2*b2+b4+a2*c2-2*b2*c2+c4);
   let v3 = c2*(-2*a4+a2*b2+b4+a2*c2-2*b2*c2+c4)*(a4+a2*b2-2*b4-2*a2*c2+b2*c2+c4);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X75(orbit) {
   /* begin vars */

   /* end vars */
   let v1 = b*c;
   let v2 = a*c;
   let v3 = a*b;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X76(orbit) {
   /* begin vars */
   let a2=a*a;
   let c2=c*c;
   let b2=b*b;
   /* end vars */
   let v1 = b2*c2;
   let v2 = a2*c2;
   let v3 = a2*b2;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X77(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   /* end vars */
   let v1 = a*(a+b-c)*(a-b+c)*(a2-b2-c2);
   let v2 = b*(a+b-c)*(-a+b+c)*(-a2+b2-c2);
   let v3 = c*(a-b+c)*(-a+b+c)*(-a2-b2+c2);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X78(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   /* end vars */
   let v1 = a*(a-b-c)*(a2-b2-c2);
   let v2 = b*(-a+b-c)*(-a2+b2-c2);
   let v3 = c*(-a-b+c)*(-a2-b2+c2);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X79(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   /* end vars */
   let v1 = (a2+a*b+b2-c2)*(a2-b2+a*c+c2);
   let v2 = (a2+a*b+b2-c2)*(-a2+b2+b*c+c2);
   let v3 = (a2-b2+a*c+c2)*(-a2+b2+b*c+c2);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X80(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   /* end vars */
   let v1 = (a2-a*b+b2-c2)*(a2-b2-a*c+c2);
   let v2 = (a2-a*b+b2-c2)*(-a2+b2-b*c+c2);
   let v3 = (a2-b2-a*c+c2)*(-a2+b2-b*c+c2);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X81(orbit) {
   /* begin vars */

   /* end vars */
   let v1 = a*(a+b)*(a+c);
   let v2 = b*(a+b)*(b+c);
   let v3 = c*(a+c)*(b+c);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X82(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   /* end vars */
   let v1 = a*(a2+b2)*(a2+c2);
   let v2 = b*(a2+b2)*(b2+c2);
   let v3 = c*(a2+c2)*(b2+c2);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X83(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   /* end vars */
   let v1 = (a2+b2)*(a2+c2);
   let v2 = (a2+b2)*(b2+c2);
   let v3 = (a2+c2)*(b2+c2);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X84(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let c3=c2*c;
   let b3=b2*b;
   let a3=a2*a;
   /* end vars */
   let v1 = a*(a3-a2*b-a*b2+b3+a2*c+2*a*b*c+b2*c-a*c2-b*c2-c3)*(a3+a2*b-a*b2-b3-a2*c+2*a*b*c-b2*c-a*c2+b*c2+c3);
   let v2 = b*(a3-a2*b-a*b2+b3+a2*c+2*a*b*c+b2*c-a*c2-b*c2-c3)*(-a3-a2*b+a*b2+b3-a2*c+2*a*b*c-b2*c+a*c2-b*c2+c3);
   let v3 = c*(-a3-a2*b+a*b2+b3-a2*c+2*a*b*c-b2*c+a*c2-b*c2+c3)*(a3+a2*b-a*b2-b3-a2*c+2*a*b*c-b2*c-a*c2+b*c2+c3);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X85(orbit) {
   /* begin vars */

   /* end vars */
   let v1 = b*(-a+b-c)*(a+b-c)*c;
   let v2 = a*c*(-a-b+c)*(-a+b+c);
   let v3 = a*b*(a-b-c)*(a-b+c);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X86(orbit) {
   /* begin vars */

   /* end vars */
   let v1 = (a+b)*(a+c);
   let v2 = (a+b)*(b+c);
   let v3 = (a+c)*(b+c);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X87(orbit) {
   /* begin vars */

   /* end vars */
   let v1 = a*(a*b-a*c-b*c)*(a*b-a*c+b*c);
   let v2 = b*(-(a*b)-a*c+b*c)*(-(a*b)+a*c+b*c);
   let v3 = c*(-(a*b)+a*c-b*c)*(a*b+a*c-b*c);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X88(orbit) {
   /* begin vars */

   /* end vars */
   let v1 = a*(a+b-2*c)*(a-2*b+c);
   let v2 = b*(a+b-2*c)*(-2*a+b+c);
   let v3 = c*(a-2*b+c)*(-2*a+b+c);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X89(orbit) {
   /* begin vars */

   /* end vars */
   let v1 = a*(2*a+2*b-c)*(2*a-b+2*c);
   let v2 = b*(2*a+2*b-c)*(-a+2*b+2*c);
   let v3 = c*(2*a-b+2*c)*(-a+2*b+2*c);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X90(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let c3=c2*c;
   let b3=b2*b;
   let a3=a2*a;
   /* end vars */
   let v1 = a*(a3-a2*b-a*b2+b3+a2*c+b2*c-a*c2-b*c2-c3)*(a3+a2*b-a*b2-b3-a2*c-b2*c-a*c2+b*c2+c3);
   let v2 = b*(a3-a2*b-a*b2+b3+a2*c+b2*c-a*c2-b*c2-c3)*(-a3-a2*b+a*b2+b3-a2*c-b2*c+a*c2-b*c2+c3);
   let v3 = c*(-a3-a2*b+a*b2+b3-a2*c-b2*c+a*c2-b*c2+c3)*(a3+a2*b-a*b2-b3-a2*c-b2*c-a*c2+b*c2+c3);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X91(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let c4=c2*c2;
   let b4=b2*b2;
   let a4=a2*a2;
   /* end vars */
   let v1 = b*c*(a4-2*a2*b2+b4-2*b2*c2+c4)*(a4+b4-2*a2*c2-2*b2*c2+c4);
   let v2 = a*c*(a4-2*a2*b2+b4-2*a2*c2+c4)*(a4+b4-2*a2*c2-2*b2*c2+c4);
   let v3 = a*b*(a4-2*a2*b2+b4-2*a2*c2+c4)*(a4-2*a2*b2+b4-2*b2*c2+c4);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X92(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   /* end vars */
   let v1 = b*c*(-a2+b2-c2)*(a2+b2-c2);
   let v2 = a*c*(-a2-b2+c2)*(-a2+b2+c2);
   let v3 = a*b*(a2-b2-c2)*(a2-b2+c2);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X93(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let c4=c2*c2;
   let b4=b2*b2;
   let a4=a2*a2;
   /* end vars */
   let v1 = b2*(-a2+b2-c2)*(a2+b2-c2)*c2*(a4-a2*b2+b4-2*a2*c2-2*b2*c2+c4)*(a4-2*a2*b2+b4-a2*c2-2*b2*c2+c4);
   let v2 = a2*c2*(-a2-b2+c2)*(-a2+b2+c2)*(a4-a2*b2+b4-2*a2*c2-2*b2*c2+c4)*(a4-2*a2*b2+b4-2*a2*c2-b2*c2+c4);
   let v3 = a2*b2*(a2-b2-c2)*(a2-b2+c2)*(a4-2*a2*b2+b4-a2*c2-2*b2*c2+c4)*(a4-2*a2*b2+b4-2*a2*c2-b2*c2+c4);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X94(orbit) {
   /* begin vars */
   let c2=c*c;
   let a2=a*a;
   let b2=b*b;
   /* end vars */
   let v1 = b2*(a2-a*b+b2-c2)*(a2+a*b+b2-c2)*(-a2+b2-a*c-c2)*(-a2+b2+a*c-c2)*c2;
   let v2 = a2*c2*(-a2-a*b-b2+c2)*(-a2+a*b-b2+c2)*(-a2+b2-b*c+c2)*(-a2+b2+b*c+c2);
   let v3 = a2*b2*(a2-b2-b*c-c2)*(a2-b2+b*c-c2)*(a2-b2-a*c+c2)*(a2-b2+a*c+c2);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X95(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let c4=c2*c2;
   let b4=b2*b2;
   let a4=a2*a2;
   /* end vars */
   let v1 = (a4-2*a2*b2+b4-a2*c2-b2*c2)*(a4-a2*b2-2*a2*c2-b2*c2+c4);
   let v2 = (a4-2*a2*b2+b4-a2*c2-b2*c2)*(-(a2*b2)+b4-a2*c2-2*b2*c2+c4);
   let v3 = (-(a2*b2)+b4-a2*c2-2*b2*c2+c4)*(a4-a2*b2-2*a2*c2-b2*c2+c4);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X96(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let c4=c2*c2;
   let b4=b2*b2;
   let a4=a2*a2;
   /* end vars */
   let v1 = (a4-2*a2*b2+b4-a2*c2-b2*c2)*(a4-2*a2*b2+b4-2*b2*c2+c4)*(a4+b4-2*a2*c2-2*b2*c2+c4)*(a4-a2*b2-2*a2*c2-b2*c2+c4);
   let v2 = (a4-2*a2*b2+b4-a2*c2-b2*c2)*(a4-2*a2*b2+b4-2*a2*c2+c4)*(a4+b4-2*a2*c2-2*b2*c2+c4)*(-(a2*b2)+b4-a2*c2-2*b2*c2+c4);
   let v3 = (a4-2*a2*b2+b4-2*a2*c2+c4)*(a4-2*a2*b2+b4-2*b2*c2+c4)*(-(a2*b2)+b4-a2*c2-2*b2*c2+c4)*(a4-a2*b2-2*a2*c2-b2*c2+c4);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X97(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let c4=c2*c2;
   let b4=b2*b2;
   let a4=a2*a2;
   /* end vars */
   let v1 = a2*(a2-b2-c2)*(a4-2*a2*b2+b4-a2*c2-b2*c2)*(a4-a2*b2-2*a2*c2-b2*c2+c4);
   let v2 = b2*(-a2+b2-c2)*(a4-2*a2*b2+b4-a2*c2-b2*c2)*(-(a2*b2)+b4-a2*c2-2*b2*c2+c4);
   let v3 = c2*(-a2-b2+c2)*(-(a2*b2)+b4-a2*c2-2*b2*c2+c4)*(a4-a2*b2-2*a2*c2-b2*c2+c4);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X98(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let c4=c2*c2;
   let b4=b2*b2;
   let a4=a2*a2;
   /* end vars */
   let v1 = (a4+b4-a2*c2-b2*c2)*(a4-a2*b2-b2*c2+c4);
   let v2 = (a4+b4-a2*c2-b2*c2)*(-(a2*b2)+b4-a2*c2+c4);
   let v3 = (-(a2*b2)+b4-a2*c2+c4)*(a4-a2*b2-b2*c2+c4);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X99(orbit) {
   /* begin vars */

   /* end vars */
   let v1 = (a-b)*(a+b)*(a-c)*(a+c);
   let v2 = (-a+b)*(a+b)*(b-c)*(b+c);
   let v3 = (-a+c)*(a+c)*(-b+c)*(b+c);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X100(orbit) {
   /* begin vars */

   /* end vars */
   let v1 = a*(a-b)*(a-c);
   let v2 = b*(-a+b)*(b-c);
   let v3 = c*(-a+c)*(-b+c);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X101(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   /* end vars */
   let v1 = a2*(a-b)*(a-c);
   let v2 = (-a+b)*b2*(b-c);
   let v3 = (-a+c)*(-b+c)*c2;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X102(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let c4=c2*c2;
   let c3=c2*c;
   let b3=b2*b;
   let a3=a2*a;
   let b4=b2*b2;
   let a4=a2*a2;
   /* end vars */
   let v1 = a2*(a4-2*a2*b2+b4-a3*c+a2*b*c+a*b2*c-b3*c+a2*c2-2*a*b*c2+b2*c2+a*c3+b*c3-2*c4)*(a4-a3*b+a2*b2+a*b3-2*b4+a2*b*c-2*a*b2*c+b3*c-2*a2*c2+a*b*c2+b2*c2-b*c3+c4);
   let v2 = b2*(a4-2*a2*b2+b4-a3*c+a2*b*c+a*b2*c-b3*c+a2*c2-2*a*b*c2+b2*c2+a*c3+b*c3-2*c4)*(-2*a4+a3*b+a2*b2-a*b3+b4+a3*c-2*a2*b*c+a*b2*c+a2*c2+a*b*c2-2*b2*c2-a*c3+c4);
   let v3 = c2*(-2*a4+a3*b+a2*b2-a*b3+b4+a3*c-2*a2*b*c+a*b2*c+a2*c2+a*b*c2-2*b2*c2-a*c3+c4)*(a4-a3*b+a2*b2+a*b3-2*b4+a2*b*c-2*a*b2*c+b3*c-2*a2*c2+a*b*c2+b2*c2-b*c3+c4);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X103(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let c3=c2*c;
   let b3=b2*b;
   let a3=a2*a;
   /* end vars */
   let v1 = a2*(a3-a2*b-a*b2+b3+a*c2+b*c2-2*c3)*(a3+a*b2-2*b3-a2*c+b2*c-a*c2+c3);
   let v2 = b2*(a3-a2*b-a*b2+b3+a*c2+b*c2-2*c3)*(-2*a3+a2*b+b3+a2*c-b2*c-b*c2+c3);
   let v3 = c2*(a3+a*b2-2*b3-a2*c+b2*c-a*c2+c3)*(-2*a3+a2*b+b3+a2*c-b2*c-b*c2+c3);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X104(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let c3=c2*c;
   let b3=b2*b;
   let a3=a2*a;
   /* end vars */
   let v1 = a*(a3-a2*b-a*b2+b3+2*a*b*c-a*c2-b*c2)*(a3-a*b2-a2*c+2*a*b*c-b2*c-a*c2+c3);
   let v2 = b*(a3-a2*b-a*b2+b3+2*a*b*c-a*c2-b*c2)*(-(a2*b)+b3-a2*c+2*a*b*c-b2*c-b*c2+c3);
   let v3 = c*(a3-a*b2-a2*c+2*a*b*c-b2*c-a*c2+c3)*(-(a2*b)+b3-a2*c+2*a*b*c-b2*c-b*c2+c3);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X105(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   /* end vars */
   let v1 = a*(a2+b2-a*c-b*c)*(a2-a*b-b*c+c2);
   let v2 = b*(a2+b2-a*c-b*c)*(-(a*b)+b2-a*c+c2);
   let v3 = c*(-(a*b)+b2-a*c+c2)*(a2-a*b-b*c+c2);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X106(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   /* end vars */
   let v1 = a2*(a+b-2*c)*(a-2*b+c);
   let v2 = b2*(a+b-2*c)*(-2*a+b+c);
   let v3 = (a-2*b+c)*(-2*a+b+c)*c2;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X107(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   /* end vars */
   let v1 = (a-b)*(a+b)*(a-c)*(a+c)*(a2+b2-c2)*(a2+b2-c2)*(a2-b2+c2)*(a2-b2+c2);
   let v2 = (-a+b)*(a+b)*(b-c)*(b+c)*(a2+b2-c2)*(a2+b2-c2)*(-a2+b2+c2)*(-a2+b2+c2);
   let v3 = (-a+c)*(a+c)*(-b+c)*(b+c)*(a2-b2+c2)*(a2-b2+c2)*(-a2+b2+c2)*(-a2+b2+c2);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X108(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   /* end vars */
   let v1 = a*(a-b)*(a-c)*(a+b-c)*(a-b+c)*(a2+b2-c2)*(a2-b2+c2);
   let v2 = b*(-a+b)*(b-c)*(a+b-c)*(-a+b+c)*(a2+b2-c2)*(-a2+b2+c2);
   let v3 = c*(-a+c)*(-b+c)*(a-b+c)*(-a+b+c)*(a2-b2+c2)*(-a2+b2+c2);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X109(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   /* end vars */
   let v1 = a2*(a-b)*(a-c)*(a+b-c)*(a-b+c);
   let v2 = (-a+b)*b2*(b-c)*(a+b-c)*(-a+b+c);
   let v3 = (-a+c)*(-b+c)*(a-b+c)*(-a+b+c)*c2;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X110(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   /* end vars */
   let v1 = a2*(a-b)*(a+b)*(a-c)*(a+c);
   let v2 = (-a+b)*(a+b)*b2*(b-c)*(b+c);
   let v3 = (-a+c)*(a+c)*(-b+c)*(b+c)*c2;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X111(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   /* end vars */
   let v1 = a2*(a2+b2-2*c2)*(a2-2*b2+c2);
   let v2 = b2*(a2+b2-2*c2)*(-2*a2+b2+c2);
   let v3 = c2*(a2-2*b2+c2)*(-2*a2+b2+c2);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X112(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   /* end vars */
   let v1 = a2*(a-b)*(a+b)*(a-c)*(a+c)*(a2+b2-c2)*(a2-b2+c2);
   let v2 = (-a+b)*(a+b)*b2*(b-c)*(b+c)*(a2+b2-c2)*(-a2+b2+c2);
   let v3 = (-a+c)*(a+c)*(-b+c)*(b+c)*c2*(a2-b2+c2)*(-a2+b2+c2);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X113(orbit) {
   /* begin vars */
   let a2=a*a;
   let a4=a2*a2;
   let c2=c*c;
   let c4=c2*c2;
   let b2=b*b;
   let b4=b2*b2;
   let a6=a2*a4;
   let c6=c2*c4;
   let b6=b2*b4;
   /* end vars */
   let v1 = (2*a4-a2*b2-b4-a2*c2+2*b2*c2-c4)*(a4*b2-2*a2*b4+b6+a4*c2+2*a2*b2*c2-b4*c2-2*a2*c4-b2*c4+c6);
   let v2 = (-a4-a2*b2+2*b4+2*a2*c2-b2*c2-c4)*(a6-2*a4*b2+a2*b4-a4*c2+2*a2*b2*c2+b4*c2-a2*c4-2*b2*c4+c6);
   let v3 = (-a4+2*a2*b2-b4-a2*c2-b2*c2+2*c4)*(a6-a4*b2-a2*b4+b6-2*a4*c2+2*a2*b2*c2-2*b4*c2+a2*c4+b2*c4);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X114(orbit) {
   /* begin vars */
   let a2=a*a;
   let c2=c*c;
   let b2=b*b;
   let a4=a2*a2;
   let c4=c2*c2;
   let b4=b2*b2;
   /* end vars */
   let v1 = (a2*b2-b4+a2*c2-c4)*(2*a4-a2*b2+b4-a2*c2-2*b2*c2+c4);
   let v2 = (-a4+a2*b2+b2*c2-c4)*(a4-a2*b2+2*b4-2*a2*c2-b2*c2+c4);
   let v3 = (-a4-b4+a2*c2+b2*c2)*(a4-2*a2*b2+b4-a2*c2-b2*c2+2*c4);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X115(orbit) {
   /* begin vars */

   /* end vars */
   let v1 = (b-c)*(b-c)*(b+c)*(b+c);
   let v2 = (-a+c)*(-a+c)*(a+c)*(a+c);
   let v3 = (a-b)*(a-b)*(a+b)*(a+b);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X116(orbit) {
   /* begin vars */
   let a2=a*a;
   let c2=c*c;
   let b2=b*b;
   /* end vars */
   let v1 = (b-c)*(b-c)*(-(a*b)+b2-a*c+b*c+c2);
   let v2 = (-a+c)*(-a+c)*(a2-a*b+a*c-b*c+c2);
   let v3 = (a-b)*(a-b)*(a2+a*b+b2-a*c-b*c);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X117(orbit) {
   /* begin vars */
   let a2=a*a;
   let a3=a2*a;
   let a4=a2*a2;
   let c2=c*c;
   let c4=c2*c2;
   let c3=c2*c;
   let b2=b*b;
   let b3=b2*b;
   let b4=b2*b2;
   let a5=a2*a3;
   let a6=a2*a4;
   let c6=c2*c4;
   let c5=c2*c3;
   let b5=b2*b3;
   let b6=b2*b4;
   /* end vars */
   let v1 = (2*a4-a3*b-a2*b2+a*b3-b4-a3*c+2*a2*b*c-a*b2*c-a2*c2-a*b*c2+2*b2*c2+a*c3-c4)*(a4*b2-2*a2*b4+b6-a3*b2*c+a2*b3*c+a*b4*c-b5*c+a4*c2-a3*b*c2+2*a2*b2*c2-a*b3*c2-b4*c2+a2*b*c3-a*b2*c3+2*b3*c3-2*a2*c4+a*b*c4-b2*c4-b*c5+c6);
   let v2 = (-a4+a3*b-a2*b2-a*b3+2*b4-a2*b*c+2*a*b2*c-b3*c+2*a2*c2-a*b*c2-b2*c2+b*c3-c4)*(a6-2*a4*b2+a2*b4-a5*c+a4*b*c+a3*b2*c-a2*b3*c-a4*c2-a3*b*c2+2*a2*b2*c2-a*b3*c2+b4*c2+2*a3*c3-a2*b*c3+a*b2*c3-a2*c4+a*b*c4-2*b2*c4-a*c5+c6);
   let v3 = (-a4+2*a2*b2-b4+a3*c-a2*b*c-a*b2*c+b3*c-a2*c2+2*a*b*c2-b2*c2-a*c3-b*c3+2*c4)*(a6-a5*b-a4*b2+2*a3*b3-a2*b4-a*b5+b6+a4*b*c-a3*b2*c-a2*b3*c+a*b4*c-2*a4*c2+a3*b*c2+2*a2*b2*c2+a*b3*c2-2*b4*c2-a2*b*c3-a*b2*c3+a2*c4+b2*c4);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X118(orbit) {
   /* begin vars */
   let a2=a*a;
   let a3=a2*a;
   let c2=c*c;
   let c3=c2*c;
   let b2=b*b;
   let b3=b2*b;
   let a4=a2*a2;
   let a5=a2*a3;
   let c5=c2*c3;
   let c4=c2*c2;
   let b5=b2*b3;
   let b4=b2*b2;
   /* end vars */
   let v1 = (2*a3-a2*b-b3-a2*c+b2*c+b*c2-c3)*(a3*b2-a2*b3-a*b4+b5+a3*c2+2*a*b2*c2-b3*c2-a2*c3-b2*c3-a*c4+c5);
   let v2 = (-a3-a*b2+2*b3+a2*c-b2*c+a*c2-c3)*(a5-a4*b-a3*b2+a2*b3-a3*c2+2*a2*b*c2+b3*c2-a2*c3-b2*c3-b*c4+c5);
   let v3 = (-a3+a2*b+a*b2-b3-a*c2-b*c2+2*c3)*(a5-a3*b2-a2*b3+b5-a4*c+2*a2*b2*c-b4*c-a3*c2-b3*c2+a2*c3+b2*c3);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X119(orbit) {
   /* begin vars */
   let a2=a*a;
   let c2=c*c;
   let b2=b*b;
   let a4=a2*a2;
   let c4=c2*c2;
   let b4=b2*b2;
   let a3=a2*a;
   let c3=c2*c;
   let b3=b2*b;
   /* end vars */
   let v1 = (a2*b-b3+a2*c-2*a*b*c+b2*c+b*c2-c3)*(a3*b-a2*b2-a*b3+b4+a3*c+a*b2*c-a2*c2+a*b*c2-2*b2*c2-a*c3+c4);
   let v2 = (-a3+a*b2+a2*c-2*a*b*c+b2*c+a*c2-c3)*(a4-a3*b-a2*b2+a*b3+a2*b*c+b3*c-2*a2*c2+a*b*c2-b2*c2-b*c3+c4);
   let v3 = (-a3+a2*b+a*b2-b3-2*a*b*c+a*c2+b*c2)*(a4-2*a2*b2+b4-a3*c+a2*b*c+a*b2*c-b3*c-a2*c2-b2*c2+a*c3+b*c3);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X120(orbit) {
   /* begin vars */
   let a2=a*a;
   let c2=c*c;
   let b2=b*b;
   let a3=a2*a;
   let c3=c2*c;
   let b3=b2*b;
   /* end vars */
   let v1 = (a*b-b2+a*c-c2)*(a2*b+b3+a2*c-2*a*b*c-b2*c-b*c2+c3);
   let v2 = (-a2+a*b+b*c-c2)*(a3+a*b2-a2*c-2*a*b*c+b2*c-a*c2+c3);
   let v3 = (-a2-b2+a*c+b*c)*(a3-a2*b-a*b2+b3-2*a*b*c+a*c2+b*c2);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X121(orbit) {
   /* begin vars */
   let a2=a*a;
   let c2=c*c;
   let b2=b*b;
   let a3=a2*a;
   let c3=c2*c;
   let b3=b2*b;
   /* end vars */
   let v1 = (2*a-b-c)*(a*b2+b3-2*b2*c+a*c2-2*b*c2+c3);
   let v2 = (-a+2*b-c)*(a3+a2*b-2*a2*c-2*a*c2+b*c2+c3);
   let v3 = (-a-b+2*c)*(a3-2*a2*b-2*a*b2+b3+a2*c+b2*c);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X122(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let c4=c2*c2;
   let b4=b2*b2;
   let a4=a2*a2;
   /* end vars */
   let v1 = (b-c)*(b-c)*(b+c)*(b+c)*(-a2+b2+c2)*(-a2+b2+c2)*(-3*a4+2*a2*b2+b4+2*a2*c2-2*b2*c2+c4);
   let v2 = (-a+c)*(-a+c)*(a+c)*(a+c)*(a2-b2+c2)*(a2-b2+c2)*(a4+2*a2*b2-3*b4-2*a2*c2+2*b2*c2+c4);
   let v3 = (a-b)*(a-b)*(a+b)*(a+b)*(a2+b2-c2)*(a2+b2-c2)*(a4-2*a2*b2+b4+2*a2*c2+2*b2*c2-3*c4);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X123(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let c4=c2*c2;
   let b4=b2*b2;
   let a4=a2*a2;
   /* end vars */
   let v1 = (b-c)*(b-c)*(-a+b+c)*(-a2+b2+c2)*(-a4+b4-2*a2*b*c+2*a*b2*c+2*a*b*c2-2*b2*c2+c4);
   let v2 = (-a+c)*(-a+c)*(a-b+c)*(a2-b2+c2)*(a4-b4+2*a2*b*c-2*a*b2*c-2*a2*c2+2*a*b*c2+c4);
   let v3 = (a-b)*(a-b)*(a+b-c)*(a2+b2-c2)*(a4-2*a2*b2+b4+2*a2*b*c+2*a*b2*c-2*a*b*c2-c4);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X124(orbit) {
   /* begin vars */
   let a2=a*a;
   let c2=c*c;
   let b2=b*b;
   let a3=a2*a;
   let c3=c2*c;
   let b3=b2*b;
   /* end vars */
   let v1 = (b-c)*(b-c)*(-a+b+c)*(-(a2*b)+b3-a2*c+a*b*c+c3);
   let v2 = (-a+c)*(-a+c)*(a-b+c)*(a3-a*b2+a*b*c-b2*c+c3);
   let v3 = (a-b)*(a-b)*(a+b-c)*(a3+b3+a*b*c-a*c2-b*c2);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X125(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   /* end vars */
   let v1 = (b-c)*(b-c)*(b+c)*(b+c)*(-a2+b2+c2);
   let v2 = (-a+c)*(-a+c)*(a+c)*(a+c)*(a2-b2+c2);
   let v3 = (a-b)*(a-b)*(a+b)*(a+b)*(a2+b2-c2);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X126(orbit) {
   /* begin vars */
   let a2=a*a;
   let c2=c*c;
   let b2=b*b;
   let a4=a2*a2;
   let c4=c2*c2;
   let b4=b2*b2;
   /* end vars */
   let v1 = (2*a2-b2-c2)*(a2*b2+b4+a2*c2-4*b2*c2+c4);
   let v2 = (-a2+2*b2-c2)*(a4+a2*b2-4*a2*c2+b2*c2+c4);
   let v3 = (-a2-b2+2*c2)*(a4-4*a2*b2+b4+a2*c2+b2*c2);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X127(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let c4=c2*c2;
   let b4=b2*b2;
   let a4=a2*a2;
   /* end vars */
   let v1 = (b-c)*(b-c)*(b+c)*(b+c)*(-a2+b2+c2)*(-a4+b4+c4);
   let v2 = (-a+c)*(-a+c)*(a+c)*(a+c)*(a2-b2+c2)*(a4-b4+c4);
   let v3 = (a-b)*(a-b)*(a+b)*(a+b)*(a2+b2-c2)*(a4+b4-c4);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X128(orbit) {
   /* begin vars */
   let c2=c*c;
   let c4=c2*c2;
   let c6=c2*c4;
   let b2=b*b;
   let b4=b2*b2;
   let b6=b2*b4;
   let a2=a*a;
   let a4=a2*a2;
   let a6=a2*a4;
   let c8=c2*c6;
   let b8=b2*b6;
   let a8=a2*a6;
   /* end vars */
   let v1 = (a2-b2-b*c-c2)*(a2-b2+b*c-c2)*(a2*b2-b4+a2*c2+2*b2*c2-c4)*(2*a8-4*a6*b2+3*a4*b4-2*a2*b6+b8-4*a6*c2+2*a2*b4*c2-4*b6*c2+3*a4*c4+2*a2*b2*c4+6*b4*c4-2*a2*c6-4*b2*c6+c8);
   let v2 = (-a2+b2-a*c-c2)*(-a2+b2+a*c-c2)*(-a4+a2*b2+2*a2*c2+b2*c2-c4)*(a8-2*a6*b2+3*a4*b4-4*a2*b6+2*b8-4*a6*c2+2*a4*b2*c2-4*b6*c2+6*a4*c4+2*a2*b2*c4+3*b4*c4-4*a2*c6-2*b2*c6+c8);
   let v3 = (-a2-a*b-b2+c2)*(-a2+a*b-b2+c2)*(-a4+2*a2*b2-b4+a2*c2+b2*c2)*(a8-4*a6*b2+6*a4*b4-4*a2*b6+b8-2*a6*c2+2*a4*b2*c2+2*a2*b4*c2-2*b6*c2+3*a4*c4+3*b4*c4-4*a2*c6-4*b2*c6+2*c8);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X129(orbit) {
   /* begin vars */
   let c2=c*c;
   let c4=c2*c2;
   let c6=c2*c4;
   let b2=b*b;
   let b4=b2*b2;
   let b6=b2*b4;
   let a2=a*a;
   let a4=a2*a2;
   let a6=a2*a4;
   let c8=c2*c6;
   let b8=b2*b6;
   let a8=a2*a6;
   /* end vars */
   let v1 = a2*(a2*b2-b4+a2*c2+2*b2*c2-c4)*(a8-2*a6*b2+a4*b4-2*a6*c2+a4*b2*c2+b6*c2+a4*c4-2*b4*c4+b2*c6)*(-4*a2*Math.pow(b,10)+Math.pow(b,12)+a8*b4-4*a6*b6+6*a4*b8-4*a2*Math.pow(c,10)-2*b2*Math.pow(c,10)+Math.pow(c,12)-2*Math.pow(b,10)*c2-2*a4*b6*c2+4*a2*b8*c2+a8*c4-2*a4*b4*c4+b8*c4-4*a6*c6-2*a4*b2*c6+6*a4*c8+4*a2*b2*c8+b4*c8);
   let v2 = b2*(-a4+a2*b2+2*a2*c2+b2*c2-c4)*(a4*b4-2*a2*b6+b8+a6*c2+a2*b4*c2-2*b6*c2-2*a4*c4+b4*c4+a2*c6)*(Math.pow(a,12)-4*Math.pow(a,10)*b2+6*a8*b4-4*a6*b6+a4*b8-2*a2*Math.pow(c,10)-4*b2*Math.pow(c,10)+Math.pow(c,12)-2*Math.pow(a,10)*c2+4*a8*b2*c2-2*a6*b4*c2+a8*c4-2*a4*b4*c4+b8*c4-2*a2*b4*c6-4*b6*c6+a4*c8+4*a2*b2*c8+6*b4*c8);
   let v3 = c2*(-a4+2*a2*b2-b4+a2*c2+b2*c2)*(a6*b2-2*a4*b4+a2*b6+a4*c4+a2*b2*c4+b4*c4-2*a2*c6-2*b2*c6+c8)*(Math.pow(a,12)-2*a2*Math.pow(b,10)+Math.pow(b,12)-2*Math.pow(a,10)*b2+a8*b4+a4*b8-4*Math.pow(a,10)*c2-4*Math.pow(b,10)*c2+4*a8*b2*c2+4*a2*b8*c2+6*a8*c4-2*a6*b2*c4-2*a4*b4*c4-2*a2*b6*c4+6*b8*c4-4*a6*c6-4*b6*c6+a4*c8+b4*c8);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X130(orbit) {
   /* begin vars */
   let c2=c*c;
   let c4=c2*c2;
   let c6=c2*c4;
   let b2=b*b;
   let b4=b2*b2;
   let b6=b2*b4;
   let a2=a*a;
   let a4=a2*a2;
   let a6=a2*a4;
   let c8=c2*c6;
   let b8=b2*b6;
   let a8=a2*a6;
   /* end vars */
   let v1 = a2*(b-c)*(b-c)*(b+c)*(b+c)*(a2-b2-c2)*(a2-b2-c2)*(a2*b2-b4+a2*c2+2*b2*c2-c4)*(a8-2*a6*b2+a4*b4-2*a6*c2+3*a4*b2*c2-b6*c2+a4*c4+2*b4*c4-b2*c6);
   let v2 = b2*(-a+c)*(-a+c)*(a+c)*(a+c)*(-a2+b2-c2)*(-a2+b2-c2)*(-a4+a2*b2+2*a2*c2+b2*c2-c4)*(a4*b4-2*a2*b6+b8-a6*c2+3*a2*b4*c2-2*b6*c2+2*a4*c4+b4*c4-a2*c6);
   let v3 = (a-b)*(a-b)*(a+b)*(a+b)*c2*(-a2-b2+c2)*(-a2-b2+c2)*(-a4+2*a2*b2-b4+a2*c2+b2*c2)*(-(a6*b2)+2*a4*b4-a2*b6+a4*c4+3*a2*b2*c4+b4*c4-2*a2*c6-2*b2*c6+c8);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X131(orbit) {
   /* begin vars */
   let c2=c*c;
   let c4=c2*c2;
   let c6=c2*c4;
   let b2=b*b;
   let b4=b2*b2;
   let b6=b2*b4;
   let a2=a*a;
   let a4=a2*a2;
   let a6=a2*a4;
   let c8=c2*c6;
   let b8=b2*b6;
   let a8=a2*a6;
   /* end vars */
   let v1 = (a2-b2-c2)*(a4*b2-2*a2*b4+b6+a4*c2+2*a2*b2*c2-b4*c2-2*a2*c4-b2*c4+c6)*(2*a8-3*a6*b2+a4*b4-a2*b6+b8-3*a6*c2+2*a4*b2*c2+a2*b4*c2-4*b6*c2+a4*c4+a2*b2*c4+6*b4*c4-a2*c6-4*b2*c6+c8);
   let v2 = (-a2+b2-c2)*(a6-2*a4*b2+a2*b4-a4*c2+2*a2*b2*c2+b4*c2-a2*c4-2*b2*c4+c6)*(a8-a6*b2+a4*b4-3*a2*b6+2*b8-4*a6*c2+a4*b2*c2+2*a2*b4*c2-3*b6*c2+6*a4*c4+a2*b2*c4+b4*c4-4*a2*c6-b2*c6+c8);
   let v3 = (-a2-b2+c2)*(a6-a4*b2-a2*b4+b6-2*a4*c2+2*a2*b2*c2-2*b4*c2+a2*c4+b2*c4)*(a8-4*a6*b2+6*a4*b4-4*a2*b6+b8-a6*c2+a4*b2*c2+a2*b4*c2-b6*c2+a4*c4+2*a2*b2*c4+b4*c4-3*a2*c6-3*b2*c6+2*c8);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X132(orbit) {
   /* begin vars */
   let c2=c*c;
   let c4=c2*c2;
   let b2=b*b;
   let b4=b2*b2;
   let a2=a*a;
   let a4=a2*a2;
   let c6=c2*c4;
   let b6=b2*b4;
   let a6=a2*a4;
   /* end vars */
   let v1 = (a2+b2-c2)*(a2-b2+c2)*(a2*b2-b4+a2*c2-c4)*(2*a6-a4*b2-b6-a4*c2+b4*c2+b2*c4-c6);
   let v2 = (a2+b2-c2)*(-a2+b2+c2)*(-a4+a2*b2+b2*c2-c4)*(-a6-a2*b4+2*b6+a4*c2-b4*c2+a2*c4-c6);
   let v3 = (a2-b2+c2)*(-a2+b2+c2)*(-a4-b4+a2*c2+b2*c2)*(-a6+a4*b2+a2*b4-b6-a2*c4-b2*c4+2*c6);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X133(orbit) {
   /* begin vars */
   let a2=a*a;
   let a4=a2*a2;
   let a6=a2*a4;
   let c2=c*c;
   let c4=c2*c2;
   let c6=c2*c4;
   let b2=b*b;
   let b4=b2*b2;
   let b6=b2*b4;
   let a8=a2*a6;
   let c8=c2*c6;
   let b8=b2*b6;
   /* end vars */
   let v1 = (a2+b2-c2)*(a2-b2+c2)*(2*a4-a2*b2-b4-a2*c2+2*b2*c2-c4)*(a6*b2-3*a4*b4+3*a2*b6-b8+a6*c2+4*a4*b2*c2-3*a2*b4*c2-2*b6*c2-3*a4*c4-3*a2*b2*c4+6*b4*c4+3*a2*c6-2*b2*c6-c8);
   let v2 = (a2+b2-c2)*(-a2+b2+c2)*(-a4-a2*b2+2*b4+2*a2*c2-b2*c2-c4)*(-a8+3*a6*b2-3*a4*b4+a2*b6-2*a6*c2-3*a4*b2*c2+4*a2*b4*c2+b6*c2+6*a4*c4-3*a2*b2*c4-3*b4*c4-2*a2*c6+3*b2*c6-c8);
   let v3 = (a2-b2+c2)*(-a2+b2+c2)*(-a4+2*a2*b2-b4-a2*c2-b2*c2+2*c4)*(-a8-2*a6*b2+6*a4*b4-2*a2*b6-b8+3*a6*c2-3*a4*b2*c2-3*a2*b4*c2+3*b6*c2-3*a4*c4+4*a2*b2*c4-3*b4*c4+a2*c6+b2*c6);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X134(orbit) {
   /* begin vars */
   let c2=c*c;
   let c4=c2*c2;
   let c6=c2*c4;
   let b2=b*b;
   let b4=b2*b2;
   let b6=b2*b4;
   let a2=a*a;
   let a4=a2*a2;
   let a6=a2*a4;
   let c8=c2*c6;
   let b8=b2*b6;
   let a8=a2*a6;
   /* end vars */
   let v1 = a2*(b-c)*(b-c)*(b+c)*(b+c)*(a2*b2-b4+a2*c2+2*b2*c2-c4)*(a4-2*a2*b2+b4-2*a2*c2+c4)*(a4-2*a2*b2+b4-2*a2*c2+c4)*(Math.pow(a,12)-4*Math.pow(a,10)*b2+6*a8*b4-4*a6*b6+a4*b8-b2*Math.pow(c,10)-4*Math.pow(a,10)*c2-Math.pow(b,10)*c2+5*a8*b2*c2+2*a6*b4*c2-4*a4*b6*c2+2*a2*b8*c2+6*a8*c4+2*a6*b2*c4+2*a4*b4*c4-2*a2*b6*c4+4*b8*c4-4*a6*c6-4*a4*b2*c6-2*a2*b4*c6-6*b6*c6+a4*c8+2*a2*b2*c8+4*b4*c8);
   let v2 = b2*(-a+c)*(-a+c)*(a+c)*(a+c)*(-a4+a2*b2+2*a2*c2+b2*c2-c4)*(a4-2*a2*b2+b4-2*b2*c2+c4)*(a4-2*a2*b2+b4-2*b2*c2+c4)*(-4*a2*Math.pow(b,10)+Math.pow(b,12)+a8*b4-4*a6*b6+6*a4*b8-a2*Math.pow(c,10)-Math.pow(a,10)*c2-4*Math.pow(b,10)*c2+2*a8*b2*c2-4*a6*b4*c2+2*a4*b6*c2+5*a2*b8*c2+4*a8*c4-2*a6*b2*c4+2*a4*b4*c4+2*a2*b6*c4+6*b8*c4-6*a6*c6-2*a4*b2*c6-4*a2*b4*c6-4*b6*c6+4*a4*c8+2*a2*b2*c8+b4*c8);
   let v3 = (a-b)*(a-b)*(a+b)*(a+b)*c2*(-a4+2*a2*b2-b4+a2*c2+b2*c2)*(a4+b4-2*a2*c2-2*b2*c2+c4)*(a4+b4-2*a2*c2-2*b2*c2+c4)*(-(a2*Math.pow(b,10))-Math.pow(a,10)*b2+4*a8*b4-6*a6*b6+4*a4*b8-4*a2*Math.pow(c,10)-4*b2*Math.pow(c,10)+Math.pow(c,12)+2*a8*b2*c2-2*a6*b4*c2-2*a4*b6*c2+2*a2*b8*c2+a8*c4-4*a6*b2*c4+2*a4*b4*c4-4*a2*b6*c4+b8*c4-4*a6*c6+2*a4*b2*c6+2*a2*b4*c6-4*b6*c6+6*a4*c8+5*a2*b2*c8+6*b4*c8);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X135(orbit) {
   /* begin vars */
   let c2=c*c;
   let c4=c2*c2;
   let b2=b*b;
   let b4=b2*b2;
   let a2=a*a;
   let a4=a2*a2;
   let c6=c2*c4;
   let b6=b2*b4;
   let a6=a2*a4;
   /* end vars */
   let v1 = (b-c)*(b-c)*(b+c)*(b+c)*(-a2+b2-c2)*(a2+b2-c2)*(a4-2*a2*b2+b4-2*a2*c2+c4)*(-a6+3*a4*b2-3*a2*b4+b6+3*a4*c2+2*a2*b2*c2-b4*c2-3*a2*c4-b2*c4+c6);
   let v2 = (-a+c)*(-a+c)*(a+c)*(a+c)*(-a2-b2+c2)*(-a2+b2+c2)*(a4-2*a2*b2+b4-2*b2*c2+c4)*(a6-3*a4*b2+3*a2*b4-b6-a4*c2+2*a2*b2*c2+3*b4*c2-a2*c4-3*b2*c4+c6);
   let v3 = (a-b)*(a-b)*(a+b)*(a+b)*(a2-b2-c2)*(a2-b2+c2)*(a4+b4-2*a2*c2-2*b2*c2+c4)*(a6-a4*b2-a2*b4+b6-3*a4*c2+2*a2*b2*c2-3*b4*c2+3*a2*c4+3*b2*c4-c6);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X136(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let c4=c2*c2;
   let b4=b2*b2;
   let a4=a2*a2;
   /* end vars */
   let v1 = (b-c)*(b-c)*(b+c)*(b+c)*(-a2+b2-c2)*(a2+b2-c2)*(a4-2*a2*b2+b4-2*a2*c2+c4);
   let v2 = (-a+c)*(-a+c)*(a+c)*(a+c)*(-a2-b2+c2)*(-a2+b2+c2)*(a4-2*a2*b2+b4-2*b2*c2+c4);
   let v3 = (a-b)*(a-b)*(a+b)*(a+b)*(a2-b2-c2)*(a2-b2+c2)*(a4+b4-2*a2*c2-2*b2*c2+c4);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X137(orbit) {
   /* begin vars */
   let a2=a*a;
   let c2=c*c;
   let b2=b*b;
   let a4=a2*a2;
   let c4=c2*c2;
   let b4=b2*b2;
   /* end vars */
   let v1 = (b-c)*(b-c)*(b+c)*(b+c)*(-(a2*b2)+b4-a2*c2-2*b2*c2+c4)*(a4-2*a2*b2+b4-2*a2*c2-b2*c2+c4);
   let v2 = (-a+c)*(-a+c)*(a+c)*(a+c)*(a4-2*a2*b2+b4-a2*c2-2*b2*c2+c4)*(a4-a2*b2-2*a2*c2-b2*c2+c4);
   let v3 = (a-b)*(a-b)*(a+b)*(a+b)*(a4-2*a2*b2+b4-a2*c2-b2*c2)*(a4-a2*b2+b4-2*a2*c2-2*b2*c2+c4);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X138(orbit) {
   /* begin vars */
   let c2=c*c;
   let c4=c2*c2;
   let c6=c2*c4;
   let b2=b*b;
   let b4=b2*b2;
   let b6=b2*b4;
   let a2=a*a;
   let a4=a2*a2;
   let a6=a2*a4;
   let c8=c2*c6;
   let b8=b2*b6;
   let a8=a2*a6;
   /* end vars */
   let v1 = (a2+b2-c2)*(a2-b2+c2)*(a2*b2-b4+a2*c2+2*b2*c2-c4)*(2*a8-4*a6*b2+a4*b4+2*a2*b6-b8-4*a6*c2+4*a4*b2*c2-2*a2*b4*c2+2*b6*c2+a4*c4-2*a2*b2*c4-2*b4*c4+2*a2*c6+2*b2*c6-c8)*(a8-2*a6*b2+2*a4*b4-2*a2*b6+b8-2*a6*c2-a4*b2*c2+2*a2*b4*c2+b6*c2+2*a4*c4+2*a2*b2*c4-4*b4*c4-2*a2*c6+b2*c6+c8);
   let v2 = (a2+b2-c2)*(-a2+b2+c2)*(-a4+a2*b2+2*a2*c2+b2*c2-c4)*(-a8+2*a6*b2+a4*b4-4*a2*b6+2*b8+2*a6*c2-2*a4*b2*c2+4*a2*b4*c2-4*b6*c2-2*a4*c4-2*a2*b2*c4+b4*c4+2*a2*c6+2*b2*c6-c8)*(a8-2*a6*b2+2*a4*b4-2*a2*b6+b8+a6*c2+2*a4*b2*c2-a2*b4*c2-2*b6*c2-4*a4*c4+2*a2*b2*c4+2*b4*c4+a2*c6-2*b2*c6+c8);
   let v3 = (a2-b2+c2)*(-a2+b2+c2)*(-a4+2*a2*b2-b4+a2*c2+b2*c2)*(a8+a6*b2-4*a4*b4+a2*b6+b8-2*a6*c2+2*a4*b2*c2+2*a2*b4*c2-2*b6*c2+2*a4*c4-a2*b2*c4+2*b4*c4-2*a2*c6-2*b2*c6+c8)*(-a8+2*a6*b2-2*a4*b4+2*a2*b6-b8+2*a6*c2-2*a4*b2*c2-2*a2*b4*c2+2*b6*c2+a4*c4+4*a2*b2*c4+b4*c4-4*a2*c6-4*b2*c6+2*c8);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X139(orbit) {
   /* begin vars */
   let c2=c*c;
   let c4=c2*c2;
   let c6=c2*c4;
   let b2=b*b;
   let b4=b2*b2;
   let b6=b2*b4;
   let a2=a*a;
   let a4=a2*a2;
   let a6=a2*a4;
   let c8=c2*c6;
   let b8=b2*b6;
   let a8=a2*a6;
   /* end vars */
   let v1 = (b-c)*(b-c)*(b+c)*(b+c)*(-a2+b2-c2)*(a2+b2-c2)*(a4-2*a2*b2+b4-2*a2*c2+c4)*(-(a2*b2)+b4-a2*c2-2*b2*c2+c4)*(Math.pow(a,12)-4*a2*Math.pow(b,10)+Math.pow(b,12)-4*Math.pow(a,10)*b2+7*a8*b4-8*a6*b6+7*a4*b8-4*a2*Math.pow(c,10)-3*b2*Math.pow(c,10)+Math.pow(c,12)-4*Math.pow(a,10)*c2-3*Math.pow(b,10)*c2+11*a8*b2*c2-10*a6*b4*c2+6*a2*b8*c2+7*a8*c4-10*a6*b2*c4+2*a4*b4*c4-2*a2*b6*c4+3*b8*c4-8*a6*c6-2*a2*b4*c6-2*b6*c6+7*a4*c8+6*a2*b2*c8+3*b4*c8);
   let v2 = (-a+c)*(-a+c)*(a+c)*(a+c)*(-a2-b2+c2)*(-a2+b2+c2)*(a4-2*a2*b2+b4-2*b2*c2+c4)*(a4-a2*b2-2*a2*c2-b2*c2+c4)*(Math.pow(a,12)-4*a2*Math.pow(b,10)+Math.pow(b,12)-4*Math.pow(a,10)*b2+7*a8*b4-8*a6*b6+7*a4*b8-3*a2*Math.pow(c,10)-4*b2*Math.pow(c,10)+Math.pow(c,12)-3*Math.pow(a,10)*c2-4*Math.pow(b,10)*c2+6*a8*b2*c2-10*a4*b6*c2+11*a2*b8*c2+3*a8*c4-2*a6*b2*c4+2*a4*b4*c4-10*a2*b6*c4+7*b8*c4-2*a6*c6-2*a4*b2*c6-8*b6*c6+3*a4*c8+6*a2*b2*c8+7*b4*c8);
   let v3 = (a-b)*(a-b)*(a+b)*(a+b)*(a2-b2-c2)*(a2-b2+c2)*(a4-2*a2*b2+b4-a2*c2-b2*c2)*(a4+b4-2*a2*c2-2*b2*c2+c4)*(Math.pow(a,12)-3*a2*Math.pow(b,10)+Math.pow(b,12)-3*Math.pow(a,10)*b2+3*a8*b4-2*a6*b6+3*a4*b8-4*a2*Math.pow(c,10)-4*b2*Math.pow(c,10)+Math.pow(c,12)-4*Math.pow(a,10)*c2-4*Math.pow(b,10)*c2+6*a8*b2*c2-2*a6*b4*c2-2*a4*b6*c2+6*a2*b8*c2+7*a8*c4+2*a4*b4*c4+7*b8*c4-8*a6*c6-10*a4*b2*c6-10*a2*b4*c6-8*b6*c6+7*a4*c8+11*a2*b2*c8+7*b4*c8);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X140(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let c4=c2*c2;
   let b4=b2*b2;
   let a4=a2*a2;
   /* end vars */
   let v1 = 2*a4-3*a2*b2+b4-3*a2*c2-2*b2*c2+c4;
   let v2 = a4-3*a2*b2+2*b4-2*a2*c2-3*b2*c2+c4;
   let v3 = a4-2*a2*b2+b4-3*a2*c2-3*b2*c2+2*c4;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X141(orbit) {
   /* begin vars */
   let a2=a*a;
   let c2=c*c;
   let b2=b*b;
   /* end vars */
   let v1 = b2+c2;
   let v2 = a2+c2;
   let v3 = a2+b2;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X142(orbit) {
   /* begin vars */
   let a2=a*a;
   let c2=c*c;
   let b2=b*b;
   /* end vars */
   let v1 = -(a*b)+b2-a*c-2*b*c+c2;
   let v2 = a2-a*b-2*a*c-b*c+c2;
   let v3 = a2-2*a*b+b2-a*c-b*c;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X143(orbit) {
   /* begin vars */
   let a2=a*a;
   let c2=c*c;
   let b2=b*b;
   let a4=a2*a2;
   let c4=c2*c2;
   let b4=b2*b2;
   /* end vars */
   let v1 = a2*(a2*b2-b4+a2*c2+2*b2*c2-c4)*(a4-2*a2*b2+b4-2*a2*c2-b2*c2+c4);
   let v2 = b2*(-a4+a2*b2+2*a2*c2+b2*c2-c4)*(a4-2*a2*b2+b4-a2*c2-2*b2*c2+c4);
   let v3 = c2*(-a4+2*a2*b2-b4+a2*c2+b2*c2)*(a4-a2*b2+b4-2*a2*c2-2*b2*c2+c4);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X144(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   /* end vars */
   let v1 = -3*a2+2*a*b+b2+2*a*c-2*b*c+c2;
   let v2 = a2+2*a*b-3*b2-2*a*c+2*b*c+c2;
   let v3 = a2-2*a*b+b2+2*a*c+2*b*c-3*c2;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X145(orbit) {
   /* begin vars */

   /* end vars */
   let v1 = -3*a+b+c;
   let v2 = a-3*b+c;
   let v3 = a+b-3*c;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X146(orbit) {
   /* begin vars */
   let c2=c*c;
   let c4=c2*c2;
   let c6=c2*c4;
   let b2=b*b;
   let b4=b2*b2;
   let b6=b2*b4;
   let a2=a*a;
   let a4=a2*a2;
   let a6=a2*a4;
   let c8=c2*c6;
   let b8=b2*b6;
   let a8=a2*a6;
   /* end vars */
   let v1 = -Math.pow(a,10)+Math.pow(b,10)-a8*b2+8*a6*b4-8*a4*b6+a2*b8+Math.pow(c,10)-a8*c2-9*a6*b2*c2+6*a4*b4*c2+7*a2*b6*c2-3*b8*c2+8*a6*c4+6*a4*b2*c4-16*a2*b4*c4+2*b6*c4-8*a4*c6+7*a2*b2*c6+2*b4*c6+a2*c8-3*b2*c8;
   let v2 = Math.pow(a,10)-Math.pow(b,10)+a8*b2-8*a6*b4+8*a4*b6-a2*b8+Math.pow(c,10)-3*a8*c2+7*a6*b2*c2+6*a4*b4*c2-9*a2*b6*c2-b8*c2+2*a6*c4-16*a4*b2*c4+6*a2*b4*c4+8*b6*c4+2*a4*c6+7*a2*b2*c6-8*b4*c6-3*a2*c8+b2*c8;
   let v3 = Math.pow(a,10)+Math.pow(b,10)-3*a8*b2+2*a6*b4+2*a4*b6-3*a2*b8-Math.pow(c,10)+a8*c2+7*a6*b2*c2-16*a4*b4*c2+7*a2*b6*c2+b8*c2-8*a6*c4+6*a4*b2*c4+6*a2*b4*c4-8*b6*c4+8*a4*c6-9*a2*b2*c6+8*b4*c6-a2*c8-b2*c8;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X147(orbit) {
   /* begin vars */
   let c2=c*c;
   let c4=c2*c2;
   let c6=c2*c4;
   let b2=b*b;
   let b4=b2*b2;
   let b6=b2*b4;
   let a2=a*a;
   let a4=a2*a2;
   let a6=a2*a4;
   let c8=c2*c6;
   let b8=b2*b6;
   let a8=a2*a6;
   /* end vars */
   let v1 = a8+a6*b2-2*a4*b4+a2*b6-b8+a6*c2-3*a4*b2*c2+a2*b4*c2+b6*c2-2*a4*c4+a2*b2*c4+a2*c6+b2*c6-c8;
   let v2 = -a8+a6*b2-2*a4*b4+a2*b6+b8+a6*c2+a4*b2*c2-3*a2*b4*c2+b6*c2+a2*b2*c4-2*b4*c4+a2*c6+b2*c6-c8;
   let v3 = -a8+a6*b2+a2*b6-b8+a6*c2+a4*b2*c2+a2*b4*c2+b6*c2-2*a4*c4-3*a2*b2*c4-2*b4*c4+a2*c6+b2*c6+c8;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X148(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let c4=c2*c2;
   let b4=b2*b2;
   let a4=a2*a2;
   /* end vars */
   let v1 = -a4+a2*b2+b4+a2*c2-3*b2*c2+c4;
   let v2 = a4+a2*b2-b4-3*a2*c2+b2*c2+c4;
   let v3 = a4-3*a2*b2+b4+a2*c2+b2*c2-c4;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X149(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let c3=c2*c;
   let b3=b2*b;
   let a3=a2*a;
   /* end vars */
   let v1 = -a3+a2*b-a*b2+b3+a2*c+a*b*c-b2*c-a*c2-b*c2+c3;
   let v2 = a3-a2*b+a*b2-b3-a2*c+a*b*c+b2*c-a*c2-b*c2+c3;
   let v3 = a3-a2*b-a*b2+b3-a2*c+a*b*c-b2*c+a*c2+b*c2-c3;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X150(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let c4=c2*c2;
   let c3=c2*c;
   let b4=b2*b2;
   let b3=b2*b;
   let a3=a2*a;
   let a4=a2*a2;
   /* end vars */
   let v1 = -a4+a3*b-a*b3+b4+a3*c-a2*b*c+a*b2*c-b3*c+a*b*c2-a*c3-b*c3+c4;
   let v2 = a4-a3*b+a*b3-b4-a3*c+a2*b*c-a*b2*c+b3*c+a*b*c2-a*c3-b*c3+c4;
   let v3 = a4-a3*b-a*b3+b4-a3*c+a2*b*c+a*b2*c-b3*c-a*b*c2+a*c3+b*c3-c4;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X151(orbit) {
   /* begin vars */
   let c2=c*c;
   let c4=c2*c2;
   let c6=c2*c4;
   let c3=c2*c;
   let b2=b*b;
   let b3=b2*b;
   let a2=a*a;
   let a3=a2*a;
   let b4=b2*b2;
   let b6=b2*b4;
   let a4=a2*a2;
   let a6=a2*a4;
   let c8=c2*c6;
   let c7=c*c6;
   let c5=c2*c3;
   let b5=b2*b3;
   let a5=a2*a3;
   let b8=b2*b6;
   let b7=b*b6;
   let a7=a*a6;
   let a8=a2*a6;
   /* end vars */
   let v1 = Math.pow(a,10)-Math.pow(a,9)*b+a*Math.pow(b,9)-Math.pow(b,10)+a8*b2+2*a7*b3-8*a6*b4+8*a4*b6-2*a3*b7-a2*b8-Math.pow(a,9)*c+3*a8*b*c+Math.pow(b,9)*c-6*a7*b2*c+2*a6*b3*c+12*a5*b4*c-12*a4*b5*c-2*a3*b6*c+6*a2*b7*c-3*a*b8*c+a*Math.pow(c,9)+b*Math.pow(c,9)-Math.pow(c,10)+a8*c2-6*a7*b*c2+16*a6*b2*c2-12*a5*b3*c2-8*a4*b4*c2+18*a3*b5*c2-12*a2*b6*c2+3*b8*c2+2*a7*c3+2*a6*b*c3-12*a5*b2*c3+24*a4*b3*c3-14*a3*b4*c3-6*a2*b5*c3+8*a*b6*c3-4*b7*c3-8*a6*c4+12*a5*b*c4-8*a4*b2*c4-14*a3*b3*c4+26*a2*b4*c4-6*a*b5*c4-2*b6*c4-12*a4*b*c5+18*a3*b2*c5-6*a2*b3*c5-6*a*b4*c5+6*b5*c5+8*a4*c6-2*a3*b*c6-12*a2*b2*c6+8*a*b3*c6-2*b4*c6-2*a3*c7+6*a2*b*c7-4*b3*c7-a2*c8-3*a*b*c8+3*b2*c8;
   let v2 = -Math.pow(a,10)+Math.pow(a,9)*b-a*Math.pow(b,9)+Math.pow(b,10)-a8*b2-2*a7*b3+8*a6*b4-8*a4*b6+2*a3*b7+a2*b8+Math.pow(a,9)*c-3*a8*b*c-Math.pow(b,9)*c+6*a7*b2*c-2*a6*b3*c-12*a5*b4*c+12*a4*b5*c+2*a3*b6*c-6*a2*b7*c+3*a*b8*c+a*Math.pow(c,9)+b*Math.pow(c,9)-Math.pow(c,10)+3*a8*c2-12*a6*b2*c2+18*a5*b3*c2-8*a4*b4*c2-12*a3*b5*c2+16*a2*b6*c2-6*a*b7*c2+b8*c2-4*a7*c3+8*a6*b*c3-6*a5*b2*c3-14*a4*b3*c3+24*a3*b4*c3-12*a2*b5*c3+2*a*b6*c3+2*b7*c3-2*a6*c4-6*a5*b*c4+26*a4*b2*c4-14*a3*b3*c4-8*a2*b4*c4+12*a*b5*c4-8*b6*c4+6*a5*c5-6*a4*b*c5-6*a3*b2*c5+18*a2*b3*c5-12*a*b4*c5-2*a4*c6+8*a3*b*c6-12*a2*b2*c6-2*a*b3*c6+8*b4*c6-4*a3*c7+6*a*b2*c7-2*b3*c7+3*a2*c8-3*a*b*c8-b2*c8;
   let v3 = -Math.pow(a,10)+Math.pow(a,9)*b+a*Math.pow(b,9)-Math.pow(b,10)+3*a8*b2-4*a7*b3-2*a6*b4+6*a5*b5-2*a4*b6-4*a3*b7+3*a2*b8+Math.pow(a,9)*c-3*a8*b*c+Math.pow(b,9)*c+8*a6*b3*c-6*a5*b4*c-6*a4*b5*c+8*a3*b6*c-3*a*b8*c-a*Math.pow(c,9)-b*Math.pow(c,9)+Math.pow(c,10)-a8*c2+6*a7*b*c2-12*a6*b2*c2-6*a5*b3*c2+26*a4*b4*c2-6*a3*b5*c2-12*a2*b6*c2+6*a*b7*c2-b8*c2-2*a7*c3-2*a6*b*c3+18*a5*b2*c3-14*a4*b3*c3-14*a3*b4*c3+18*a2*b5*c3-2*a*b6*c3-2*b7*c3+8*a6*c4-12*a5*b*c4-8*a4*b2*c4+24*a3*b3*c4-8*a2*b4*c4-12*a*b5*c4+8*b6*c4+12*a4*b*c5-12*a3*b2*c5-12*a2*b3*c5+12*a*b4*c5-8*a4*c6+2*a3*b*c6+16*a2*b2*c6+2*a*b3*c6-8*b4*c6+2*a3*c7-6*a2*b*c7-6*a*b2*c7+2*b3*c7+a2*c8+3*a*b*c8+b2*c8;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X152(orbit) {
   /* begin vars */
   let c2=c*c;
   let c4=c2*c2;
   let c6=c2*c4;
   let c3=c2*c;
   let b2=b*b;
   let a2=a*a;
   let b4=b2*b2;
   let b6=b2*b4;
   let b3=b2*b;
   let a3=a2*a;
   let a4=a2*a2;
   let a6=a2*a4;
   let c8=c2*c6;
   let c7=c*c6;
   let c5=c2*c3;
   let b8=b2*b6;
   let b7=b*b6;
   let b5=b2*b3;
   let a5=a2*a3;
   let a7=a*a6;
   let a8=a2*a6;
   /* end vars */
   let v1 = -a8+a7*b-2*a6*b2+5*a5*b3-5*a3*b5+2*a2*b6-a*b7+b8+a7*c-a6*b*c-a5*b2*c+a4*b3*c-a3*b4*c+a2*b5*c+a*b6*c-b7*c-2*a6*c2-a5*b*c2-6*a4*b2*c2+6*a3*b3*c2+2*a2*b4*c2+3*a*b5*c2-2*b6*c2+5*a5*c3+a4*b*c3+6*a3*b2*c3-10*a2*b3*c3-3*a*b4*c3+b5*c3-a3*b*c4+2*a2*b2*c4-3*a*b3*c4+2*b4*c4-5*a3*c5+a2*b*c5+3*a*b2*c5+b3*c5+2*a2*c6+a*b*c6-2*b2*c6-a*c7-b*c7+c8;
   let v2 = a8-a7*b+2*a6*b2-5*a5*b3+5*a3*b5-2*a2*b6+a*b7-b8-a7*c+a6*b*c+a5*b2*c-a4*b3*c+a3*b4*c-a2*b5*c-a*b6*c+b7*c-2*a6*c2+3*a5*b*c2+2*a4*b2*c2+6*a3*b3*c2-6*a2*b4*c2-a*b5*c2-2*b6*c2+a5*c3-3*a4*b*c3-10*a3*b2*c3+6*a2*b3*c3+a*b4*c3+5*b5*c3+2*a4*c4-3*a3*b*c4+2*a2*b2*c4-a*b3*c4+a3*c5+3*a2*b*c5+a*b2*c5-5*b3*c5-2*a2*c6+a*b*c6+2*b2*c6-a*c7-b*c7+c8;
   let v3 = a8-a7*b-2*a6*b2+a5*b3+2*a4*b4+a3*b5-2*a2*b6-a*b7+b8-a7*c+a6*b*c+3*a5*b2*c-3*a4*b3*c-3*a3*b4*c+3*a2*b5*c+a*b6*c-b7*c+2*a6*c2+a5*b*c2+2*a4*b2*c2-10*a3*b3*c2+2*a2*b4*c2+a*b5*c2+2*b6*c2-5*a5*c3-a4*b*c3+6*a3*b2*c3+6*a2*b3*c3-a*b4*c3-5*b5*c3+a3*b*c4-6*a2*b2*c4+a*b3*c4+5*a3*c5-a2*b*c5-a*b2*c5+5*b3*c5-2*a2*c6-a*b*c6-2*b2*c6+a*c7+b*c7-c8;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X153(orbit) {
   /* begin vars */
   let c2=c*c;
   let c4=c2*c2;
   let c6=c2*c4;
   let c3=c2*c;
   let b2=b*b;
   let b4=b2*b2;
   let b6=b2*b4;
   let b3=b2*b;
   let a2=a*a;
   let a3=a2*a;
   let a4=a2*a2;
   let a6=a2*a4;
   let c7=c*c6;
   let c5=c2*c3;
   let b7=b*b6;
   let b5=b2*b3;
   let a5=a2*a3;
   let a7=a*a6;
   /* end vars */
   let v1 = -a7+a6*b+a5*b2-a4*b3+a3*b4-a2*b5-a*b6+b7+a6*c-7*a5*b*c+5*a4*b2*c+2*a3*b3*c-5*a2*b4*c+5*a*b5*c-b6*c+a5*c2+5*a4*b*c2-10*a3*b2*c2+6*a2*b3*c2+a*b4*c2-3*b5*c2-a4*c3+2*a3*b*c3+6*a2*b2*c3-10*a*b3*c3+3*b4*c3+a3*c4-5*a2*b*c4+a*b2*c4+3*b3*c4-a2*c5+5*a*b*c5-3*b2*c5-a*c6-b*c6+c7;
   let v2 = a7-a6*b-a5*b2+a4*b3-a3*b4+a2*b5+a*b6-b7-a6*c+5*a5*b*c-5*a4*b2*c+2*a3*b3*c+5*a2*b4*c-7*a*b5*c+b6*c-3*a5*c2+a4*b*c2+6*a3*b2*c2-10*a2*b3*c2+5*a*b4*c2+b5*c2+3*a4*c3-10*a3*b*c3+6*a2*b2*c3+2*a*b3*c3-b4*c3+3*a3*c4+a2*b*c4-5*a*b2*c4+b3*c4-3*a2*c5+5*a*b*c5-b2*c5-a*c6-b*c6+c7;
   let v3 = a7-a6*b-3*a5*b2+3*a4*b3+3*a3*b4-3*a2*b5-a*b6+b7-a6*c+5*a5*b*c+a4*b2*c-10*a3*b3*c+a2*b4*c+5*a*b5*c-b6*c-a5*c2-5*a4*b*c2+6*a3*b2*c2+6*a2*b3*c2-5*a*b4*c2-b5*c2+a4*c3+2*a3*b*c3-10*a2*b2*c3+2*a*b3*c3+b4*c3-a3*c4+5*a2*b*c4+5*a*b2*c4-b3*c4+a2*c5-7*a*b*c5+b2*c5+a*c6+b*c6-c7;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X154(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let c4=c2*c2;
   let b4=b2*b2;
   let a4=a2*a2;
   /* end vars */
   let v1 = a2*(3*a4-2*a2*b2-b4-2*a2*c2+2*b2*c2-c4);
   let v2 = b2*(-a4-2*a2*b2+3*b4+2*a2*c2-2*b2*c2-c4);
   let v3 = c2*(-a4+2*a2*b2-b4-2*a2*c2-2*b2*c2+3*c4);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X155(orbit) {
   /* begin vars */
   let c2=c*c;
   let c4=c2*c2;
   let b2=b*b;
   let b4=b2*b2;
   let a2=a*a;
   let a4=a2*a2;
   let c6=c2*c4;
   let b6=b2*b4;
   let a6=a2*a4;
   /* end vars */
   let v1 = a2*(a2-b2-c2)*(a6-3*a4*b2+3*a2*b4-b6-3*a4*c2-2*a2*b2*c2+b4*c2+3*a2*c4+b2*c4-c6);
   let v2 = b2*(-a2+b2-c2)*(-a6+3*a4*b2-3*a2*b4+b6+a4*c2-2*a2*b2*c2-3*b4*c2+a2*c4+3*b2*c4-c6);
   let v3 = c2*(-a2-b2+c2)*(-a6+a4*b2+a2*b4-b6+3*a4*c2-2*a2*b2*c2+3*b4*c2-3*a2*c4-3*b2*c4+c6);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X156(orbit) {
   /* begin vars */
   let c2=c*c;
   let c4=c2*c2;
   let c6=c2*c4;
   let b2=b*b;
   let b4=b2*b2;
   let b6=b2*b4;
   let a2=a*a;
   let a4=a2*a2;
   let a6=a2*a4;
   let c8=c2*c6;
   let b8=b2*b6;
   let a8=a2*a6;
   /* end vars */
   let v1 = a2*(a8-3*a6*b2+3*a4*b4-a2*b6-3*a6*c2+2*a4*b2*c2+b6*c2+3*a4*c4-2*b4*c4-a2*c6+b2*c6);
   let v2 = b2*(-(a6*b2)+3*a4*b4-3*a2*b6+b8+a6*c2+2*a2*b4*c2-3*b6*c2-2*a4*c4+3*b4*c4+a2*c6-b2*c6);
   let v3 = c2*(a6*b2-2*a4*b4+a2*b6-a6*c2-b6*c2+3*a4*c4+2*a2*b2*c4+3*b4*c4-3*a2*c6-3*b2*c6+c8);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X157(orbit) {
   /* begin vars */
   let c2=c*c;
   let c4=c2*c2;
   let b2=b*b;
   let b4=b2*b2;
   let a2=a*a;
   let a4=a2*a2;
   let c6=c2*c4;
   let b6=b2*b4;
   let a6=a2*a4;
   /* end vars */
   let v1 = a2*(a6-a4*b2+a2*b4-b6-a4*c2+b4*c2+a2*c4+b2*c4-c6);
   let v2 = b2*(-a6+a4*b2-a2*b4+b6+a4*c2-b4*c2+a2*c4+b2*c4-c6);
   let v3 = c2*(-a6+a4*b2+a2*b4-b6+a4*c2+b4*c2-a2*c4-b2*c4+c6);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X158(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   /* end vars */
   let v1 = b*c*(-a2+b2-c2)*(-a2+b2-c2)*(a2+b2-c2)*(a2+b2-c2);
   let v2 = a*c*(-a2-b2+c2)*(-a2-b2+c2)*(-a2+b2+c2)*(-a2+b2+c2);
   let v3 = a*b*(a2-b2-c2)*(a2-b2-c2)*(a2-b2+c2)*(a2-b2+c2);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X159(orbit) {
   /* begin vars */
   let c2=c*c;
   let c4=c2*c2;
   let b2=b*b;
   let b4=b2*b2;
   let a2=a*a;
   let a4=a2*a2;
   let c6=c2*c4;
   let b6=b2*b4;
   let a6=a2*a4;
   /* end vars */
   let v1 = a2*(a6+a4*b2-a2*b4-b6+a4*c2-2*a2*b2*c2+b4*c2-a2*c4+b2*c4-c6);
   let v2 = b2*(-a6-a4*b2+a2*b4+b6+a4*c2-2*a2*b2*c2+b4*c2+a2*c4-b2*c4-c6);
   let v3 = c2*(-a6+a4*b2+a2*b4-b6-a4*c2-2*a2*b2*c2-b4*c2+a2*c4+b2*c4+c6);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X160(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let c4=c2*c2;
   let b4=b2*b2;
   let a4=a2*a2;
   /* end vars */
   let v1 = a4*(a2*b2-b4+a2*c2-b2*c2-c4);
   let v2 = b4*(-a4+a2*b2-a2*c2+b2*c2-c4);
   let v3 = (-a4-a2*b2-b4+a2*c2+b2*c2)*c4;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X161(orbit) {
   /* begin vars */
   let c2=c*c;
   let c4=c2*c2;
   let c6=c2*c4;
   let b2=b*b;
   let b4=b2*b2;
   let b6=b2*b4;
   let a2=a*a;
   let a4=a2*a2;
   let a6=a2*a4;
   let c8=c2*c6;
   let b8=b2*b6;
   let a8=a2*a6;
   /* end vars */
   let v1 = a2*(Math.pow(a,10)-Math.pow(b,10)-a8*b2-2*a6*b4+2*a4*b6+a2*b8-Math.pow(c,10)-a8*c2-2*a6*b2*c2+2*a4*b4*c2-2*a2*b6*c2+3*b8*c2-2*a6*c4+2*a4*b2*c4+2*a2*b4*c4-2*b6*c4+2*a4*c6-2*a2*b2*c6-2*b4*c6+a2*c8+3*b2*c8);
   let v2 = b2*(-Math.pow(a,10)+Math.pow(b,10)+a8*b2+2*a6*b4-2*a4*b6-a2*b8-Math.pow(c,10)+3*a8*c2-2*a6*b2*c2+2*a4*b4*c2-2*a2*b6*c2-b8*c2-2*a6*c4+2*a4*b2*c4+2*a2*b4*c4-2*b6*c4-2*a4*c6-2*a2*b2*c6+2*b4*c6+3*a2*c8+b2*c8);
   let v3 = c2*(-Math.pow(a,10)-Math.pow(b,10)+3*a8*b2-2*a6*b4-2*a4*b6+3*a2*b8+Math.pow(c,10)+a8*c2-2*a6*b2*c2+2*a4*b4*c2-2*a2*b6*c2+b8*c2+2*a6*c4+2*a4*b2*c4+2*a2*b4*c4+2*b6*c4-2*a4*c6-2*a2*b2*c6-2*b4*c6-a2*c8-b2*c8);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X162(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   /* end vars */
   let v1 = a*(a-b)*(a+b)*(a-c)*(a+c)*(a2+b2-c2)*(a2-b2+c2);
   let v2 = b*(-a+b)*(a+b)*(b-c)*(b+c)*(a2+b2-c2)*(-a2+b2+c2);
   let v3 = c*(-a+c)*(a+c)*(-b+c)*(b+c)*(a2-b2+c2)*(-a2+b2+c2);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X163(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let c3=c2*c;
   let b3=b2*b;
   let a3=a2*a;
   /* end vars */
   let v1 = a3*(a-b)*(a+b)*(a-c)*(a+c);
   let v2 = (-a+b)*(a+b)*b3*(b-c)*(b+c);
   let v3 = (-a+c)*(a+c)*(-b+c)*(b+c)*c3;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X164(orbit) {
   /* begin vars */
   let Sqrt=Math.sqrt;
   /* end vars */
   let v1 = a*(-Sqrt(a*(a+b-c)*(a-b+c))+Sqrt(b*(a+b-c)*(-a+b+c))+Sqrt(c*(a-b+c)*(-a+b+c)));
   let v2 = b*(Sqrt(a*(a+b-c)*(a-b+c))-Sqrt(b*(a+b-c)*(-a+b+c))+Sqrt(c*(a-b+c)*(-a+b+c)));
   let v3 = c*(Sqrt(a*(a+b-c)*(a-b+c))+Sqrt(b*(a+b-c)*(-a+b+c))-Sqrt(c*(a-b+c)*(-a+b+c)));
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X165(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   /* end vars */
   let v1 = a*(3*a2-2*a*b-b2-2*a*c+2*b*c-c2);
   let v2 = b*(-a2-2*a*b+3*b2+2*a*c-2*b*c-c2);
   let v3 = c*(-a2+2*a*b-b2-2*a*c-2*b*c+3*c2);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X166(orbit) {
   /* begin vars */
   let sb=(c+a-b)/2;
   let sa=(b+c-a)/2;
   let s=(a+b+c)/2;
   let Sqrt=Math.sqrt;
   let sc=(a+b-c)/2;
   /* end vars */
   let v1 = a*(1/(sc*(Sqrt(a*s*sa)+Sqrt(b*s*sb)-Sqrt(c*s*sc)))+1/(sb*(Sqrt(a*s*sa)-Sqrt(b*s*sb)+Sqrt(c*s*sc)))-1/(sa*(-Sqrt(a*s*sa)+Sqrt(b*s*sb)+Sqrt(c*s*sc))));
   let v2 = b*(1/(sc*(Sqrt(a*s*sa)+Sqrt(b*s*sb)-Sqrt(c*s*sc)))-1/(sb*(Sqrt(a*s*sa)-Sqrt(b*s*sb)+Sqrt(c*s*sc)))+1/(sa*(-Sqrt(a*s*sa)+Sqrt(b*s*sb)+Sqrt(c*s*sc))));
   let v3 = c*(-(1/(sc*(Sqrt(a*s*sa)+Sqrt(b*s*sb)-Sqrt(c*s*sc))))+1/(sb*(Sqrt(a*s*sa)-Sqrt(b*s*sb)+Sqrt(c*s*sc)))+1/(sa*(-Sqrt(a*s*sa)+Sqrt(b*s*sb)+Sqrt(c*s*sc))));
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X167(orbit) {
   /* begin vars */
   let b2=b*b;
   let a2=a*a;
   let c2=c*c;
   let SB=(c2+a2-b2)/2;
   let SA=(b2+c2-a2)/2;
   let SC=(a2+b2-c2)/2;
   let sb=(c+a-b)/2;
   let sc=(a+b-c)/2;
   let sa=(b+c-a)/2;
   let Sqrt=Math.sqrt;
   /* end vars */
   let v1 = a*(-(c*(a*Sqrt((sa*sc)/(a*c))-b*Sqrt((sb*sc)/(b*c)))*(2*(-a+c)*sb*Sqrt((sa*sc)/(a*c))+2*sa*sc*(-Sqrt((sa*sb)/(a*b))+Sqrt((sb*sc)/(b*c)))))+b*(-(a*Sqrt((sa*sb)/(a*b)))+c*Sqrt((sb*sc)/(b*c)))*(2*(a-b)*Sqrt((sa*sb)/(a*b))*sc+(-Sqrt((sa*sc)/(a*c))+Sqrt((sb*sc)/(b*c)))*(-(a*b)+SC)));
   let v2 = b*(c*(a*Sqrt((sa*sc)/(a*c))-b*Sqrt((sb*sc)/(b*c)))*(2*(b-c)*sa*Sqrt((sb*sc)/(b*c))+(-(b*c)+SA)*(-Sqrt((sa*sb)/(a*b))+Sqrt((sa*sc)/(a*c))))-a*(b*Sqrt((sa*sb)/(a*b))-c*Sqrt((sa*sc)/(a*c)))*(2*(a-b)*Sqrt((sa*sb)/(a*b))*sc+2*sa*sb*(Sqrt((sa*sc)/(a*c))-Sqrt((sb*sc)/(b*c)))));
   let v3 = c*(-(b*(-(a*Sqrt((sa*sb)/(a*b)))+c*Sqrt((sb*sc)/(b*c)))*(2*(b-c)*sa*Sqrt((sb*sc)/(b*c))+2*sb*sc*(Sqrt((sa*sb)/(a*b))-Sqrt((sa*sc)/(a*c)))))+a*(b*Sqrt((sa*sb)/(a*b))-c*Sqrt((sa*sc)/(a*c)))*(2*(-a+c)*sb*Sqrt((sa*sc)/(a*c))+(-(a*c)+SB)*(Sqrt((sa*sb)/(a*b))-Sqrt((sb*sc)/(b*c)))));
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X168(orbit) {
   /* begin vars */
   let sc=(a+b-c)/2;
   let sb=(c+a-b)/2;
   let s=(a+b+c)/2;
   let Sqrt=Math.sqrt;
   let sa=(b+c-a)/2;
   let area=triAreaHeron(a,b,c);
   /* end vars */
   let v1 = a*(2*area*sa-a*(c*Sqrt((s*sb)/(a*c))*(sa-sc)+b*(c*Sqrt((s*sa)/(b*c))+(sa-sb)*Sqrt((s*sc)/(a*b)))));
   let v2 = b*(2*area*sb-b*(c*(a*Sqrt((s*sb)/(a*c))+Sqrt((s*sa)/(b*c))*(sb-sc))+a*(-sa+sb)*Sqrt((s*sc)/(a*b))));
   let v3 = c*(2*area*sc-c*(b*Sqrt((s*sa)/(b*c))*(-sb+sc)+a*(b*Sqrt((s*sc)/(a*b))+Sqrt((s*sb)/(a*c))*(-sa+sc))));
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X169(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let c3=c2*c;
   let b3=b2*b;
   let a3=a2*a;
   /* end vars */
   let v1 = a*(a3-a2*b+a*b2-b3-a2*c+b2*c+a*c2+b*c2-c3);
   let v2 = b*(-a3+a2*b-a*b2+b3+a2*c-b2*c+a*c2+b*c2-c3);
   let v3 = c*(-a3+a2*b+a*b2-b3+a2*c+b2*c-a*c2-b*c2+c3);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X170(orbit) {
   /* begin vars */
   let c2=c*c;
   let c3=c2*c;
   let b2=b*b;
   let b3=b2*b;
   let a2=a*a;
   let a3=a2*a;
   let c5=c2*c3;
   let c4=c2*c2;
   let b5=b2*b3;
   let b4=b2*b2;
   let a4=a2*a2;
   let a5=a2*a3;
   /* end vars */
   let v1 = a*(a5*b-4*a4*b2+6*a3*b3-4*a2*b4+a*b5+a5*c-a4*b*c-2*a3*b2*c+2*a2*b3*c+a*b4*c-b5*c-4*a4*c2-2*a3*b*c2+4*a2*b2*c2-2*a*b3*c2+4*b4*c2+6*a3*c3+2*a2*b*c3-2*a*b2*c3-6*b3*c3-4*a2*c4+a*b*c4+4*b2*c4+a*c5-b*c5);
   let v2 = b*(a5*b-4*a4*b2+6*a3*b3-4*a2*b4+a*b5-a5*c+a4*b*c+2*a3*b2*c-2*a2*b3*c-a*b4*c+b5*c+4*a4*c2-2*a3*b*c2+4*a2*b2*c2-2*a*b3*c2-4*b4*c2-6*a3*c3-2*a2*b*c3+2*a*b2*c3+6*b3*c3+4*a2*c4+a*b*c4-4*b2*c4-a*c5+b*c5);
   let v3 = c*(-(a5*b)+4*a4*b2-6*a3*b3+4*a2*b4-a*b5+a5*c+a4*b*c-2*a3*b2*c-2*a2*b3*c+a*b4*c+b5*c-4*a4*c2+2*a3*b*c2+4*a2*b2*c2+2*a*b3*c2-4*b4*c2+6*a3*c3-2*a2*b*c3-2*a*b2*c3+6*b3*c3-4*a2*c4-a*b*c4-4*b2*c4+a*c5+b*c5);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X171(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   /* end vars */
   let v1 = a*(a2+b*c);
   let v2 = b*(b2+a*c);
   let v3 = c*(a*b+c2);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X172(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   /* end vars */
   let v1 = a2*(a2+b*c);
   let v2 = b2*(b2+a*c);
   let v3 = c2*(a*b+c2);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X173(orbit) {
   /* begin vars */
   let sc=(a+b-c)/2;
   let sb=(c+a-b)/2;
   let sa=(b+c-a)/2;
   let s=(a+b+c)/2;
   let Sqrt=Math.sqrt;
   /* end vars */
   let v1 = a*(-Sqrt((s*sa)/(b*c))+Sqrt((s*sb)/(a*c))+Sqrt((s*sc)/(a*b)));
   let v2 = b*(Sqrt((s*sa)/(b*c))-Sqrt((s*sb)/(a*c))+Sqrt((s*sc)/(a*b)));
   let v3 = c*(Sqrt((s*sa)/(b*c))+Sqrt((s*sb)/(a*c))-Sqrt((s*sc)/(a*b)));
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X174(orbit) {
   /* begin vars */
   let sa=(b+c-a)/2;
   let sc=(a+b-c)/2;
   let sb=(c+a-b)/2;
   let Sqrt=Math.sqrt;
   /* end vars */
   let v1 = Sqrt((sb*sc)/(b*c));
   let v2 = Sqrt((sa*sc)/(a*c));
   let v3 = Sqrt((sa*sb)/(a*b));
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X175(orbit) {
   /* begin vars */
   let area=triAreaHeron(a,b,c);
   let S=2*area;
   /* end vars */
   let v1 = (a+b-c)*(a-b+c)*(a*(a-b-c)+S);
   let v2 = (a+b-c)*(-a+b+c)*(b*(-a+b-c)+S);
   let v3 = (a-b+c)*(-a+b+c)*(c*(-a-b+c)+S);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X176(orbit) {
   /* begin vars */
   let area=triAreaHeron(a,b,c);
   let S=2*area;
   /* end vars */
   let v1 = (a+b-c)*(a-b+c)*(a*(a-b-c)-S);
   let v2 = (a+b-c)*(-a+b+c)*(b*(-a+b-c)-S);
   let v3 = (a-b+c)*(-a+b+c)*(c*(-a-b+c)-S);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X177(orbit) {
   /* begin vars */
   let sqrtc=sqrt(c);
   let sqrtb=sqrt(b);
   let sqrta=sqrt(a);
   let Sqrt=Math.sqrt;
   /* end vars */
   let v1 = (Sqrt(b*(a+b-c))*(a-b+c)+(a+b-c)*Sqrt(c*(a-b+c)))*sqrta;
   let v2 = (Sqrt(a*(a+b-c))*(-a+b+c)+(a+b-c)*Sqrt(c*(-a+b+c)))*sqrtb;
   let v3 = (Sqrt(a*(a-b+c))*(-a+b+c)+(a-b+c)*Sqrt(b*(-a+b+c)))*sqrtc;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X178(orbit) {
   /* begin vars */
   let Sqrt=Math.sqrt;
   /* end vars */
   let v1 = Sqrt((a+b-c)*c)+Sqrt(b*(a-b+c));
   let v2 = Sqrt((a+b-c)*c)+Sqrt(a*(-a+b+c));
   let v3 = Sqrt(b*(a-b+c))+Sqrt(a*(-a+b+c));
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X179(orbit) {
   /* begin vars */
   let cosC=lawOfCosines(c,a,b);
   let cosHalfC=cosHalfAngle(cosC);
   let cosQuarterC=cosHalfAngle(cosHalfC);
   let cosB=lawOfCosines(b,a,c);
   let cosHalfB=cosHalfAngle(cosB);
   let cosQuarterB=cosHalfAngle(cosHalfB);
   let cosA=lawOfCosines(a,b,c);
   let cosHalfA=cosHalfAngle(cosA);
   let cosQuarterA=cosHalfAngle(cosHalfA);
   let sinC=getSin(cosC);
   let secQuarterC=1/cosQuarterC;
   let sinB=getSin(cosB);
   let secQuarterB=1/cosQuarterB;
   let sinA=getSin(cosA);
   let secQuarterA=1/cosQuarterA;
   /* end vars */
   let v1 = Math.pow(secQuarterA,4)*sinA;
   let v2 = Math.pow(secQuarterB,4)*sinB;
   let v3 = Math.pow(secQuarterC,4)*sinC;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X180(orbit) {
   /* begin vars */
   let cosC=lawOfCosines(c,a,b);
   let cosB=lawOfCosines(b,a,c);
   let cosA=lawOfCosines(a,b,c);
   let cosHalfC=cosHalfAngle(cosC);
   let cosQuarterC=cosHalfAngle(cosHalfC);
   let cosHalfB=cosHalfAngle(cosB);
   let cosQuarterB=cosHalfAngle(cosHalfB);
   let cosHalfA=cosHalfAngle(cosA);
   let cosQuarterA=cosHalfAngle(cosHalfA);
   let sinC=getSin(cosC);
   let sinB=getSin(cosB);
   let sinA=getSin(cosA);
   let secQuarterC=1/cosQuarterC;
   let secQuarterB=1/cosQuarterB;
   let secQuarterA=1/cosQuarterA;
   /* end vars */
   let v1 = (1/(-1-2*cosQuarterB*cosQuarterB*cosQuarterC*cosQuarterC*secQuarterA*secQuarterA)+1/(1+2*cosQuarterA*cosQuarterA*cosQuarterC*cosQuarterC*secQuarterB*secQuarterB)+1/(1+2*cosQuarterA*cosQuarterA*cosQuarterB*cosQuarterB*secQuarterC*secQuarterC))*sinA;
   let v2 = (1/(1+2*cosQuarterB*cosQuarterB*cosQuarterC*cosQuarterC*secQuarterA*secQuarterA)+1/(-1-2*cosQuarterA*cosQuarterA*cosQuarterC*cosQuarterC*secQuarterB*secQuarterB)+1/(1+2*cosQuarterA*cosQuarterA*cosQuarterB*cosQuarterB*secQuarterC*secQuarterC))*sinB;
   let v3 = (1/(1+2*cosQuarterB*cosQuarterB*cosQuarterC*cosQuarterC*secQuarterA*secQuarterA)+1/(1+2*cosQuarterA*cosQuarterA*cosQuarterC*cosQuarterC*secQuarterB*secQuarterB)+1/(-1-2*cosQuarterA*cosQuarterA*cosQuarterB*cosQuarterB*secQuarterC*secQuarterC))*sinC;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X181(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   /* end vars */
   let v1 = a2*(a+b-c)*(a-b+c)*(b+c)*(b+c);
   let v2 = b2*(a+b-c)*(a+c)*(a+c)*(-a+b+c);
   let v3 = (a+b)*(a+b)*(a-b+c)*(-a+b+c)*c2;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X182(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let c4=c2*c2;
   let b4=b2*b2;
   let a4=a2*a2;
   /* end vars */
   let v1 = a2*(a4-a2*b2-a2*c2-2*b2*c2);
   let v2 = b2*(-(a2*b2)+b4-2*a2*c2-b2*c2);
   let v3 = c2*(-2*a2*b2-a2*c2-b2*c2+c4);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X183(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let c4=c2*c2;
   let b4=b2*b2;
   let a4=a2*a2;
   /* end vars */
   let v1 = a4-a2*b2-a2*c2-2*b2*c2;
   let v2 = -(a2*b2)+b4-2*a2*c2-b2*c2;
   let v3 = -2*a2*b2-a2*c2-b2*c2+c4;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X184(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let c4=c2*c2;
   let b4=b2*b2;
   let a4=a2*a2;
   /* end vars */
   let v1 = a4*(a2-b2-c2);
   let v2 = b4*(-a2+b2-c2);
   let v3 = (-a2-b2+c2)*c4;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X185(orbit) {
   /* begin vars */
   let a2=a*a;
   let a4=a2*a2;
   let c2=c*c;
   let c4=c2*c2;
   let b2=b*b;
   let b4=b2*b2;
   let a6=a2*a4;
   let c6=c2*c4;
   let b6=b2*b4;
   /* end vars */
   let v1 = a2*(a2-b2-c2)*(a4*b2-2*a2*b4+b6+a4*c2+4*a2*b2*c2-b4*c2-2*a2*c4-b2*c4+c6);
   let v2 = b2*(-a2+b2-c2)*(a6-2*a4*b2+a2*b4-a4*c2+4*a2*b2*c2+b4*c2-a2*c4-2*b2*c4+c6);
   let v3 = c2*(-a2-b2+c2)*(a6-a4*b2-a2*b4+b6-2*a4*c2+4*a2*b2*c2-2*b4*c2+a2*c4+b2*c4);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X186(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   /* end vars */
   let v1 = a2*(a2+b2-c2)*(a2-b2-b*c-c2)*(a2-b2+b*c-c2)*(a2-b2+c2);
   let v2 = b2*(a2+b2-c2)*(-a2+b2-a*c-c2)*(-a2+b2+a*c-c2)*(-a2+b2+c2);
   let v3 = c2*(a2-b2+c2)*(-a2-a*b-b2+c2)*(-a2+a*b-b2+c2)*(-a2+b2+c2);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X187(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   /* end vars */
   let v1 = a2*(2*a2-b2-c2);
   let v2 = b2*(-a2+2*b2-c2);
   let v3 = c2*(-a2-b2+2*c2);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X188(orbit) {
   /* begin vars */
   let Sqrt=Math.sqrt;
   /* end vars */
   let v1 = a*Sqrt(b*c*(-a+b+c));
   let v2 = b*Sqrt(a*c*(a-b+c));
   let v3 = Sqrt(a*b*(a+b-c))*c;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X189(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let c3=c2*c;
   let b3=b2*b;
   let a3=a2*a;
   /* end vars */
   let v1 = (a3-a2*b-a*b2+b3+a2*c+2*a*b*c+b2*c-a*c2-b*c2-c3)*(a3+a2*b-a*b2-b3-a2*c+2*a*b*c-b2*c-a*c2+b*c2+c3);
   let v2 = (a3-a2*b-a*b2+b3+a2*c+2*a*b*c+b2*c-a*c2-b*c2-c3)*(-a3-a2*b+a*b2+b3-a2*c+2*a*b*c-b2*c+a*c2-b*c2+c3);
   let v3 = (-a3-a2*b+a*b2+b3-a2*c+2*a*b*c-b2*c+a*c2-b*c2+c3)*(a3+a2*b-a*b2-b3-a2*c+2*a*b*c-b2*c-a*c2+b*c2+c3);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X190(orbit) {
   /* begin vars */

   /* end vars */
   let v1 = (a-b)*(a-c);
   let v2 = (-a+b)*(b-c);
   let v3 = (-a+c)*(-b+c);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X191(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let c3=c2*c;
   let b3=b2*b;
   let a3=a2*a;
   /* end vars */
   let v1 = a*(a3+a2*b-a*b2-b3+a2*c-a*b*c-b2*c-a*c2-b*c2-c3);
   let v2 = b*(-a3-a2*b+a*b2+b3-a2*c-a*b*c+b2*c-a*c2-b*c2-c3);
   let v3 = c*(-a3-a2*b-a*b2-b3-a2*c-a*b*c-b2*c+a*c2+b*c2+c3);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X192(orbit) {
   /* begin vars */

   /* end vars */
   let v1 = a*b+a*c-b*c;
   let v2 = a*b-a*c+b*c;
   let v3 = -(a*b)+a*c+b*c;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X193(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   /* end vars */
   let v1 = -3*a2+b2+c2;
   let v2 = a2-3*b2+c2;
   let v3 = a2+b2-3*c2;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X194(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   /* end vars */
   let v1 = a2*b2+a2*c2-b2*c2;
   let v2 = a2*b2-a2*c2+b2*c2;
   let v3 = -(a2*b2)+a2*c2+b2*c2;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X195(orbit) {
   /* begin vars */
   let c2=c*c;
   let c4=c2*c2;
   let c6=c2*c4;
   let b2=b*b;
   let b4=b2*b2;
   let b6=b2*b4;
   let a2=a*a;
   let a4=a2*a2;
   let a6=a2*a4;
   let c8=c2*c6;
   let b8=b2*b6;
   let a8=a2*a6;
   /* end vars */
   let v1 = a2*(a8-4*a6*b2+6*a4*b4-4*a2*b6+b8-4*a6*c2+5*a4*b2*c2+a2*b4*c2-2*b6*c2+6*a4*c4+a2*b2*c4+2*b4*c4-4*a2*c6-2*b2*c6+c8);
   let v2 = b2*(a8-4*a6*b2+6*a4*b4-4*a2*b6+b8-2*a6*c2+a4*b2*c2+5*a2*b4*c2-4*b6*c2+2*a4*c4+a2*b2*c4+6*b4*c4-2*a2*c6-4*b2*c6+c8);
   let v3 = c2*(a8-2*a6*b2+2*a4*b4-2*a2*b6+b8-4*a6*c2+a4*b2*c2+a2*b4*c2-4*b6*c2+6*a4*c4+5*a2*b2*c4+6*b4*c4-4*a2*c6-4*b2*c6+c8);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X196(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let c3=c2*c;
   let b3=b2*b;
   let a3=a2*a;
   /* end vars */
   let v1 = (a+b-c)*(a-b+c)*(a2+b2-c2)*(a2-b2+c2)*(a3+a2*b-a*b2-b3+a2*c-2*a*b*c+b2*c-a*c2+b*c2-c3);
   let v2 = (a+b-c)*(-a+b+c)*(a2+b2-c2)*(-a2+b2+c2)*(-a3-a2*b+a*b2+b3+a2*c-2*a*b*c+b2*c+a*c2-b*c2-c3);
   let v3 = (a-b+c)*(-a+b+c)*(a2-b2+c2)*(-a2+b2+c2)*(-a3+a2*b+a*b2-b3-a2*c-2*a*b*c-b2*c+a*c2+b*c2+c3);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X197(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let c4=c2*c2;
   let b4=b2*b2;
   let a4=a2*a2;
   /* end vars */
   let v1 = a2*(a4-b4+2*a2*b*c-2*a*b2*c-2*a*b*c2+2*b2*c2-c4);
   let v2 = b2*(-a4+b4-2*a2*b*c+2*a*b2*c+2*a2*c2-2*a*b*c2-c4);
   let v3 = c2*(-a4+2*a2*b2-b4-2*a2*b*c-2*a*b2*c+2*a*b*c2+c4);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X198(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let c3=c2*c;
   let b3=b2*b;
   let a3=a2*a;
   /* end vars */
   let v1 = a2*(a3+a2*b-a*b2-b3+a2*c-2*a*b*c+b2*c-a*c2+b*c2-c3);
   let v2 = b2*(-a3-a2*b+a*b2+b3+a2*c-2*a*b*c+b2*c+a*c2-b*c2-c3);
   let v3 = c2*(-a3+a2*b+a*b2-b3-a2*c-2*a*b*c-b2*c+a*c2+b*c2+c3);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X199(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let c4=c2*c2;
   let c3=c2*c;
   let b4=b2*b2;
   let b3=b2*b;
   let a3=a2*a;
   let a4=a2*a2;
   /* end vars */
   let v1 = a2*(a4+a3*b-a*b3-b4+a3*c+a2*b*c-a*b2*c-b3*c-a*b*c2-a*c3-b*c3-c4);
   let v2 = b2*(-a4-a3*b+a*b3+b4-a3*c-a2*b*c+a*b2*c+b3*c-a*b*c2-a*c3-b*c3-c4);
   let v3 = c2*(-a4-a3*b-a*b3-b4-a3*c-a2*b*c-a*b2*c-b3*c+a*b*c2+a*c3+b*c3+c4);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X200(orbit) {
   /* begin vars */

   /* end vars */
   let v1 = a*(a-b-c)*(a-b-c);
   let v2 = b*(-a+b-c)*(-a+b-c);
   let v3 = c*(-a-b+c)*(-a-b+c);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X201(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   /* end vars */
   let v1 = a*(a+b-c)*(a-b+c)*(b+c)*(b+c)*(a2-b2-c2);
   let v2 = b*(a+b-c)*(a+c)*(a+c)*(-a+b+c)*(-a2+b2-c2);
   let v3 = (a+b)*(a+b)*c*(a-b+c)*(-a+b+c)*(-a2-b2+c2);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X202(orbit) {
   /* begin vars */
   let area=triAreaHeron(a,b,c);
   let S=2*area;
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   /* end vars */
   let v1 = a2*(a2-b2+4*b*c-c2+2*S*sqrt3);
   let v2 = b2*(-a2+b2+4*a*c-c2+2*S*sqrt3);
   let v3 = c2*(-a2+4*a*b-b2+c2+2*S*sqrt3);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X203(orbit) {
   /* begin vars */
   let area=triAreaHeron(a,b,c);
   let S=2*area;
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   /* end vars */
   let v1 = a2*(a2-b2+4*b*c-c2-2*S*sqrt3);
   let v2 = b2*(-a2+b2+4*a*c-c2-2*S*sqrt3);
   let v3 = c2*(-a2+4*a*b-b2+c2-2*S*sqrt3);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X204(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let c4=c2*c2;
   let b4=b2*b2;
   let a4=a2*a2;
   /* end vars */
   let v1 = a*(a2+b2-c2)*(a2-b2+c2)*(3*a4-2*a2*b2-b4-2*a2*c2+2*b2*c2-c4);
   let v2 = b*(a2+b2-c2)*(-a2+b2+c2)*(-a4-2*a2*b2+3*b4+2*a2*c2-2*b2*c2-c4);
   let v3 = c*(a2-b2+c2)*(-a2+b2+c2)*(-a4+2*a2*b2-b4-2*a2*c2-2*b2*c2+3*c4);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X205(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let c3=c2*c;
   let b3=b2*b;
   let c4=c2*c2;
   let b4=b2*b2;
   let a4=a2*a2;
   let a3=a2*a;
   /* end vars */
   let v1 = a3*(a4-b4+2*a2*b*c-2*a*b2*c-2*a*b*c2+2*b2*c2-c4);
   let v2 = b3*(-a4+b4-2*a2*b*c+2*a*b2*c+2*a2*c2-2*a*b*c2-c4);
   let v3 = c3*(-a4+2*a2*b2-b4-2*a2*b*c-2*a*b2*c+2*a*b*c2+c4);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X206(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let c4=c2*c2;
   let b4=b2*b2;
   let a4=a2*a2;
   /* end vars */
   let v1 = a4*(a4-b4-c4);
   let v2 = b4*(-a4+b4-c4);
   let v3 = c4*(-a4-b4+c4);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X207(orbit) {
   /* begin vars */
   let c2=c*c;
   let c4=c2*c2;
   let c3=c2*c;
   let b2=b*b;
   let b4=b2*b2;
   let b3=b2*b;
   let a2=a*a;
   let a3=a2*a;
   let a4=a2*a2;
   let c6=c2*c4;
   let c5=c2*c3;
   let b6=b2*b4;
   let b5=b2*b3;
   let a5=a2*a3;
   let a6=a2*a4;
   /* end vars */
   let v1 = a*(a+b-c)*(a-b+c)*(a2+b2-c2)*(a2-b2+c2)*(a6-2*a5*b-a4*b2+4*a3*b3-a2*b4-2*a*b5+b6-2*a5*c-2*a4*b*c+2*a*b4*c+2*b5*c-a4*c2+2*a2*b2*c2-b4*c2+4*a3*c3-4*b3*c3-a2*c4+2*a*b*c4-b2*c4-2*a*c5+2*b*c5+c6);
   let v2 = b*(a+b-c)*(-a+b+c)*(a2+b2-c2)*(-a2+b2+c2)*(a6-2*a5*b-a4*b2+4*a3*b3-a2*b4-2*a*b5+b6+2*a5*c+2*a4*b*c-2*a*b4*c-2*b5*c-a4*c2+2*a2*b2*c2-b4*c2-4*a3*c3+4*b3*c3-a2*c4+2*a*b*c4-b2*c4+2*a*c5-2*b*c5+c6);
   let v3 = c*(a-b+c)*(-a+b+c)*(a2-b2+c2)*(-a2+b2+c2)*(a6+2*a5*b-a4*b2-4*a3*b3-a2*b4+2*a*b5+b6-2*a5*c+2*a4*b*c+2*a*b4*c-2*b5*c-a4*c2+2*a2*b2*c2-b4*c2+4*a3*c3+4*b3*c3-a2*c4-2*a*b*c4-b2*c4-2*a*c5-2*b*c5+c6);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X208(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let c3=c2*c;
   let b3=b2*b;
   let a3=a2*a;
   /* end vars */
   let v1 = a*(a+b-c)*(a-b+c)*(a2+b2-c2)*(a2-b2+c2)*(a3+a2*b-a*b2-b3+a2*c-2*a*b*c+b2*c-a*c2+b*c2-c3);
   let v2 = b*(a+b-c)*(-a+b+c)*(a2+b2-c2)*(-a2+b2+c2)*(-a3-a2*b+a*b2+b3+a2*c-2*a*b*c+b2*c+a*c2-b*c2-c3);
   let v3 = c*(a-b+c)*(-a+b+c)*(a2-b2+c2)*(-a2+b2+c2)*(-a3+a2*b+a*b2-b3-a2*c-2*a*b*c-b2*c+a*c2+b*c2+c3);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X209(orbit) {
   /* begin vars */
   let a2=a*a;
   let c2=c*c;
   let b2=b*b;
   let a3=a2*a;
   let c3=c2*c;
   let b3=b2*b;
   /* end vars */
   let v1 = a2*(b+c)*(a2*b-b3+a2*c+a*b*c-c3);
   let v2 = b2*(a+c)*(-a3+a*b2+a*b*c+b2*c-c3);
   let v3 = (a+b)*c2*(-a3-b3+a*b*c+a*c2+b*c2);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X210(orbit) {
   /* begin vars */

   /* end vars */
   let v1 = a*(a-b-c)*(b+c);
   let v2 = b*(-a+b-c)*(a+c);
   let v3 = (a+b)*c*(-a-b+c);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X211(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let c4=c2*c2;
   let b4=b2*b2;
   let a4=a2*a2;
   /* end vars */
   let v1 = a4*(b2+c2)*(a2*b2-b4+a2*c2+b2*c2-c4);
   let v2 = b4*(a2+c2)*(-a4+a2*b2+a2*c2+b2*c2-c4);
   let v3 = (a2+b2)*(-a4+a2*b2-b4+a2*c2+b2*c2)*c4;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X212(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let c3=c2*c;
   let b3=b2*b;
   let a3=a2*a;
   /* end vars */
   let v1 = a3*(a-b-c)*(a2-b2-c2);
   let v2 = b3*(-a+b-c)*(-a2+b2-c2);
   let v3 = (-a-b+c)*(-a2-b2+c2)*c3;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X213(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let c3=c2*c;
   let b3=b2*b;
   let a3=a2*a;
   /* end vars */
   let v1 = a3*(b+c);
   let v2 = b3*(a+c);
   let v3 = (a+b)*c3;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X214(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   /* end vars */
   let v1 = a*(2*a-b-c)*(a2-b2+b*c-c2);
   let v2 = b*(-a+2*b-c)*(-a2+b2+a*c-c2);
   let v3 = c*(-a-b+2*c)*(-a2+a*b-b2+c2);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X215(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let c4=c2*c2;
   let b4=b2*b2;
   let a4=a2*a2;
   /* end vars */
   let v1 = a4*(a-b-c)*(a2-b2+b*c-c2)*(a2-b2+b*c-c2);
   let v2 = b4*(-a+b-c)*(-a2+b2+a*c-c2)*(-a2+b2+a*c-c2);
   let v3 = (-a-b+c)*(-a2+a*b-b2+c2)*(-a2+a*b-b2+c2)*c4;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X216(orbit) {
   /* begin vars */
   let a2=a*a;
   let c2=c*c;
   let b2=b*b;
   let a4=a2*a2;
   let c4=c2*c2;
   let b4=b2*b2;
   /* end vars */
   let v1 = a2*(a2-b2-c2)*(a2*b2-b4+a2*c2+2*b2*c2-c4);
   let v2 = b2*(-a2+b2-c2)*(-a4+a2*b2+2*a2*c2+b2*c2-c4);
   let v3 = c2*(-a2-b2+c2)*(-a4+2*a2*b2-b4+a2*c2+b2*c2);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X217(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let c4=c2*c2;
   let b4=b2*b2;
   let a4=a2*a2;
   /* end vars */
   let v1 = a4*(a2-b2-c2)*(a2*b2-b4+a2*c2+2*b2*c2-c4);
   let v2 = b4*(-a2+b2-c2)*(-a4+a2*b2+2*a2*c2+b2*c2-c4);
   let v3 = (-a2-b2+c2)*(-a4+2*a2*b2-b4+a2*c2+b2*c2)*c4;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X218(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   /* end vars */
   let v1 = a2*(a2-2*a*b+b2-2*a*c+c2);
   let v2 = b2*(a2-2*a*b+b2-2*b*c+c2);
   let v3 = c2*(a2+b2-2*a*c-2*b*c+c2);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X219(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   /* end vars */
   let v1 = a2*(a-b-c)*(a2-b2-c2);
   let v2 = b2*(-a+b-c)*(-a2+b2-c2);
   let v3 = (-a-b+c)*c2*(-a2-b2+c2);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X220(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   /* end vars */
   let v1 = a2*(a-b-c)*(a-b-c);
   let v2 = b2*(-a+b-c)*(-a+b-c);
   let v3 = (-a-b+c)*(-a-b+c)*c2;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X221(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let c3=c2*c;
   let b3=b2*b;
   let a3=a2*a;
   /* end vars */
   let v1 = a2*(a+b-c)*(a-b+c)*(a3+a2*b-a*b2-b3+a2*c-2*a*b*c+b2*c-a*c2+b*c2-c3);
   let v2 = b2*(a+b-c)*(-a+b+c)*(-a3-a2*b+a*b2+b3+a2*c-2*a*b*c+b2*c+a*c2-b*c2-c3);
   let v3 = (a-b+c)*(-a+b+c)*c2*(-a3+a2*b+a*b2-b3-a2*c-2*a*b*c-b2*c+a*c2+b*c2+c3);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X222(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   /* end vars */
   let v1 = a2*(a+b-c)*(a-b+c)*(a2-b2-c2);
   let v2 = b2*(a+b-c)*(-a+b+c)*(-a2+b2-c2);
   let v3 = (a-b+c)*(-a+b+c)*c2*(-a2-b2+c2);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X223(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let c3=c2*c;
   let b3=b2*b;
   let a3=a2*a;
   /* end vars */
   let v1 = a*(a+b-c)*(a-b+c)*(a3+a2*b-a*b2-b3+a2*c-2*a*b*c+b2*c-a*c2+b*c2-c3);
   let v2 = b*(a+b-c)*(-a+b+c)*(-a3-a2*b+a*b2+b3+a2*c-2*a*b*c+b2*c+a*c2-b*c2-c3);
   let v3 = c*(a-b+c)*(-a+b+c)*(-a3+a2*b+a*b2-b3-a2*c-2*a*b*c-b2*c+a*c2+b*c2+c3);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X224(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let c4=c2*c2;
   let c3=c2*c;
   let b4=b2*b2;
   let b3=b2*b;
   let a3=a2*a;
   let a4=a2*a2;
   /* end vars */
   let v1 = a*(a2-b2-c2)*(a4-2*a3*b+2*a*b3-b4-2*a3*c-2*a2*b*c+2*b2*c2+2*a*c3-c4);
   let v2 = b*(-a2+b2-c2)*(-a4+2*a3*b-2*a*b3+b4-2*a*b2*c-2*b3*c+2*a2*c2+2*b*c3-c4);
   let v3 = c*(-a2-b2+c2)*(-a4+2*a2*b2-b4+2*a3*c+2*b3*c-2*a*b*c2-2*a*c3-2*b*c3+c4);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X225(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   /* end vars */
   let v1 = (a+b-c)*(a-b+c)*(b+c)*(a2+b2-c2)*(a2-b2+c2);
   let v2 = (a+b-c)*(a+c)*(-a+b+c)*(a2+b2-c2)*(-a2+b2+c2);
   let v3 = (a+b)*(a-b+c)*(-a+b+c)*(a2-b2+c2)*(-a2+b2+c2);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X226(orbit) {
   /* begin vars */

   /* end vars */
   let v1 = (a+b-c)*(a-b+c)*(b+c);
   let v2 = (a+b-c)*(a+c)*(-a+b+c);
   let v3 = (a+b)*(a-b+c)*(-a+b+c);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X227(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let c3=c2*c;
   let b3=b2*b;
   let a3=a2*a;
   /* end vars */
   let v1 = a*(a+b-c)*(a-b+c)*(b+c)*(a3+a2*b-a*b2-b3+a2*c-2*a*b*c+b2*c-a*c2+b*c2-c3);
   let v2 = b*(a+b-c)*(a+c)*(-a+b+c)*(-a3-a2*b+a*b2+b3+a2*c-2*a*b*c+b2*c+a*c2-b*c2-c3);
   let v3 = (a+b)*c*(a-b+c)*(-a+b+c)*(-a3+a2*b+a*b2-b3-a2*c-2*a*b*c-b2*c+a*c2+b*c2+c3);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X228(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let c3=c2*c;
   let b3=b2*b;
   let a3=a2*a;
   /* end vars */
   let v1 = a3*(b+c)*(a2-b2-c2);
   let v2 = b3*(a+c)*(-a2+b2-c2);
   let v3 = (a+b)*(-a2-b2+c2)*c3;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X229(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let c4=c2*c2;
   let b4=b2*b2;
   let a4=a2*a2;
   /* end vars */
   let v1 = a*(a+b)*(a+c)*(a4-b4+a2*b*c+a*b2*c+a*b*c2+2*b2*c2-c4);
   let v2 = b*(a+b)*(b+c)*(-a4+b4+a2*b*c+a*b2*c+2*a2*c2+a*b*c2-c4);
   let v3 = c*(a+c)*(b+c)*(-a4+2*a2*b2-b4+a2*b*c+a*b2*c+a*b*c2+c4);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X230(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let c4=c2*c2;
   let b4=b2*b2;
   let a4=a2*a2;
   /* end vars */
   let v1 = 2*a4-a2*b2+b4-a2*c2-2*b2*c2+c4;
   let v2 = a4-a2*b2+2*b4-2*a2*c2-b2*c2+c4;
   let v3 = a4-2*a2*b2+b4-a2*c2-b2*c2+2*c4;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X231(orbit) {
   /* begin vars */
   let c2=c*c;
   let c4=c2*c2;
   let c6=c2*c4;
   let b2=b*b;
   let b4=b2*b2;
   let b6=b2*b4;
   let a2=a*a;
   let a4=a2*a2;
   let a6=a2*a4;
   let c8=c2*c6;
   let b8=b2*b6;
   let a8=a2*a6;
   /* end vars */
   let v1 = 2*a8-4*a6*b2+3*a4*b4-2*a2*b6+b8-4*a6*c2+2*a2*b4*c2-4*b6*c2+3*a4*c4+2*a2*b2*c4+6*b4*c4-2*a2*c6-4*b2*c6+c8;
   let v2 = a8-2*a6*b2+3*a4*b4-4*a2*b6+2*b8-4*a6*c2+2*a4*b2*c2-4*b6*c2+6*a4*c4+2*a2*b2*c4+3*b4*c4-4*a2*c6-2*b2*c6+c8;
   let v3 = a8-4*a6*b2+6*a4*b4-4*a2*b6+b8-2*a6*c2+2*a4*b2*c2+2*a2*b4*c2-2*b6*c2+3*a4*c4+3*b4*c4-4*a2*c6-4*b2*c6+2*c8;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X232(orbit) {
   /* begin vars */
   let a2=a*a;
   let c2=c*c;
   let b2=b*b;
   let a4=a2*a2;
   let c4=c2*c2;
   let b4=b2*b2;
   /* end vars */
   let v1 = a2*(a2+b2-c2)*(a2-b2+c2)*(a2*b2-b4+a2*c2-c4);
   let v2 = b2*(a2+b2-c2)*(-a2+b2+c2)*(-a4+a2*b2+b2*c2-c4);
   let v3 = c2*(a2-b2+c2)*(-a2+b2+c2)*(-a4-b4+a2*c2+b2*c2);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X233(orbit) {
   /* begin vars */
   let a2=a*a;
   let c2=c*c;
   let b2=b*b;
   let a4=a2*a2;
   let c4=c2*c2;
   let b4=b2*b2;
   /* end vars */
   let v1 = (a2*b2-b4+a2*c2+2*b2*c2-c4)*(2*a4-3*a2*b2+b4-3*a2*c2-2*b2*c2+c4);
   let v2 = (-a4+a2*b2+2*a2*c2+b2*c2-c4)*(a4-3*a2*b2+2*b4-2*a2*c2-3*b2*c2+c4);
   let v3 = (-a4+2*a2*b2-b4+a2*c2+b2*c2)*(a4-2*a2*b2+b4-3*a2*c2-3*b2*c2+2*c4);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X234(orbit) {
   /* begin vars */
   let s=(a+b+c)/2;
   let sa=(b+c-a)/2;
   let Sqrt=Math.sqrt;
   let sc=(a+b-c)/2;
   let sb=(c+a-b)/2;
   let s2=s*s;
   /* end vars */
   let v1 = s2*sb*sc*(Sqrt((s*sb)/(a*c))+Sqrt((s*sc)/(a*b)));
   let v2 = s2*sa*sc*(Sqrt((s*sa)/(b*c))+Sqrt((s*sc)/(a*b)));
   let v3 = s2*sa*sb*(Sqrt((s*sa)/(b*c))+Sqrt((s*sb)/(a*c)));
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X235(orbit) {
   /* begin vars */
   let a2=a*a;
   let a4=a2*a2;
   let c2=c*c;
   let c4=c2*c2;
   let b2=b*b;
   let b4=b2*b2;
   let a6=a2*a4;
   let c6=c2*c4;
   let b6=b2*b4;
   /* end vars */
   let v1 = (a2+b2-c2)*(a2-b2+c2)*(a4*b2-2*a2*b4+b6+a4*c2+4*a2*b2*c2-b4*c2-2*a2*c4-b2*c4+c6);
   let v2 = (a2+b2-c2)*(-a2+b2+c2)*(a6-2*a4*b2+a2*b4-a4*c2+4*a2*b2*c2+b4*c2-a2*c4-2*b2*c4+c6);
   let v3 = (a2-b2+c2)*(-a2+b2+c2)*(a6-a4*b2-a2*b4+b6-2*a4*c2+4*a2*b2*c2-2*b4*c2+a2*c4+b2*c4);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X236(orbit) {
   /* begin vars */
   let sc=(a+b-c)/2;
   let sb=(c+a-b)/2;
   let sa=(b+c-a)/2;
   let s=(a+b+c)/2;
   let Sqrt=Math.sqrt;
   /* end vars */
   let v1 = Sqrt((s*sa)/(b*c))*(Sqrt((s*sa)/(b*c))-Sqrt((s*sb)/(a*c))-Sqrt((s*sc)/(a*b)));
   let v2 = Sqrt((s*sb)/(a*c))*(-Sqrt((s*sa)/(b*c))+Sqrt((s*sb)/(a*c))-Sqrt((s*sc)/(a*b)));
   let v3 = Sqrt((s*sc)/(a*b))*(-Sqrt((s*sa)/(b*c))-Sqrt((s*sb)/(a*c))+Sqrt((s*sc)/(a*b)));
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X237(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let c4=c2*c2;
   let b4=b2*b2;
   let a4=a2*a2;
   /* end vars */
   let v1 = a4*(a2*b2-b4+a2*c2-c4);
   let v2 = b4*(-a4+a2*b2+b2*c2-c4);
   let v3 = (-a4-b4+a2*c2+b2*c2)*c4;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X238(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   /* end vars */
   let v1 = a*(a2-b*c);
   let v2 = b*(b2-a*c);
   let v3 = c*(-(a*b)+c2);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X239(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   /* end vars */
   let v1 = a2-b*c;
   let v2 = b2-a*c;
   let v3 = -(a*b)+c2;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X240(orbit) {
   /* begin vars */
   let a2=a*a;
   let c2=c*c;
   let b2=b*b;
   let a4=a2*a2;
   let c4=c2*c2;
   let b4=b2*b2;
   /* end vars */
   let v1 = a*(a2+b2-c2)*(a2-b2+c2)*(a2*b2-b4+a2*c2-c4);
   let v2 = b*(a2+b2-c2)*(-a2+b2+c2)*(-a4+a2*b2+b2*c2-c4);
   let v3 = c*(a2-b2+c2)*(-a2+b2+c2)*(-a4-b4+a2*c2+b2*c2);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X241(orbit) {
   /* begin vars */
   let a2=a*a;
   let c2=c*c;
   let b2=b*b;
   /* end vars */
   let v1 = a*(a+b-c)*(a-b+c)*(a*b-b2+a*c-c2);
   let v2 = b*(a+b-c)*(-a+b+c)*(-a2+a*b+b*c-c2);
   let v3 = c*(a-b+c)*(-a+b+c)*(-a2-b2+a*c+b*c);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X242(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   /* end vars */
   let v1 = (a2-b*c)*(a2+b2-c2)*(a2-b2+c2);
   let v2 = (b2-a*c)*(a2+b2-c2)*(-a2+b2+c2);
   let v3 = (-(a*b)+c2)*(a2-b2+c2)*(-a2+b2+c2);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X243(orbit) {
   /* begin vars */
   let c2=c*c;
   let a2=a*a;
   let b2=b*b;
   let c4=c2*c2;
   let a3=a2*a;
   let b4=b2*b2;
   let c3=c2*c;
   let b3=b2*b;
   let a4=a2*a2;
   /* end vars */
   let v1 = (a-b-c)*(a2+b2-c2)*(a2-b2+c2)*(a4-a2*b2+a2*b*c-b3*c-a2*c2+2*b2*c2-b*c3);
   let v2 = (-a+b-c)*(a2+b2-c2)*(-a2+b2+c2)*(-(a2*b2)+b4-a3*c+a*b2*c+2*a2*c2-b2*c2-a*c3);
   let v3 = (-a-b+c)*(a2-b2+c2)*(-a2+b2+c2)*(-(a3*b)+2*a2*b2-a*b3-a2*c2+a*b*c2-b2*c2+c4);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X244(orbit) {
   /* begin vars */

   /* end vars */
   let v1 = a*(b-c)*(b-c);
   let v2 = b*(-a+c)*(-a+c);
   let v3 = (a-b)*(a-b)*c;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X245(orbit) {
   /* begin vars */
   let a2=a*a;
   let a4=a2*a2;
   let c2=c*c;
   let c4=c2*c2;
   let c6=c2*c4;
   let c3=c2*c;
   let b2=b*b;
   let b4=b2*b2;
   let b6=b2*b4;
   let b3=b2*b;
   let a3=a2*a;
   let a6=a2*a4;
   let c7=c*c6;
   let c5=c2*c3;
   let b7=b*b6;
   let b5=b2*b3;
   let a5=a2*a3;
   let a7=a*a6;
   /* end vars */
   let v1 = a*(b-c)*(b-c)*(b+c)*(a7-2*a5*b2+a4*b3+a3*b4-2*a2*b5+b7-a5*b*c+a4*b2*c+a3*b3*c-2*a2*b4*c+b6*c-2*a5*c2+a4*b*c2+3*a3*b2*c2-a*b4*c2+a4*c3+a3*b*c3-a*b3*c3+a3*c4-2*a2*b*c4-a*b2*c4-2*a2*c5+b*c6+c7);
   let v2 = b*(-a+c)*(-a+c)*(a+c)*(a7-2*a5*b2+a4*b3+a3*b4-2*a2*b5+b7+a6*c-2*a4*b2*c+a3*b3*c+a2*b4*c-a*b5*c-a4*b*c2+3*a2*b3*c2+a*b4*c2-2*b5*c2-a3*b*c3+a*b3*c3+b4*c3-a2*b*c4-2*a*b2*c4+b3*c4-2*b2*c5+a*c6+c7);
   let v3 = (a-b)*(a-b)*(a+b)*c*(a7+a6*b+a*b6+b7-a4*b2*c-a3*b3*c-a2*b4*c-2*a5*c2-2*a4*b*c2-2*a*b4*c2-2*b5*c2+a4*c3+a3*b*c3+3*a2*b2*c3+a*b3*c3+b4*c3+a3*c4+a2*b*c4+a*b2*c4+b3*c4-2*a2*c5-a*b*c5-2*b2*c5+c7);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X246(orbit) {
   /* begin vars */
   let c2=c*c;
   let c4=c2*c2;
   let c6=c2*c4;
   let b2=b*b;
   let b4=b2*b2;
   let b6=b2*b4;
   let a2=a*a;
   let a4=a2*a2;
   let a6=a2*a4;
   let c8=c2*c6;
   let b8=b2*b6;
   let a8=a2*a6;
   /* end vars */
   let v1 = a2*(b-c)*(b-c)*(b+c)*(b+c)*(a8-4*a6*b2+6*a4*b4-4*a2*b6+b8-4*a6*c2+7*a4*b2*c2-5*a2*b4*c2+2*b6*c2+6*a4*c4-5*a2*b2*c4-4*a2*c6+2*b2*c6+c8);
   let v2 = b2*(-a+c)*(-a+c)*(a+c)*(a+c)*(a8-4*a6*b2+6*a4*b4-4*a2*b6+b8+2*a6*c2-5*a4*b2*c2+7*a2*b4*c2-4*b6*c2-5*a2*b2*c4+6*b4*c4+2*a2*c6-4*b2*c6+c8);
   let v3 = (a-b)*(a-b)*(a+b)*(a+b)*c2*(a8+2*a6*b2+2*a2*b6+b8-4*a6*c2-5*a4*b2*c2-5*a2*b4*c2-4*b6*c2+6*a4*c4+7*a2*b2*c4+6*b4*c4-4*a2*c6-4*b2*c6+c8);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X247(orbit) {
   /* begin vars */
   let c2=c*c;
   let c4=c2*c2;
   let c6=c2*c4;
   let b2=b*b;
   let b4=b2*b2;
   let b6=b2*b4;
   let a2=a*a;
   let a4=a2*a2;
   let a6=a2*a4;
   let c8=c2*c6;
   let b8=b2*b6;
   let a8=a2*a6;
   /* end vars */
   let v1 = (b-c)*(b-c)*(b+c)*(b+c)*(Math.pow(a,10)+Math.pow(b,10)-4*a8*b2+5*a6*b4-a4*b6-2*a2*b8+Math.pow(c,10)-4*a8*c2+7*a6*b2*c2-6*a4*b4*c2+5*a2*b6*c2-2*b8*c2+5*a6*c4-6*a4*b2*c4-2*a2*b4*c4+b6*c4-a4*c6+5*a2*b2*c6+b4*c6-2*a2*c8-2*b2*c8);
   let v2 = (-a+c)*(-a+c)*(a+c)*(a+c)*(Math.pow(a,10)+Math.pow(b,10)-2*a8*b2-a6*b4+5*a4*b6-4*a2*b8+Math.pow(c,10)-2*a8*c2+5*a6*b2*c2-6*a4*b4*c2+7*a2*b6*c2-4*b8*c2+a6*c4-2*a4*b2*c4-6*a2*b4*c4+5*b6*c4+a4*c6+5*a2*b2*c6-b4*c6-2*a2*c8-2*b2*c8);
   let v3 = (a-b)*(a-b)*(a+b)*(a+b)*(Math.pow(a,10)+Math.pow(b,10)-2*a8*b2+a6*b4+a4*b6-2*a2*b8+Math.pow(c,10)-2*a8*c2+5*a6*b2*c2-2*a4*b4*c2+5*a2*b6*c2-2*b8*c2-a6*c4-6*a4*b2*c4-6*a2*b4*c4-b6*c4+5*a4*c6+7*a2*b2*c6+5*b4*c6-4*a2*c8-4*b2*c8);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X248(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let c4=c2*c2;
   let b4=b2*b2;
   let a4=a2*a2;
   /* end vars */
   let v1 = a2*(a2-b2-c2)*(a4+b4-a2*c2-b2*c2)*(a4-a2*b2-b2*c2+c4);
   let v2 = b2*(-a2+b2-c2)*(a4+b4-a2*c2-b2*c2)*(-(a2*b2)+b4-a2*c2+c4);
   let v3 = c2*(-a2-b2+c2)*(-(a2*b2)+b4-a2*c2+c4)*(a4-a2*b2-b2*c2+c4);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X249(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   /* end vars */
   let v1 = a2*(a-b)*(a-b)*(a+b)*(a+b)*(a-c)*(a-c)*(a+c)*(a+c);
   let v2 = (-a+b)*(-a+b)*(a+b)*(a+b)*b2*(b-c)*(b-c)*(b+c)*(b+c);
   let v3 = (-a+c)*(-a+c)*(a+c)*(a+c)*(-b+c)*(-b+c)*(b+c)*(b+c)*c2;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X250(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   /* end vars */
   let v1 = a2*(a-b)*(a-b)*(a+b)*(a+b)*(a-c)*(a-c)*(a+c)*(a+c)*(a2+b2-c2)*(a2-b2+c2);
   let v2 = (-a+b)*(-a+b)*(a+b)*(a+b)*b2*(b-c)*(b-c)*(b+c)*(b+c)*(a2+b2-c2)*(-a2+b2+c2);
   let v3 = (-a+c)*(-a+c)*(a+c)*(a+c)*(-b+c)*(-b+c)*(b+c)*(b+c)*c2*(a2-b2+c2)*(-a2+b2+c2);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X251(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   /* end vars */
   let v1 = a2*(a2+b2)*(a2+c2);
   let v2 = b2*(a2+b2)*(b2+c2);
   let v3 = c2*(a2+c2)*(b2+c2);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X252(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let c4=c2*c2;
   let b4=b2*b2;
   let a4=a2*a2;
   /* end vars */
   let v1 = (a4-2*a2*b2+b4-a2*c2-b2*c2)*(a4-a2*b2+b4-2*a2*c2-2*b2*c2+c4)*(a4-2*a2*b2+b4-a2*c2-2*b2*c2+c4)*(a4-a2*b2-2*a2*c2-b2*c2+c4);
   let v2 = (a4-2*a2*b2+b4-a2*c2-b2*c2)*(a4-a2*b2+b4-2*a2*c2-2*b2*c2+c4)*(-(a2*b2)+b4-a2*c2-2*b2*c2+c4)*(a4-2*a2*b2+b4-2*a2*c2-b2*c2+c4);
   let v3 = (a4-2*a2*b2+b4-a2*c2-2*b2*c2+c4)*(-(a2*b2)+b4-a2*c2-2*b2*c2+c4)*(a4-a2*b2-2*a2*c2-b2*c2+c4)*(a4-2*a2*b2+b4-2*a2*c2-b2*c2+c4);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X253(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let c4=c2*c2;
   let b4=b2*b2;
   let a4=a2*a2;
   /* end vars */
   let v1 = (a4-2*a2*b2+b4+2*a2*c2+2*b2*c2-3*c4)*(a4+2*a2*b2-3*b4-2*a2*c2+2*b2*c2+c4);
   let v2 = (a4-2*a2*b2+b4+2*a2*c2+2*b2*c2-3*c4)*(-3*a4+2*a2*b2+b4+2*a2*c2-2*b2*c2+c4);
   let v3 = (-3*a4+2*a2*b2+b4+2*a2*c2-2*b2*c2+c4)*(a4+2*a2*b2-3*b4-2*a2*c2+2*b2*c2+c4);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X254(orbit) {
   /* begin vars */
   let c2=c*c;
   let c4=c2*c2;
   let b2=b*b;
   let b4=b2*b2;
   let a2=a*a;
   let a4=a2*a2;
   let c6=c2*c4;
   let b6=b2*b4;
   let a6=a2*a4;
   /* end vars */
   let v1 = (a2+b2-c2)*(a2-b2+c2)*(a6-a4*b2-a2*b4+b6-3*a4*c2+2*a2*b2*c2-3*b4*c2+3*a2*c4+3*b2*c4-c6)*(a6-3*a4*b2+3*a2*b4-b6-a4*c2+2*a2*b2*c2+3*b4*c2-a2*c4-3*b2*c4+c6);
   let v2 = (a2+b2-c2)*(-a2+b2+c2)*(a6-a4*b2-a2*b4+b6-3*a4*c2+2*a2*b2*c2-3*b4*c2+3*a2*c4+3*b2*c4-c6)*(-a6+3*a4*b2-3*a2*b4+b6+3*a4*c2+2*a2*b2*c2-b4*c2-3*a2*c4-b2*c4+c6);
   let v3 = (a2-b2+c2)*(-a2+b2+c2)*(a6-3*a4*b2+3*a2*b4-b6-a4*c2+2*a2*b2*c2+3*b4*c2-a2*c4-3*b2*c4+c6)*(-a6+3*a4*b2-3*a2*b4+b6+3*a4*c2+2*a2*b2*c2-b4*c2-3*a2*c4-b2*c4+c6);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X255(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let c3=c2*c;
   let b3=b2*b;
   let a3=a2*a;
   /* end vars */
   let v1 = a3*(a2-b2-c2)*(a2-b2-c2);
   let v2 = b3*(-a2+b2-c2)*(-a2+b2-c2);
   let v3 = (-a2-b2+c2)*(-a2-b2+c2)*c3;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X256(orbit) {
   /* begin vars */
   let a2=a*a;
   let c2=c*c;
   let b2=b*b;
   /* end vars */
   let v1 = a*(b2+a*c)*(a*b+c2);
   let v2 = b*(a2+b*c)*(a*b+c2);
   let v3 = c*(b2+a*c)*(a2+b*c);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X257(orbit) {
   /* begin vars */
   let a2=a*a;
   let c2=c*c;
   let b2=b*b;
   /* end vars */
   let v1 = (b2+a*c)*(a*b+c2);
   let v2 = (a2+b*c)*(a*b+c2);
   let v3 = (b2+a*c)*(a2+b*c);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X258(orbit) {
   /* begin vars */
   let sc=(a+b-c)/2;
   let sb=(c+a-b)/2;
   let sa=(b+c-a)/2;
   let s=(a+b+c)/2;
   let Sqrt=Math.sqrt;
   /* end vars */
   let v1 = a/(Sqrt((s*sa)/(b*c))-Sqrt((s*sb)/(a*c))-Sqrt((s*sc)/(a*b)));
   let v2 = b/(-Sqrt((s*sa)/(b*c))+Sqrt((s*sb)/(a*c))-Sqrt((s*sc)/(a*b)));
   let v3 = c/(-Sqrt((s*sa)/(b*c))-Sqrt((s*sb)/(a*c))+Sqrt((s*sc)/(a*b)));
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X259(orbit) {
   /* begin vars */
   let c2=c*c;
   let sa=(b+c-a)/2;
   let b2=b*b;
   let sc=(a+b-c)/2;
   let sb=(c+a-b)/2;
   let Sqrt=Math.sqrt;
   let a2=a*a;
   /* end vars */
   let v1 = a2*Sqrt((b*c)/(sb*sc));
   let v2 = b2*Sqrt((a*c)/(sa*sc));
   let v3 = c2*Sqrt((a*b)/(sa*sb));
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X260(orbit) {
   /* begin vars */
   let c2=c*c;
   let sa=(b+c-a)/2;
   let b2=b*b;
   let s=(a+b+c)/2;
   let sc=(a+b-c)/2;
   let sb=(c+a-b)/2;
   let Sqrt=Math.sqrt;
   let a2=a*a;
   /* end vars */
   let v1 = (a2*Sqrt((b*c)/(sb*sc)))/(Sqrt((s*sb)/(a*c))+Sqrt((s*sc)/(a*b)));
   let v2 = (b2*Sqrt((a*c)/(sa*sc)))/(Sqrt((s*sa)/(b*c))+Sqrt((s*sc)/(a*b)));
   let v3 = (c2*Sqrt((a*b)/(sa*sb)))/(Sqrt((s*sa)/(b*c))+Sqrt((s*sb)/(a*c)));
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X261(orbit) {
   /* begin vars */

   /* end vars */
   let v1 = (a+b)*(a+b)*(a-b-c)*(a+c)*(a+c);
   let v2 = (a+b)*(a+b)*(-a+b-c)*(b+c)*(b+c);
   let v3 = (a+c)*(a+c)*(-a-b+c)*(b+c)*(b+c);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X262(orbit) {
   /* begin vars */
   let a2=a*a;
   let c2=c*c;
   let b2=b*b;
   let a4=a2*a2;
   let c4=c2*c2;
   let b4=b2*b2;
   /* end vars */
   let v1 = (a2*b2-b4+2*a2*c2+b2*c2)*(2*a2*b2+a2*c2+b2*c2-c4);
   let v2 = (-a4+a2*b2+a2*c2+2*b2*c2)*(2*a2*b2+a2*c2+b2*c2-c4);
   let v3 = (a2*b2-b4+2*a2*c2+b2*c2)*(-a4+a2*b2+a2*c2+2*b2*c2);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X263(orbit) {
   /* begin vars */
   let a2=a*a;
   let c2=c*c;
   let b2=b*b;
   let a4=a2*a2;
   let c4=c2*c2;
   let b4=b2*b2;
   /* end vars */
   let v1 = a2*(a2*b2-b4+2*a2*c2+b2*c2)*(2*a2*b2+a2*c2+b2*c2-c4);
   let v2 = b2*(-a4+a2*b2+a2*c2+2*b2*c2)*(2*a2*b2+a2*c2+b2*c2-c4);
   let v3 = c2*(a2*b2-b4+2*a2*c2+b2*c2)*(-a4+a2*b2+a2*c2+2*b2*c2);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X264(orbit) {
   /* begin vars */
   let c2=c*c;
   let a2=a*a;
   let b2=b*b;
   /* end vars */
   let v1 = b2*(-a2+b2-c2)*(a2+b2-c2)*c2;
   let v2 = a2*c2*(-a2-b2+c2)*(-a2+b2+c2);
   let v3 = a2*b2*(a2-b2-c2)*(a2-b2+c2);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X265(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   /* end vars */
   let v1 = (a2-b2-c2)*(a2-a*b+b2-c2)*(a2+a*b+b2-c2)*(a2-b2-a*c+c2)*(a2-b2+a*c+c2);
   let v2 = (-a2+b2-c2)*(a2-a*b+b2-c2)*(a2+a*b+b2-c2)*(-a2+b2-b*c+c2)*(-a2+b2+b*c+c2);
   let v3 = (-a2-b2+c2)*(a2-b2-a*c+c2)*(a2-b2+a*c+c2)*(-a2+b2-b*c+c2)*(-a2+b2+b*c+c2);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X266(orbit) {
   /* begin vars */
   let sa=(b+c-a)/2;
   let sc=(a+b-c)/2;
   let sb=(c+a-b)/2;
   let Sqrt=Math.sqrt;
   /* end vars */
   let v1 = a*Sqrt((sb*sc)/(b*c));
   let v2 = b*Sqrt((sa*sc)/(a*c));
   let v3 = c*Sqrt((sa*sb)/(a*b));
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X267(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let c3=c2*c;
   let b3=b2*b;
   let a3=a2*a;
   /* end vars */
   let v1 = a*(a3+a2*b+a*b2+b3+a2*c+a*b*c+b2*c-a*c2-b*c2-c3)*(a3+a2*b-a*b2-b3+a2*c+a*b*c-b2*c+a*c2+b*c2+c3);
   let v2 = b*(a3+a2*b+a*b2+b3+a2*c+a*b*c+b2*c-a*c2-b*c2-c3)*(-a3-a2*b+a*b2+b3-a2*c+a*b*c+b2*c+a*c2+b*c2+c3);
   let v3 = c*(a3+a2*b-a*b2-b3+a2*c+a*b*c-b2*c+a*c2+b*c2+c3)*(-a3-a2*b+a*b2+b3-a2*c+a*b*c+b2*c+a*c2+b*c2+c3);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X268(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let c3=c2*c;
   let b3=b2*b;
   let a3=a2*a;
   /* end vars */
   let v1 = a2*(a-b-c)*(a2-b2-c2)*(a3-a2*b-a*b2+b3+a2*c+2*a*b*c+b2*c-a*c2-b*c2-c3)*(a3+a2*b-a*b2-b3-a2*c+2*a*b*c-b2*c-a*c2+b*c2+c3);
   let v2 = b2*(-a+b-c)*(-a2+b2-c2)*(a3-a2*b-a*b2+b3+a2*c+2*a*b*c+b2*c-a*c2-b*c2-c3)*(-a3-a2*b+a*b2+b3-a2*c+2*a*b*c-b2*c+a*c2-b*c2+c3);
   let v3 = (-a-b+c)*c2*(-a2-b2+c2)*(-a3-a2*b+a*b2+b3-a2*c+2*a*b*c-b2*c+a*c2-b*c2+c3)*(a3+a2*b-a*b2-b3-a2*c+2*a*b*c-b2*c-a*c2+b*c2+c3);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X269(orbit) {
   /* begin vars */

   /* end vars */
   let v1 = a*(a+b-c)*(a+b-c)*(a-b+c)*(a-b+c);
   let v2 = b*(a+b-c)*(a+b-c)*(-a+b+c)*(-a+b+c);
   let v3 = c*(a-b+c)*(a-b+c)*(-a+b+c)*(-a+b+c);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X270(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   /* end vars */
   let v1 = a*(a+b)*(a+b)*(a-b-c)*(a+c)*(a+c)*(a2+b2-c2)*(a2-b2+c2);
   let v2 = b*(a+b)*(a+b)*(-a+b-c)*(b+c)*(b+c)*(a2+b2-c2)*(-a2+b2+c2);
   let v3 = c*(a+c)*(a+c)*(-a-b+c)*(b+c)*(b+c)*(a2-b2+c2)*(-a2+b2+c2);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X271(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let c3=c2*c;
   let b3=b2*b;
   let a3=a2*a;
   /* end vars */
   let v1 = a*(a-b-c)*(a2-b2-c2)*(a3-a2*b-a*b2+b3+a2*c+2*a*b*c+b2*c-a*c2-b*c2-c3)*(a3+a2*b-a*b2-b3-a2*c+2*a*b*c-b2*c-a*c2+b*c2+c3);
   let v2 = b*(-a+b-c)*(-a2+b2-c2)*(a3-a2*b-a*b2+b3+a2*c+2*a*b*c+b2*c-a*c2-b*c2-c3)*(-a3-a2*b+a*b2+b3-a2*c+2*a*b*c-b2*c+a*c2-b*c2+c3);
   let v3 = c*(-a-b+c)*(-a2-b2+c2)*(-a3-a2*b+a*b2+b3-a2*c+2*a*b*c-b2*c+a*c2-b*c2+c3)*(a3+a2*b-a*b2-b3-a2*c+2*a*b*c-b2*c-a*c2+b*c2+c3);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X272(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let c3=c2*c;
   let b3=b2*b;
   let a3=a2*a;
   /* end vars */
   let v1 = (a+b)*(a+c)*(a3+b3-a*b*c-a*c2-b*c2)*(a3-a*b2-a*b*c-b2*c+c3);
   let v2 = (a+b)*(b+c)*(a3+b3-a*b*c-a*c2-b*c2)*(-(a2*b)+b3-a2*c-a*b*c+c3);
   let v3 = (a+c)*(b+c)*(-(a2*b)+b3-a2*c-a*b*c+c3)*(a3-a*b2-a*b*c-b2*c+c3);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X273(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   /* end vars */
   let v1 = b*(-a+b-c)*(a+b-c)*c*(-a2+b2-c2)*(a2+b2-c2);
   let v2 = a*c*(-a-b+c)*(-a+b+c)*(-a2-b2+c2)*(-a2+b2+c2);
   let v3 = a*b*(a-b-c)*(a-b+c)*(a2-b2-c2)*(a2-b2+c2);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X274(orbit) {
   /* begin vars */

   /* end vars */
   let v1 = b*(a+b)*c*(a+c);
   let v2 = a*(a+b)*c*(b+c);
   let v3 = a*b*(a+c)*(b+c);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X275(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let c4=c2*c2;
   let b4=b2*b2;
   let a4=a2*a2;
   /* end vars */
   let v1 = (a2+b2-c2)*(a2-b2+c2)*(a4-2*a2*b2+b4-a2*c2-b2*c2)*(a4-a2*b2-2*a2*c2-b2*c2+c4);
   let v2 = (a2+b2-c2)*(-a2+b2+c2)*(a4-2*a2*b2+b4-a2*c2-b2*c2)*(-(a2*b2)+b4-a2*c2-2*b2*c2+c4);
   let v3 = (a2-b2+c2)*(-a2+b2+c2)*(-(a2*b2)+b4-a2*c2-2*b2*c2+c4)*(a4-a2*b2-2*a2*c2-b2*c2+c4);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X276(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let c4=c2*c2;
   let b4=b2*b2;
   let a4=a2*a2;
   /* end vars */
   let v1 = b2*(-a2+b2-c2)*(a2+b2-c2)*c2*(a4-2*a2*b2+b4-a2*c2-b2*c2)*(-a4+a2*b2+2*a2*c2+b2*c2-c4);
   let v2 = a2*c2*(-a2-b2+c2)*(-a2+b2+c2)*(-a4+2*a2*b2-b4+a2*c2+b2*c2)*(-(a2*b2)+b4-a2*c2-2*b2*c2+c4);
   let v3 = a2*b2*(a2-b2-c2)*(a2-b2+c2)*(a2*b2-b4+a2*c2+2*b2*c2-c4)*(a4-a2*b2-2*a2*c2-b2*c2+c4);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X277(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   /* end vars */
   let v1 = (a2-2*a*b+b2-2*b*c+c2)*(a2+b2-2*a*c-2*b*c+c2);
   let v2 = (a2-2*a*b+b2-2*a*c+c2)*(a2+b2-2*a*c-2*b*c+c2);
   let v3 = (a2-2*a*b+b2-2*a*c+c2)*(a2-2*a*b+b2-2*b*c+c2);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X278(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   /* end vars */
   let v1 = (a+b-c)*(a-b+c)*(a2+b2-c2)*(a2-b2+c2);
   let v2 = (a+b-c)*(-a+b+c)*(a2+b2-c2)*(-a2+b2+c2);
   let v3 = (a-b+c)*(-a+b+c)*(a2-b2+c2)*(-a2+b2+c2);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X279(orbit) {
   /* begin vars */

   /* end vars */
   let v1 = (a+b-c)*(a+b-c)*(a-b+c)*(a-b+c);
   let v2 = (a+b-c)*(a+b-c)*(-a+b+c)*(-a+b+c);
   let v3 = (a-b+c)*(a-b+c)*(-a+b+c)*(-a+b+c);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X280(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let c3=c2*c;
   let b3=b2*b;
   let a3=a2*a;
   /* end vars */
   let v1 = (a-b-c)*(a3-a2*b-a*b2+b3+a2*c+2*a*b*c+b2*c-a*c2-b*c2-c3)*(a3+a2*b-a*b2-b3-a2*c+2*a*b*c-b2*c-a*c2+b*c2+c3);
   let v2 = (-a+b-c)*(a3-a2*b-a*b2+b3+a2*c+2*a*b*c+b2*c-a*c2-b*c2-c3)*(-a3-a2*b+a*b2+b3-a2*c+2*a*b*c-b2*c+a*c2-b*c2+c3);
   let v3 = (-a-b+c)*(-a3-a2*b+a*b2+b3-a2*c+2*a*b*c-b2*c+a*c2-b*c2+c3)*(a3+a2*b-a*b2-b3-a2*c+2*a*b*c-b2*c-a*c2+b*c2+c3);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X281(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   /* end vars */
   let v1 = (a-b-c)*(a2+b2-c2)*(a2-b2+c2);
   let v2 = (-a+b-c)*(a2+b2-c2)*(-a2+b2+c2);
   let v3 = (-a-b+c)*(a2-b2+c2)*(-a2+b2+c2);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X282(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let c3=c2*c;
   let b3=b2*b;
   let a3=a2*a;
   /* end vars */
   let v1 = a*(a-b-c)*(a3-a2*b-a*b2+b3+a2*c+2*a*b*c+b2*c-a*c2-b*c2-c3)*(a3+a2*b-a*b2-b3-a2*c+2*a*b*c-b2*c-a*c2+b*c2+c3);
   let v2 = b*(-a+b-c)*(a3-a2*b-a*b2+b3+a2*c+2*a*b*c+b2*c-a*c2-b*c2-c3)*(-a3-a2*b+a*b2+b3-a2*c+2*a*b*c-b2*c+a*c2-b*c2+c3);
   let v3 = c*(-a-b+c)*(-a3-a2*b+a*b2+b3-a2*c+2*a*b*c-b2*c+a*c2-b*c2+c3)*(a3+a2*b-a*b2-b3-a2*c+2*a*b*c-b2*c-a*c2+b*c2+c3);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X283(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   /* end vars */
   let v1 = a2*(a+b)*(a-b-c)*(a+c)*(a2-b2-c2);
   let v2 = (a+b)*b2*(-a+b-c)*(b+c)*(-a2+b2-c2);
   let v3 = (a+c)*(-a-b+c)*(b+c)*c2*(-a2-b2+c2);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X284(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   /* end vars */
   let v1 = a2*(a+b)*(a-b-c)*(a+c);
   let v2 = (a+b)*b2*(-a+b-c)*(b+c);
   let v3 = (a+c)*(-a-b+c)*(b+c)*c2;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X285(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let c3=c2*c;
   let b3=b2*b;
   let a3=a2*a;
   /* end vars */
   let v1 = a*(a+b)*(a-b-c)*(a+c)*(a3-a2*b-a*b2+b3+a2*c+2*a*b*c+b2*c-a*c2-b*c2-c3)*(a3+a2*b-a*b2-b3-a2*c+2*a*b*c-b2*c-a*c2+b*c2+c3);
   let v2 = b*(a+b)*(-a+b-c)*(b+c)*(a3-a2*b-a*b2+b3+a2*c+2*a*b*c+b2*c-a*c2-b*c2-c3)*(-a3-a2*b+a*b2+b3-a2*c+2*a*b*c-b2*c+a*c2-b*c2+c3);
   let v3 = c*(a+c)*(-a-b+c)*(b+c)*(-a3-a2*b+a*b2+b3-a2*c+2*a*b*c-b2*c+a*c2-b*c2+c3)*(a3+a2*b-a*b2-b3-a2*c+2*a*b*c-b2*c-a*c2+b*c2+c3);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X286(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   /* end vars */
   let v1 = b*(a+b)*c*(a+c)*(-a2+b2-c2)*(a2+b2-c2);
   let v2 = a*(a+b)*c*(b+c)*(-a2-b2+c2)*(-a2+b2+c2);
   let v3 = a*b*(a+c)*(b+c)*(a2-b2-c2)*(a2-b2+c2);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X287(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let c4=c2*c2;
   let b4=b2*b2;
   let a4=a2*a2;
   /* end vars */
   let v1 = (a2-b2-c2)*(a4+b4-a2*c2-b2*c2)*(a4-a2*b2-b2*c2+c4);
   let v2 = (-a2+b2-c2)*(a4+b4-a2*c2-b2*c2)*(-(a2*b2)+b4-a2*c2+c4);
   let v3 = (-a2-b2+c2)*(-(a2*b2)+b4-a2*c2+c4)*(a4-a2*b2-b2*c2+c4);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X288(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let c4=c2*c2;
   let b4=b2*b2;
   let a4=a2*a2;
   /* end vars */
   let v1 = a2*(a4-2*a2*b2+b4-a2*c2-b2*c2)*(a4-3*a2*b2+2*b4-2*a2*c2-3*b2*c2+c4)*(a4-a2*b2-2*a2*c2-b2*c2+c4)*(a4-2*a2*b2+b4-3*a2*c2-3*b2*c2+2*c4);
   let v2 = b2*(a4-2*a2*b2+b4-a2*c2-b2*c2)*(2*a4-3*a2*b2+b4-3*a2*c2-2*b2*c2+c4)*(-(a2*b2)+b4-a2*c2-2*b2*c2+c4)*(a4-2*a2*b2+b4-3*a2*c2-3*b2*c2+2*c4);
   let v3 = c2*(a4-3*a2*b2+2*b4-2*a2*c2-3*b2*c2+c4)*(2*a4-3*a2*b2+b4-3*a2*c2-2*b2*c2+c4)*(-(a2*b2)+b4-a2*c2-2*b2*c2+c4)*(a4-a2*b2-2*a2*c2-b2*c2+c4);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X289(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let sc=(a+b-c)/2;
   let sb=(c+a-b)/2;
   let sa=(b+c-a)/2;
   let s=(a+b+c)/2;
   let Sqrt=Math.sqrt;
   let a2=a*a;
   /* end vars */
   let v1 = a2/(Sqrt((s*sa)/(b*c))*(Sqrt((s*sa)/(b*c))-Sqrt((s*sb)/(a*c))-Sqrt((s*sc)/(a*b))));
   let v2 = b2/(Sqrt((s*sb)/(a*c))*(-Sqrt((s*sa)/(b*c))+Sqrt((s*sb)/(a*c))-Sqrt((s*sc)/(a*b))));
   let v3 = c2/(Sqrt((s*sc)/(a*b))*(-Sqrt((s*sa)/(b*c))-Sqrt((s*sb)/(a*c))+Sqrt((s*sc)/(a*b))));
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X290(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let c4=c2*c2;
   let b4=b2*b2;
   let a4=a2*a2;
   /* end vars */
   let v1 = b2*c2*(a4+b4-a2*c2-b2*c2)*(-a4+a2*b2+b2*c2-c4);
   let v2 = a2*c2*(-a4-b4+a2*c2+b2*c2)*(-(a2*b2)+b4-a2*c2+c4);
   let v3 = a2*b2*(a2*b2-b4+a2*c2-c4)*(a4-a2*b2-b2*c2+c4);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X291(orbit) {
   /* begin vars */
   let a2=a*a;
   let c2=c*c;
   let b2=b*b;
   /* end vars */
   let v1 = a*(-b2+a*c)*(a*b-c2);
   let v2 = b*(-a2+b*c)*(a*b-c2);
   let v3 = c*(-b2+a*c)*(-a2+b*c);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X292(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   /* end vars */
   let v1 = a2*(-b2+a*c)*(a*b-c2);
   let v2 = b2*(-a2+b*c)*(a*b-c2);
   let v3 = (-b2+a*c)*(-a2+b*c)*c2;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X293(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let c4=c2*c2;
   let b4=b2*b2;
   let a4=a2*a2;
   /* end vars */
   let v1 = a*(a2-b2-c2)*(a4+b4-a2*c2-b2*c2)*(a4-a2*b2-b2*c2+c4);
   let v2 = b*(-a2+b2-c2)*(a4+b4-a2*c2-b2*c2)*(-(a2*b2)+b4-a2*c2+c4);
   let v3 = c*(-a2-b2+c2)*(-(a2*b2)+b4-a2*c2+c4)*(a4-a2*b2-b2*c2+c4);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X294(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   /* end vars */
   let v1 = a*(a-b-c)*(a2+b2-a*c-b*c)*(a2-a*b-b*c+c2);
   let v2 = b*(-a+b-c)*(a2+b2-a*c-b*c)*(-(a*b)+b2-a*c+c2);
   let v3 = c*(-a-b+c)*(-(a*b)+b2-a*c+c2)*(a2-a*b-b*c+c2);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X295(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   /* end vars */
   let v1 = a2*(-b2+a*c)*(a*b-c2)*(a2-b2-c2);
   let v2 = b2*(-a2+b*c)*(a*b-c2)*(-a2+b2-c2);
   let v3 = (-b2+a*c)*(-a2+b*c)*c2*(-a2-b2+c2);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X296(orbit) {
   /* begin vars */
   let a2=a*a;
   let c2=c*c;
   let b2=b*b;
   let a4=a2*a2;
   let c4=c2*c2;
   let b3=b2*b;
   let c3=c2*c;
   let a3=a2*a;
   let b4=b2*b2;
   /* end vars */
   let v1 = a2*(a+b-c)*(a-b+c)*(a2-b2-c2)*(a2*b2-b4+a3*c-a*b2*c-2*a2*c2+b2*c2+a*c3)*(a3*b-2*a2*b2+a*b3+a2*c2-a*b*c2+b2*c2-c4);
   let v2 = b2*(a+b-c)*(-a+b+c)*(-a2+b2-c2)*(-a4+a2*b2-a2*b*c+b3*c+a2*c2-2*b2*c2+b*c3)*(a3*b-2*a2*b2+a*b3+a2*c2-a*b*c2+b2*c2-c4);
   let v3 = (a-b+c)*(-a+b+c)*c2*(-a2-b2+c2)*(a2*b2-b4+a3*c-a*b2*c-2*a2*c2+b2*c2+a*c3)*(-a4+a2*b2-a2*b*c+b3*c+a2*c2-2*b2*c2+b*c3);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X297(orbit) {
   /* begin vars */
   let a2=a*a;
   let c2=c*c;
   let b2=b*b;
   let a4=a2*a2;
   let c4=c2*c2;
   let b4=b2*b2;
   /* end vars */
   let v1 = (a2+b2-c2)*(a2-b2+c2)*(a2*b2-b4+a2*c2-c4);
   let v2 = (a2+b2-c2)*(-a2+b2+c2)*(-a4+a2*b2+b2*c2-c4);
   let v3 = (a2-b2+c2)*(-a2+b2+c2)*(-a4-b4+a2*c2+b2*c2);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X298(orbit) {
   /* begin vars */
   let area=triAreaHeron(a,b,c);
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let S=2*area;
   let c4=c2*c2;
   let b4=b2*b2;
   let a4=a2*a2;
   /* end vars */
   let v1 = (2*a4-a2*b2-b4-4*a2*c2-b2*c2+2*c4-2*b2*S*sqrt3)*(2*a4-4*a2*b2+2*b4-a2*c2-b2*c2-c4-2*c2*S*sqrt3);
   let v2 = (-a4-a2*b2+2*b4-a2*c2-4*b2*c2+2*c4-2*a2*S*sqrt3)*(2*a4-4*a2*b2+2*b4-a2*c2-b2*c2-c4-2*c2*S*sqrt3);
   let v3 = (-a4-a2*b2+2*b4-a2*c2-4*b2*c2+2*c4-2*a2*S*sqrt3)*(2*a4-a2*b2-b4-4*a2*c2-b2*c2+2*c4-2*b2*S*sqrt3);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X299(orbit) {
   /* begin vars */
   let area=triAreaHeron(a,b,c);
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let S=2*area;
   let c4=c2*c2;
   let b4=b2*b2;
   let a4=a2*a2;
   /* end vars */
   let v1 = (2*a4-a2*b2-b4-4*a2*c2-b2*c2+2*c4+2*b2*S*sqrt3)*(2*a4-4*a2*b2+2*b4-a2*c2-b2*c2-c4+2*c2*S*sqrt3);
   let v2 = (-a4-a2*b2+2*b4-a2*c2-4*b2*c2+2*c4+2*a2*S*sqrt3)*(2*a4-4*a2*b2+2*b4-a2*c2-b2*c2-c4+2*c2*S*sqrt3);
   let v3 = (-a4-a2*b2+2*b4-a2*c2-4*b2*c2+2*c4+2*a2*S*sqrt3)*(2*a4-a2*b2-b4-4*a2*c2-b2*c2+2*c4+2*b2*S*sqrt3);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X300(orbit) {
   /* begin vars */
   let area=triAreaHeron(a,b,c);
   let a2=a*a;
   let S=2*area;
   let c2=c*c;
   let b2=b*b;
   /* end vars */
   let v1 = b2*c2*(-2*S+(-a2+b2-c2)*sqrt3)*(-2*S+(-a2-b2+c2)*sqrt3);
   let v2 = a2*c2*(-2*S+(a2-b2-c2)*sqrt3)*(-2*S+(-a2-b2+c2)*sqrt3);
   let v3 = a2*b2*(-2*S+(a2-b2-c2)*sqrt3)*(-2*S+(-a2+b2-c2)*sqrt3);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X301(orbit) {
   /* begin vars */
   let area=triAreaHeron(a,b,c);
   let a2=a*a;
   let S=2*area;
   let c2=c*c;
   let b2=b*b;
   /* end vars */
   let v1 = b2*c2*(2*S+(-a2+b2-c2)*sqrt3)*(2*S+(-a2-b2+c2)*sqrt3);
   let v2 = a2*c2*(2*S+(a2-b2-c2)*sqrt3)*(2*S+(-a2-b2+c2)*sqrt3);
   let v3 = a2*b2*(2*S+(a2-b2-c2)*sqrt3)*(2*S+(-a2+b2-c2)*sqrt3);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X302(orbit) {
   /* begin vars */
   let area=triAreaHeron(a,b,c);
   let S=2*area;
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   /* end vars */
   let v1 = -a2+b2+c2+2*S*sqrt3;
   let v2 = a2-b2+c2+2*S*sqrt3;
   let v3 = a2+b2-c2+2*S*sqrt3;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X303(orbit) {
   /* begin vars */
   let area=triAreaHeron(a,b,c);
   let S=2*area;
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   /* end vars */
   let v1 = -a2+b2+c2-2*S*sqrt3;
   let v2 = a2-b2+c2-2*S*sqrt3;
   let v3 = a2+b2-c2-2*S*sqrt3;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X304(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   /* end vars */
   let v1 = b*c*(-a2+b2+c2);
   let v2 = a*c*(a2-b2+c2);
   let v3 = a*b*(a2+b2-c2);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X305(orbit) {
   /* begin vars */
   let a2=a*a;
   let c2=c*c;
   let b2=b*b;
   /* end vars */
   let v1 = b2*c2*(-a2+b2+c2);
   let v2 = a2*c2*(a2-b2+c2);
   let v3 = a2*b2*(a2+b2-c2);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X306(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   /* end vars */
   let v1 = (b+c)*(-a2+b2+c2);
   let v2 = (a+c)*(a2-b2+c2);
   let v3 = (a+b)*(a2+b2-c2);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X307(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   /* end vars */
   let v1 = (a+b-c)*(a-b+c)*(b+c)*(a2-b2-c2);
   let v2 = (a+b-c)*(a+c)*(-a+b+c)*(-a2+b2-c2);
   let v3 = (a+b)*(a-b+c)*(-a+b+c)*(-a2-b2+c2);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X308(orbit) {
   /* begin vars */
   let c2=c*c;
   let a2=a*a;
   let b2=b*b;
   /* end vars */
   let v1 = b2*(a2+b2)*c2*(a2+c2);
   let v2 = a2*(a2+b2)*c2*(b2+c2);
   let v3 = a2*b2*(a2+c2)*(b2+c2);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X309(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let c3=c2*c;
   let b3=b2*b;
   let a3=a2*a;
   /* end vars */
   let v1 = b*c*(a3-a2*b-a*b2+b3+a2*c+2*a*b*c+b2*c-a*c2-b*c2-c3)*(-a3-a2*b+a*b2+b3+a2*c-2*a*b*c+b2*c+a*c2-b*c2-c3);
   let v2 = a*c*(-a3-a2*b+a*b2+b3-a2*c+2*a*b*c-b2*c+a*c2-b*c2+c3)*(-a3+a2*b+a*b2-b3-a2*c-2*a*b*c-b2*c+a*c2+b*c2+c3);
   let v3 = a*b*(a3+a2*b-a*b2-b3+a2*c-2*a*b*c+b2*c-a*c2+b*c2-c3)*(a3+a2*b-a*b2-b3-a2*c+2*a*b*c-b2*c-a*c2+b*c2+c3);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X310(orbit) {
   /* begin vars */
   let a2=a*a;
   let c2=c*c;
   let b2=b*b;
   /* end vars */
   let v1 = (a+b)*b2*(a+c)*c2;
   let v2 = a2*(a+b)*(b+c)*c2;
   let v3 = a2*b2*(a+c)*(b+c);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X311(orbit) {
   /* begin vars */
   let a2=a*a;
   let c2=c*c;
   let b2=b*b;
   let a4=a2*a2;
   let c4=c2*c2;
   let b4=b2*b2;
   /* end vars */
   let v1 = b2*c2*(-(a2*b2)+b4-a2*c2-2*b2*c2+c4);
   let v2 = a2*c2*(a4-a2*b2-2*a2*c2-b2*c2+c4);
   let v3 = a2*b2*(a4-2*a2*b2+b4-a2*c2-b2*c2);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X312(orbit) {
   /* begin vars */

   /* end vars */
   let v1 = b*c*(-a+b+c);
   let v2 = a*c*(a-b+c);
   let v3 = a*b*(a+b-c);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X313(orbit) {
   /* begin vars */
   let a2=a*a;
   let c2=c*c;
   let b2=b*b;
   /* end vars */
   let v1 = b2*(b+c)*c2;
   let v2 = a2*(a+c)*c2;
   let v3 = a2*(a+b)*b2;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X314(orbit) {
   /* begin vars */

   /* end vars */
   let v1 = b*(a+b)*c*(a+c)*(-a+b+c);
   let v2 = a*(a+b)*c*(a-b+c)*(b+c);
   let v3 = a*b*(a+b-c)*(a+c)*(b+c);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X315(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let c4=c2*c2;
   let b4=b2*b2;
   let a4=a2*a2;
   /* end vars */
   let v1 = -a4+b4+c4;
   let v2 = a4-b4+c4;
   let v3 = a4+b4-c4;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X316(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let c4=c2*c2;
   let b4=b2*b2;
   let a4=a2*a2;
   /* end vars */
   let v1 = -a4+b4-b2*c2+c4;
   let v2 = a4-b4-a2*c2+c4;
   let v3 = a4-a2*b2+b4-c4;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X317(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let c4=c2*c2;
   let b4=b2*b2;
   let a4=a2*a2;
   /* end vars */
   let v1 = (a2+b2-c2)*(a2-b2+c2)*(a4-2*a2*b2+b4-2*a2*c2+c4);
   let v2 = (a2+b2-c2)*(-a2+b2+c2)*(a4-2*a2*b2+b4-2*b2*c2+c4);
   let v3 = (a2-b2+c2)*(-a2+b2+c2)*(a4+b4-2*a2*c2-2*b2*c2+c4);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X318(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   /* end vars */
   let v1 = b*c*(-a+b+c)*(-a2+b2-c2)*(a2+b2-c2);
   let v2 = a*c*(a-b+c)*(-a2-b2+c2)*(-a2+b2+c2);
   let v3 = a*b*(a+b-c)*(a2-b2-c2)*(a2-b2+c2);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X319(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   /* end vars */
   let v1 = -a2+b2+b*c+c2;
   let v2 = a2-b2+a*c+c2;
   let v3 = a2+a*b+b2-c2;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X320(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   /* end vars */
   let v1 = -a2+b2-b*c+c2;
   let v2 = a2-b2-a*c+c2;
   let v3 = a2-a*b+b2-c2;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X321(orbit) {
   /* begin vars */

   /* end vars */
   let v1 = b*c*(b+c);
   let v2 = a*c*(a+c);
   let v3 = a*b*(a+b);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X322(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let c3=c2*c;
   let b3=b2*b;
   let a3=a2*a;
   /* end vars */
   let v1 = b*c*(-a3-a2*b+a*b2+b3-a2*c+2*a*b*c-b2*c+a*c2-b*c2+c3);
   let v2 = a*c*(a3+a2*b-a*b2-b3-a2*c+2*a*b*c-b2*c-a*c2+b*c2+c3);
   let v3 = a*b*(a3-a2*b-a*b2+b3+a2*c+2*a*b*c+b2*c-a*c2-b*c2-c3);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X323(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   /* end vars */
   let v1 = a2*(a2-b2-b*c-c2)*(a2-b2+b*c-c2);
   let v2 = b2*(-a2+b2-a*c-c2)*(-a2+b2+a*c-c2);
   let v3 = c2*(-a2-a*b-b2+c2)*(-a2+a*b-b2+c2);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X324(orbit) {
   /* begin vars */
   let a2=a*a;
   let c2=c*c;
   let b2=b*b;
   let a4=a2*a2;
   let c4=c2*c2;
   let b4=b2*b2;
   /* end vars */
   let v1 = b2*(-a2+b2-c2)*(a2+b2-c2)*c2*(-(a2*b2)+b4-a2*c2-2*b2*c2+c4);
   let v2 = a2*c2*(-a2-b2+c2)*(-a2+b2+c2)*(a4-a2*b2-2*a2*c2-b2*c2+c4);
   let v3 = a2*b2*(a2-b2-c2)*(a2-b2+c2)*(a4-2*a2*b2+b4-a2*c2-b2*c2);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X325(orbit) {
   /* begin vars */
   let a2=a*a;
   let c2=c*c;
   let b2=b*b;
   let a4=a2*a2;
   let c4=c2*c2;
   let b4=b2*b2;
   /* end vars */
   let v1 = -(a2*b2)+b4-a2*c2+c4;
   let v2 = a4-a2*b2-b2*c2+c4;
   let v3 = a4+b4-a2*c2-b2*c2;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X326(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   /* end vars */
   let v1 = a*(a2-b2-c2)*(a2-b2-c2);
   let v2 = b*(-a2+b2-c2)*(-a2+b2-c2);
   let v3 = c*(-a2-b2+c2)*(-a2-b2+c2);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X327(orbit) {
   /* begin vars */
   let a2=a*a;
   let c2=c*c;
   let b2=b*b;
   let a4=a2*a2;
   let c4=c2*c2;
   let b4=b2*b2;
   /* end vars */
   let v1 = b2*c2*(-(a2*b2)+b4-2*a2*c2-b2*c2)*(2*a2*b2+a2*c2+b2*c2-c4);
   let v2 = a2*c2*(-a4+a2*b2+a2*c2+2*b2*c2)*(-2*a2*b2-a2*c2-b2*c2+c4);
   let v3 = a2*b2*(a4-a2*b2-a2*c2-2*b2*c2)*(a2*b2-b4+2*a2*c2+b2*c2);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X328(orbit) {
   /* begin vars */
   let c2=c*c;
   let a2=a*a;
   let b2=b*b;
   /* end vars */
   let v1 = b2*(a2-a*b+b2-c2)*(a2+a*b+b2-c2)*(-a2+b2-a*c-c2)*(-a2+b2+a*c-c2)*c2*(-a2+b2+c2);
   let v2 = a2*c2*(a2-b2+c2)*(-a2-a*b-b2+c2)*(-a2+a*b-b2+c2)*(-a2+b2-b*c+c2)*(-a2+b2+b*c+c2);
   let v3 = a2*b2*(a2+b2-c2)*(a2-b2-b*c-c2)*(a2-b2+b*c-c2)*(a2-b2-a*c+c2)*(a2-b2+a*c+c2);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X329(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let c3=c2*c;
   let b3=b2*b;
   let a3=a2*a;
   /* end vars */
   let v1 = a3+a2*b-a*b2-b3+a2*c-2*a*b*c+b2*c-a*c2+b*c2-c3;
   let v2 = -a3-a2*b+a*b2+b3+a2*c-2*a*b*c+b2*c+a*c2-b*c2-c3;
   let v3 = -a3+a2*b+a*b2-b3-a2*c-2*a*b*c-b2*c+a*c2+b*c2+c3;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X330(orbit) {
   /* begin vars */

   /* end vars */
   let v1 = (a*b-a*c-b*c)*(a*b-a*c+b*c);
   let v2 = (-(a*b)-a*c+b*c)*(-(a*b)+a*c+b*c);
   let v3 = (-(a*b)+a*c-b*c)*(a*b+a*c-b*c);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X331(orbit) {
   /* begin vars */
   let c2=c*c;
   let a2=a*a;
   let b2=b*b;
   /* end vars */
   let v1 = b2*(-a+b-c)*(a+b-c)*(-a2+b2-c2)*(a2+b2-c2)*c2;
   let v2 = a2*(-a-b+c)*(-a+b+c)*c2*(-a2-b2+c2)*(-a2+b2+c2);
   let v3 = a2*b2*(a-b-c)*(a-b+c)*(a2-b2-c2)*(a2-b2+c2);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X332(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   /* end vars */
   let v1 = (a+b)*(a-b-c)*(a+c)*(a2-b2-c2);
   let v2 = (a+b)*(-a+b-c)*(b+c)*(-a2+b2-c2);
   let v3 = (a+c)*(-a-b+c)*(b+c)*(-a2-b2+c2);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X333(orbit) {
   /* begin vars */

   /* end vars */
   let v1 = (a+b)*(a-b-c)*(a+c);
   let v2 = (a+b)*(-a+b-c)*(b+c);
   let v3 = (a+c)*(-a-b+c)*(b+c);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X334(orbit) {
   /* begin vars */
   let a2=a*a;
   let c2=c*c;
   let b2=b*b;
   /* end vars */
   let v1 = b*c*(b2-a*c)*(a*b-c2);
   let v2 = a*c*(-a2+b*c)*(-(a*b)+c2);
   let v3 = a*b*(-b2+a*c)*(a2-b*c);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X335(orbit) {
   /* begin vars */
   let a2=a*a;
   let c2=c*c;
   let b2=b*b;
   /* end vars */
   let v1 = (b2-a*c)*(a*b-c2);
   let v2 = (-a2+b*c)*(-(a*b)+c2);
   let v3 = (-b2+a*c)*(a2-b*c);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X336(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let c4=c2*c2;
   let b4=b2*b2;
   let a4=a2*a2;
   /* end vars */
   let v1 = b*c*(-a2+b2+c2)*(a4+b4-a2*c2-b2*c2)*(-a4+a2*b2+b2*c2-c4);
   let v2 = a*c*(a2-b2+c2)*(-a4-b4+a2*c2+b2*c2)*(-(a2*b2)+b4-a2*c2+c4);
   let v3 = a*b*(a2+b2-c2)*(a2*b2-b4+a2*c2-c4)*(a4-a2*b2-b2*c2+c4);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X337(orbit) {
   /* begin vars */
   let a2=a*a;
   let c2=c*c;
   let b2=b*b;
   /* end vars */
   let v1 = (b2-a*c)*(a*b-c2)*(-a2+b2+c2);
   let v2 = (-a2+b*c)*(-(a*b)+c2)*(a2-b2+c2);
   let v3 = (-b2+a*c)*(a2-b*c)*(a2+b2-c2);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X338(orbit) {
   /* begin vars */
   let a2=a*a;
   let c2=c*c;
   let b2=b*b;
   /* end vars */
   let v1 = b2*(b-c)*(b-c)*(b+c)*(b+c)*c2;
   let v2 = a2*(-a+c)*(-a+c)*(a+c)*(a+c)*c2;
   let v3 = a2*(a-b)*(a-b)*(a+b)*(a+b)*b2;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X339(orbit) {
   /* begin vars */
   let a2=a*a;
   let c2=c*c;
   let b2=b*b;
   /* end vars */
   let v1 = b2*(b-c)*(b-c)*(b+c)*(b+c)*c2*(-a2+b2+c2);
   let v2 = a2*(-a+c)*(-a+c)*(a+c)*(a+c)*c2*(a2-b2+c2);
   let v3 = a2*(a-b)*(a-b)*(a+b)*(a+b)*b2*(a2+b2-c2);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X340(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   /* end vars */
   let v1 = (a2+b2-c2)*(a2-b2-b*c-c2)*(a2-b2+b*c-c2)*(a2-b2+c2);
   let v2 = (a2+b2-c2)*(-a2+b2-a*c-c2)*(-a2+b2+a*c-c2)*(-a2+b2+c2);
   let v3 = (a2-b2+c2)*(-a2-a*b-b2+c2)*(-a2+a*b-b2+c2)*(-a2+b2+c2);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X341(orbit) {
   /* begin vars */

   /* end vars */
   let v1 = b*c*(-a+b+c)*(-a+b+c);
   let v2 = a*c*(a-b+c)*(a-b+c);
   let v3 = a*b*(a+b-c)*(a+b-c);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X342(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let c3=c2*c;
   let b3=b2*b;
   let a3=a2*a;
   /* end vars */
   let v1 = b*(-a+b-c)*(a+b-c)*c*(-a2+b2-c2)*(a2+b2-c2)*(-a3-a2*b+a*b2+b3-a2*c+2*a*b*c-b2*c+a*c2-b*c2+c3);
   let v2 = a*c*(-a-b+c)*(-a+b+c)*(-a2-b2+c2)*(-a2+b2+c2)*(a3+a2*b-a*b2-b3-a2*c+2*a*b*c-b2*c-a*c2+b*c2+c3);
   let v3 = a*b*(a-b-c)*(a-b+c)*(a2-b2-c2)*(a2-b2+c2)*(a3-a2*b-a*b2+b3+a2*c+2*a*b*c+b2*c-a*c2-b*c2-c3);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X343(orbit) {
   /* begin vars */
   let a2=a*a;
   let c2=c*c;
   let b2=b*b;
   let a4=a2*a2;
   let c4=c2*c2;
   let b4=b2*b2;
   /* end vars */
   let v1 = (a2-b2-c2)*(a2*b2-b4+a2*c2+2*b2*c2-c4);
   let v2 = (-a2+b2-c2)*(-a4+a2*b2+2*a2*c2+b2*c2-c4);
   let v3 = (-a2-b2+c2)*(-a4+2*a2*b2-b4+a2*c2+b2*c2);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X344(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   /* end vars */
   let v1 = a2-2*a*b+b2-2*a*c+c2;
   let v2 = a2-2*a*b+b2-2*b*c+c2;
   let v3 = a2+b2-2*a*c-2*b*c+c2;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X345(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   /* end vars */
   let v1 = (a-b-c)*(a2-b2-c2);
   let v2 = (-a+b-c)*(-a2+b2-c2);
   let v3 = (-a-b+c)*(-a2-b2+c2);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X346(orbit) {
   /* begin vars */

   /* end vars */
   let v1 = (a-b-c)*(a-b-c);
   let v2 = (-a+b-c)*(-a+b-c);
   let v3 = (-a-b+c)*(-a-b+c);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X347(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let c3=c2*c;
   let b3=b2*b;
   let a3=a2*a;
   /* end vars */
   let v1 = (a+b-c)*(a-b+c)*(a3+a2*b-a*b2-b3+a2*c-2*a*b*c+b2*c-a*c2+b*c2-c3);
   let v2 = (a+b-c)*(-a+b+c)*(-a3-a2*b+a*b2+b3+a2*c-2*a*b*c+b2*c+a*c2-b*c2-c3);
   let v3 = (a-b+c)*(-a+b+c)*(-a3+a2*b+a*b2-b3-a2*c-2*a*b*c-b2*c+a*c2+b*c2+c3);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X348(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   /* end vars */
   let v1 = (a+b-c)*(a-b+c)*(a2-b2-c2);
   let v2 = (a+b-c)*(-a+b+c)*(-a2+b2-c2);
   let v3 = (a-b+c)*(-a+b+c)*(-a2-b2+c2);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X349(orbit) {
   /* begin vars */
   let a2=a*a;
   let c2=c*c;
   let b2=b*b;
   /* end vars */
   let v1 = b2*(-a+b-c)*(a+b-c)*(b+c)*c2;
   let v2 = a2*(a+c)*(-a-b+c)*(-a+b+c)*c2;
   let v3 = a2*(a+b)*b2*(a-b-c)*(a-b+c);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X350(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   /* end vars */
   let v1 = b*c*(-a2+b*c);
   let v2 = a*c*(-b2+a*c);
   let v3 = a*b*(a*b-c2);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X351(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   /* end vars */
   let v1 = a2*(b2-c2)*(-2*a2+b2+c2);
   let v2 = b2*(-a2+c2)*(a2-2*b2+c2);
   let v3 = (a2-b2)*(a2+b2-2*c2)*c2;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X352(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let c4=c2*c2;
   let b4=b2*b2;
   let a4=a2*a2;
   /* end vars */
   let v1 = a2*(a4-4*a2*b2+b4-4*a2*c2+5*b2*c2+c4);
   let v2 = b2*(a4-4*a2*b2+b4+5*a2*c2-4*b2*c2+c4);
   let v3 = c2*(a4+5*a2*b2+b4-4*a2*c2-4*b2*c2+c4);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X353(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let c4=c2*c2;
   let b4=b2*b2;
   let a4=a2*a2;
   /* end vars */
   let v1 = a2*(4*a4-4*a2*b2-2*b4-4*a2*c2-b2*c2-2*c4);
   let v2 = b2*(-2*a4-4*a2*b2+4*b4-a2*c2-4*b2*c2-2*c4);
   let v3 = c2*(-2*a4-a2*b2-2*b4-4*a2*c2-4*b2*c2+4*c4);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X354(orbit) {
   /* begin vars */
   let a2=a*a;
   let c2=c*c;
   let b2=b*b;
   /* end vars */
   let v1 = a*(a*b-b2+a*c+2*b*c-c2);
   let v2 = b*(-a2+a*b+2*a*c+b*c-c2);
   let v3 = c*(-a2+2*a*b-b2+a*c+b*c);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X355(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let c4=c2*c2;
   let c3=c2*c;
   let b4=b2*b2;
   let b3=b2*b;
   let a3=a2*a;
   let a4=a2*a2;
   /* end vars */
   let v1 = a4-a3*b+a*b3-b4-a3*c+2*a2*b*c-a*b2*c-a*b*c2+2*b2*c2+a*c3-c4;
   let v2 = -a4+a3*b-a*b3+b4-a2*b*c+2*a*b2*c-b3*c+2*a2*c2-a*b*c2+b*c3-c4;
   let v3 = -a4+2*a2*b2-b4+a3*c-a2*b*c-a*b2*c+b3*c+2*a*b*c2-a*c3-b*c3+c4;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X356(orbit) {
   /* begin vars */
   let cosC=lawOfCosines(c,a,b);
   let cosB=lawOfCosines(b,a,c);
   let cosA=lawOfCosines(a,b,c);
   let sinC=getSin(cosC);
   let sinB=getSin(cosB);
   let sinA=getSin(cosA);
   let cosThirdC=cosThirdAngle(cosC);
   let cosThirdB=cosThirdAngle(cosB);
   let cosThirdA=cosThirdAngle(cosA);
   /* end vars */
   let v1 = (cosThirdA+2*cosThirdB*cosThirdC)*sinA;
   let v2 = (cosThirdB+2*cosThirdA*cosThirdC)*sinB;
   let v3 = (2*cosThirdA*cosThirdB+cosThirdC)*sinC;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X357(orbit) {
   /* begin vars */
   let cosC=lawOfCosines(c,a,b);
   let cosThirdC=cosThirdAngle(cosC);
   let cosB=lawOfCosines(b,a,c);
   let cosThirdB=cosThirdAngle(cosB);
   let cosA=lawOfCosines(a,b,c);
   let cosThirdA=cosThirdAngle(cosA);
   let sinC=getSin(cosC);
   let secThirdC=1/cosThirdC;
   let sinB=getSin(cosB);
   let secThirdB=1/cosThirdB;
   let sinA=getSin(cosA);
   let secThirdA=1/cosThirdA;
   /* end vars */
   let v1 = secThirdA*sinA;
   let v2 = secThirdB*sinB;
   let v3 = secThirdC*sinC;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X358(orbit) {
   /* begin vars */
   let cosC=lawOfCosines(c,a,b);
   let cosB=lawOfCosines(b,a,c);
   let cosA=lawOfCosines(a,b,c);
   let sinC=getSin(cosC);
   let cosThirdC=cosThirdAngle(cosC);
   let sinB=getSin(cosB);
   let cosThirdB=cosThirdAngle(cosB);
   let sinA=getSin(cosA);
   let cosThirdA=cosThirdAngle(cosA);
   /* end vars */
   let v1 = cosThirdA*sinA;
   let v2 = cosThirdB*sinB;
   let v3 = cosThirdC*sinC;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X359(orbit) {
   /* begin vars */
   let cosC=lawOfCosines(c,a,b);
   let cosB=lawOfCosines(b,a,c);
   let cosA=lawOfCosines(a,b,c);
   let C=Math.acos(cosC);
   let c2=c*c;
   let B=Math.acos(cosB);
   let b2=b*b;
   let A=Math.acos(cosA);
   let a2=a*a;
   /* end vars */
   let v1 = a2/A;
   let v2 = b2/B;
   let v3 = c2/C;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X360(orbit) {
   /* begin vars */
   let cosC=lawOfCosines(c,a,b);
   let cosB=lawOfCosines(b,a,c);
   let cosA=lawOfCosines(a,b,c);
   let C=Math.acos(cosC);
   let B=Math.acos(cosB);
   let A=Math.acos(cosA);
   /* end vars */
   let v1 = A;
   let v2 = B;
   let v3 = C;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X361(orbit) {
   /* begin vars */
   let sc=(a+b-c)/2;
   let sb=(c+a-b)/2;
   let sa=(b+c-a)/2;
   let Sqrt=Math.sqrt;
   /* end vars */
   let v1 = a*(Sqrt((a*b)/(sa*sb))+Sqrt((a*c)/(sa*sc))-Sqrt((b*c)/(sb*sc)));
   let v2 = b*(Sqrt((a*b)/(sa*sb))-Sqrt((a*c)/(sa*sc))+Sqrt((b*c)/(sb*sc)));
   let v3 = c*(-Sqrt((a*b)/(sa*sb))+Sqrt((a*c)/(sa*sc))+Sqrt((b*c)/(sb*sc)));
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X362(orbit) {
   /* begin vars */
   let sc=(a+b-c)/2;
   let sb=(c+a-b)/2;
   let sa=(b+c-a)/2;
   let s=(a+b+c)/2;
   let Sqrt=Math.sqrt;
   /* end vars */
   let v1 = a*(-(a*Sqrt((s*sa)/(b*c)))+b*Sqrt((s*sb)/(a*c))+c*Sqrt((s*sc)/(a*b)));
   let v2 = b*(a*Sqrt((s*sa)/(b*c))-b*Sqrt((s*sb)/(a*c))+c*Sqrt((s*sc)/(a*b)));
   let v3 = c*(a*Sqrt((s*sa)/(b*c))+b*Sqrt((s*sb)/(a*c))-c*Sqrt((s*sc)/(a*b)));
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X363(orbit) {
   /* begin vars */
   let sc=(a+b-c)/2;
   let sb=(c+a-b)/2;
   let sa=(b+c-a)/2;
   let Sqrt=Math.sqrt;
   /* end vars */
   let v1 = a*(c/(1+Sqrt((sa*sb)/(a*b)))+b/(1+Sqrt((sa*sc)/(a*c)))-a/(1+Sqrt((sb*sc)/(b*c))));
   let v2 = b*(c/(1+Sqrt((sa*sb)/(a*b)))-b/(1+Sqrt((sa*sc)/(a*c)))+a/(1+Sqrt((sb*sc)/(b*c))));
   let v3 = c*(-(c/(1+Sqrt((sa*sb)/(a*b))))+b/(1+Sqrt((sa*sc)/(a*c)))+a/(1+Sqrt((sb*sc)/(b*c))));
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X364(orbit) {
   /* begin vars */
   let sqrtc=sqrt(c);
   let sqrtb=sqrt(b);
   let sqrta=sqrt(a);
   /* end vars */
   let v1 = a*(-sqrta+sqrtb+sqrtc);
   let v2 = b*(sqrta-sqrtb+sqrtc);
   let v3 = c*(sqrta+sqrtb-sqrtc);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X365(orbit) {
   /* begin vars */

   /* end vars */
   let v1 = Math.pow(a,1.5);
   let v2 = Math.pow(b,1.5);
   let v3 = Math.pow(c,1.5);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X366(orbit) {
   /* begin vars */
   let sqrtc=sqrt(c);
   let sqrtb=sqrt(b);
   let sqrta=sqrt(a);
   /* end vars */
   let v1 = sqrta;
   let v2 = sqrtb;
   let v3 = sqrtc;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X367(orbit) {
   /* begin vars */
   let sqrta=sqrt(a);
   let sqrtc=sqrt(c);
   let sqrtb=sqrt(b);
   /* end vars */
   let v1 = a*(sqrtb+sqrtc);
   let v2 = b*(sqrta+sqrtc);
   let v3 = c*(sqrta+sqrtb);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X368(orbit) {
   /* begin vars */
   let a2=a*a;
   let c2=c*c;
   let b2=b*b;
   let area=triAreaHeron(a,b,c);
   let S=2*area;
   let SA=(b2+c2-a2)/2;
   let SC=(a2+b2-c2)/2;
   let SB=(c2+a2-b2)/2;
   let S2=S*S;
   /* end vars */
   let v1 = S2-3*SB*SC;
   let v2 = S2-3*SA*SC;
   let v3 = S2-3*SA*SB;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X369(orbit) {
   /* begin vars */
   let Sqrt=Math.sqrt;
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   /* end vars */
   let v1 = a2-2*a*b-b2-3*a*c-2*b*c-2*c2+(a+2*c)*((a+b+c)/2.+(a2+b2-10*(a*b+a*c+b*c)+c2)/(Math.pow(2,0.6666666666666666)*Math.pow(216*(a-b)*(b-c)*(-a+c)+Sqrt(46656*(a-b)*(a-b)*(b-c)*(b-c)*(-a+c)*(-a+c)-108*Math.pow(a2+b2-10*(a*b+a*c+b*c)+c2,3)),0.3333333333333333))+Math.pow(216*(a-b)*(b-c)*(-a+c)+Sqrt(46656*(a-b)*(a-b)*(b-c)*(b-c)*(-a+c)*(-a+c)-108*Math.pow(a2+b2-10*(a*b+a*c+b*c)+c2,3)),0.3333333333333333)/(6.*Math.pow(2,0.3333333333333333)))-Math.pow((a+b+c)/2.+(a2+b2-10*(a*b+a*c+b*c)+c2)/(Math.pow(2,0.6666666666666666)*Math.pow(216*(a-b)*(b-c)*(-a+c)+Sqrt(46656*(a-b)*(a-b)*(b-c)*(b-c)*(-a+c)*(-a+c)-108*Math.pow(a2+b2-10*(a*b+a*c+b*c)+c2,3)),0.3333333333333333))+Math.pow(216*(a-b)*(b-c)*(-a+c)+Sqrt(46656*(a-b)*(a-b)*(b-c)*(b-c)*(-a+c)*(-a+c)-108*Math.pow(a2+b2-10*(a*b+a*c+b*c)+c2,3)),0.3333333333333333)/(6.*Math.pow(2,0.3333333333333333)),2);
   let v2 = -2*a2-3*a*b+b2-2*a*c-2*b*c-c2+(2*a+b)*((a+b+c)/2.+(a2+b2-10*(a*b+a*c+b*c)+c2)/(Math.pow(2,0.6666666666666666)*Math.pow(216*(a-b)*(b-c)*(-a+c)+Sqrt(46656*(a-b)*(a-b)*(b-c)*(b-c)*(-a+c)*(-a+c)-108*Math.pow(a2+b2-10*(a*b+a*c+b*c)+c2,3)),0.3333333333333333))+Math.pow(216*(a-b)*(b-c)*(-a+c)+Sqrt(46656*(a-b)*(a-b)*(b-c)*(b-c)*(-a+c)*(-a+c)-108*Math.pow(a2+b2-10*(a*b+a*c+b*c)+c2,3)),0.3333333333333333)/(6.*Math.pow(2,0.3333333333333333)))-Math.pow((a+b+c)/2.+(a2+b2-10*(a*b+a*c+b*c)+c2)/(Math.pow(2,0.6666666666666666)*Math.pow(216*(a-b)*(b-c)*(-a+c)+Sqrt(46656*(a-b)*(a-b)*(b-c)*(b-c)*(-a+c)*(-a+c)-108*Math.pow(a2+b2-10*(a*b+a*c+b*c)+c2,3)),0.3333333333333333))+Math.pow(216*(a-b)*(b-c)*(-a+c)+Sqrt(46656*(a-b)*(a-b)*(b-c)*(b-c)*(-a+c)*(-a+c)-108*Math.pow(a2+b2-10*(a*b+a*c+b*c)+c2,3)),0.3333333333333333)/(6.*Math.pow(2,0.3333333333333333)),2);
   let v3 = -a2-2*a*b-2*b2-2*a*c-3*b*c+c2+(2*b+c)*((a+b+c)/2.+(a2+b2-10*(a*b+a*c+b*c)+c2)/(Math.pow(2,0.6666666666666666)*Math.pow(216*(a-b)*(b-c)*(-a+c)+Sqrt(46656*(a-b)*(a-b)*(b-c)*(b-c)*(-a+c)*(-a+c)-108*Math.pow(a2+b2-10*(a*b+a*c+b*c)+c2,3)),0.3333333333333333))+Math.pow(216*(a-b)*(b-c)*(-a+c)+Sqrt(46656*(a-b)*(a-b)*(b-c)*(b-c)*(-a+c)*(-a+c)-108*Math.pow(a2+b2-10*(a*b+a*c+b*c)+c2,3)),0.3333333333333333)/(6.*Math.pow(2,0.3333333333333333)))-Math.pow((a+b+c)/2.+(a2+b2-10*(a*b+a*c+b*c)+c2)/(Math.pow(2,0.6666666666666666)*Math.pow(216*(a-b)*(b-c)*(-a+c)+Sqrt(46656*(a-b)*(a-b)*(b-c)*(b-c)*(-a+c)*(-a+c)-108*Math.pow(a2+b2-10*(a*b+a*c+b*c)+c2,3)),0.3333333333333333))+Math.pow(216*(a-b)*(b-c)*(-a+c)+Sqrt(46656*(a-b)*(a-b)*(b-c)*(b-c)*(-a+c)*(-a+c)-108*Math.pow(a2+b2-10*(a*b+a*c+b*c)+c2,3)),0.3333333333333333)/(6.*Math.pow(2,0.3333333333333333)),2);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X370(orbit) {
   /* begin vars */
   let a2=a*a;
   let c2=c*c;
   let b2=b*b;
   let area=triAreaHeron(a,b,c);
   let S=2*area;
   let SA=(b2+c2-a2)/2;
   let SC=(a2+b2-c2)/2;
   let SB=(c2+a2-b2)/2;
   let S2=S*S;
   /* end vars */
   let v1 = S2-3*SB*SC;
   let v2 = S2-3*SA*SC;
   let v3 = S2-3*SA*SB;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X371(orbit) {
   /* begin vars */
   let area=triAreaHeron(a,b,c);
   let S=2*area;
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   /* end vars */
   let v1 = a2*(a2-b2-c2-2*S);
   let v2 = b2*(-a2+b2-c2-2*S);
   let v3 = c2*(-a2-b2+c2-2*S);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X372(orbit) {
   /* begin vars */
   let area=triAreaHeron(a,b,c);
   let S=2*area;
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   /* end vars */
   let v1 = a2*(a2-b2-c2+2*S);
   let v2 = b2*(-a2+b2-c2+2*S);
   let v3 = c2*(-a2-b2+c2+2*S);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X373(orbit) {
   /* begin vars */
   let a2=a*a;
   let c2=c*c;
   let b2=b*b;
   let a4=a2*a2;
   let c4=c2*c2;
   let b4=b2*b2;
   /* end vars */
   let v1 = a2*(a2*b2-b4+a2*c2+6*b2*c2-c4);
   let v2 = b2*(-a4+a2*b2+6*a2*c2+b2*c2-c4);
   let v3 = c2*(-a4+6*a2*b2-b4+a2*c2+b2*c2);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X374(orbit) {
   /* begin vars */
   let a2=a*a;
   let c2=c*c;
   let b2=b*b;
   let a4=a2*a2;
   let c4=c2*c2;
   let c3=c2*c;
   let b4=b2*b2;
   let b3=b2*b;
   let a3=a2*a;
   /* end vars */
   let v1 = a*(a3*b+a2*b2-a*b3-b4+a3*c-6*a2*b*c+5*a*b2*c+a2*c2+5*a*b*c2+2*b2*c2-a*c3-c4);
   let v2 = b*(-a4-a3*b+a2*b2+a*b3+5*a2*b*c-6*a*b2*c+b3*c+2*a2*c2+5*a*b*c2+b2*c2-b*c3-c4);
   let v3 = c*(-a4+2*a2*b2-b4-a3*c+5*a2*b*c+5*a*b2*c-b3*c+a2*c2-6*a*b*c2+b2*c2+a*c3+b*c3);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X375(orbit) {
   /* begin vars */
   let a2=a*a;
   let c2=c*c;
   let b2=b*b;
   let a3=a2*a;
   let a4=a2*a2;
   let c4=c2*c2;
   let c3=c2*c;
   let b3=b2*b;
   let b4=b2*b2;
   /* end vars */
   let v1 = a2*(a2*b2-b4-a*b2*c+b3*c+a2*c2-a*b*c2+4*b2*c2+b*c3-c4);
   let v2 = b2*(-a4+a2*b2+a3*c-a2*b*c+4*a2*c2-a*b*c2+b2*c2+a*c3-c4);
   let v3 = c2*(-a4+a3*b+4*a2*b2+a*b3-b4-a2*b*c-a*b2*c+a2*c2+b2*c2);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X376(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let SC=(a2+b2-c2)/2;
   let SB=(c2+a2-b2)/2;
   let SA=(b2+c2-a2)/2;
   /* end vars */
   let v1 = 2*a2*SA-SB*SC;
   let v2 = 2*b2*SB-SA*SC;
   let v3 = -(SA*SB)+2*c2*SC;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X377(orbit) {
   /* begin vars */
   let a2=a*a;
   let c2=c*c;
   let b2=b*b;
   let SA=(b2+c2-a2)/2;
   let SC=(a2+b2-c2)/2;
   let SB=(c2+a2-b2)/2;
   /* end vars */
   let v1 = a*b*c*(a+b+c)+2*SB*SC;
   let v2 = a*b*c*(a+b+c)+2*SA*SC;
   let v3 = a*b*c*(a+b+c)+2*SA*SB;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X378(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let SC=(a2+b2-c2)/2;
   let SB=(c2+a2-b2)/2;
   let SA=(b2+c2-a2)/2;
   let SC2=SC*SC;
   let SB2=SB*SB;
   let SA2=SA*SA;
   let area=triAreaHeron(a,b,c);
   /* end vars */
   let v1 = a2*(4*area*area+3*SA2)*SB*SC;
   let v2 = b2*SA*(4*area*area+3*SB2)*SC;
   let v3 = c2*SA*SB*(4*area*area+3*SC2);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X379(orbit) {
   /* begin vars */
   let c2=c*c;
   let c3=c2*c;
   let a2=a*a;
   let b2=b*b;
   let b3=b2*b;
   let a3=a2*a;
   let c5=c2*c3;
   let b5=b2*b3;
   let a4=a2*a2;
   let c4=c2*c2;
   let b4=b2*b2;
   let a5=a2*a3;
   /* end vars */
   let v1 = a5-a*b4+a2*b2*c-b4*c+a2*b*c2+2*a*b2*c2+b3*c2+b2*c3-a*c4-b*c4;
   let v2 = -(a4*b)+b5-a4*c+a2*b2*c+a3*c2+2*a2*b*c2+a*b2*c2+a2*c3-a*c4-b*c4;
   let v3 = -(a4*b)+a3*b2+a2*b3-a*b4-a4*c+2*a2*b2*c-b4*c+a2*b*c2+a*b2*c2+c5;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X380(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let c4=c2*c2;
   let b4=b2*b2;
   let SC=(a2+b2-c2)/2;
   let SB=(c2+a2-b2)/2;
   let SA=(b2+c2-a2)/2;
   let a4=a2*a2;
   /* end vars */
   let v1 = a*(a4-4*a*b*c*(b+c)-4*a2*(b*c+SA)-(SB-SC)*(SB-SC));
   let v2 = b*(b4-4*a*b*c*(a+c)-4*b2*(a*c+SB)-(-SA+SC)*(-SA+SC));
   let v3 = c*(-4*a*b*(a+b)*c+c4-(SA-SB)*(SA-SB)-4*c2*(a*b+SC));
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X381(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let SC=(a2+b2-c2)/2;
   let SB=(c2+a2-b2)/2;
   let SA=(b2+c2-a2)/2;
   /* end vars */
   let v1 = a2*SA+4*SB*SC;
   let v2 = b2*SB+4*SA*SC;
   let v3 = 4*SA*SB+c2*SC;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X382(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let SC=(a2+b2-c2)/2;
   let SB=(c2+a2-b2)/2;
   let SA=(b2+c2-a2)/2;
   /* end vars */
   let v1 = a2*SA-4*SB*SC;
   let v2 = b2*SB-4*SA*SC;
   let v3 = -4*SA*SB+c2*SC;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X383(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let SC=(a2+b2-c2)/2;
   let SB=(c2+a2-b2)/2;
   let SA=(b2+c2-a2)/2;
   let SC2=SC*SC;
   let SB2=SB*SB;
   let SA2=SA*SA;
   let area=triAreaHeron(a,b,c);
   /* end vars */
   let v1 = (c2+SC)*(a2*SA+SB*SC)*(a2*SA+4*SB*SC)-2*area*(SA2*SB2+c2*SA*SC*(3*SB+SC)+SB*SC*(c2*(SB+3*SC)+SC2))*sqrt3;
   let v2 = (a2+SA)*(b2*SB+SA*SC)*(b2*SB+4*SA*SC)-2*area*(a2*SA*SB*(SA+3*SC)+SA*SC*(SA2+a2*(3*SA+SC))+SB2*SC2)*sqrt3;
   let v3 = (b2+SB)*(SA*SB+c2*SC)*(4*SA*SB+c2*SC)-2*area*(SA*SB*(b2*(SA+3*SB)+SB2)+b2*SB*(3*SA+SB)*SC+SA2*SC2)*sqrt3;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X384(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let c4=c2*c2;
   let b4=b2*b2;
   let a4=a2*a2;
   /* end vars */
   let v1 = a4+b2*c2;
   let v2 = b4+a2*c2;
   let v3 = a2*b2+c4;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X385(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let c4=c2*c2;
   let b4=b2*b2;
   let a4=a2*a2;
   /* end vars */
   let v1 = a4-b2*c2;
   let v2 = b4-a2*c2;
   let v3 = -(a2*b2)+c4;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X386(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let SC=(a2+b2-c2)/2;
   let SB=(c2+a2-b2)/2;
   let SA=(b2+c2-a2)/2;
   /* end vars */
   let v1 = a2*((a+b)*(a+c)+2*SA);
   let v2 = b2*((a+b)*(b+c)+2*SB);
   let v3 = c2*((a+c)*(b+c)+2*SC);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X387(orbit) {
   /* begin vars */
   let b2=b*b;
   let a2=a*a;
   let c2=c*c;
   let SB=(c2+a2-b2)/2;
   let SA=(b2+c2-a2)/2;
   let SC=(a2+b2-c2)/2;
   let SB2=SB*SB;
   let SA2=SA*SA;
   let SC2=SC*SC;
   /* end vars */
   let v1 = a2*(a+c)*(b+c)+SC2;
   let v2 = (a+b)*b2*(a+c)+SA2;
   let v3 = (a+b)*(b+c)*c2+SB2;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X388(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   /* end vars */
   let v1 = (a2+b2+2*b*c+c2)/(a-b-c);
   let v2 = (a2+b2+2*a*c+c2)/(-a+b-c);
   let v3 = (a2+2*a*b+b2+c2)/(-a-b+c);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X389(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let SC=(a2+b2-c2)/2;
   let SB=(c2+a2-b2)/2;
   let SA=(b2+c2-a2)/2;
   let SC2=SC*SC;
   let c4=c2*c2;
   let SB2=SB*SB;
   let b4=b2*b2;
   let SA2=SA*SA;
   let area=triAreaHeron(a,b,c);
   let a4=a2*a2;
   /* end vars */
   let v1 = a4*b2*c2*SA+4*a2*area*area*(-SA2+SB*SC);
   let v2 = a2*b4*c2*SB+4*area*area*b2*(-SB2+SA*SC);
   let v3 = a2*b2*c4*SC+4*area*area*c2*(SA*SB-SC2);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X390(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   /* end vars */
   let v1 = (a-b-c)*(3*a2+b2-2*b*c+c2);
   let v2 = (-a+b-c)*(a2+3*b2-2*a*c+c2);
   let v3 = (-a-b+c)*(a2-2*a*b+b2+3*c2);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X391(orbit) {
   /* begin vars */

   /* end vars */
   let v1 = (a-b-c)*(3*a+b+c);
   let v2 = (-a+b-c)*(a+3*b+c);
   let v3 = (-a-b+c)*(a+b+3*c);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X392(orbit) {
   /* begin vars */
   let a2=a*a;
   let c2=c*c;
   let b2=b*b;
   let a3=a2*a;
   let c3=c2*c;
   let b3=b2*b;
   /* end vars */
   let v1 = a*(a2*b-b3+a2*c-4*a*b*c-b2*c-b*c2-c3);
   let v2 = b*(-a3+a*b2-a2*c-4*a*b*c+b2*c-a*c2-c3);
   let v3 = c*(-a3-a2*b-a*b2-b3-4*a*b*c+a*c2+b*c2);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X393(orbit) {
   /* begin vars */
   let a2=a*a;
   let c2=c*c;
   let b2=b*b;
   let SA=(b2+c2-a2)/2;
   let SC=(a2+b2-c2)/2;
   let SB=(c2+a2-b2)/2;
   let SA2=SA*SA;
   let SC2=SC*SC;
   let SB2=SB*SB;
   /* end vars */
   let v1 = SB2*SC2;
   let v2 = SA2*SC2;
   let v3 = SA2*SB2;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X394(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let SC=(a2+b2-c2)/2;
   let SB=(c2+a2-b2)/2;
   let SA=(b2+c2-a2)/2;
   let SC2=SC*SC;
   let SB2=SB*SB;
   let SA2=SA*SA;
   /* end vars */
   let v1 = a2*SA2;
   let v2 = b2*SB2;
   let v3 = c2*SC2;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X395(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let area=triAreaHeron(a,b,c);
   /* end vars */
   let v1 = -4*area+a2*sqrt3;
   let v2 = -4*area+b2*sqrt3;
   let v3 = -4*area+c2*sqrt3;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X396(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let area=triAreaHeron(a,b,c);
   /* end vars */
   let v1 = 4*area+a2*sqrt3;
   let v2 = 4*area+b2*sqrt3;
   let v3 = 4*area+c2*sqrt3;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X397(orbit) {
   /* begin vars */
   let a2=a*a;
   let c2=c*c;
   let b2=b*b;
   let SA=(b2+c2-a2)/2;
   let area=triAreaHeron(a,b,c);
   let SC=(a2+b2-c2)/2;
   let SB=(c2+a2-b2)/2;
   /* end vars */
   let v1 = SB*SC+a2*area*sqrt3;
   let v2 = SA*SC+area*b2*sqrt3;
   let v3 = SA*SB+area*c2*sqrt3;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X398(orbit) {
   /* begin vars */
   let a2=a*a;
   let c2=c*c;
   let b2=b*b;
   let SA=(b2+c2-a2)/2;
   let area=triAreaHeron(a,b,c);
   let SC=(a2+b2-c2)/2;
   let SB=(c2+a2-b2)/2;
   /* end vars */
   let v1 = SB*SC-a2*area*sqrt3;
   let v2 = SA*SC-area*b2*sqrt3;
   let v3 = SA*SB-area*c2*sqrt3;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X399(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let SC=(a2+b2-c2)/2;
   let SB=(c2+a2-b2)/2;
   let SA=(b2+c2-a2)/2;
   let SC2=SC*SC;
   let SB2=SB*SB;
   let SA2=SA*SA;
   let area=triAreaHeron(a,b,c);
   /* end vars */
   let v1 = a2*(-32*area*area*SA2+b2*c2*(5*a2*SA-4*SB*SC));
   let v2 = b2*(-32*area*area*SB2+a2*c2*(5*b2*SB-4*SA*SC));
   let v3 = c2*(a2*b2*(-4*SA*SB+5*c2*SC)-32*area*area*SC2);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X400(orbit) {
   /* begin vars */
   let cosC=lawOfCosines(c,a,b);
   let cosHalfC=cosHalfAngle(cosC);
   let sinQuarterC=sinHalfAngle(cosHalfC);
   let cosB=lawOfCosines(b,a,c);
   let cosHalfB=cosHalfAngle(cosB);
   let sinQuarterB=sinHalfAngle(cosHalfB);
   let cosA=lawOfCosines(a,b,c);
   let cosHalfA=cosHalfAngle(cosA);
   let sinQuarterA=sinHalfAngle(cosHalfA);
   let sinC=getSin(cosC);
   let cscQuarterC=1/sinQuarterC;
   let sinB=getSin(cosB);
   let cscQuarterB=1/sinQuarterB;
   let sinA=getSin(cosA);
   let cscQuarterA=1/sinQuarterA;
   /* end vars */
   let v1 = Math.pow(cscQuarterA,4)*sinA;
   let v2 = Math.pow(cscQuarterB,4)*sinB;
   let v3 = Math.pow(cscQuarterC,4)*sinC;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X401(orbit) {
   /* begin vars */
   let a2=a*a;
   let c2=c*c;
   let b2=b*b;
   let SA=(b2+c2-a2)/2;
   let SC=(a2+b2-c2)/2;
   let SB=(c2+a2-b2)/2;
   let SA2=SA*SA;
   let SC2=SC*SC;
   let SB2=SB*SB;
   /* end vars */
   let v1 = a2*SA*SB*SC+SB2*SC2-SA2*(SB2+SB*SC+SC2);
   let v2 = b2*SA*SB*SC+SA2*SC2-SB2*(SA2+SA*SC+SC2);
   let v3 = SA2*SB2+c2*SA*SB*SC-(SA2+SA*SB+SB2)*SC2;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X402(orbit) {
   /* begin vars */
   let a2=a*a;
   let c2=c*c;
   let b2=b*b;
   let SA=(b2+c2-a2)/2;
   let SC=(a2+b2-c2)/2;
   let SB=(c2+a2-b2)/2;
   let SA2=SA*SA;
   let SC2=SC*SC;
   let SB2=SB*SB;
   /* end vars */
   let v1 = (a2*SA-2*SB*SC)*(a2*SA*SB*SC-SB2*SC2+SA2*(SB2-3*SB*SC+SC2));
   let v2 = (b2*SB-2*SA*SC)*(b2*SA*SB*SC-SA2*SC2+SB2*(SA2-3*SA*SC+SC2));
   let v3 = (-2*SA*SB+c2*SC)*(-(SA2*SB2)+c2*SA*SB*SC+(SA2-3*SA*SB+SB2)*SC2);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X403(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let SC=(a2+b2-c2)/2;
   let SB=(c2+a2-b2)/2;
   let SA=(b2+c2-a2)/2;
   let SC2=SC*SC;
   let SB2=SB*SB;
   let SA2=SA*SA;
   /* end vars */
   let v1 = SB*SC*(SA*(SB-SC)*(SB-SC)+a2*(-SA2+SB*SC));
   let v2 = SA*SC*(SB*(-SA+SC)*(-SA+SC)+b2*(-SB2+SA*SC));
   let v3 = SA*SB*((SA-SB)*(SA-SB)*SC+c2*(SA*SB-SC2));
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X404(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let SC=(a2+b2-c2)/2;
   let SB=(c2+a2-b2)/2;
   let SA=(b2+c2-a2)/2;
   /* end vars */
   let v1 = a*(b*c*(a+b+c)-2*a*SA);
   let v2 = b*(a*c*(a+b+c)-2*b*SB);
   let v3 = c*(a*b*(a+b+c)-2*c*SC);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X405(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let SC=(a2+b2-c2)/2;
   let SB=(c2+a2-b2)/2;
   let SA=(b2+c2-a2)/2;
   /* end vars */
   let v1 = a*(b*c*(a+b+c)+a*SA);
   let v2 = b*(a*c*(a+b+c)+b*SB);
   let v3 = c*(a*b*(a+b+c)+c*SC);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X406(orbit) {
   /* begin vars */
   let a2=a*a;
   let c2=c*c;
   let b2=b*b;
   let a3=a2*a;
   let SC=(a2+b2-c2)/2;
   let SA=(b2+c2-a2)/2;
   let c3=c2*c;
   let b3=b2*b;
   let SB=(c2+a2-b2)/2;
   /* end vars */
   let v1 = SB*(b3+c3+2*a*(b*c+SA)+c*(SA-SB)+b*(SA-SC))*SC;
   let v2 = SA*(a3+c3+2*b*(a*c+SB)+c*(-SA+SB)+a*(SB-SC))*SC;
   let v3 = SA*SB*(a3+b3+2*c*(a*b+SC)+b*(-SA+SC)+a*(-SB+SC));
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X407(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let SC=(a2+b2-c2)/2;
   let SB=(c2+a2-b2)/2;
   let SA=(b2+c2-a2)/2;
   /* end vars */
   let v1 = (b+c)*(a2+2*b*c+a*(b+c)-2*SA)*SB*SC;
   let v2 = (a+c)*SA*(b2+2*a*c+b*(a+c)-2*SB)*SC;
   let v3 = (a+b)*SA*SB*(2*a*b+(a+b)*c+c2-2*SC);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X408(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let c4=c2*c2;
   let c3=c2*c;
   let b4=b2*b2;
   let b3=b2*b;
   let a4=a2*a2;
   let a3=a2*a;
   /* end vars */
   let v1 = a3*(b+c)*(-a2+b2+c2)*(-a2+b2+c2)*(-(a4*b)-a3*b2+a2*b3+a*b4-a4*c-a2*b2*c+2*b4*c-a3*c2-a2*b*c2-2*a*b2*c2-2*b3*c2+a2*c3-2*b2*c3+a*c4+2*b*c4);
   let v2 = b3*(a+c)*(a2-b2+c2)*(a2-b2+c2)*(a4*b+a3*b2-a2*b3-a*b4+2*a4*c-a2*b2*c-b4*c-2*a3*c2-2*a2*b*c2-a*b2*c2-b3*c2-2*a2*c3+b2*c3+2*a*c4+b*c4);
   let v3 = (a+b)*(a2+b2-c2)*(a2+b2-c2)*c3*(2*a4*b-2*a3*b2-2*a2*b3+2*a*b4+a4*c-2*a2*b2*c+b4*c+a3*c2-a2*b*c2-a*b2*c2+b3*c2-a2*c3-b2*c3-a*c4-b*c4);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X409(orbit) {
   /* begin vars */

   /* end vars */
   let v1 = a*(a+b)*(a+c)*((a+b-c)*(a-b+c)*(b+c)*(b+c)+(a+b)*(a+c)*(-a+b+c)*(-a+b+c));
   let v2 = b*(a+b)*(b+c)*((a+b)*(a-b+c)*(a-b+c)*(b+c)+(a+b-c)*(a+c)*(a+c)*(-a+b+c));
   let v3 = c*(a+c)*(b+c)*((a+b-c)*(a+b-c)*(a+c)*(b+c)+(a+b)*(a+b)*(a-b+c)*(-a+b+c));
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X410(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let SC=(a2+b2-c2)/2;
   let SB=(c2+a2-b2)/2;
   let SA=(b2+c2-a2)/2;
   let SC2=SC*SC;
   let SB2=SB*SB;
   let SA2=SA*SA;
   /* end vars */
   let v1 = (a+b)*(a+c)*SB*SC*(a2*(a+b-c)*(a-b+c)*(b+c)*(b+c)*SA2+b*(a+b)*c*(a+c)*(-a+b+c)*(-a+b+c)*SB*SC);
   let v2 = (a+b)*(b+c)*SA*SC*(b2*(a+b-c)*(a+c)*(a+c)*(-a+b+c)*SB2+a*(a+b)*c*(a-b+c)*(a-b+c)*(b+c)*SA*SC);
   let v3 = (a+c)*(b+c)*SA*SB*(a*b*(a+b-c)*(a+b-c)*(a+c)*(b+c)*SA*SB+(a+b)*(a+b)*(a-b+c)*(-a+b+c)*c2*SC2);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X411(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let SC=(a2+b2-c2)/2;
   let c3=c2*c;
   let SB=(c2+a2-b2)/2;
   let b3=b2*b;
   let SA=(b2+c2-a2)/2;
   let a3=a2*a;
   let SC2=SC*SC;
   let c5=c2*c3;
   let c4=c2*c2;
   let SB2=SB*SB;
   let b5=b2*b3;
   let b4=b2*b2;
   let SA2=SA*SA;
   let a5=a2*a3;
   let a4=a2*a2;
   /* end vars */
   let v1 = a*(a4*b*c+a5*(b+c)-b*c*(b2-c2)*(b2-c2)+a*(b-c)*(b-c)*(b+c)*(b2+c2)-2*a3*(b3+c3)-4*a2*SA2);
   let v2 = b*(a*b4*c+b5*(a+c)-a*c*(-a2+c2)*(-a2+c2)+b*(-a+c)*(-a+c)*(a+c)*(a2+c2)-2*b3*(a3+c3)-4*b2*SB2);
   let v3 = c*(-(a*b*(a2-b2)*(a2-b2))+(a-b)*(a-b)*(a+b)*(a2+b2)*c-2*(a3+b3)*c3+a*b*c4+(a+b)*c5-4*c2*SC2);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X412(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let SC=(a2+b2-c2)/2;
   let SB=(c2+a2-b2)/2;
   let SA=(b2+c2-a2)/2;
   let c3=c2*c;
   let b3=b2*b;
   let a3=a2*a;
   let SC2=SC*SC;
   let SB2=SB*SB;
   let SA2=SA*SA;
   let c5=c2*c3;
   let c4=c2*c2;
   let b4=b2*b2;
   let b5=b2*b3;
   let a4=a2*a2;
   let a5=a2*a3;
   /* end vars */
   let v1 = (-(a5*(b+c))+b*c*(a4-(b-c)*(b-c)*(b+c)*(b+c))+2*a3*(b3+c3)+a*(-b5+b4*c+b*c4-c5)-4*a2*SA2)*SB*SC;
   let v2 = SA*(-(b5*(a+c))+a*c*(b4-(-a+c)*(-a+c)*(a+c)*(a+c))+2*b3*(a3+c3)+b*(-a5+a4*c+a*c4-c5)-4*b2*SB2)*SC;
   let v3 = SA*SB*((-a5+a4*b+a*b4-b5)*c+2*(a3+b3)*c3+a*b*(-((a-b)*(a-b)*(a+b)*(a+b))+c4)-(a+b)*c5-4*c2*SC2);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X413(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let c4=c2*c2;
   let c3=c2*c;
   let b3=b2*b;
   let b4=b2*b2;
   let a3=a2*a;
   let a4=a2*a2;
   /* end vars */
   let v1 = a*(a+b)*Math.pow(a-b-c,3)*(a+c)*(a4+a2*b*c+a3*(b+c)+a*(b-c)*(b-c)*(b+c)+(b-c)*(b-c)*(b2+b*c+c2));
   let v2 = b*(a+b)*Math.pow(-a+b-c,3)*(b+c)*(b4+a*b2*c+b3*(a+c)+b*(-a+c)*(-a+c)*(a+c)+(-a+c)*(-a+c)*(a2+a*c+c2));
   let v3 = c*(a+c)*Math.pow(-a-b+c,3)*(b+c)*((a-b)*(a-b)*(a2+a*b+b2)+(a-b)*(a-b)*(a+b)*c+a*b*c2+(a+b)*c3+c4);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X414(orbit) {
   /* begin vars */
   let a2=a*a;
   let c2=c*c;
   let b2=b*b;
   let SA=(b2+c2-a2)/2;
   let SC=(a2+b2-c2)/2;
   let SB=(c2+a2-b2)/2;
   let SA2=SA*SA;
   let SC2=SC*SC;
   let sb=(c+a-b)/2;
   let sc=(a+b-c)/2;
   let SB2=SB*SB;
   let sa=(b+c-a)/2;
   /* end vars */
   let v1 = (a+b)*(a+c)*Math.pow(sa,3)*SB*SC*(b2*(a+c)*(a+c)*SB2*sc*sc-b*(a+b)*c*(a+c)*sb*SB*sc*SC+(a+b)*(a+b)*c2*sb*sb*SC2);
   let v2 = (a+b)*(b+c)*SA*Math.pow(sb,3)*SC*(a2*(b+c)*(b+c)*SA2*sc*sc-a*(a+b)*c*(b+c)*sa*SA*sc*SC+(a+b)*(a+b)*c2*sa*sa*SC2);
   let v3 = (a+c)*(b+c)*SA*SB*(a2*(b+c)*(b+c)*SA2*sb*sb-a*b*(a+c)*(b+c)*sa*SA*sb*SB+b2*(a+c)*(a+c)*sa*sa*SB2)*Math.pow(sc,3);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X415(orbit) {
   /* begin vars */
   let a2=a*a;
   let c2=c*c;
   let b2=b*b;
   let SA=(b2+c2-a2)/2;
   let SC=(a2+b2-c2)/2;
   let SB=(c2+a2-b2)/2;
   let c3=c2*c;
   let b3=b2*b;
   let a3=a2*a;
   /* end vars */
   let v1 = (a+b)*(a+c)*(a3-2*a2*b+b3-2*a2*c+a*b*c+c3)*SB*SC;
   let v2 = (a+b)*(b+c)*(a3-2*a*b2+b3+a*b*c-2*b2*c+c3)*SA*SC;
   let v3 = (a+c)*(b+c)*(a3+b3+a*b*c-2*a*c2-2*b*c2+c3)*SA*SB;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X416(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let SC=(a2+b2-c2)/2;
   let SB=(c2+a2-b2)/2;
   let SA=(b2+c2-a2)/2;
   let SC2=SC*SC;
   let SB2=SB*SB;
   let sa=(b+c-a)/2;
   let sc=(a+b-c)/2;
   let sb=(c+a-b)/2;
   let SA2=SA*SA;
   /* end vars */
   let v1 = a*(a+b)*(a+c)*(-(a2*(b+c)*(b+c)*SA2*sb*sc)+b*(a+b)*c*(a+c)*sa*sa*SB*SC);
   let v2 = b*(a+b)*(b+c)*(-(b2*(a+c)*(a+c)*sa*SB2*sc)+a*(a+b)*c*(b+c)*SA*sb*sb*SC);
   let v3 = c*(a+c)*(b+c)*(a*b*(a+c)*(b+c)*SA*SB*sc*sc-(a+b)*(a+b)*c2*sa*sb*SC2);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X417(orbit) {
   /* begin vars */
   let c2=c*c;
   let a2=a*a;
   let b2=b*b;
   let SA=(b2+c2-a2)/2;
   let SC=(a2+b2-c2)/2;
   let SB=(c2+a2-b2)/2;
   let c4=c2*c2;
   let SA2=SA*SA;
   let b4=b2*b2;
   let SC2=SC*SC;
   let SB2=SB*SB;
   let a4=a2*a2;
   /* end vars */
   let v1 = a4*Math.pow(SA,3)*(a2*SB*SC+SA*(SB2+SC2));
   let v2 = b4*Math.pow(SB,3)*(b2*SA*SC+SB*(SA2+SC2));
   let v3 = c4*Math.pow(SC,3)*(c2*SA*SB+(SA2+SB2)*SC);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X418(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let SC=(a2+b2-c2)/2;
   let SB=(c2+a2-b2)/2;
   let SA=(b2+c2-a2)/2;
   let SC2=SC*SC;
   let c4=c2*c2;
   let SB2=SB*SB;
   let b4=b2*b2;
   let SA2=SA*SA;
   let a4=a2*a2;
   /* end vars */
   let v1 = a4*SA2*(a2*SA+2*SB*SC);
   let v2 = b4*SB2*(b2*SB+2*SA*SC);
   let v3 = c4*(2*SA*SB+c2*SC)*SC2;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X419(orbit) {
   /* begin vars */
   let c2=c*c;
   let a2=a*a;
   let b2=b*b;
   let c4=c2*c2;
   let SA=(b2+c2-a2)/2;
   let b4=b2*b2;
   let SC=(a2+b2-c2)/2;
   let SB=(c2+a2-b2)/2;
   let a4=a2*a2;
   /* end vars */
   let v1 = (a4-b2*c2)*SB*SC;
   let v2 = (b4-a2*c2)*SA*SC;
   let v3 = (-(a2*b2)+c4)*SA*SB;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X420(orbit) {
   /* begin vars */
   let a2=a*a;
   let c2=c*c;
   let b2=b*b;
   let SA=(b2+c2-a2)/2;
   let SC=(a2+b2-c2)/2;
   let SB=(c2+a2-b2)/2;
   let SA2=SA*SA;
   let SC2=SC*SC;
   let SB2=SB*SB;
   /* end vars */
   let v1 = SB*SC*(-(SA*(a2+3*SA))+SB2+3*SB*SC+SC2);
   let v2 = SA*SC*(SA2-SB*(b2+3*SB)+3*SA*SC+SC2);
   let v3 = SA*SB*(SA2+3*SA*SB+SB2-SC*(c2+3*SC));
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X421(orbit) {
   /* begin vars */
   let a2=a*a;
   let c2=c*c;
   let b2=b*b;
   let SA=(b2+c2-a2)/2;
   let SC=(a2+b2-c2)/2;
   let SB=(c2+a2-b2)/2;
   let SA2=SA*SA;
   let SC2=SC*SC;
   let SB2=SB*SB;
   /* end vars */
   let v1 = SB*SC*(a2*SA*(SB2+SC2)+SB*SC*(SB2+SC2)-SA2*(SB2+4*SB*SC+SC2));
   let v2 = SA*SC*(b2*SB*(SA2+SC2)+SA*SC*(SA2+SC2)-SB2*(SA2+4*SA*SC+SC2));
   let v3 = SA*SB*(SA*SB*(SA2+SB2)+c2*(SA2+SB2)*SC-(SA2+4*SA*SB+SB2)*SC2);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X422(orbit) {
   /* begin vars */
   let c2=c*c;
   let a2=a*a;
   let b2=b*b;
   let c3=c2*c;
   let SA=(b2+c2-a2)/2;
   let b3=b2*b;
   let SC=(a2+b2-c2)/2;
   let SB=(c2+a2-b2)/2;
   let a3=a2*a;
   /* end vars */
   let v1 = (a+b)*(a+c)*(a3+a*b*c-b*c*(b+c))*SB*SC;
   let v2 = (a+b)*(b+c)*(b3+a*b*c-a*c*(a+c))*SA*SC;
   let v3 = (a+c)*(b+c)*(-(a*b*(a+b))+a*b*c+c3)*SA*SB;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X423(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let SC=(a2+b2-c2)/2;
   let SB=(c2+a2-b2)/2;
   let SA=(b2+c2-a2)/2;
   /* end vars */
   let v1 = (a+b)*(a+c)*(a*b+a*c-b*c-2*SA)*SB*SC;
   let v2 = (a+b)*(b+c)*SA*(a*b-a*c+b*c-2*SB)*SC;
   let v3 = (a+c)*(b+c)*SA*SB*(-(a*b)+a*c+b*c-2*SC);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X424(orbit) {
   /* begin vars */
   let a2=a*a;
   let c2=c*c;
   let b2=b*b;
   let SA=(b2+c2-a2)/2;
   let a4=a2*a2;
   let SC=(a2+b2-c2)/2;
   let SB=(c2+a2-b2)/2;
   let c4=c2*c2;
   let c3=c2*c;
   let b4=b2*b2;
   let b3=b2*b;
   let a3=a2*a;
   /* end vars */
   let v1 = (b+c)*(a3*b+a2*b2-a*b3-b4+a3*c+a2*c2-a*c3-c4)*SB*SC;
   let v2 = (a+c)*(-a4-a3*b+a2*b2+a*b3+b3*c+b2*c2-b*c3-c4)*SA*SC;
   let v3 = (a+b)*(-a4-b4-a3*c-b3*c+a2*c2+b2*c2+a*c3+b*c3)*SA*SB;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X425(orbit) {
   /* begin vars */
   let c2=c*c;
   let c4=c2*c2;
   let a2=a*a;
   let b2=b*b;
   let b4=b2*b2;
   let c3=c2*c;
   let b3=b2*b;
   let a3=a2*a;
   let a4=a2*a2;
   let c6=c2*c4;
   let SA=(b2+c2-a2)/2;
   let b6=b2*b4;
   let SC=(a2+b2-c2)/2;
   let SB=(c2+a2-b2)/2;
   let c5=c2*c3;
   let b5=b2*b3;
   let a5=a2*a3;
   let a6=a2*a4;
   /* end vars */
   let v1 = (a+b)*(a+c)*(-a6+a5*b+a4*b2-a3*b3+a5*c-a4*b*c-a3*b2*c+b5*c+a4*c2-a3*b*c2-2*a2*b2*c2+2*a*b3*c2-a3*c3+2*a*b2*c3-2*b3*c3+b*c5)*SB*SC;
   let v2 = (a+b)*(b+c)*(-(a3*b3)+a2*b4+a*b5-b6+a5*c-a2*b3*c-a*b4*c+b5*c+2*a3*b*c2-2*a2*b2*c2-a*b3*c2+b4*c2-2*a3*c3+2*a2*b*c3-b3*c3+a*c5)*SA*SC;
   let v3 = (a+c)*(b+c)*(a5*b-2*a3*b3+a*b5+2*a3*b2*c+2*a2*b3*c-2*a2*b2*c2-a3*c3-a2*b*c3-a*b2*c3-b3*c3+a2*c4-a*b*c4+b2*c4+a*c5+b*c5-c6)*SA*SB;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X426(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let SA=(b2+c2-a2)/2;
   let SC=(a2+b2-c2)/2;
   let SB=(c2+a2-b2)/2;
   let SA2=SA*SA;
   let SC2=SC*SC;
   let SB2=SB*SB;
   /* end vars */
   let v1 = a2*Math.pow(SA,3)*(SB2+SC2);
   let v2 = b2*Math.pow(SB,3)*(SA2+SC2);
   let v3 = c2*(SA2+SB2)*Math.pow(SC,3);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X427(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let SC=(a2+b2-c2)/2;
   let SB=(c2+a2-b2)/2;
   let SA=(b2+c2-a2)/2;
   /* end vars */
   let v1 = (a2+2*SA)*SB*SC;
   let v2 = SA*(b2+2*SB)*SC;
   let v3 = SA*SB*(c2+2*SC);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X428(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let SC=(a2+b2-c2)/2;
   let SB=(c2+a2-b2)/2;
   let SA=(b2+c2-a2)/2;
   /* end vars */
   let v1 = (3*a2+2*SA)*SB*SC;
   let v2 = SA*(3*b2+2*SB)*SC;
   let v3 = SA*SB*(3*c2+2*SC);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X429(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let SC=(a2+b2-c2)/2;
   let SB=(c2+a2-b2)/2;
   let SA=(b2+c2-a2)/2;
   /* end vars */
   let v1 = (b+c)*(a*(a+b+c)+2*SA)*SB*SC;
   let v2 = (a+c)*SA*(b*(a+b+c)+2*SB)*SC;
   let v3 = (a+b)*SA*SB*(c*(a+b+c)+2*SC);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X430(orbit) {
   /* begin vars */
   let a2=a*a;
   let c2=c*c;
   let b2=b*b;
   let SA=(b2+c2-a2)/2;
   let SC=(a2+b2-c2)/2;
   let SB=(c2+a2-b2)/2;
   /* end vars */
   let v1 = (b+c)*(2*a+b+c)*SB*SC;
   let v2 = (a+c)*(a+2*b+c)*SA*SC;
   let v3 = (a+b)*(a+b+2*c)*SA*SB;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X431(orbit) {
   /* begin vars */
   let a2=a*a;
   let a3=a2*a;
   let c2=c*c;
   let c3=c2*c;
   let b2=b*b;
   let b3=b2*b;
   let a5=a2*a3;
   let c5=c2*c3;
   let c4=c2*c2;
   let b4=b2*b2;
   let b5=b2*b3;
   let a4=a2*a2;
   /* end vars */
   let v1 = (a4*b-2*a2*b3+b5+a4*c+2*a3*b*c+2*a2*b2*c-b4*c+2*a2*b*c2-2*a2*c3-b*c4+c5)/((a+b)*(a+c)*(a2-b2-c2));
   let v2 = (a5-2*a3*b2+a*b4-a4*c+2*a2*b2*c+2*a*b3*c+b4*c+2*a*b2*c2-2*b2*c3-a*c4+c5)/((a+b)*(b+c)*(-a2+b2-c2));
   let v3 = (a5-a4*b-a*b4+b5-2*a3*c2+2*a2*b*c2+2*a*b2*c2-2*b3*c2+2*a*b*c3+a*c4+b*c4)/((a+c)*(b+c)*(-a2-b2+c2));
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X432(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let SC=(a2+b2-c2)/2;
   let SA=(b2+c2-a2)/2;
   let SB=(c2+a2-b2)/2;
   let SC2=SC*SC;
   let SA2=SA*SA;
   let SB2=SB*SB;
   /* end vars */
   let v1 = SB*SC*(b2*SB2*Math.pow(b2*(SB2-SA*SC)-SB*(SA2+SC2),2)+c2*SC2*Math.pow((-SA2-SB2)*SC+c2*(-(SA*SB)+SC2),2));
   let v2 = SA*SC*(a2*SA2*Math.pow(a2*(SA2-SB*SC)+SA*(-SB2-SC2),2)+c2*SC2*Math.pow(-((SA2+SB2)*SC)+c2*(-(SA*SB)+SC2),2));
   let v3 = SA*SB*(b2*SB2*Math.pow(b2*(SB2-SA*SC)+SB*(-SA2-SC2),2)+a2*SA2*Math.pow(a2*(SA2-SB*SC)-SA*(SB2+SC2),2));
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X433(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let SC=(a2+b2-c2)/2;
   let SB=(c2+a2-b2)/2;
   let SA=(b2+c2-a2)/2;
   let SC2=SC*SC;
   let SB2=SB*SB;
   let SA2=SA*SA;
   /* end vars */
   let v1 = SB*SC*(c2*Math.pow(-((SA2+SB2)*SC)+c2*(SA*SB-SC2),2)+b2*Math.pow(b2*(-SB2+SA*SC)-SB*(SA2+SC2),2));
   let v2 = SA*SC*(c2*Math.pow(-((SA2+SB2)*SC)+c2*(SA*SB-SC2),2)+a2*Math.pow(a2*(-SA2+SB*SC)-SA*(SB2+SC2),2));
   let v3 = SA*SB*(b2*Math.pow(b2*(-SB2+SA*SC)-SB*(SA2+SC2),2)+a2*Math.pow(a2*(-SA2+SB*SC)-SA*(SB2+SC2),2));
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X434(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let SC=(a2+b2-c2)/2;
   let SB=(c2+a2-b2)/2;
   let SA=(b2+c2-a2)/2;
   let SC2=SC*SC;
   let SB2=SB*SB;
   let SA2=SA*SA;
   /* end vars */
   let v1 = SB*SC*(c2*Math.pow(4*SA2*SB2+9*SA2*SB*SC+9*SA*SB2*SC-3*SA*Math.pow(SC,3)-3*SB*Math.pow(SC,3)+5*SA2*SC2+6*SA*SB*SC2+5*SB2*SC2,2)+b2*Math.pow(-3*SA*Math.pow(SB,3)+5*SA2*SB2+9*SA2*SB*SC-3*Math.pow(SB,3)*SC+6*SA*SB2*SC+4*SA2*SC2+9*SA*SB*SC2+5*SB2*SC2,2));
   let v2 = SA*SC*(a2*Math.pow(-3*Math.pow(SA,3)*SB+5*SA2*SB2-3*Math.pow(SA,3)*SC+6*SA2*SB*SC+9*SA*SB2*SC+5*SA2*SC2+9*SA*SB*SC2+4*SB2*SC2,2)+c2*Math.pow(4*SA2*SB2+9*SA2*SB*SC+9*SA*SB2*SC-3*SA*Math.pow(SC,3)-3*SB*Math.pow(SC,3)+5*SA2*SC2+6*SA*SB*SC2+5*SB2*SC2,2));
   let v3 = SA*SB*(a2*Math.pow(-3*Math.pow(SA,3)*SB+5*SA2*SB2-3*Math.pow(SA,3)*SC+6*SA2*SB*SC+9*SA*SB2*SC+5*SA2*SC2+9*SA*SB*SC2+4*SB2*SC2,2)+b2*Math.pow(-3*SA*Math.pow(SB,3)+5*SA2*SB2+9*SA2*SB*SC-3*Math.pow(SB,3)*SC+6*SA*SB2*SC+4*SA2*SC2+9*SA*SB*SC2+5*SB2*SC2,2));
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X435(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let SC=(a2+b2-c2)/2;
   let SB=(c2+a2-b2)/2;
   let SA=(b2+c2-a2)/2;
   let SC2=SC*SC;
   let SB2=SB*SB;
   let SA2=SA*SA;
   /* end vars */
   let v1 = SB*SC*(b2*Math.pow(-(a2*SA2*(5*SB-4*SC))+a2*SA*SB*(3*SB-SC)+SB2*(3*SB-5*SC)*SC,2)+c2*Math.pow(b2*SB*SC*(-SA+3*SC)-b2*SB2*(-4*SA+5*SC)+SA*(-5*SA+3*SC)*SC2,2));
   let v2 = SA*SC*(a2*Math.pow(SA2*(3*SA-5*SB)*SB+c2*SA*(3*SA-SB)*SC-c2*(5*SA-4*SB)*SC2,2)+c2*Math.pow(b2*SB*SC*(-SA+3*SC)-b2*SB2*(-4*SA+5*SC)+SA*(-5*SA+3*SC)*SC2,2));
   let v3 = SA*SB*(b2*Math.pow(-(a2*SA2*(5*SB-4*SC))+a2*SA*SB*(3*SB-SC)+SB2*(3*SB-5*SC)*SC,2)+a2*Math.pow(SA2*(3*SA-5*SB)*SB+c2*SA*(3*SA-SB)*SC-c2*(5*SA-4*SB)*SC2,2));
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X436(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let c4=c2*c2;
   let b4=b2*b2;
   let a4=a2*a2;
   /* end vars */
   let v1 = (a2+b2-c2)*(a2-b2+c2)*(a4*(a2-b2-c2)*(a2-b2-c2)+b2*(-a2+b2-c2)*c2*(-a2-b2+c2));
   let v2 = (a2+b2-c2)*(-a2+b2+c2)*(b4*(-a2+b2-c2)*(-a2+b2-c2)+a2*(a2-b2-c2)*c2*(-a2-b2+c2));
   let v3 = (a2-b2+c2)*(-a2+b2+c2)*(a2*b2*(a2-b2-c2)*(-a2+b2-c2)+(-a2-b2+c2)*(-a2-b2+c2)*c4);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X437(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   /* end vars */
   let v1 = (a2+b2-c2)*(a2-b2+c2)*((2*a-b-c)*(2*a-b-c)*(a2-b2+b*c-c2)*(a2-b2+b*c-c2)+(-a+2*b-c)*(-a-b+2*c)*(-a2+b2+a*c-c2)*(-a2+a*b-b2+c2));
   let v2 = (a2+b2-c2)*(-a2+b2+c2)*((-a+2*b-c)*(-a+2*b-c)*(-a2+b2+a*c-c2)*(-a2+b2+a*c-c2)+(2*a-b-c)*(-a-b+2*c)*(a2-b2+b*c-c2)*(-a2+a*b-b2+c2));
   let v3 = (a2-b2+c2)*(-a2+b2+c2)*((2*a-b-c)*(-a+2*b-c)*(-a2+b2+a*c-c2)*(a2-b2+b*c-c2)+(-a-b+2*c)*(-a-b+2*c)*(-a2+a*b-b2+c2)*(-a2+a*b-b2+c2));
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X438(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let c4=c2*c2;
   let b4=b2*b2;
   let a4=a2*a2;
   /* end vars */
   let v1 = (a2+b2-c2)*(a2+b2-c2)*(a2-b2+c2)*(a2-b2+c2)*(3*a4-2*a2*b2-b4-2*a2*c2+2*b2*c2-c4)*(3*a4-2*a2*b2-b4-2*a2*c2+2*b2*c2-c4)+(a2+b2-c2)*(a2-b2+c2)*(-a2+b2+c2)*(-a2+b2+c2)*(-a4-2*a2*b2+3*b4+2*a2*c2-2*b2*c2-c4)*(-a4+2*a2*b2-b4-2*a2*c2-2*b2*c2+3*c4);
   let v2 = (a2+b2-c2)*(a2+b2-c2)*(-a2+b2+c2)*(-a2+b2+c2)*(-a4-2*a2*b2+3*b4+2*a2*c2-2*b2*c2-c4)*(-a4-2*a2*b2+3*b4+2*a2*c2-2*b2*c2-c4)+(a2+b2-c2)*(a2-b2+c2)*(a2-b2+c2)*(-a2+b2+c2)*(3*a4-2*a2*b2-b4-2*a2*c2+2*b2*c2-c4)*(-a4+2*a2*b2-b4-2*a2*c2-2*b2*c2+3*c4);
   let v3 = (a2+b2-c2)*(a2+b2-c2)*(a2-b2+c2)*(-a2+b2+c2)*(-a4-2*a2*b2+3*b4+2*a2*c2-2*b2*c2-c4)*(3*a4-2*a2*b2-b4-2*a2*c2+2*b2*c2-c4)+(a2-b2+c2)*(a2-b2+c2)*(-a2+b2+c2)*(-a2+b2+c2)*(-a4+2*a2*b2-b4-2*a2*c2-2*b2*c2+3*c4)*(-a4+2*a2*b2-b4-2*a2*c2-2*b2*c2+3*c4);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X439(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   /* end vars */
   let v1 = (3*a2-b2-c2)*(3*a2-b2-c2);
   let v2 = (-a2+3*b2-c2)*(-a2+3*b2-c2);
   let v3 = (-a2-b2+3*c2)*(-a2-b2+3*c2);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X440(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let SC=(a2+b2-c2)/2;
   let SB=(c2+a2-b2)/2;
   let SA=(b2+c2-a2)/2;
   let c3=c2*c;
   let b3=b2*b;
   let a3=a2*a;
   /* end vars */
   let v1 = (b+c)*(2*a3+a2*b+b3+a2*c-b2*c-b*c2+c3)*SA;
   let v2 = (a+c)*(a3+a*b2+2*b3-a2*c+b2*c-a*c2+c3)*SB;
   let v3 = (a+b)*(a3-a2*b-a*b2+b3+a*c2+b*c2+2*c3)*SC;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X441(orbit) {
   /* begin vars */
   let a2=a*a;
   let c2=c*c;
   let b2=b*b;
   let SA=(b2+c2-a2)/2;
   let SC=(a2+b2-c2)/2;
   let SB=(c2+a2-b2)/2;
   let SA2=SA*SA;
   let SC2=SC*SC;
   let SB2=SB*SB;
   /* end vars */
   let v1 = SA*(-(a2*SB*SC)+SA*(SB2+SC2));
   let v2 = SB*(-(b2*SA*SC)+SB*(SA2+SC2));
   let v3 = SC*(-(c2*SA*SB)+(SA2+SB2)*SC);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X442(orbit) {
   /* begin vars */
   let a2=a*a;
   let c2=c*c;
   let b2=b*b;
   let a3=a2*a;
   let c3=c2*c;
   let b3=b2*b;
   /* end vars */
   let v1 = (b+c)*(a2*b-b3+a2*c+2*a*b*c+b2*c+b*c2-c3);
   let v2 = (a+c)*(-a3+a*b2+a2*c+2*a*b*c+b2*c+a*c2-c3);
   let v3 = (a+b)*(-a3+a2*b+a*b2-b3+2*a*b*c+a*c2+b*c2);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X443(orbit) {
   /* begin vars */
   let a2=a*a;
   let c2=c*c;
   let b2=b*b;
   let SA=(b2+c2-a2)/2;
   let SC=(a2+b2-c2)/2;
   let SB=(c2+a2-b2)/2;
   /* end vars */
   let v1 = a*b*c*(a+b+c)+SB*SC;
   let v2 = a*b*c*(a+b+c)+SA*SC;
   let v3 = a*b*c*(a+b+c)+SA*SB;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X444(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let SC=(a2+b2-c2)/2;
   let SB=(c2+a2-b2)/2;
   let SA=(b2+c2-a2)/2;
   /* end vars */
   let v1 = a*(a2+b*c)*(a*(a+b+c)+2*SA)*SB*SC;
   let v2 = b*(b2+a*c)*SA*(b*(a+b+c)+2*SB)*SC;
   let v3 = c*(a*b+c2)*SA*SB*(c*(a+b+c)+2*SC);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X445(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let SC=(a2+b2-c2)/2;
   let SB=(c2+a2-b2)/2;
   let SA=(b2+c2-a2)/2;
   /* end vars */
   let v1 = (2*a*b*c+a2*(b+c)-(b-c)*(b-c)*(b+c))*(b*c+2*SA)*SB*SC;
   let v2 = (2*a*b*c+b2*(a+c)-(-a+c)*(-a+c)*(a+c))*SA*(a*c+2*SB)*SC;
   let v3 = (-((a-b)*(a-b)*(a+b))+2*a*b*c+(a+b)*c2)*SA*SB*(a*b+2*SC);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X446(orbit) {
   /* begin vars */
   let c2=c*c;
   let c4=c2*c2;
   let c6=c2*c4;
   let b2=b*b;
   let b4=b2*b2;
   let b6=b2*b4;
   let a2=a*a;
   let a4=a2*a2;
   let a6=a2*a4;
   let c8=c2*c6;
   let b8=b2*b6;
   let a8=a2*a6;
   /* end vars */
   let v1 = a2*(a2*b2-b4+a2*c2-c4)*(a8*b2-2*a6*b4+a4*b6+a8*c2+2*a2*b6*c2+b8*c2-2*a6*c4-4*a2*b4*c4-b6*c4+a4*c6+2*a2*b2*c6-b4*c6+b2*c8);
   let v2 = b2*(-a4+a2*b2+b2*c2-c4)*(a6*b4-2*a4*b6+a2*b8+a8*c2+2*a6*b2*c2+b8*c2-a6*c4-4*a4*b2*c4-2*b6*c4-a4*c6+2*a2*b2*c6+b4*c6+a2*c8);
   let v3 = c2*(-a4-b4+a2*c2+b2*c2)*(a8*b2-a6*b4-a4*b6+a2*b8+2*a6*b2*c2-4*a4*b4*c2+2*a2*b6*c2+a6*c4+b6*c4-2*a4*c6-2*b4*c6+a2*c8+b2*c8);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X447(orbit) {
   /* begin vars */
   let a2=a*a;
   let c2=c*c;
   let b2=b*b;
   let SA=(b2+c2-a2)/2;
   let SC=(a2+b2-c2)/2;
   let SB=(c2+a2-b2)/2;
   let c4=c2*c2;
   let c3=c2*c;
   let b4=b2*b2;
   let b3=b2*b;
   let a3=a2*a;
   let a4=a2*a2;
   /* end vars */
   let v1 = (a+b)*(a+c)*(a4-a3*b+a*b3-b4-a3*c+a2*b*c+a*b2*c-b3*c+a*b*c2+a*c3-b*c3-c4)*SB*SC;
   let v2 = (a+b)*(b+c)*(-a4+a3*b-a*b3+b4-a3*c+a2*b*c+a*b2*c-b3*c+a*b*c2-a*c3+b*c3-c4)*SA*SC;
   let v3 = (a+c)*(b+c)*(-a4-a3*b-a*b3-b4+a3*c+a2*b*c+a*b2*c+b3*c+a*b*c2-a*c3-b*c3+c4)*SA*SB;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X448(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   /* end vars */
   let v1 = -(b*(a+b)*(a+b-c)*c*(a+c)*(a-b+c)*(b+c)*(b+c))+a2*(a+b)*(a+b)*(a+c)*(a+c)*(-a+b+c)*(-a+b+c);
   let v2 = (a+b)*(a+b)*b2*(a-b+c)*(a-b+c)*(b+c)*(b+c)-a*(a+b)*(a+b-c)*c*(a+c)*(a+c)*(b+c)*(-a+b+c);
   let v3 = -(a*b*(a+b)*(a+b)*(a+c)*(a-b+c)*(b+c)*(-a+b+c))+(a+b-c)*(a+b-c)*(a+c)*(a+c)*(b+c)*(b+c)*c2;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X449(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let c3=c2*c;
   let b3=b2*b;
   let a3=a2*a;
   /* end vars */
   let v1 = (a-b-c)*(a-b-c)*(3*a3+3*a2*b+a*b2+b3+3*a2*c+2*a*b*c-b2*c+a*c2-b*c2+c3)*(3*a3+3*a2*b+a*b2+b3+3*a2*c+2*a*b*c-b2*c+a*c2-b*c2+c3)-(-a+b-c)*(-a-b+c)*(a3+a2*b+3*a*b2+3*b3-a2*c+2*a*b*c+3*b2*c-a*c2+b*c2+c3)*(a3-a2*b-a*b2+b3+a2*c+2*a*b*c+b2*c+3*a*c2+3*b*c2+3*c3);
   let v2 = (-a+b-c)*(-a+b-c)*(a3+a2*b+3*a*b2+3*b3-a2*c+2*a*b*c+3*b2*c-a*c2+b*c2+c3)*(a3+a2*b+3*a*b2+3*b3-a2*c+2*a*b*c+3*b2*c-a*c2+b*c2+c3)-(a-b-c)*(-a-b+c)*(3*a3+3*a2*b+a*b2+b3+3*a2*c+2*a*b*c-b2*c+a*c2-b*c2+c3)*(a3-a2*b-a*b2+b3+a2*c+2*a*b*c+b2*c+3*a*c2+3*b*c2+3*c3);
   let v3 = -((a-b-c)*(-a+b-c)*(3*a3+3*a2*b+a*b2+b3+3*a2*c+2*a*b*c-b2*c+a*c2-b*c2+c3)*(a3+a2*b+3*a*b2+3*b3-a2*c+2*a*b*c+3*b2*c-a*c2+b*c2+c3))+(-a-b+c)*(-a-b+c)*(a3-a2*b-a*b2+b3+a2*c+2*a*b*c+b2*c+3*a*c2+3*b*c2+3*c3)*(a3-a2*b-a*b2+b3+a2*c+2*a*b*c+b2*c+3*a*c2+3*b*c2+3*c3);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X450(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let c4=c2*c2;
   let b4=b2*b2;
   let SC=(a2+b2-c2)/2;
   let SB=(c2+a2-b2)/2;
   let SA=(b2+c2-a2)/2;
   let a4=a2*a2;
   /* end vars */
   let v1 = -(a4*Math.pow(SA,4)*SB*SC)+b2*c2*Math.pow(SB,3)*Math.pow(SC,3);
   let v2 = -(b4*SA*Math.pow(SB,4)*SC)+a2*c2*Math.pow(SA,3)*Math.pow(SC,3);
   let v3 = a2*b2*Math.pow(SA,3)*Math.pow(SB,3)-c4*SA*SB*Math.pow(SC,4);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X451(orbit) {
   /* begin vars */
   let a2=a*a;
   let c2=c*c;
   let b2=b*b;
   let SA=(b2+c2-a2)/2;
   let SC=(a2+b2-c2)/2;
   let SB=(c2+a2-b2)/2;
   let c3=c2*c;
   let b3=b2*b;
   let a3=a2*a;
   /* end vars */
   let v1 = (-a3-a2*b+a*b2+b3-a2*c+a*b*c+b2*c+a*c2+b*c2+c3)*SB*SC;
   let v2 = (a3+a2*b-a*b2-b3+a2*c+a*b*c-b2*c+a*c2+b*c2+c3)*SA*SC;
   let v3 = (a3+a2*b+a*b2+b3+a2*c+a*b*c+b2*c-a*c2-b*c2-c3)*SA*SB;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X452(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let c3=c2*c;
   let b3=b2*b;
   let a3=a2*a;
   /* end vars */
   let v1 = (a-b-c)*(3*a3+3*a2*b+a*b2+b3+3*a2*c+2*a*b*c-b2*c+a*c2-b*c2+c3);
   let v2 = (-a+b-c)*(a3+a2*b+3*a*b2+3*b3-a2*c+2*a*b*c+3*b2*c-a*c2+b*c2+c3);
   let v3 = (-a-b+c)*(a3-a2*b-a*b2+b3+a2*c+2*a*b*c+b2*c+3*a*c2+3*b*c2+3*c3);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X453(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let SC=(a2+b2-c2)/2;
   let SB=(c2+a2-b2)/2;
   let SA=(b2+c2-a2)/2;
   /* end vars */
   let v1 = a*(a+b)*(a+c)*(-a+b+c)*Math.pow(-(a*SA)+b*SB+c*SC,2);
   let v2 = b*(a+b)*(a-b+c)*(b+c)*(a*SA-b*SB+c*SC)*(a*SA-b*SB+c*SC);
   let v3 = (a+b-c)*c*(a+c)*(b+c)*(a*SA+b*SB-c*SC)*(a*SA+b*SB-c*SC);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X454(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let SC=(a2+b2-c2)/2;
   let SB=(c2+a2-b2)/2;
   let SA=(b2+c2-a2)/2;
   let SC2=SC*SC;
   let SB2=SB*SB;
   let SA2=SA*SA;
   /* end vars */
   let v1 = a2*SA*Math.pow(a2*(SA2-SB*SC)-SA*(SB2+SC2),2);
   let v2 = b2*SB*Math.pow(b2*(SB2-SA*SC)-SB*(SA2+SC2),2);
   let v3 = c2*SC*Math.pow(-((SA2+SB2)*SC)+c2*(-(SA*SB)+SC2),2);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X455(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let SC=(a2+b2-c2)/2;
   let SB=(c2+a2-b2)/2;
   let SA=(b2+c2-a2)/2;
   let SC2=SC*SC;
   let SB2=SB*SB;
   let SA2=SA*SA;
   /* end vars */
   let v1 = a2*SB*SC*Math.pow(a2*(SA2-SB*SC)+SA*(SB2+SC2),2);
   let v2 = b2*SA*SC*Math.pow(b2*(SB2-SA*SC)+SB*(SA2+SC2),2);
   let v3 = c2*SA*SB*Math.pow((SA2+SB2)*SC+c2*(-(SA*SB)+SC2),2);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X456(orbit) {
   /* begin vars */
   let b2=b*b;
   let a2=a*a;
   let c2=c*c;
   let SB=(c2+a2-b2)/2;
   let SC=(a2+b2-c2)/2;
   let SA=(b2+c2-a2)/2;
   let SB2=SB*SB;
   let SC2=SC*SC;
   let SA2=SA*SA;
   /* end vars */
   let v1 = a2*SB*SC*Math.pow(SA2*(3*SA-5*SB)*SB+3*c2*SA*(SA-3*SB)*SC-c2*(5*SA+4*SB)*SC2,2);
   let v2 = b2*SA*SC*Math.pow(3*a2*SA*SB*(SB-3*SC)+SB2*(3*SB-5*SC)*SC-a2*SA2*(5*SB+4*SC),2);
   let v3 = c2*SA*SB*Math.pow(3*b2*SB*SC*(-3*SA+SC)-b2*SB2*(4*SA+5*SC)+SA*(-5*SA+3*SC)*SC2,2);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X457(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let SC=(a2+b2-c2)/2;
   let SB=(c2+a2-b2)/2;
   let SA=(b2+c2-a2)/2;
   let SC2=SC*SC;
   let SB2=SB*SB;
   let SA2=SA*SA;
   /* end vars */
   let v1 = a2*SB*SC*Math.pow(-3*Math.pow(SA,3)*SB+c2*SA*SC*(SB+5*SC)+SA2*(5*SB2-3*c2*SC)-4*c2*SB*SC2,2);
   let v2 = b2*SA*SC*Math.pow(-4*a2*SA2*SC-3*Math.pow(SB,3)*SC+a2*SA*SB*(5*SA+SC)+SB2*(-3*a2*SA+5*SC2),2);
   let v3 = c2*SA*SB*Math.pow(-4*b2*SA*SB2+b2*SB*(SA+5*SB)*SC-3*SA*Math.pow(SC,3)+(5*SA2-3*b2*SB)*SC2,2);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X458(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let SC=(a2+b2-c2)/2;
   let SB=(c2+a2-b2)/2;
   let SA=(b2+c2-a2)/2;
   let SC2=SC*SC;
   let SB2=SB*SB;
   let SA2=SA*SA;
   /* end vars */
   let v1 = SB*SC*(2*a2*SA+SA2+SB*SC);
   let v2 = SA*SC*(2*b2*SB+SB2+SA*SC);
   let v3 = SA*SB*(SA*SB+2*c2*SC+SC2);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X459(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let area=triAreaHeron(a,b,c);
   let S=2*area;
   let SC=(a2+b2-c2)/2;
   let SB=(c2+a2-b2)/2;
   let S2=S*S;
   let SA=(b2+c2-a2)/2;
   /* end vars */
   let v1 = 1/(SA*(S2-2*SB*SC));
   let v2 = 1/(SB*(S2-2*SA*SC));
   let v3 = 1/((S2-2*SA*SB)*SC);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X460(orbit) {
   /* begin vars */
   let a2=a*a;
   let c2=c*c;
   let b2=b*b;
   let SA=(b2+c2-a2)/2;
   let SC=(a2+b2-c2)/2;
   let SB=(c2+a2-b2)/2;
   let SA2=SA*SA;
   let SC2=SC*SC;
   let SB2=SB*SB;
   /* end vars */
   let v1 = SB*SC*(-(a2*SA)+SB2+SC2);
   let v2 = SA*SC*(SA2-b2*SB+SC2);
   let v3 = SA*SB*(SA2+SB2-c2*SC);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X461(orbit) {
   /* begin vars */
   let a2=a*a;
   let c2=c*c;
   let b2=b*b;
   let SA=(b2+c2-a2)/2;
   let SC=(a2+b2-c2)/2;
   let SB=(c2+a2-b2)/2;
   /* end vars */
   let v1 = (-a+b+c)*(3*a+b+c)*SB*SC;
   let v2 = (a-b+c)*(a+3*b+c)*SA*SC;
   let v3 = (a+b-c)*(a+b+3*c)*SA*SB;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X462(orbit) {
   /* begin vars */
   let a2=a*a;
   let c2=c*c;
   let b2=b*b;
   let SA=(b2+c2-a2)/2;
   let area=triAreaHeron(a,b,c);
   let SC=(a2+b2-c2)/2;
   let SB=(c2+a2-b2)/2;
   /* end vars */
   let v1 = SB*SC*(-4*area+a2*sqrt3);
   let v2 = SA*SC*(-4*area+b2*sqrt3);
   let v3 = SA*SB*(-4*area+c2*sqrt3);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X463(orbit) {
   /* begin vars */
   let a2=a*a;
   let c2=c*c;
   let b2=b*b;
   let SA=(b2+c2-a2)/2;
   let area=triAreaHeron(a,b,c);
   let SC=(a2+b2-c2)/2;
   let SB=(c2+a2-b2)/2;
   /* end vars */
   let v1 = SB*SC*(4*area+a2*sqrt3);
   let v2 = SA*SC*(4*area+b2*sqrt3);
   let v3 = SA*SB*(4*area+c2*sqrt3);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X464(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let SB=(c2+a2-b2)/2;
   let SA=(b2+c2-a2)/2;
   let SC=(a2+b2-c2)/2;
   let SB2=SB*SB;
   let SA2=SA*SA;
   let SC2=SC*SC;
   /* end vars */
   let v1 = SA*(a2*(a+c)*(b+c)+SC2);
   let v2 = ((a+b)*b2*(a+c)+SA2)*SB;
   let v3 = ((a+b)*(b+c)*c2+SB2)*SC;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X465(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let area=triAreaHeron(a,b,c);
   let SC=(a2+b2-c2)/2;
   let SB=(c2+a2-b2)/2;
   let SA=(b2+c2-a2)/2;
   /* end vars */
   let v1 = SA*(SB*SC+a2*area*sqrt3);
   let v2 = SB*(SA*SC+area*b2*sqrt3);
   let v3 = SC*(SA*SB+area*c2*sqrt3);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X466(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let area=triAreaHeron(a,b,c);
   let SC=(a2+b2-c2)/2;
   let SB=(c2+a2-b2)/2;
   let SA=(b2+c2-a2)/2;
   /* end vars */
   let v1 = SA*(SB*SC-a2*area*sqrt3);
   let v2 = SB*(SA*SC-area*b2*sqrt3);
   let v3 = SC*(SA*SB-area*c2*sqrt3);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X467(orbit) {
   /* begin vars */
   let a2=a*a;
   let c2=c*c;
   let b2=b*b;
   let SA=(b2+c2-a2)/2;
   let SC=(a2+b2-c2)/2;
   let SB=(c2+a2-b2)/2;
   /* end vars */
   let v1 = SB*SC*((a2-SA)*SA+SB*SC)*(a2*SA+2*SB*SC);
   let v2 = SA*SC*((b2-SB)*SB+SA*SC)*(b2*SB+2*SA*SC);
   let v3 = SA*SB*(2*SA*SB+c2*SC)*(SA*SB+(c2-SC)*SC);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X468(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   /* end vars */
   let v1 = (2*a2-b2-c2)*(a2+b2-c2)*(a2-b2+c2);
   let v2 = (a2+b2-c2)*(-a2+2*b2-c2)*(-a2+b2+c2);
   let v3 = (a2-b2+c2)*(-a2+b2+c2)*(-a2-b2+2*c2);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X469(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let SC=(a2+b2-c2)/2;
   let SB=(c2+a2-b2)/2;
   let SA=(b2+c2-a2)/2;
   /* end vars */
   let v1 = (a2+a*b+a*c+b*c+2*SA)*SB*SC;
   let v2 = SA*(a*b+b2+a*c+b*c+2*SB)*SC;
   let v3 = SA*SB*(a*b+a*c+b*c+c2+2*SC);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X470(orbit) {
   /* begin vars */
   let a2=a*a;
   let c2=c*c;
   let b2=b*b;
   let SA=(b2+c2-a2)/2;
   let area=triAreaHeron(a,b,c);
   let SC=(a2+b2-c2)/2;
   let SB=(c2+a2-b2)/2;
   /* end vars */
   let v1 = SB*SC*(2*area+SA*sqrt3);
   let v2 = SA*SC*(2*area+SB*sqrt3);
   let v3 = SA*SB*(2*area+SC*sqrt3);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X471(orbit) {
   /* begin vars */
   let a2=a*a;
   let c2=c*c;
   let b2=b*b;
   let SA=(b2+c2-a2)/2;
   let area=triAreaHeron(a,b,c);
   let SC=(a2+b2-c2)/2;
   let SB=(c2+a2-b2)/2;
   /* end vars */
   let v1 = SB*SC*(-2*area+SA*sqrt3);
   let v2 = SA*SC*(-2*area+SB*sqrt3);
   let v3 = SA*SB*(-2*area+SC*sqrt3);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X472(orbit) {
   /* begin vars */
   let a2=a*a;
   let c2=c*c;
   let b2=b*b;
   let area=triAreaHeron(a,b,c);
   let SA=(b2+c2-a2)/2;
   let SC=(a2+b2-c2)/2;
   let SB=(c2+a2-b2)/2;
   /* end vars */
   let v1 = SB*SC*(-SA+2*area*sqrt3);
   let v2 = SA*SC*(-SB+2*area*sqrt3);
   let v3 = SA*SB*(-SC+2*area*sqrt3);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X473(orbit) {
   /* begin vars */
   let a2=a*a;
   let c2=c*c;
   let b2=b*b;
   let area=triAreaHeron(a,b,c);
   let SA=(b2+c2-a2)/2;
   let SC=(a2+b2-c2)/2;
   let SB=(c2+a2-b2)/2;
   /* end vars */
   let v1 = SB*SC*(SA+2*area*sqrt3);
   let v2 = SA*SC*(SB+2*area*sqrt3);
   let v3 = SA*SB*(SC+2*area*sqrt3);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X474(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let SC=(a2+b2-c2)/2;
   let SB=(c2+a2-b2)/2;
   let SA=(b2+c2-a2)/2;
   /* end vars */
   let v1 = a*(b*c*(a+b+c)-a*SA);
   let v2 = b*(a*c*(a+b+c)-b*SB);
   let v3 = c*(a*b*(a+b+c)-c*SC);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X475(orbit) {
   /* begin vars */
   let a2=a*a;
   let c2=c*c;
   let b2=b*b;
   let SA=(b2+c2-a2)/2;
   let SC=(a2+b2-c2)/2;
   let SB=(c2+a2-b2)/2;
   let c3=c2*c;
   let b3=b2*b;
   let a3=a2*a;
   /* end vars */
   let v1 = (a3+a2*b-a*b2-b3+a2*c+2*a*b*c-b2*c-a*c2-b*c2-c3)*SB*SC;
   let v2 = (-a3-a2*b+a*b2+b3-a2*c+2*a*b*c+b2*c-a*c2-b*c2-c3)*SA*SC;
   let v3 = (-a3-a2*b-a*b2-b3-a2*c+2*a*b*c-b2*c+a*c2+b*c2+c3)*SA*SB;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X476(orbit) {
   /* begin vars */
   let a2=a*a;
   let b2=b*b;
   let b4=b2*b2;
   let c2=c*c;
   let a4=a2*a2;
   let c4=c2*c2;
   let a3=a2*a;
   let b6=b2*b4;
   let c3=c2*c;
   let a6=a2*a4;
   let c6=c2*c4;
   let b3=b2*b;
   /* end vars */
   let v1 = 1/(Math.pow(-(a2*b)+b3,2)-a4*c2+2*a2*c4-c6);
   let v2 = 1/(-a6+2*a4*b2-a2*b4+Math.pow(-(b2*c)+c3,2));
   let v3 = 1/(-b6+2*b4*c2+(a3-a*c2)*(a3-a*c2)-b2*c4);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X477(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let SC=(a2+b2-c2)/2;
   let SB=(c2+a2-b2)/2;
   let SA=(b2+c2-a2)/2;
   let SC2=SC*SC;
   let SB2=SB*SB;
   let SA2=SA*SA;
   /* end vars */
   let v1 = 1/(a2*SA*(SA2-SB*SC)+SA2*(-3*SB2+4*SB*SC-3*SC2)+2*SB2*SC2);
   let v2 = 1/(b2*SB*(SB2-SA*SC)+SB2*(-3*SA2+4*SA*SC-3*SC2)+2*SA2*SC2);
   let v3 = 1/(2*SA2*SB2+(-3*SA2+4*SA*SB-3*SB2)*SC2+c2*SC*(-(SA*SB)+SC2));
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X478(orbit) {
   /* begin vars */
   let a2=a*a;
   let c2=c*c;
   let b2=b*b;
   let SA=(b2+c2-a2)/2;
   let SC=(a2+b2-c2)/2;
   let SB=(c2+a2-b2)/2;
   /* end vars */
   let v1 = a2*(a+b-c)*(a-b+c)*(a*b*(a-b-c)*c+2*SB*SC);
   let v2 = b2*(a+b-c)*(-a+b+c)*(a*b*(-a+b-c)*c+2*SA*SC);
   let v3 = (a-b+c)*(-a+b+c)*c2*(a*b*c*(-a-b+c)+2*SA*SB);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X479(orbit) {
   /* begin vars */

   /* end vars */
   let v1 = Math.pow(-a+b+c,-3);
   let v2 = Math.pow(a-b+c,-3);
   let v3 = Math.pow(a+b-c,-3);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X480(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   /* end vars */
   let v1 = a2*Math.pow(-a+b+c,3);
   let v2 = b2*Math.pow(a-b+c,3);
   let v3 = Math.pow(a+b-c,3)*c2;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X481(orbit) {
   /* begin vars */
   let area=triAreaHeron(a,b,c);
   /* end vars */
   let v1 = a-(4*area)/(-a+b+c);
   let v2 = b-(4*area)/(a-b+c);
   let v3 = (-4*area)/(a+b-c)+c;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X482(orbit) {
   /* begin vars */
   let area=triAreaHeron(a,b,c);
   /* end vars */
   let v1 = a+(4*area)/(-a+b+c);
   let v2 = b+(4*area)/(a-b+c);
   let v3 = (4*area)/(a+b-c)+c;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X483(orbit) {
   /* begin vars */
   let cosC=lawOfCosines(c,a,b);
   let cosHalfC=cosHalfAngle(cosC);
   let cosQuarterC=cosHalfAngle(cosHalfC);
   let cosB=lawOfCosines(b,a,c);
   let cosHalfB=cosHalfAngle(cosB);
   let cosQuarterB=cosHalfAngle(cosHalfB);
   let cosA=lawOfCosines(a,b,c);
   let cosHalfA=cosHalfAngle(cosA);
   let cosQuarterA=cosHalfAngle(cosHalfA);
   let sinC=getSin(cosC);
   let secQuarterC=1/cosQuarterC;
   let sinB=getSin(cosB);
   let secQuarterB=1/cosQuarterB;
   let sinA=getSin(cosA);
   let secQuarterA=1/cosQuarterA;
   /* end vars */
   let v1 = secQuarterA*secQuarterA*sinA;
   let v2 = secQuarterB*secQuarterB*sinB;
   let v3 = secQuarterC*secQuarterC*sinC;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X484(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let c3=c2*c;
   let b3=b2*b;
   let a3=a2*a;
   /* end vars */
   let v1 = a*(a3+a2*b-a*b2-b3+a2*c-a*b*c+b2*c-a*c2+b*c2-c3);
   let v2 = b*(-a3-a2*b+a*b2+b3+a2*c-a*b*c+b2*c+a*c2-b*c2-c3);
   let v3 = c*(-a3+a2*b+a*b2-b3-a2*c-a*b*c-b2*c+a*c2+b*c2+c3);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X485(orbit) {
   /* begin vars */
   let a2=a*a;
   let c2=c*c;
   let b2=b*b;
   let SA=(b2+c2-a2)/2;
   let SC=(a2+b2-c2)/2;
   let SB=(c2+a2-b2)/2;
   let area=triAreaHeron(a,b,c);
   /* end vars */
   let v1 = (2*area+SB)*(2*area+SC);
   let v2 = (2*area+SA)*(2*area+SC);
   let v3 = (2*area+SA)*(2*area+SB);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X486(orbit) {
   /* begin vars */
   let a2=a*a;
   let c2=c*c;
   let b2=b*b;
   let SA=(b2+c2-a2)/2;
   let SC=(a2+b2-c2)/2;
   let SB=(c2+a2-b2)/2;
   let area=triAreaHeron(a,b,c);
   /* end vars */
   let v1 = (-2*area+SB)*(-2*area+SC);
   let v2 = (-2*area+SA)*(-2*area+SC);
   let v3 = (-2*area+SA)*(-2*area+SB);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X487(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let SC=(a2+b2-c2)/2;
   let SB=(c2+a2-b2)/2;
   let area=triAreaHeron(a,b,c);
   let SA=(b2+c2-a2)/2;
   /* end vars */
   let v1 = a2*SA-4*area*(-area+SA)-SB*SC;
   let v2 = b2*SB-4*area*(-area+SB)-SA*SC;
   let v3 = -(SA*SB)+c2*SC-4*area*(-area+SC);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X488(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let SC=(a2+b2-c2)/2;
   let SB=(c2+a2-b2)/2;
   let area=triAreaHeron(a,b,c);
   let SA=(b2+c2-a2)/2;
   /* end vars */
   let v1 = a2*SA+4*area*(area+SA)-SB*SC;
   let v2 = b2*SB+4*area*(area+SB)-SA*SC;
   let v3 = -(SA*SB)+c2*SC+4*area*(area+SC);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X489(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let SC=(a2+b2-c2)/2;
   let SB=(c2+a2-b2)/2;
   let area=triAreaHeron(a,b,c);
   let SA=(b2+c2-a2)/2;
   /* end vars */
   let v1 = a2*SA-2*area*SA-SB*SC;
   let v2 = -2*area*SB+b2*SB-SA*SC;
   let v3 = -(SA*SB)-2*area*SC+c2*SC;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X490(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let SC=(a2+b2-c2)/2;
   let SB=(c2+a2-b2)/2;
   let area=triAreaHeron(a,b,c);
   let SA=(b2+c2-a2)/2;
   /* end vars */
   let v1 = a2*SA+2*area*SA-SB*SC;
   let v2 = 2*area*SB+b2*SB-SA*SC;
   let v3 = -(SA*SB)+2*area*SC+c2*SC;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X491(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let SC=(a2+b2-c2)/2;
   let SB=(c2+a2-b2)/2;
   let SA=(b2+c2-a2)/2;
   let area=triAreaHeron(a,b,c);
   /* end vars */
   let v1 = -2*area+SA;
   let v2 = -2*area+SB;
   let v3 = -2*area+SC;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X492(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let SC=(a2+b2-c2)/2;
   let SB=(c2+a2-b2)/2;
   let SA=(b2+c2-a2)/2;
   let area=triAreaHeron(a,b,c);
   /* end vars */
   let v1 = 2*area+SA;
   let v2 = 2*area+SB;
   let v3 = 2*area+SC;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X493(orbit) {
   /* begin vars */
   let area=triAreaHeron(a,b,c);
   let c2=c*c;
   let S=2*area;
   let b2=b*b;
   let a2=a*a;
   /* end vars */
   let v1 = a2*(b2+S)*(c2+S);
   let v2 = b2*(a2+S)*(c2+S);
   let v3 = c2*(a2+S)*(b2+S);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X494(orbit) {
   /* begin vars */
   let area=triAreaHeron(a,b,c);
   let c2=c*c;
   let S=2*area;
   let b2=b*b;
   let a2=a*a;
   /* end vars */
   let v1 = a2*(b2-S)*(c2-S);
   let v2 = b2*(a2-S)*(c2-S);
   let v3 = c2*(a2-S)*(b2-S);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X495(orbit) {
   /* begin vars */
   let a2=a*a;
   let c2=c*c;
   let b2=b*b;
   let a4=a2*a2;
   let c4=c2*c2;
   let b4=b2*b2;
   /* end vars */
   let v1 = a2*b2-b4+4*a2*b*c+a2*c2+2*b2*c2-c4;
   let v2 = -a4+a2*b2+4*a*b2*c+2*a2*c2+b2*c2-c4;
   let v3 = -a4+2*a2*b2-b4+a2*c2+4*a*b*c2+b2*c2;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X496(orbit) {
   /* begin vars */
   let a2=a*a;
   let c2=c*c;
   let b2=b*b;
   let a4=a2*a2;
   let c4=c2*c2;
   let b4=b2*b2;
   /* end vars */
   let v1 = a2*b2-b4-4*a2*b*c+a2*c2+2*b2*c2-c4;
   let v2 = -a4+a2*b2-4*a*b2*c+2*a2*c2+b2*c2-c4;
   let v3 = -a4+2*a2*b2-b4+a2*c2-4*a*b*c2+b2*c2;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X497(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   /* end vars */
   let v1 = (a-b-c)*(a2+b2-2*b*c+c2);
   let v2 = (-a+b-c)*(a2+b2-2*a*c+c2);
   let v3 = (-a-b+c)*(a2-2*a*b+b2+c2);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X498(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let c4=c2*c2;
   let b4=b2*b2;
   let a4=a2*a2;
   /* end vars */
   let v1 = a4-2*a2*b2+b4-2*a2*b*c-2*a2*c2-2*b2*c2+c4;
   let v2 = a4-2*a2*b2+b4-2*a*b2*c-2*a2*c2-2*b2*c2+c4;
   let v3 = a4-2*a2*b2+b4-2*a2*c2-2*a*b*c2-2*b2*c2+c4;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X499(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let c4=c2*c2;
   let b4=b2*b2;
   let a4=a2*a2;
   /* end vars */
   let v1 = a4-2*a2*b2+b4+2*a2*b*c-2*a2*c2-2*b2*c2+c4;
   let v2 = a4-2*a2*b2+b4+2*a*b2*c-2*a2*c2-2*b2*c2+c4;
   let v3 = a4-2*a2*b2+b4-2*a2*c2+2*a*b*c2-2*b2*c2+c4;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X500(orbit) {
   /* begin vars */
   let a2=a*a;
   let c2=c*c;
   let b2=b*b;
   let a3=a2*a;
   let c3=c2*c;
   let b3=b2*b;
   /* end vars */
   let v1 = a2*(a2-b2-b*c-c2)*(a2*b-b3+a2*c+2*a*b*c+b2*c+b*c2-c3);
   let v2 = b2*(-a2+b2-a*c-c2)*(-a3+a*b2+a2*c+2*a*b*c+b2*c+a*c2-c3);
   let v3 = c2*(-a2-a*b-b2+c2)*(-a3+a2*b+a*b2-b3+2*a*b*c+a*c2+b*c2);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X501(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let c3=c2*c;
   let b3=b2*b;
   let a3=a2*a;
   /* end vars */
   let v1 = a2*(a+b)*(a+c)*(a3+a2*b-a*b2-b3+a2*c-a*b*c-b2*c-a*c2-b*c2-c3);
   let v2 = (a+b)*b2*(b+c)*(-a3-a2*b+a*b2+b3-a2*c-a*b*c+b2*c-a*c2-b*c2-c3);
   let v3 = (a+c)*(b+c)*c2*(-a3-a2*b-a*b2-b3-a2*c-a*b*c-b2*c+a*c2+b*c2+c3);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X502(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let c3=c2*c;
   let b3=b2*b;
   let a3=a2*a;
   /* end vars */
   let v1 = 1/((a+b)*(a+c)*(a3+a2*b-a*b2-b3+a2*c-a*b*c-b2*c-a*c2-b*c2-c3));
   let v2 = 1/((a+b)*(b+c)*(-a3-a2*b+a*b2+b3-a2*c-a*b*c+b2*c-a*c2-b*c2-c3));
   let v3 = 1/((a+c)*(b+c)*(-a3-a2*b-a*b2-b3-a2*c-a*b*c-b2*c+a*c2+b*c2+c3));
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X503(orbit) {
   /* begin vars */
   let sc=(a+b-c)/2;
   let sb=(c+a-b)/2;
   let sa=(b+c-a)/2;
   let s=(a+b+c)/2;
   let Sqrt=Math.sqrt;
   /* end vars */
   let v1 = a*(-Sqrt((b*c)/(s*sa))+Sqrt((a*c)/(s*sb))+Sqrt((a*b)/(s*sc)));
   let v2 = b*(Sqrt((b*c)/(s*sa))-Sqrt((a*c)/(s*sb))+Sqrt((a*b)/(s*sc)));
   let v3 = c*(Sqrt((b*c)/(s*sa))+Sqrt((a*c)/(s*sb))-Sqrt((a*b)/(s*sc)));
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X504(orbit) {
   /* begin vars */
   let sc=(a+b-c)/2;
   let sb=(c+a-b)/2;
   let sa=(b+c-a)/2;
   let Sqrt=Math.sqrt;
   /* end vars */
   let v1 = a*(c*Sqrt((sa*sb)/(a*b))+b*Sqrt((sa*sc)/(a*c))-a*Sqrt((sb*sc)/(b*c)));
   let v2 = b*(c*Sqrt((sa*sb)/(a*b))-b*Sqrt((sa*sc)/(a*c))+a*Sqrt((sb*sc)/(b*c)));
   let v3 = c*(-(c*Sqrt((sa*sb)/(a*b)))+b*Sqrt((sa*sc)/(a*c))+a*Sqrt((sb*sc)/(b*c)));
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X505(orbit) {
   /* begin vars */
   let sc=(a+b-c)/2;
   let sb=(c+a-b)/2;
   let sa=(b+c-a)/2;
   let Sqrt=Math.sqrt;
   /* end vars */
   let v1 = a/(Sqrt((sa*sb)/(a*b))+Sqrt((sa*sc)/(a*c))-Sqrt((sb*sc)/(b*c)));
   let v2 = b/(Sqrt((sa*sb)/(a*b))-Sqrt((sa*sc)/(a*c))+Sqrt((sb*sc)/(b*c)));
   let v3 = c/(-Sqrt((sa*sb)/(a*b))+Sqrt((sa*sc)/(a*c))+Sqrt((sb*sc)/(b*c)));
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X506(orbit) {
   /* begin vars */
   let sc=(a+b-c)/2;
   let sb=(c+a-b)/2;
   let sa=(b+c-a)/2;
   let s=(a+b+c)/2;
   /* end vars */
   let v1 = Math.pow(a,0.6666666666666666)/(Math.pow(s,0.3333333333333333)*Math.pow(sa,0.3333333333333333));
   let v2 = Math.pow(b,0.6666666666666666)/(Math.pow(s,0.3333333333333333)*Math.pow(sb,0.3333333333333333));
   let v3 = Math.pow(c,0.6666666666666666)/(Math.pow(s,0.3333333333333333)*Math.pow(sc,0.3333333333333333));
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X507(orbit) {
   /* begin vars */
   let sc=(a+b-c)/2;
   let sb=(c+a-b)/2;
   let sa=(b+c-a)/2;
   let s=(a+b+c)/2;
   /* end vars */
   let v1 = Math.pow(a,0.75)/(Math.pow(s,0.25)*Math.pow(sa,0.25));
   let v2 = Math.pow(b,0.75)/(Math.pow(s,0.25)*Math.pow(sb,0.25));
   let v3 = Math.pow(c,0.75)/(Math.pow(s,0.25)*Math.pow(sc,0.25));
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X508(orbit) {
   /* begin vars */
   let sc=(a+b-c)/2;
   let sb=(c+a-b)/2;
   let sa=(b+c-a)/2;
   let s=(a+b+c)/2;
   let Sqrt=Math.sqrt;
   /* end vars */
   let v1 = Sqrt(1/(s*sa));
   let v2 = Sqrt(1/(s*sb));
   let v3 = Sqrt(1/(s*sc));
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X509(orbit) {
   /* begin vars */
   let sc=(a+b-c)/2;
   let sb=(c+a-b)/2;
   let sa=(b+c-a)/2;
   let s=(a+b+c)/2;
   let Sqrt=Math.sqrt;
   /* end vars */
   let v1 = a*Sqrt(1/(s*sa));
   let v2 = b*Sqrt(1/(s*sb));
   let v3 = c*Sqrt(1/(s*sc));
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X510(orbit) {
   /* begin vars */

   /* end vars */
   let v1 = a*(-Math.pow(a,1.5)+Math.pow(b,1.5)+Math.pow(c,1.5));
   let v2 = b*(Math.pow(a,1.5)-Math.pow(b,1.5)+Math.pow(c,1.5));
   let v3 = c*(Math.pow(a,1.5)+Math.pow(b,1.5)-Math.pow(c,1.5));
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X511(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let SC=(a2+b2-c2)/2;
   let SB=(c2+a2-b2)/2;
   let SA=(b2+c2-a2)/2;
   let SC2=SC*SC;
   let SB2=SB*SB;
   let SA2=SA*SA;
   /* end vars */
   let v1 = a2*(-SA2+SB*SC);
   let v2 = b2*(-SB2+SA*SC);
   let v3 = c2*(SA*SB-SC2);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X512(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   /* end vars */
   let v1 = a2*(b2-c2);
   let v2 = b2*(-a2+c2);
   let v3 = (a2-b2)*c2;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X513(orbit) {
   /* begin vars */

   /* end vars */
   let v1 = a*(b-c);
   let v2 = b*(-a+c);
   let v3 = (a-b)*c;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X514(orbit) {
   /* begin vars */

   /* end vars */
   let v1 = b-c;
   let v2 = -a+c;
   let v3 = a-b;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X515(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let c3=c2*c;
   let b3=b2*b;
   let SC=(a2+b2-c2)/2;
   let SB=(c2+a2-b2)/2;
   let SA=(b2+c2-a2)/2;
   let a3=a2*a;
   /* end vars */
   let v1 = a3*SA-(b+c)*SB*SC;
   let v2 = b3*SB-(a+c)*SA*SC;
   let v3 = -((a+b)*SA*SB)+c3*SC;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X516(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let c3=c2*c;
   let b3=b2*b;
   let a3=a2*a;
   /* end vars */
   let v1 = 2*a3-a2*(b+c)-(b-c)*(b-c)*(b+c);
   let v2 = 2*b3-b2*(a+c)-(-a+c)*(-a+c)*(a+c);
   let v3 = -((a-b)*(a-b)*(a+b))-(a+b)*c2+2*c3;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X517(orbit) {
   /* begin vars */
   let a2=a*a;
   let c2=c*c;
   let b2=b*b;
   let SA=(b2+c2-a2)/2;
   let SC=(a2+b2-c2)/2;
   let SB=(c2+a2-b2)/2;
   /* end vars */
   let v1 = a*(a*b*c-b*SB-c*SC);
   let v2 = b*(a*b*c-a*SA-c*SC);
   let v3 = c*(a*b*c-a*SA-b*SB);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X518(orbit) {
   /* begin vars */
   let a2=a*a;
   let c2=c*c;
   let b2=b*b;
   /* end vars */
   let v1 = a*(a*b-b2+a*c-c2);
   let v2 = b*(-a2+a*b+b*c-c2);
   let v3 = c*(-a2-b2+a*c+b*c);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X519(orbit) {
   /* begin vars */

   /* end vars */
   let v1 = 2*a-b-c;
   let v2 = -a+2*b-c;
   let v3 = -a-b+2*c;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X520(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let SC=(a2+b2-c2)/2;
   let SB=(c2+a2-b2)/2;
   let SA=(b2+c2-a2)/2;
   let SC2=SC*SC;
   let SB2=SB*SB;
   let SA2=SA*SA;
   /* end vars */
   let v1 = a2*SA2*(-SB+SC);
   let v2 = b2*SB2*(SA-SC);
   let v3 = c2*(-SA+SB)*SC2;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X521(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let SC=(a2+b2-c2)/2;
   let SB=(c2+a2-b2)/2;
   let SA=(b2+c2-a2)/2;
   /* end vars */
   let v1 = a*(a-b-c)*(b-c)*SA;
   let v2 = b*(-a+b-c)*(-a+c)*SB;
   let v3 = (a-b)*c*(-a-b+c)*SC;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X522(orbit) {
   /* begin vars */

   /* end vars */
   let v1 = (a-b-c)*(b-c);
   let v2 = (-a+b-c)*(-a+c);
   let v3 = (a-b)*(-a-b+c);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X523(orbit) {
   /* begin vars */
   let a2=a*a;
   let c2=c*c;
   let b2=b*b;
   /* end vars */
   let v1 = b2-c2;
   let v2 = -a2+c2;
   let v3 = a2-b2;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X524(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   /* end vars */
   let v1 = 2*a2-b2-c2;
   let v2 = -a2+2*b2-c2;
   let v3 = -a2-b2+2*c2;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X525(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let SC=(a2+b2-c2)/2;
   let SB=(c2+a2-b2)/2;
   let SA=(b2+c2-a2)/2;
   /* end vars */
   let v1 = (b2-c2)*SA;
   let v2 = (-a2+c2)*SB;
   let v3 = (a2-b2)*SC;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X526(orbit) {
   /* begin vars */
   let a2=a*a;
   let b2=b*b;
   let b4=b2*b2;
   let c2=c*c;
   let a4=a2*a2;
   let c4=c2*c2;
   let a3=a2*a;
   let b6=b2*b4;
   let c3=c2*c;
   let a6=a2*a4;
   let c6=c2*c4;
   let b3=b2*b;
   /* end vars */
   let v1 = a2*(Math.pow(-(a2*b)+b3,2)-a4*c2+2*a2*c4-c6);
   let v2 = b2*(-a6+2*a4*b2-a2*b4+Math.pow(-(b2*c)+c3,2));
   let v3 = c2*(-b6+2*b4*c2+(a3-a*c2)*(a3-a*c2)-b2*c4);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X527(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let SC=(a2+b2-c2)/2;
   let SB=(c2+a2-b2)/2;
   let SA=(b2+c2-a2)/2;
   /* end vars */
   let v1 = a2+2*b*c-a*(b+c)-2*SA;
   let v2 = b2+2*a*c-b*(a+c)-2*SB;
   let v3 = 2*a*b-(a+b)*c+c2-2*SC;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X528(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let c3=c2*c;
   let b3=b2*b;
   let a3=a2*a;
   /* end vars */
   let v1 = 2*a3-2*a2*b+a*b2-b3-2*a2*c+b2*c+a*c2+b*c2-c3;
   let v2 = -a3+a2*b-2*a*b2+2*b3+a2*c-2*b2*c+a*c2+b*c2-c3;
   let v3 = -a3+a2*b+a*b2-b3+a2*c+b2*c-2*a*c2-2*b*c2+2*c3;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X529(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let SC=(a2+b2-c2)/2;
   let SB=(c2+a2-b2)/2;
   let SA=(b2+c2-a2)/2;
   /* end vars */
   let v1 = a*(-(b*c*(-2*a+b+c))-a*SA)+2*SB*SC;
   let v2 = b*(-(a*c*(a-2*b+c))-b*SB)+2*SA*SC;
   let v3 = 2*SA*SB+c*(-(a*b*(a+b-2*c))-c*SC);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X530(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let area=triAreaHeron(a,b,c);
   let SC=(a2+b2-c2)/2;
   let SB=(c2+a2-b2)/2;
   let SA=(b2+c2-a2)/2;
   /* end vars */
   let v1 = 3*(a2*SA-2*SB*SC)+2*area*(-a2+2*SA)*sqrt3;
   let v2 = 3*(b2*SB-2*SA*SC)+2*area*(-b2+2*SB)*sqrt3;
   let v3 = 3*(-2*SA*SB+c2*SC)+2*area*(-c2+2*SC)*sqrt3;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X531(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let area=triAreaHeron(a,b,c);
   let SC=(a2+b2-c2)/2;
   let SB=(c2+a2-b2)/2;
   let SA=(b2+c2-a2)/2;
   /* end vars */
   let v1 = 3*(a2*SA-2*SB*SC)-2*area*(-a2+2*SA)*sqrt3;
   let v2 = 3*(b2*SB-2*SA*SC)-2*area*(-b2+2*SB)*sqrt3;
   let v3 = 3*(-2*SA*SB+c2*SC)-2*area*(-c2+2*SC)*sqrt3;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X532(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let area=triAreaHeron(a,b,c);
   let SC=(a2+b2-c2)/2;
   let SB=(c2+a2-b2)/2;
   let SA=(b2+c2-a2)/2;
   /* end vars */
   let v1 = a2*SA-2*SB*SC+2*area*(-a2+2*SA)*sqrt3;
   let v2 = b2*SB-2*SA*SC+2*area*(-b2+2*SB)*sqrt3;
   let v3 = -2*SA*SB+c2*SC+2*area*(-c2+2*SC)*sqrt3;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X533(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let area=triAreaHeron(a,b,c);
   let SC=(a2+b2-c2)/2;
   let SB=(c2+a2-b2)/2;
   let SA=(b2+c2-a2)/2;
   /* end vars */
   let v1 = a2*SA-2*SB*SC-2*area*(-a2+2*SA)*sqrt3;
   let v2 = b2*SB-2*SA*SC-2*area*(-b2+2*SB)*sqrt3;
   let v3 = -2*SA*SB+c2*SC-2*area*(-c2+2*SC)*sqrt3;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X534(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let c3=c2*c;
   let b3=b2*b;
   let a3=a2*a;
   let SC=(a2+b2-c2)/2;
   let SB=(c2+a2-b2)/2;
   let SA=(b2+c2-a2)/2;
   let c5=c2*c3;
   let b5=b2*b3;
   let a5=a2*a3;
   /* end vars */
   let v1 = -2*a5+b5+c5-c*(SA-SB)*(SA-SB)-b*(SA-SC)*(SA-SC)+2*a*(SB-SC)*(SB-SC);
   let v2 = a5-2*b5+c5-c*(-SA+SB)*(-SA+SB)-a*(SB-SC)*(SB-SC)+2*b*(-SA+SC)*(-SA+SC);
   let v3 = a5+b5-2*c5+2*c*(SA-SB)*(SA-SB)-b*(-SA+SC)*(-SA+SC)-a*(-SB+SC)*(-SB+SC);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X535(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let SC=(a2+b2-c2)/2;
   let SB=(c2+a2-b2)/2;
   let SA=(b2+c2-a2)/2;
   /* end vars */
   let v1 = a*(b*c*(-2*a+b+c)+2*a*SA)-4*SB*SC;
   let v2 = b*(a*c*(a-2*b+c)+2*b*SB)-4*SA*SC;
   let v3 = -4*SA*SB+c*(a*b*(a+b-2*c)+2*c*SC);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X536(orbit) {
   /* begin vars */

   /* end vars */
   let v1 = a*b+a*c-2*b*c;
   let v2 = a*b-2*a*c+b*c;
   let v3 = -2*a*b+a*c+b*c;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X537(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let SC=(a2+b2-c2)/2;
   let SB=(c2+a2-b2)/2;
   let SA=(b2+c2-a2)/2;
   let c3=c2*c;
   let b3=b2*b;
   let a3=a2*a;
   /* end vars */
   let v1 = 2*a3-b3-c3+4*a*SA-2*b*SB-2*c*SC;
   let v2 = -a3+2*b3-c3-2*a*SA+4*b*SB-2*c*SC;
   let v3 = -a3-b3+2*c3-2*a*SA-2*b*SB+4*c*SC;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X538(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let SC=(a2+b2-c2)/2;
   let SB=(c2+a2-b2)/2;
   let SA=(b2+c2-a2)/2;
   let SC2=SC*SC;
   let SB2=SB*SB;
   let SA2=SA*SA;
   /* end vars */
   let v1 = -2*SA2+SB2+SC2;
   let v2 = SA2-2*SB2+SC2;
   let v3 = SA2+SB2-2*SC2;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X539(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let SC=(a2+b2-c2)/2;
   let SB=(c2+a2-b2)/2;
   let SA=(b2+c2-a2)/2;
   let c4=c2*c2;
   let b4=b2*b2;
   let SC2=SC*SC;
   let SB2=SB*SB;
   let SA2=SA*SA;
   let a4=a2*a2;
   /* end vars */
   let v1 = SA*(3*a4*SA2+SB*SC*(-3*SB2+2*SB*SC-3*SC2)-a2*SA*(3*SB2-2*SB*SC+3*SC2));
   let v2 = SB*(3*b4*SB2+SA*SC*(-3*SA2+2*SA*SC-3*SC2)-b2*SB*(3*SA2-2*SA*SC+3*SC2));
   let v3 = SC*(SA*SB*(-3*SA2+2*SA*SB-3*SB2)-c2*(3*SA2-2*SA*SB+3*SB2)*SC+3*c4*SC2);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X540(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let c4=c2*c2;
   let c3=c2*c;
   let b4=b2*b2;
   let b3=b2*b;
   let a3=a2*a;
   let a4=a2*a2;
   /* end vars */
   let v1 = -2*a4-2*a3*b+a*b3+b4-2*a3*c-2*a2*b*c+a*b2*c+b3*c+a*b*c2+a*c3+b*c3+c4;
   let v2 = a4+a3*b-2*a*b3-2*b4+a3*c+a2*b*c-2*a*b2*c-2*b3*c+a*b*c2+a*c3+b*c3+c4;
   let v3 = a4+a3*b+a*b3+b4+a3*c+a2*b*c+a*b2*c+b3*c-2*a*b*c2-2*a*c3-2*b*c3-2*c4;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X541(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let SA=(b2+c2-a2)/2;
   let SC=(a2+b2-c2)/2;
   let SB=(c2+a2-b2)/2;
   let c4=c2*c2;
   let b4=b2*b2;
   let SA2=SA*SA;
   let SC2=SC*SC;
   let SB2=SB*SB;
   let a4=a2*a2;
   /* end vars */
   let v1 = -(a4*Math.pow(SA,3))+SA*SB*SC*(SB2+10*SB*SC+SC2)+a2*(5*SA2*(SB-SC)*(SB-SC)-4*SB2*SC2);
   let v2 = -(b4*Math.pow(SB,3))+SA*SB*SC*(SA2+10*SA*SC+SC2)+b2*(5*SB2*(-SA+SC)*(-SA+SC)-4*SA2*SC2);
   let v3 = SA*SB*(SA2+10*SA*SB+SB2)*SC-c4*Math.pow(SC,3)+c2*(-4*SA2*SB2+5*(SA-SB)*(SA-SB)*SC2);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X542(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let SC=(a2+b2-c2)/2;
   let SB=(c2+a2-b2)/2;
   let SA=(b2+c2-a2)/2;
   let SC2=SC*SC;
   let SB2=SB*SB;
   let SA2=SA*SA;
   /* end vars */
   let v1 = -(a2*(SA2+SB*SC))+2*SA*(SB2+SC2);
   let v2 = -(b2*(SB2+SA*SC))+2*SB*(SA2+SC2);
   let v3 = 2*(SA2+SB2)*SC-c2*(SA*SB+SC2);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X543(orbit) {
   /* begin vars */
   let a2=a*a;
   let c2=c*c;
   let b2=b*b;
   let SA=(b2+c2-a2)/2;
   let SC=(a2+b2-c2)/2;
   let SB=(c2+a2-b2)/2;
   let SA2=SA*SA;
   let SC2=SC*SC;
   let SB2=SB*SB;
   /* end vars */
   let v1 = 2*(a2-SA)*SA+SB2-4*SB*SC+SC2;
   let v2 = SA2+2*(b2-SB)*SB-4*SA*SC+SC2;
   let v3 = SA2-4*SA*SB+SB2+2*(c2-SC)*SC;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X544(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let c4=c2*c2;
   let c3=c2*c;
   let b4=b2*b2;
   let b3=b2*b;
   let a3=a2*a;
   let a4=a2*a2;
   /* end vars */
   let v1 = 2*a4-2*a3*b+a*b3-b4-2*a3*c+2*a2*b*c-a*b2*c+b3*c-a*b*c2+a*c3+b*c3-c4;
   let v2 = -a4+a3*b-2*a*b3+2*b4+a3*c-a2*b*c+2*a*b2*c-2*b3*c-a*b*c2+a*c3+b*c3-c4;
   let v3 = -a4+a3*b+a*b3-b4+a3*c-a2*b*c-a*b2*c+b3*c+2*a*b*c2-2*a*c3-2*b*c3+2*c4;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X545(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let SC=(a2+b2-c2)/2;
   let SB=(c2+a2-b2)/2;
   let SA=(b2+c2-a2)/2;
   /* end vars */
   let v1 = (a-2*b)*(a-2*c)-2*SA;
   let v2 = (-2*a+b)*(b-2*c)-2*SB;
   let v3 = (-2*a+c)*(-2*b+c)-2*SC;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X546(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let SC=(a2+b2-c2)/2;
   let SB=(c2+a2-b2)/2;
   let SA=(b2+c2-a2)/2;
   /* end vars */
   let v1 = a2*SA+6*SB*SC;
   let v2 = b2*SB+6*SA*SC;
   let v3 = 6*SA*SB+c2*SC;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X547(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let SC=(a2+b2-c2)/2;
   let SB=(c2+a2-b2)/2;
   let SA=(b2+c2-a2)/2;
   /* end vars */
   let v1 = 7*a2*SA+10*SB*SC;
   let v2 = 7*b2*SB+10*SA*SC;
   let v3 = 10*SA*SB+7*c2*SC;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X548(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let SC=(a2+b2-c2)/2;
   let SB=(c2+a2-b2)/2;
   let SA=(b2+c2-a2)/2;
   /* end vars */
   let v1 = 5*a2*SA-2*SB*SC;
   let v2 = 5*b2*SB-2*SA*SC;
   let v3 = -2*SA*SB+5*c2*SC;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X549(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let SC=(a2+b2-c2)/2;
   let SB=(c2+a2-b2)/2;
   let SA=(b2+c2-a2)/2;
   /* end vars */
   let v1 = 5*a2*SA+2*SB*SC;
   let v2 = 5*b2*SB+2*SA*SC;
   let v3 = 2*SA*SB+5*c2*SC;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X550(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let SC=(a2+b2-c2)/2;
   let SB=(c2+a2-b2)/2;
   let SA=(b2+c2-a2)/2;
   /* end vars */
   let v1 = 3*a2*SA-2*SB*SC;
   let v2 = 3*b2*SB-2*SA*SC;
   let v3 = -2*SA*SB+3*c2*SC;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X551(orbit) {
   /* begin vars */

   /* end vars */
   let v1 = 4*a+b+c;
   let v2 = a+4*b+c;
   let v3 = a+b+4*c;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X552(orbit) {
   /* begin vars */

   /* end vars */
   let v1 = 1/((b+c)*(b+c)*(-a+b+c));
   let v2 = 1/((a+c)*(a+c)*(a-b+c));
   let v3 = 1/((a+b)*(a+b)*(a+b-c));
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X553(orbit) {
   /* begin vars */

   /* end vars */
   let v1 = (2*a+b+c)/(-a+b+c);
   let v2 = (a+2*b+c)/(a-b+c);
   let v3 = (a+b+2*c)/(a+b-c);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X554(orbit) {
   /* begin vars */
   let area=triAreaHeron(a,b,c);
   let sc=(a+b-c)/2;
   let sb=(c+a-b)/2;
   let sa=(b+c-a)/2;
   /* end vars */
   let v1 = 1/(sa*(sb*sc+area*sqrt3));
   let v2 = 1/(sb*(sa*sc+area*sqrt3));
   let v3 = 1/(sc*(sa*sb+area*sqrt3));
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X555(orbit) {
   /* begin vars */
   let Sqrt=Math.sqrt;
   /* end vars */
   let v1 = (a+b-c)*(a-b+c)*Sqrt(b*(a+b-c)*c*(a-b+c));
   let v2 = (a+b-c)*(-a+b+c)*Sqrt(a*(a+b-c)*c*(-a+b+c));
   let v3 = (a-b+c)*(-a+b+c)*Sqrt(a*b*(a-b+c)*(-a+b+c));
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X556(orbit) {
   /* begin vars */
   let sa=(b+c-a)/2;
   let sc=(a+b-c)/2;
   let sb=(c+a-b)/2;
   let Sqrt=Math.sqrt;
   /* end vars */
   let v1 = Sqrt((b*c)/(sb*sc));
   let v2 = Sqrt((a*c)/(sa*sc));
   let v3 = Sqrt((a*b)/(sa*sb));
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X557(orbit) {
   /* begin vars */
   let cosC=lawOfCosines(c,a,b);
   let cosHalfC=cosHalfAngle(cosC);
   let cosB=lawOfCosines(b,a,c);
   let cosHalfB=cosHalfAngle(cosB);
   let cosA=lawOfCosines(a,b,c);
   let cosHalfA=cosHalfAngle(cosA);
   let cosQuarterC=cosHalfAngle(cosHalfC);
   let cosQuarterB=cosHalfAngle(cosHalfB);
   let cosQuarterA=cosHalfAngle(cosHalfA);
   /* end vars */
   let v1 = cosQuarterA*cosQuarterA;
   let v2 = cosQuarterB*cosQuarterB;
   let v3 = cosQuarterC*cosQuarterC;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X558(orbit) {
   /* begin vars */
   let cosC=lawOfCosines(c,a,b);
   let cosHalfC=cosHalfAngle(cosC);
   let cosB=lawOfCosines(b,a,c);
   let cosHalfB=cosHalfAngle(cosB);
   let cosA=lawOfCosines(a,b,c);
   let cosHalfA=cosHalfAngle(cosA);
   let sinQuarterC=sinHalfAngle(cosHalfC);
   let sinQuarterB=sinHalfAngle(cosHalfB);
   let sinQuarterA=sinHalfAngle(cosHalfA);
   /* end vars */
   let v1 = sinQuarterA*sinQuarterA;
   let v2 = sinQuarterB*sinQuarterB;
   let v3 = sinQuarterC*sinQuarterC;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X559(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let SC=(a2+b2-c2)/2;
   let SB=(c2+a2-b2)/2;
   let area=triAreaHeron(a,b,c);
   let SA=(b2+c2-a2)/2;
   /* end vars */
   let v1 = a*(b*c-SA+2*area*sqrt3);
   let v2 = b*(a*c-SB+2*area*sqrt3);
   let v3 = c*(a*b-SC+2*area*sqrt3);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X560(orbit) {
   /* begin vars */
   let c2=c*c;
   let c3=c2*c;
   let b2=b*b;
   let b3=b2*b;
   let a2=a*a;
   let a3=a2*a;
   let c5=c2*c3;
   let b5=b2*b3;
   let a5=a2*a3;
   /* end vars */
   let v1 = a5;
   let v2 = b5;
   let v3 = c5;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X561(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let c3=c2*c;
   let b3=b2*b;
   let a3=a2*a;
   /* end vars */
   let v1 = 1/a3;
   let v2 = 1/b3;
   let v3 = 1/c3;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X562(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let SC=(a2+b2-c2)/2;
   let SB=(c2+a2-b2)/2;
   let SA=(b2+c2-a2)/2;
   let SC2=SC*SC;
   let SB2=SB*SB;
   let area=triAreaHeron(a,b,c);
   let SA2=SA*SA;
   /* end vars */
   let v1 = (b2*c2-4*SA2)/(SA*(-12*area*area+SA2));
   let v2 = (a2*c2-4*SB2)/(SB*(-12*area*area+SB2));
   let v3 = (a2*b2-4*SC2)/(SC*(-12*area*area+SC2));
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X563(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let SC=(a2+b2-c2)/2;
   let c3=c2*c;
   let SB=(c2+a2-b2)/2;
   let b3=b2*b;
   let SA=(b2+c2-a2)/2;
   let a3=a2*a;
   let SC2=SC*SC;
   let c5=c2*c3;
   let SB2=SB*SB;
   let b5=b2*b3;
   let SA2=SA*SA;
   let area=triAreaHeron(a,b,c);
   let a5=a2*a3;
   /* end vars */
   let v1 = a5*SA*(-4*area*area+SA2);
   let v2 = b5*SB*(-4*area*area+SB2);
   let v3 = c5*SC*(-4*area*area+SC2);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X564(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let SC=(a2+b2-c2)/2;
   let SB=(c2+a2-b2)/2;
   let SA=(b2+c2-a2)/2;
   let SC2=SC*SC;
   let c4=c2*c2;
   let SB2=SB*SB;
   let b4=b2*b2;
   let area=triAreaHeron(a,b,c);
   let SA2=SA*SA;
   let a4=a2*a2;
   /* end vars */
   let v1 = b*c*(a4*SA2-a2*SA*(SB-SC)*(SB-SC)+SB*(16*area*area-(SB-SC)*(SB-SC))*SC);
   let v2 = a*c*(b4*SB2-b2*SB*(-SA+SC)*(-SA+SC)+SA*SC*(16*area*area-(-SA+SC)*(-SA+SC)));
   let v3 = a*b*(SA*(16*area*area-(SA-SB)*(SA-SB))*SB-c2*(SA-SB)*(SA-SB)*SC+c4*SC2);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X565(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let SC=(a2+b2-c2)/2;
   let SB=(c2+a2-b2)/2;
   let SA=(b2+c2-a2)/2;
   let SC2=SC*SC;
   let c4=c2*c2;
   let SB2=SB*SB;
   let b4=b2*b2;
   let SA2=SA*SA;
   let a4=a2*a2;
   /* end vars */
   let v1 = b2*c2*(a2*SA+2*SB*SC)*(a4*SA2-a2*SA*(SB-3*SC)*(3*SB-SC)-SB*(SB-3*SC)*(3*SB-SC)*SC);
   let v2 = a2*c2*(b2*SB+2*SA*SC)*(b4*SB2-b2*SB*(-3*SA+SC)*(-SA+3*SC)-SA*SC*(-3*SA+SC)*(-SA+3*SC));
   let v3 = a2*b2*(2*SA*SB+c2*SC)*(-(SA*(SA-3*SB)*(3*SA-SB)*SB)-c2*(SA-3*SB)*(3*SA-SB)*SC+c4*SC2);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X566(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let SC=(a2+b2-c2)/2;
   let SB=(c2+a2-b2)/2;
   let SA=(b2+c2-a2)/2;
   let SC2=SC*SC;
   let SB2=SB*SB;
   let SA2=SA*SA;
   /* end vars */
   let v1 = a2*(a2*(5*SA2+SB*SC)+SA*(SB2+10*SB*SC+SC2));
   let v2 = b2*(b2*(5*SB2+SA*SC)+SB*(SA2+10*SA*SC+SC2));
   let v3 = c2*((SA2+10*SA*SB+SB2)*SC+c2*(SA*SB+5*SC2));
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X567(orbit) {
   /* begin vars */
   let c2=c*c;
   let c4=c2*c2;
   let c6=c2*c4;
   let b2=b*b;
   let b4=b2*b2;
   let b6=b2*b4;
   let a2=a*a;
   let a4=a2*a2;
   let a6=a2*a4;
   let c8=c2*c6;
   let b8=b2*b6;
   let a8=a2*a6;
   /* end vars */
   let v1 = a2*(a8-3*a6*b2+3*a4*b4-a2*b6-3*a6*c2+3*a4*b2*c2+2*a2*b4*c2-2*b6*c2+3*a4*c4+2*a2*b2*c4+4*b4*c4-a2*c6-2*b2*c6);
   let v2 = b2*(-(a6*b2)+3*a4*b4-3*a2*b6+b8-2*a6*c2+2*a4*b2*c2+3*a2*b4*c2-3*b6*c2+4*a4*c4+2*a2*b2*c4+3*b4*c4-2*a2*c6-b2*c6);
   let v3 = c2*(-2*a6*b2+4*a4*b4-2*a2*b6-a6*c2+2*a4*b2*c2+2*a2*b4*c2-b6*c2+3*a4*c4+3*a2*b2*c4+3*b4*c4-3*a2*c6-3*b2*c6+c8);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X568(orbit) {
   /* begin vars */
   let a2=a*a;
   let a4=a2*a2;
   let a6=a2*a4;
   let c2=c*c;
   let c4=c2*c2;
   let c6=c2*c4;
   let b2=b*b;
   let b4=b2*b2;
   let b6=b2*b4;
   let a8=a2*a6;
   let c8=c2*c6;
   let b8=b2*b6;
   /* end vars */
   let v1 = a2*(a6*b2-3*a4*b4+3*a2*b6-b8+a6*c2-a4*b2*c2-2*a2*b4*c2+2*b6*c2-3*a4*c4-2*a2*b2*c4-2*b4*c4+3*a2*c6+2*b2*c6-c8);
   let v2 = b2*(-a8+3*a6*b2-3*a4*b4+a2*b6+2*a6*c2-2*a4*b2*c2-a2*b4*c2+b6*c2-2*a4*c4-2*a2*b2*c4-3*b4*c4+2*a2*c6+3*b2*c6-c8);
   let v3 = c2*(-a8+2*a6*b2-2*a4*b4+2*a2*b6-b8+3*a6*c2-2*a4*b2*c2-2*a2*b4*c2+3*b6*c2-3*a4*c4-a2*b2*c4-3*b4*c4+a2*c6+b2*c6);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X569(orbit) {
   /* begin vars */
   let c2=c*c;
   let c4=c2*c2;
   let c6=c2*c4;
   let b2=b*b;
   let b4=b2*b2;
   let b6=b2*b4;
   let a2=a*a;
   let a4=a2*a2;
   let a6=a2*a4;
   let c8=c2*c6;
   let b8=b2*b6;
   let a8=a2*a6;
   /* end vars */
   let v1 = a2*(a8-3*a6*b2+3*a4*b4-a2*b6-3*a6*c2+2*a4*b2*c2+3*a2*b4*c2-2*b6*c2+3*a4*c4+3*a2*b2*c4+4*b4*c4-a2*c6-2*b2*c6);
   let v2 = b2*(-(a6*b2)+3*a4*b4-3*a2*b6+b8-2*a6*c2+3*a4*b2*c2+2*a2*b4*c2-3*b6*c2+4*a4*c4+3*a2*b2*c4+3*b4*c4-2*a2*c6-b2*c6);
   let v3 = c2*(-2*a6*b2+4*a4*b4-2*a2*b6-a6*c2+3*a4*b2*c2+3*a2*b4*c2-b6*c2+3*a4*c4+2*a2*b2*c4+3*b4*c4-3*a2*c6-3*b2*c6+c8);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X570(orbit) {
   /* begin vars */
   let a2=a*a;
   let a4=a2*a2;
   let c2=c*c;
   let c4=c2*c2;
   let b2=b*b;
   let b4=b2*b2;
   let a6=a2*a4;
   let c6=c2*c4;
   let b6=b2*b4;
   /* end vars */
   let v1 = a2*(a4*b2-2*a2*b4+b6+a4*c2-2*a2*b2*c2-b4*c2-2*a2*c4-b2*c4+c6);
   let v2 = b2*(a6-2*a4*b2+a2*b4-a4*c2-2*a2*b2*c2+b4*c2-a2*c4-2*b2*c4+c6);
   let v3 = c2*(a6-a4*b2-a2*b4+b6-2*a4*c2-2*a2*b2*c2-2*b4*c2+a2*c4+b2*c4);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X571(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let SC=(a2+b2-c2)/2;
   let SB=(c2+a2-b2)/2;
   let SA=(b2+c2-a2)/2;
   let SC2=SC*SC;
   let c4=c2*c2;
   let SB2=SB*SB;
   let b4=b2*b2;
   let SA2=SA*SA;
   let a4=a2*a2;
   /* end vars */
   let v1 = a4*SA*(a4*SA2-b4*SB2-c4*SC2);
   let v2 = b4*SB*(-(a4*SA2)+b4*SB2-c4*SC2);
   let v3 = c4*SC*(-(a4*SA2)-b4*SB2+c4*SC2);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X572(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let SC=(a2+b2-c2)/2;
   let SB=(c2+a2-b2)/2;
   let SA=(b2+c2-a2)/2;
   /* end vars */
   let v1 = a2*(-(b*c*(-a+b+c))-2*a*SA);
   let v2 = b2*(-(a*c*(a-b+c))-2*b*SB);
   let v3 = c2*(-(a*b*(a+b-c))-2*c*SC);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X573(orbit) {
   /* begin vars */
   let a2=a*a;
   let c2=c*c;
   let b2=b*b;
   let a3=a2*a;
   let c3=c2*c;
   let b3=b2*b;
   /* end vars */
   let v1 = a2*(a2*b-b3+a2*c-a*b*c-c3);
   let v2 = b2*(-a3+a*b2-a*b*c+b2*c-c3);
   let v3 = c2*(-a3-b3-a*b*c+a*c2+b*c2);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X574(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   /* end vars */
   let v1 = a2*(-a2+2*b2+2*c2);
   let v2 = b2*(2*a2-b2+2*c2);
   let v3 = (2*a2+2*b2-c2)*c2;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X575(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let SC=(a2+b2-c2)/2;
   let SB=(c2+a2-b2)/2;
   let SA=(b2+c2-a2)/2;
   let SC2=SC*SC;
   let SB2=SB*SB;
   let SA2=SA*SA;
   /* end vars */
   let v1 = a2*(4*a2*SA+SA2+3*SB*SC);
   let v2 = b2*(4*b2*SB+SB2+3*SA*SC);
   let v3 = c2*(3*SA*SB+4*c2*SC+SC2);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X576(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let c4=c2*c2;
   let b4=b2*b2;
   let a4=a2*a2;
   /* end vars */
   let v1 = a2*(a4-3*a2*b2+2*b4-3*a2*c2-2*b2*c2+2*c4);
   let v2 = b2*(2*a4-3*a2*b2+b4-2*a2*c2-3*b2*c2+2*c4);
   let v3 = c2*(2*a4-2*a2*b2+2*b4-3*a2*c2-3*b2*c2+c4);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X577(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let SC=(a2+b2-c2)/2;
   let SB=(c2+a2-b2)/2;
   let SA=(b2+c2-a2)/2;
   let SC2=SC*SC;
   let c4=c2*c2;
   let SB2=SB*SB;
   let b4=b2*b2;
   let SA2=SA*SA;
   let a4=a2*a2;
   /* end vars */
   let v1 = a4*SA2;
   let v2 = b4*SB2;
   let v3 = c4*SC2;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X578(orbit) {
   /* begin vars */
   let c2=c*c;
   let c4=c2*c2;
   let c6=c2*c4;
   let b2=b*b;
   let b4=b2*b2;
   let b6=b2*b4;
   let a2=a*a;
   let a4=a2*a2;
   let a6=a2*a4;
   let c8=c2*c6;
   let b8=b2*b6;
   let a8=a2*a6;
   /* end vars */
   let v1 = a2*(a8-3*a6*b2+3*a4*b4-a2*b6-3*a6*c2+4*a4*b2*c2+a2*b4*c2-2*b6*c2+3*a4*c4+a2*b2*c4+4*b4*c4-a2*c6-2*b2*c6);
   let v2 = b2*(-(a6*b2)+3*a4*b4-3*a2*b6+b8-2*a6*c2+a4*b2*c2+4*a2*b4*c2-3*b6*c2+4*a4*c4+a2*b2*c4+3*b4*c4-2*a2*c6-b2*c6);
   let v3 = c2*(-2*a6*b2+4*a4*b4-2*a2*b6-a6*c2+a4*b2*c2+a2*b4*c2-b6*c2+3*a4*c4+4*a2*b2*c4+3*b4*c4-3*a2*c6-3*b2*c6+c8);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X579(orbit) {
   /* begin vars */
   let a2=a*a;
   let c2=c*c;
   let b2=b*b;
   let a3=a2*a;
   let c3=c2*c;
   let b3=b2*b;
   /* end vars */
   let v1 = a2*(a2*b-b3+a2*c+a*b*c-c3);
   let v2 = b2*(-a3+a*b2+a*b*c+b2*c-c3);
   let v3 = c2*(-a3-b3+a*b*c+a*c2+b*c2);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X580(orbit) {
   /* begin vars */
   let c2=c*c;
   let c3=c2*c;
   let b2=b*b;
   let b3=b2*b;
   let a2=a*a;
   let a3=a2*a;
   let c5=c2*c3;
   let b5=b2*b3;
   let a4=a2*a2;
   let c4=c2*c2;
   let b4=b2*b2;
   let a5=a2*a3;
   /* end vars */
   let v1 = a2*(a5-2*a3*b2+a*b4-a3*b*c-a2*b2*c+a*b3*c+b4*c-2*a3*c2-a2*b*c2-b3*c2+a*b*c3-b2*c3+a*c4+b*c4);
   let v2 = b2*(a4*b-2*a2*b3+b5+a4*c+a3*b*c-a2*b2*c-a*b3*c-a3*c2-a*b2*c2-2*b3*c2-a2*c3+a*b*c3+a*c4+b*c4);
   let v3 = c2*(a4*b-a3*b2-a2*b3+a*b4+a4*c+a3*b*c+a*b3*c+b4*c-a2*b*c2-a*b2*c2-2*a2*c3-a*b*c3-2*b2*c3+c5);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X581(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let a3=a2*a;
   let c3=c2*c;
   let b3=b2*b;
   let c4=c2*c2;
   let b4=b2*b2;
   let a5=a2*a3;
   let c5=c2*c3;
   let b5=b2*b3;
   let a4=a2*a2;
   /* end vars */
   let v1 = a2*(a4*b-2*a2*b3+b5+a4*c+a3*b*c-a2*b2*c-a*b3*c-a2*b*c2-2*a*b2*c2-b3*c2-2*a2*c3-a*b*c3-b2*c3+c5);
   let v2 = b2*(a5-2*a3*b2+a*b4-a3*b*c-a2*b2*c+a*b3*c+b4*c-a3*c2-2*a2*b*c2-a*b2*c2-a2*c3-a*b*c3-2*b2*c3+c5);
   let v3 = c2*(a5-a3*b2-a2*b3+b5-a3*b*c-2*a2*b2*c-a*b3*c-2*a3*c2-a2*b*c2-a*b2*c2-2*b3*c2+a*b*c3+a*c4+b*c4);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X582(orbit) {
   /* begin vars */
   let c2=c*c;
   let c3=c2*c;
   let b2=b*b;
   let b3=b2*b;
   let a2=a*a;
   let a3=a2*a;
   let c5=c2*c3;
   let b5=b2*b3;
   let a4=a2*a2;
   let c4=c2*c2;
   let b4=b2*b2;
   let a5=a2*a3;
   /* end vars */
   let v1 = a2*(a5-2*a3*b2+a*b4-2*a3*b*c-a2*b2*c+2*a*b3*c+b4*c-2*a3*c2-a2*b*c2-b3*c2+2*a*b*c3-b2*c3+a*c4+b*c4);
   let v2 = b2*(a4*b-2*a2*b3+b5+a4*c+2*a3*b*c-a2*b2*c-2*a*b3*c-a3*c2-a*b2*c2-2*b3*c2-a2*c3+2*a*b*c3+a*c4+b*c4);
   let v3 = c2*(a4*b-a3*b2-a2*b3+a*b4+a4*c+2*a3*b*c+2*a*b3*c+b4*c-a2*b*c2-a*b2*c2-2*a2*c3-2*a*b*c3-2*b2*c3+c5);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X583(orbit) {
   /* begin vars */
   let a2=a*a;
   let c2=c*c;
   let b2=b*b;
   let a3=a2*a;
   let c3=c2*c;
   let b3=b2*b;
   /* end vars */
   let v1 = a2*(a2*b-b3+a2*c+2*a*b*c-c3);
   let v2 = b2*(-a3+a*b2+2*a*b*c+b2*c-c3);
   let v3 = c2*(-a3-b3+2*a*b*c+a*c2+b*c2);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X584(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let c3=c2*c;
   let b3=b2*b;
   let a3=a2*a;
   /* end vars */
   let v1 = a2*(a3-a*b2-2*a*b*c-b2*c-a*c2-b*c2);
   let v2 = b2*(-(a2*b)+b3-a2*c-2*a*b*c-a*c2-b*c2);
   let v3 = c2*(-(a2*b)-a*b2-a2*c-2*a*b*c-b2*c+c3);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X585(orbit) {
   /* begin vars */
   let area=triAreaHeron(a,b,c);
   /* end vars */
   let v1 = -a+b+2*area*(-(1/a)+1/b+1/c)+c;
   let v2 = a-b+2*area*(1/a-1/b+1/c)+c;
   let v3 = a+b+2*area*(1/a+1/b-1/c)-c;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X586(orbit) {
   /* begin vars */
   let area=triAreaHeron(a,b,c);
   /* end vars */
   let v1 = -a+b-2*area*(-(1/a)+1/b+1/c)+c;
   let v2 = a-b-2*area*(1/a-1/b+1/c)+c;
   let v3 = a+b-2*area*(1/a+1/b-1/c)-c;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X587(orbit) {
   /* begin vars */
   let a2=a*a;
   let c2=c*c;
   let b2=b*b;
   let a3=a2*a;
   let SA=(b2+c2-a2)/2;
   let c3=c2*c;
   let b3=b2*b;
   let SC=(a2+b2-c2)/2;
   let SB=(c2+a2-b2)/2;
   let area=triAreaHeron(a,b,c);
   /* end vars */
   let v1 = area*(2*a-2*(b+c))*SB*SC+SB*SC*(b3+c3+2*a*SA+b*SA+c*SA-c*SB-b*SC);
   let v2 = area*(2*b-2*(a+c))*SA*SC+SA*SC*(a3+c3-c*SA+a*SB+2*b*SB+c*SB-a*SC);
   let v3 = area*(-2*(a+b)+2*c)*SA*SB+SA*SB*(a3+b3-b*SA-a*SB+a*SC+b*SC+2*c*SC);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X588(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let area=triAreaHeron(a,b,c);
   let a2=a*a;
   /* end vars */
   let v1 = a2/(a2+4*area);
   let v2 = b2/(4*area+b2);
   let v3 = c2/(4*area+c2);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X589(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let area=triAreaHeron(a,b,c);
   let a2=a*a;
   /* end vars */
   let v1 = a2/(a2-4*area);
   let v2 = b2/(-4*area+b2);
   let v3 = c2/(-4*area+c2);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X590(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let area=triAreaHeron(a,b,c);
   let a2=a*a;
   /* end vars */
   let v1 = a2+4*area;
   let v2 = 4*area+b2;
   let v3 = 4*area+c2;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X591(orbit) {
   /* begin vars */
   let area=triAreaHeron(a,b,c);
   let S=2*area;
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   /* end vars */
   let v1 = 2*a2-b2-c2-2*S;
   let v2 = -a2+2*b2-c2-2*S;
   let v3 = -a2-b2+2*c2-2*S;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X592(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let SC=(a2+b2-c2)/2;
   let SB=(c2+a2-b2)/2;
   let SA=(b2+c2-a2)/2;
   let SC2=SC*SC;
   let SB2=SB*SB;
   let SA2=SA*SA;
   /* end vars */
   let v1 = a2/(2*SB*(a2+SB)*SC*(a2+SC)+a2*SA*(SA2+3*SB2+17*SB*SC+3*SC2)+SA2*(7*SB2+16*SB*SC+7*SC2));
   let v2 = b2/(2*SA*(b2+SA)*SC*(b2+SC)+b2*SB*(3*SA2+SB2+17*SA*SC+3*SC2)+SB2*(7*SA2+16*SA*SC+7*SC2));
   let v3 = c2/(2*SA*(c2+SA)*SB*(c2+SB)+(7*SA2+16*SA*SB+7*SB2)*SC2+c2*SC*(3*SA2+17*SA*SB+3*SB2+SC2));
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X593(orbit) {
   /* begin vars */

   /* end vars */
   let v1 = a*a/(b+c)*(b+c);
   let v2 = b*b/(a+c)*(a+c);
   let v3 = c*c/(a+b)*(a+b);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X594(orbit) {
   /* begin vars */

   /* end vars */
   let v1 = (b+c)*(b+c);
   let v2 = (a+c)*(a+c);
   let v3 = (a+b)*(a+b);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X595(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   /* end vars */
   let v1 = a2*(a2+a*b+a*c-b*c);
   let v2 = b2*(a*b+b2-a*c+b*c);
   let v3 = c2*(-(a*b)+a*c+b*c+c2);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X596(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   /* end vars */
   let v1 = 1/(a2+a*b+a*c-b*c);
   let v2 = 1/(a*b+b2-a*c+b*c);
   let v3 = 1/(-(a*b)+a*c+b*c+c2);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X597(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   /* end vars */
   let v1 = 4*a2+b2+c2;
   let v2 = a2+4*b2+c2;
   let v3 = a2+b2+4*c2;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X598(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   /* end vars */
   let v1 = 1/(a2-2*b2-2*c2);
   let v2 = 1/(-2*a2+b2-2*c2);
   let v3 = 1/(-2*a2-2*b2+c2);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X599(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   /* end vars */
   let v1 = a2-2*b2-2*c2;
   let v2 = -2*a2+b2-2*c2;
   let v3 = -2*a2-2*b2+c2;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X600(orbit) {
   /* begin vars */
   let area=triAreaHeron(a,b,c);
   let c2=c*c;
   let b2=b*b;
   let S=2*area;
   let a2=a*a;
   /* end vars */
   let v1 = a2*(b*c+2*S)*(a2*b*c-a*b2*c-a*b*c2+2*(a2-b2-c2)*S);
   let v2 = b2*(a*c+2*S)*(-(a2*b*c)+a*b2*c-a*b*c2+2*(-a2+b2-c2)*S);
   let v3 = c2*(a*b+2*S)*(-(a2*b*c)-a*b2*c+a*b*c2+2*(-a2-b2+c2)*S);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X601(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let SC=(a2+b2-c2)/2;
   let c3=c2*c;
   let SB=(c2+a2-b2)/2;
   let b3=b2*b;
   let SA=(b2+c2-a2)/2;
   let area=triAreaHeron(a,b,c);
   let a3=a2*a;
   /* end vars */
   let v1 = a3*(4*area*area+b*c*SA);
   let v2 = b3*(4*area*area+a*c*SB);
   let v3 = c3*(4*area*area+a*b*SC);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X602(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let SC=(a2+b2-c2)/2;
   let c3=c2*c;
   let SB=(c2+a2-b2)/2;
   let b3=b2*b;
   let SA=(b2+c2-a2)/2;
   let area=triAreaHeron(a,b,c);
   let a3=a2*a;
   /* end vars */
   let v1 = a3*(-4*area*area+b*c*SA);
   let v2 = b3*(-4*area*area+a*c*SB);
   let v3 = c3*(-4*area*area+a*b*SC);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X603(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let SC=(a2+b2-c2)/2;
   let c3=c2*c;
   let SB=(c2+a2-b2)/2;
   let b3=b2*b;
   let SA=(b2+c2-a2)/2;
   let a3=a2*a;
   /* end vars */
   let v1 = a3*(a+b-c)*(a-b+c)*SA;
   let v2 = b3*(a+b-c)*(-a+b+c)*SB;
   let v3 = (a-b+c)*(-a+b+c)*c3*SC;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X604(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let c3=c2*c;
   let b3=b2*b;
   let a3=a2*a;
   /* end vars */
   let v1 = a3*(a+b-c)*(a-b+c);
   let v2 = b3*(a+b-c)*(-a+b+c);
   let v3 = (a-b+c)*(-a+b+c)*c3;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X605(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let c3=c2*c;
   let b3=b2*b;
   let area=triAreaHeron(a,b,c);
   let a3=a2*a;
   /* end vars */
   let v1 = a3*(2*area+b*c);
   let v2 = b3*(2*area+a*c);
   let v3 = (2*area+a*b)*c3;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X606(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let c3=c2*c;
   let b3=b2*b;
   let area=triAreaHeron(a,b,c);
   let a3=a2*a;
   /* end vars */
   let v1 = a3*(-2*area+b*c);
   let v2 = b3*(-2*area+a*c);
   let v3 = (-2*area+a*b)*c3;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X607(orbit) {
   /* begin vars */
   let a2=a*a;
   let c2=c*c;
   let b2=b*b;
   let SA=(b2+c2-a2)/2;
   let SC=(a2+b2-c2)/2;
   let SB=(c2+a2-b2)/2;
   /* end vars */
   let v1 = a2*(a-b-c)*SB*SC;
   let v2 = b2*(-a+b-c)*SA*SC;
   let v3 = (-a-b+c)*c2*SA*SB;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X608(orbit) {
   /* begin vars */
   let a2=a*a;
   let c2=c*c;
   let b2=b*b;
   let SA=(b2+c2-a2)/2;
   let SC=(a2+b2-c2)/2;
   let SB=(c2+a2-b2)/2;
   /* end vars */
   let v1 = a2*(a+b-c)*(a-b+c)*SB*SC;
   let v2 = b2*(a+b-c)*(-a+b+c)*SA*SC;
   let v3 = (a-b+c)*(-a+b+c)*c2*SA*SB;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X609(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   /* end vars */
   let v1 = a2*(2*a2+b*c);
   let v2 = b2*(2*b2+a*c);
   let v3 = c2*(a*b+2*c2);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X610(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let SC=(a2+b2-c2)/2;
   let SB=(c2+a2-b2)/2;
   let SA=(b2+c2-a2)/2;
   /* end vars */
   let v1 = a*(-(a2*SA)+SB*SC);
   let v2 = b*(-(b2*SB)+SA*SC);
   let v3 = c*(SA*SB-c2*SC);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X611(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let area=triAreaHeron(a,b,c);
   let a2=a*a;
   /* end vars */
   let v1 = a2*(8*area*area+b*c*(a2+b2+c2));
   let v2 = b2*(8*area*area+a*c*(a2+b2+c2));
   let v3 = c2*(8*area*area+a*b*(a2+b2+c2));
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X612(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   /* end vars */
   let v1 = a*(a2+b2+2*b*c+c2);
   let v2 = b*(a2+b2+2*a*c+c2);
   let v3 = c*(a2+2*a*b+b2+c2);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X613(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let area=triAreaHeron(a,b,c);
   let a2=a*a;
   /* end vars */
   let v1 = a2*(-8*area*area+b*c*(a2+b2+c2));
   let v2 = b2*(-8*area*area+a*c*(a2+b2+c2));
   let v3 = c2*(-8*area*area+a*b*(a2+b2+c2));
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X614(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   /* end vars */
   let v1 = a*(a2+b2-2*b*c+c2);
   let v2 = b*(a2+b2-2*a*c+c2);
   let v3 = c*(a2-2*a*b+b2+c2);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X615(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let area=triAreaHeron(a,b,c);
   let a2=a*a;
   /* end vars */
   let v1 = a2-4*area;
   let v2 = -4*area+b2;
   let v3 = -4*area+c2;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X616(orbit) {
   /* begin vars */
   let area=triAreaHeron(a,b,c);
   let S=2*area;
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   /* end vars */
   let v1 = (a+b-c)*(a-b+c)*(-a+b+c)*(a+b+c)+6*a2*(a2-b2-c2)+2*(a2-b2-c2)*S*sqrt3;
   let v2 = (a+b-c)*(a-b+c)*(-a+b+c)*(a+b+c)+6*b2*(-a2+b2-c2)+2*(-a2+b2-c2)*S*sqrt3;
   let v3 = (a+b-c)*(a-b+c)*(-a+b+c)*(a+b+c)+6*c2*(-a2-b2+c2)+2*(-a2-b2+c2)*S*sqrt3;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X617(orbit) {
   /* begin vars */
   let area=triAreaHeron(a,b,c);
   let S=2*area;
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   /* end vars */
   let v1 = (a+b-c)*(a-b+c)*(-a+b+c)*(a+b+c)+6*a2*(a2-b2-c2)-2*(a2-b2-c2)*S*sqrt3;
   let v2 = (a+b-c)*(a-b+c)*(-a+b+c)*(a+b+c)+6*b2*(-a2+b2-c2)-2*(-a2+b2-c2)*S*sqrt3;
   let v3 = (a+b-c)*(a-b+c)*(-a+b+c)*(a+b+c)+6*c2*(-a2-b2+c2)-2*(-a2-b2+c2)*S*sqrt3;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X618(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let area=triAreaHeron(a,b,c);
   let SC=(a2+b2-c2)/2;
   let SB=(c2+a2-b2)/2;
   let SA=(b2+c2-a2)/2;
   /* end vars */
   let v1 = 5*a2*SA+2*SB*SC+2*area*(b2+c2)*sqrt3;
   let v2 = 5*b2*SB+2*SA*SC+2*area*(a2+c2)*sqrt3;
   let v3 = 2*SA*SB+5*c2*SC+2*area*(a2+b2)*sqrt3;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X619(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let area=triAreaHeron(a,b,c);
   let SC=(a2+b2-c2)/2;
   let SB=(c2+a2-b2)/2;
   let SA=(b2+c2-a2)/2;
   /* end vars */
   let v1 = 5*a2*SA+2*SB*SC-2*area*(b2+c2)*sqrt3;
   let v2 = 5*b2*SB+2*SA*SC-2*area*(a2+c2)*sqrt3;
   let v3 = 2*SA*SB+5*c2*SC-2*area*(a2+b2)*sqrt3;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X620(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let SC=(a2+b2-c2)/2;
   let SB=(c2+a2-b2)/2;
   let a4=a2*a2;
   let SA=(b2+c2-a2)/2;
   let c4=c2*c2;
   let b4=b2*b2;
   /* end vars */
   let v1 = b4+c4-4*a2*SA;
   let v2 = a4+c4-4*b2*SB;
   let v3 = a4+b4-4*c2*SC;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X621(orbit) {
   /* begin vars */
   let area=triAreaHeron(a,b,c);
   let S=2*area;
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   /* end vars */
   let v1 = -2*(a2-b2-c2)*S+(a2+b2-c2)*(a2-b2+c2)*sqrt3;
   let v2 = -2*(-a2+b2-c2)*S+(a2+b2-c2)*(-a2+b2+c2)*sqrt3;
   let v3 = -2*(-a2-b2+c2)*S+(a2-b2+c2)*(-a2+b2+c2)*sqrt3;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X622(orbit) {
   /* begin vars */
   let area=triAreaHeron(a,b,c);
   let S=2*area;
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   /* end vars */
   let v1 = 2*(a2-b2-c2)*S+(a2+b2-c2)*(a2-b2+c2)*sqrt3;
   let v2 = 2*(-a2+b2-c2)*S+(a2+b2-c2)*(-a2+b2+c2)*sqrt3;
   let v3 = 2*(-a2-b2+c2)*S+(a2-b2+c2)*(-a2+b2+c2)*sqrt3;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X623(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let SC=(a2+b2-c2)/2;
   let SB=(c2+a2-b2)/2;
   let SA=(b2+c2-a2)/2;
   let area=triAreaHeron(a,b,c);
   /* end vars */
   let v1 = 2*area*(b2+c2)+(a2*SA+2*SB*SC)*sqrt3;
   let v2 = 2*area*(a2+c2)+(b2*SB+2*SA*SC)*sqrt3;
   let v3 = 2*area*(a2+b2)+(2*SA*SB+c2*SC)*sqrt3;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X624(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let SC=(a2+b2-c2)/2;
   let SB=(c2+a2-b2)/2;
   let SA=(b2+c2-a2)/2;
   let area=triAreaHeron(a,b,c);
   /* end vars */
   let v1 = -2*area*(b2+c2)+(a2*SA+2*SB*SC)*sqrt3;
   let v2 = -2*area*(a2+c2)+(b2*SB+2*SA*SC)*sqrt3;
   let v3 = -2*area*(a2+b2)+(2*SA*SB+c2*SC)*sqrt3;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X625(orbit) {
   /* begin vars */
   let a2=a*a;
   let c2=c*c;
   let b2=b*b;
   let a4=a2*a2;
   let c4=c2*c2;
   let b4=b2*b2;
   /* end vars */
   let v1 = a2*(b2+c2)-2*(b4-b2*c2+c4);
   let v2 = b2*(a2+c2)-2*(a4-a2*c2+c4);
   let v3 = -2*(a4-a2*b2+b4)+(a2+b2)*c2;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X626(orbit) {
   /* begin vars */
   let a2=a*a;
   let c2=c*c;
   let b2=b*b;
   let a4=a2*a2;
   let c4=c2*c2;
   let b4=b2*b2;
   /* end vars */
   let v1 = b4+c4;
   let v2 = a4+c4;
   let v3 = a4+b4;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X627(orbit) {
   /* begin vars */
   let area=triAreaHeron(a,b,c);
   let S=2*area;
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   /* end vars */
   let v1 = (a+b-c)*(a-b+c)*(-a+b+c)*(a+b+c)-2*(a2-b2-c2)*(a2+S*sqrt3);
   let v2 = (a+b-c)*(a-b+c)*(-a+b+c)*(a+b+c)-2*(-a2+b2-c2)*(b2+S*sqrt3);
   let v3 = (a+b-c)*(a-b+c)*(-a+b+c)*(a+b+c)-2*(-a2-b2+c2)*(c2+S*sqrt3);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X628(orbit) {
   /* begin vars */
   let area=triAreaHeron(a,b,c);
   let S=2*area;
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   /* end vars */
   let v1 = (a+b-c)*(a-b+c)*(-a+b+c)*(a+b+c)-2*(a2-b2-c2)*(a2-S*sqrt3);
   let v2 = (a+b-c)*(a-b+c)*(-a+b+c)*(a+b+c)-2*(-a2+b2-c2)*(b2-S*sqrt3);
   let v3 = (a+b-c)*(a-b+c)*(-a+b+c)*(a+b+c)-2*(-a2-b2+c2)*(c2-S*sqrt3);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X629(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let area=triAreaHeron(a,b,c);
   let SC=(a2+b2-c2)/2;
   let SB=(c2+a2-b2)/2;
   let SA=(b2+c2-a2)/2;
   /* end vars */
   let v1 = 7*a2*SA+6*SB*SC+2*area*(b2+c2)*sqrt3;
   let v2 = 7*b2*SB+6*SA*SC+2*area*(a2+c2)*sqrt3;
   let v3 = 6*SA*SB+7*c2*SC+2*area*(a2+b2)*sqrt3;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X630(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let area=triAreaHeron(a,b,c);
   let SC=(a2+b2-c2)/2;
   let SB=(c2+a2-b2)/2;
   let SA=(b2+c2-a2)/2;
   /* end vars */
   let v1 = 7*a2*SA+6*SB*SC-2*area*(b2+c2)*sqrt3;
   let v2 = 7*b2*SB+6*SA*SC-2*area*(a2+c2)*sqrt3;
   let v3 = 6*SA*SB+7*c2*SC-2*area*(a2+b2)*sqrt3;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X631(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let SC=(a2+b2-c2)/2;
   let SB=(c2+a2-b2)/2;
   let SA=(b2+c2-a2)/2;
   /* end vars */
   let v1 = 2*a2*SA+SB*SC;
   let v2 = 2*b2*SB+SA*SC;
   let v3 = SA*SB+2*c2*SC;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X632(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let SC=(a2+b2-c2)/2;
   let SB=(c2+a2-b2)/2;
   let SA=(b2+c2-a2)/2;
   /* end vars */
   let v1 = 7*a2*SA+6*SB*SC;
   let v2 = 7*b2*SB+6*SA*SC;
   let v3 = 6*SA*SB+7*c2*SC;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X633(orbit) {
   /* begin vars */
   let area=triAreaHeron(a,b,c);
   let S=2*area;
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   /* end vars */
   let v1 = (a+b-c)*(a-b+c)*(-a+b+c)*(a+b+c)+2*(a2-b2-c2)*(a2-S*sqrt3);
   let v2 = (a+b-c)*(a-b+c)*(-a+b+c)*(a+b+c)+2*(-a2+b2-c2)*(b2-S*sqrt3);
   let v3 = (a+b-c)*(a-b+c)*(-a+b+c)*(a+b+c)+2*(-a2-b2+c2)*(c2-S*sqrt3);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X634(orbit) {
   /* begin vars */
   let area=triAreaHeron(a,b,c);
   let S=2*area;
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   /* end vars */
   let v1 = (a+b-c)*(a-b+c)*(-a+b+c)*(a+b+c)+2*(a2-b2-c2)*(a2+S*sqrt3);
   let v2 = (a+b-c)*(a-b+c)*(-a+b+c)*(a+b+c)+2*(-a2+b2-c2)*(b2+S*sqrt3);
   let v3 = (a+b-c)*(a-b+c)*(-a+b+c)*(a+b+c)+2*(-a2-b2+c2)*(c2+S*sqrt3);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X635(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let area=triAreaHeron(a,b,c);
   let SC=(a2+b2-c2)/2;
   let SB=(c2+a2-b2)/2;
   let SA=(b2+c2-a2)/2;
   /* end vars */
   let v1 = a2*SA+2*SB*SC+2*area*(b2+c2)*sqrt3;
   let v2 = b2*SB+2*SA*SC+2*area*(a2+c2)*sqrt3;
   let v3 = 2*SA*SB+c2*SC+2*area*(a2+b2)*sqrt3;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X636(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let area=triAreaHeron(a,b,c);
   let SC=(a2+b2-c2)/2;
   let SB=(c2+a2-b2)/2;
   let SA=(b2+c2-a2)/2;
   /* end vars */
   let v1 = a2*SA+2*SB*SC-2*area*(b2+c2)*sqrt3;
   let v2 = b2*SB+2*SA*SC-2*area*(a2+c2)*sqrt3;
   let v3 = 2*SA*SB+c2*SC-2*area*(a2+b2)*sqrt3;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X637(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let SC=(a2+b2-c2)/2;
   let SB=(c2+a2-b2)/2;
   let SA=(b2+c2-a2)/2;
   let area=triAreaHeron(a,b,c);
   /* end vars */
   let v1 = 2*area*SA+SB*SC;
   let v2 = 2*area*SB+SA*SC;
   let v3 = SA*SB+2*area*SC;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X638(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let SC=(a2+b2-c2)/2;
   let SB=(c2+a2-b2)/2;
   let SA=(b2+c2-a2)/2;
   let area=triAreaHeron(a,b,c);
   /* end vars */
   let v1 = -2*area*SA+SB*SC;
   let v2 = -2*area*SB+SA*SC;
   let v3 = SA*SB-2*area*SC;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X639(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let SC=(a2+b2-c2)/2;
   let SB=(c2+a2-b2)/2;
   let SA=(b2+c2-a2)/2;
   let area=triAreaHeron(a,b,c);
   /* end vars */
   let v1 = 2*area*(b2+c2)+a2*SA+2*SB*SC;
   let v2 = 2*area*(a2+c2)+b2*SB+2*SA*SC;
   let v3 = 2*area*(a2+b2)+2*SA*SB+c2*SC;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X640(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let SC=(a2+b2-c2)/2;
   let SB=(c2+a2-b2)/2;
   let SA=(b2+c2-a2)/2;
   let area=triAreaHeron(a,b,c);
   /* end vars */
   let v1 = -2*area*(b2+c2)+a2*SA+2*SB*SC;
   let v2 = -2*area*(a2+c2)+b2*SB+2*SA*SC;
   let v3 = -2*area*(a2+b2)+2*SA*SB+c2*SC;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X641(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let SC=(a2+b2-c2)/2;
   let SB=(c2+a2-b2)/2;
   let SA=(b2+c2-a2)/2;
   let area=triAreaHeron(a,b,c);
   /* end vars */
   let v1 = 2*area*(b2+c2)+3*a2*SA+2*SB*SC;
   let v2 = 2*area*(a2+c2)+3*b2*SB+2*SA*SC;
   let v3 = 2*area*(a2+b2)+2*SA*SB+3*c2*SC;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X642(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let SC=(a2+b2-c2)/2;
   let SB=(c2+a2-b2)/2;
   let SA=(b2+c2-a2)/2;
   let area=triAreaHeron(a,b,c);
   /* end vars */
   let v1 = -2*area*(b2+c2)+3*a2*SA+2*SB*SC;
   let v2 = -2*area*(a2+c2)+3*b2*SB+2*SA*SC;
   let v3 = -2*area*(a2+b2)+2*SA*SB+3*c2*SC;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X643(orbit) {
   /* begin vars */
   let a2=a*a;
   let c2=c*c;
   let b2=b*b;
   /* end vars */
   let v1 = (a*(-a+b+c))/(b2-c2);
   let v2 = (b*(a-b+c))/(-a2+c2);
   let v3 = ((a+b-c)*c)/(a2-b2);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X644(orbit) {
   /* begin vars */

   /* end vars */
   let v1 = (a*(-a+b+c))/(b-c);
   let v2 = (b*(a-b+c))/(-a+c);
   let v3 = ((a+b-c)*c)/(a-b);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X645(orbit) {
   /* begin vars */
   let a2=a*a;
   let c2=c*c;
   let b2=b*b;
   /* end vars */
   let v1 = (-a+b+c)/(b2-c2);
   let v2 = (a-b+c)/(-a2+c2);
   let v3 = (a+b-c)/(a2-b2);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X646(orbit) {
   /* begin vars */

   /* end vars */
   let v1 = (-a+b+c)/(a*(b-c));
   let v2 = (a-b+c)/(b*(-a+c));
   let v3 = (a+b-c)/((a-b)*c);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X647(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let SC=(a2+b2-c2)/2;
   let SB=(c2+a2-b2)/2;
   let SA=(b2+c2-a2)/2;
   /* end vars */
   let v1 = a2*(b2-c2)*SA;
   let v2 = b2*(-a2+c2)*SB;
   let v3 = (a2-b2)*c2*SC;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X648(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let SC=(a2+b2-c2)/2;
   let SA=(b2+c2-a2)/2;
   let SB=(c2+a2-b2)/2;
   /* end vars */
   let v1 = SB*(-SA+SB)*(SA-SC)*SC;
   let v2 = SA*(-SA+SB)*SC*(-SB+SC);
   let v3 = SA*SB*(SA-SC)*(-SB+SC);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X649(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   /* end vars */
   let v1 = a2*(b-c);
   let v2 = b2*(-a+c);
   let v3 = (a-b)*c2;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X650(orbit) {
   /* begin vars */

   /* end vars */
   let v1 = a*(b-c)*(-a+b+c);
   let v2 = b*(-a+c)*(a-b+c);
   let v3 = (a-b)*(a+b-c)*c;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X651(orbit) {
   /* begin vars */

   /* end vars */
   let v1 = a/((b-c)*(-a+b+c));
   let v2 = b/((-a+c)*(a-b+c));
   let v3 = c/((a-b)*(a+b-c));
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X652(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let SC=(a2+b2-c2)/2;
   let SB=(c2+a2-b2)/2;
   let SA=(b2+c2-a2)/2;
   /* end vars */
   let v1 = a2*(a-b-c)*(b-c)*SA;
   let v2 = b2*(-a+b-c)*(-a+c)*SB;
   let v3 = (a-b)*(-a-b+c)*c2*SC;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X653(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let SC=(a2+b2-c2)/2;
   let SB=(c2+a2-b2)/2;
   let SA=(b2+c2-a2)/2;
   /* end vars */
   let v1 = 1/((a-b-c)*(b-c)*SA);
   let v2 = 1/((-a+b-c)*(-a+c)*SB);
   let v3 = 1/((a-b)*(-a-b+c)*SC);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X654(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let SC=(a2+b2-c2)/2;
   let SB=(c2+a2-b2)/2;
   let SA=(b2+c2-a2)/2;
   /* end vars */
   let v1 = a2*(a-b-c)*(b-c)*(b*c-2*SA);
   let v2 = b2*(-a+b-c)*(-a+c)*(a*c-2*SB);
   let v3 = (a-b)*(-a-b+c)*c2*(a*b-2*SC);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X655(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let SC=(a2+b2-c2)/2;
   let SB=(c2+a2-b2)/2;
   let SA=(b2+c2-a2)/2;
   /* end vars */
   let v1 = 1/((a-b-c)*(b-c)*(b*c-2*SA));
   let v2 = 1/((-a+b-c)*(-a+c)*(a*c-2*SB));
   let v3 = 1/((a-b)*(-a-b+c)*(a*b-2*SC));
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X656(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let SC=(a2+b2-c2)/2;
   let SB=(c2+a2-b2)/2;
   let SA=(b2+c2-a2)/2;
   /* end vars */
   let v1 = a*SA*(SB-SC);
   let v2 = b*SB*(-SA+SC);
   let v3 = c*(SA-SB)*SC;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X657(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   /* end vars */
   let v1 = a2*(a-b-c)*(a-b-c)*(b-c);
   let v2 = b2*(-a+b-c)*(-a+b-c)*(-a+c);
   let v3 = (a-b)*(-a-b+c)*(-a-b+c)*c2;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X658(orbit) {
   /* begin vars */

   /* end vars */
   let v1 = 1/((a-b-c)*(a-b-c)*(b-c));
   let v2 = 1/((-a+b-c)*(-a+b-c)*(-a+c));
   let v3 = 1/((a-b)*(-a-b+c)*(-a-b+c));
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X659(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   /* end vars */
   let v1 = a*(b-c)*(a2-b*c);
   let v2 = b*(-a+c)*(b2-a*c);
   let v3 = (a-b)*c*(-(a*b)+c2);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X660(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   /* end vars */
   let v1 = a/((b-c)*(a2-b*c));
   let v2 = b/((-a+c)*(b2-a*c));
   let v3 = c/((a-b)*(-(a*b)+c2));
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X661(orbit) {
   /* begin vars */
   let a2=a*a;
   let c2=c*c;
   let b2=b*b;
   /* end vars */
   let v1 = a*(b2-c2);
   let v2 = b*(-a2+c2);
   let v3 = (a2-b2)*c;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X662(orbit) {
   /* begin vars */
   let a2=a*a;
   let c2=c*c;
   let b2=b*b;
   /* end vars */
   let v1 = a/(b2-c2);
   let v2 = b/(-a2+c2);
   let v3 = c/(a2-b2);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X663(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   /* end vars */
   let v1 = a2*(b-c)*(-a+b+c);
   let v2 = b2*(-a+c)*(a-b+c);
   let v3 = (a-b)*(a+b-c)*c2;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X664(orbit) {
   /* begin vars */

   /* end vars */
   let v1 = 1/((a-b-c)*(b-c));
   let v2 = 1/((-a+b-c)*(-a+c));
   let v3 = 1/((a-b)*(-a-b+c));
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X665(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   /* end vars */
   let v1 = a2*((a-b)*(a-b)*(a+b-c)-(-a+c)*(-a+c)*(a-b+c));
   let v2 = b2*(-((a-b)*(a-b)*(a+b-c))+(b-c)*(b-c)*(-a+b+c));
   let v3 = ((-a+c)*(-a+c)*(a-b+c)-(b-c)*(b-c)*(-a+b+c))*c2;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X666(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   /* end vars */
   let v1 = (a-b)*(a-c)*(a2+b2-a*c-b*c)*(a2-a*b-b*c+c2);
   let v2 = (-a+b)*(b-c)*(a2+b2-a*c-b*c)*(-(a*b)+b2-a*c+c2);
   let v3 = (-a+c)*(-b+c)*(-(a*b)+b2-a*c+c2)*(a2-a*b-b*c+c2);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X667(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let c3=c2*c;
   let b3=b2*b;
   let a3=a2*a;
   /* end vars */
   let v1 = a3*(b-c);
   let v2 = b3*(-a+c);
   let v3 = (a-b)*c3;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X668(orbit) {
   /* begin vars */

   /* end vars */
   let v1 = 1/(a*(b-c));
   let v2 = 1/(b*(-a+c));
   let v3 = 1/((a-b)*c);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X669(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let c4=c2*c2;
   let b4=b2*b2;
   let a4=a2*a2;
   /* end vars */
   let v1 = a4*(b2-c2);
   let v2 = b4*(-a2+c2);
   let v3 = (a2-b2)*c4;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X670(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   /* end vars */
   let v1 = 1/(a2*(b2-c2));
   let v2 = 1/(b2*(-a2+c2));
   let v3 = 1/((a2-b2)*c2);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X671(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   /* end vars */
   let v1 = 1/(2*a2-b2-c2);
   let v2 = 1/(-a2+2*b2-c2);
   let v3 = 1/(-a2-b2+2*c2);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X672(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   /* end vars */
   let v1 = a2*(b2-a*(b+c)+c2);
   let v2 = b2*(a2-b*(a+c)+c2);
   let v3 = (a2+b2-(a+b)*c)*c2;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X673(orbit) {
   /* begin vars */
   let a2=a*a;
   let c2=c*c;
   let b2=b*b;
   /* end vars */
   let v1 = 1/(b2-a*(b+c)+c2);
   let v2 = 1/(a2-b*(a+c)+c2);
   let v3 = 1/(a2+b2-(a+b)*c);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X674(orbit) {
   /* begin vars */
   let a2=a*a;
   let c2=c*c;
   let b2=b*b;
   let a3=a2*a;
   let c3=c2*c;
   let b3=b2*b;
   /* end vars */
   let v1 = a2*(b3-a*(b2+c2)+c3);
   let v2 = b2*(a3-b*(a2+c2)+c3);
   let v3 = (a3+b3-(a2+b2)*c)*c2;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X675(orbit) {
   /* begin vars */
   let a2=a*a;
   let c2=c*c;
   let b2=b*b;
   let a3=a2*a;
   let c3=c2*c;
   let b3=b2*b;
   /* end vars */
   let v1 = 1/(b3-a*(b2+c2)+c3);
   let v2 = 1/(a3-b*(a2+c2)+c3);
   let v3 = 1/(a3+b3-(a2+b2)*c);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X676(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let c3=c2*c;
   let b3=b2*b;
   let a3=a2*a;
   /* end vars */
   let v1 = (b-c)*(2*a3-a2*b-b3-a2*c+b2*c+b*c2-c3);
   let v2 = (-a+c)*(-a3-a*b2+2*b3+a2*c-b2*c+a*c2-c3);
   let v3 = (a-b)*(-a3+a2*b+a*b2-b3-a*c2-b*c2+2*c3);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X677(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let c3=c2*c;
   let b3=b2*b;
   let a3=a2*a;
   /* end vars */
   let v1 = a2/((b-c)*(2*a3-a2*b-b3-a2*c+b2*c+b*c2-c3));
   let v2 = b2/((-a+c)*(-a3-a*b2+2*b3+a2*c-b2*c+a*c2-c3));
   let v3 = c2/((a-b)*(-a3+a2*b+a*b2-b3-a*c2-b*c2+2*c3));
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X678(orbit) {
   /* begin vars */

   /* end vars */
   let v1 = a*(-2*a+b+c)*(-2*a+b+c);
   let v2 = b*(a-2*b+c)*(a-2*b+c);
   let v3 = (a+b-2*c)*(a+b-2*c)*c;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X679(orbit) {
   /* begin vars */

   /* end vars */
   let v1 = a/(-2*a+b+c)*(-2*a+b+c);
   let v2 = b/(a-2*b+c)*(a-2*b+c);
   let v3 = c/(a+b-2*c)*(a+b-2*c);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X680(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let SC=(a2+b2-c2)/2;
   let SB=(c2+a2-b2)/2;
   let SA=(b2+c2-a2)/2;
   let SC2=SC*SC;
   let c3=c2*c;
   let SB2=SB*SB;
   let b3=b2*b;
   let SA2=SA*SA;
   let a3=a2*a;
   /* end vars */
   let v1 = a3*SA2*(-(b3*SB2)+c3*SC2);
   let v2 = b3*SB2*(a3*SA2-c3*SC2);
   let v3 = c3*(-(a3*SA2)+b3*SB2)*SC2;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X681(orbit) {
   /* begin vars */
   let a2=a*a;
   let c2=c*c;
   let b2=b*b;
   let SC=(a2+b2-c2)/2;
   let SB=(c2+a2-b2)/2;
   let SA=(b2+c2-a2)/2;
   let a3=a2*a;
   let SC2=SC*SC;
   let c3=c2*c;
   let SB2=SB*SB;
   let b3=b2*b;
   let SA2=SA*SA;
   /* end vars */
   let v1 = 1/(a*SA2*(-(b3*SB2)+c3*SC2));
   let v2 = 1/(b*SB2*(a3*SA2-c3*SC2));
   let v3 = 1/(c*(-(a3*SA2)+b3*SB2)*SC2);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X682(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let SA=(b2+c2-a2)/2;
   let SC=(a2+b2-c2)/2;
   let SB=(c2+a2-b2)/2;
   let c4=c2*c2;
   let SA2=SA*SA;
   let b4=b2*b2;
   let SC2=SC*SC;
   let SB2=SB*SB;
   let a4=a2*a2;
   /* end vars */
   let v1 = a4*SA*(a2*SA+SB2+SC2);
   let v2 = b4*SB*(SA2+b2*SB+SC2);
   let v3 = c4*SC*(SA2+SB2+c2*SC);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X683(orbit) {
   /* begin vars */
   let a2=a*a;
   let c2=c*c;
   let b2=b*b;
   let SA=(b2+c2-a2)/2;
   let SC=(a2+b2-c2)/2;
   let SB=(c2+a2-b2)/2;
   /* end vars */
   let v1 = b2*c2*SB*SC*(SA*(SA+SB)+a2*SC)*(a2*SB+SA*(SA+SC));
   let v2 = a2*c2*SA*SC*(SB*(SA+SB)+b2*SC)*(b2*SA+SB*(SB+SC));
   let v3 = a2*b2*SA*SB*(c2*SB+SC*(SA+SC))*(c2*SA+SC*(SB+SC));
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X684(orbit) {
   /* begin vars */
   let a2=a*a;
   let c2=c*c;
   let b2=b*b;
   let a4=a2*a2;
   let c4=c2*c2;
   let b4=b2*b2;
   /* end vars */
   let v1 = a2*(b2+c2+b4/(-a2+c2)+c4/(-a2+b2));
   let v2 = b2*(a2+c2+a4/(-b2+c2)+c4/(a2-b2));
   let v3 = (a2+b2+b4/(a2-c2)+a4/(b2-c2))*c2;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X685(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let SC=(a2+b2-c2)/2;
   let SB=(c2+a2-b2)/2;
   let SA=(b2+c2-a2)/2;
   let SC2=SC*SC;
   let SB2=SB*SB;
   let SA2=SA*SA;
   /* end vars */
   let v1 = 1/((-b2+c2)*SA*(SA2-SB*SC));
   let v2 = 1/((a2-c2)*SB*(SB2-SA*SC));
   let v3 = 1/((-a2+b2)*SC*(-(SA*SB)+SC2));
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X686(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let SC=(a2+b2-c2)/2;
   let SB=(c2+a2-b2)/2;
   let SA=(b2+c2-a2)/2;
   let SC2=SC*SC;
   let SB2=SB*SB;
   let SA2=SA*SA;
   /* end vars */
   let v1 = a2*SA*(SB-SC)*(-(SA*(SB-SC)*(SB-SC))+a2*(SA2-SB*SC));
   let v2 = b2*SB*(-SA+SC)*(-(SB*(-SA+SC)*(-SA+SC))+b2*(SB2-SA*SC));
   let v3 = c2*(SA-SB)*SC*(-((SA-SB)*(SA-SB)*SC)+c2*(-(SA*SB)+SC2));
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X687(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let SC=(a2+b2-c2)/2;
   let SB=(c2+a2-b2)/2;
   let SA=(b2+c2-a2)/2;
   let SC2=SC*SC;
   let SB2=SB*SB;
   let SA2=SA*SA;
   /* end vars */
   let v1 = 1/(SA*(SB-SC)*(-(SA*(SB-SC)*(SB-SC))+a2*(SA2-SB*SC)));
   let v2 = 1/(SB*(-SA+SC)*(-(SB*(-SA+SC)*(-SA+SC))+b2*(SB2-SA*SC)));
   let v3 = 1/((SA-SB)*SC*(-((SA-SB)*(SA-SB)*SC)+c2*(-(SA*SB)+SC2)));
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X688(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let c4=c2*c2;
   let b4=b2*b2;
   let a4=a2*a2;
   /* end vars */
   let v1 = a4*(b4-c4);
   let v2 = b4*(-a4+c4);
   let v3 = (a4-b4)*c4;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X689(orbit) {
   /* begin vars */
   let a2=a*a;
   let c2=c*c;
   let b2=b*b;
   let a4=a2*a2;
   let c4=c2*c2;
   let b4=b2*b2;
   /* end vars */
   let v1 = 1/(a2*(b4-c4));
   let v2 = 1/(b2*(-a4+c4));
   let v3 = 1/((a4-b4)*c2);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X690(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   /* end vars */
   let v1 = (2*a2-b2-c2)*(b2-c2);
   let v2 = (-a2+2*b2-c2)*(-a2+c2);
   let v3 = (a2-b2)*(-a2-b2+2*c2);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X691(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   /* end vars */
   let v1 = a2/((2*a2-b2-c2)*(b2-c2));
   let v2 = b2/((-a2+2*b2-c2)*(-a2+c2));
   let v3 = c2/((a2-b2)*(-a2-b2+2*c2));
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X692(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let c3=c2*c;
   let b3=b2*b;
   let a3=a2*a;
   /* end vars */
   let v1 = a3/(b-c);
   let v2 = b3/(-a+c);
   let v3 = c3/(a-b);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X693(orbit) {
   /* begin vars */

   /* end vars */
   let v1 = (b-c)/a;
   let v2 = (-a+c)/b;
   let v3 = (a-b)/c;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X694(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let c4=c2*c2;
   let b4=b2*b2;
   let a4=a2*a2;
   /* end vars */
   let v1 = a2/(a4-b2*c2);
   let v2 = b2/(b4-a2*c2);
   let v3 = c2/(-(a2*b2)+c4);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X695(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let c4=c2*c2;
   let b4=b2*b2;
   let a4=a2*a2;
   /* end vars */
   let v1 = a2/(a4+b2*c2);
   let v2 = b2/(b4+a2*c2);
   let v3 = c2/(a2*b2+c4);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X696(orbit) {
   /* begin vars */
   let a2=a*a;
   let c2=c*c;
   let b2=b*b;
   let a4=a2*a2;
   let c4=c2*c2;
   let b4=b2*b2;
   /* end vars */
   let v1 = a*b4-b4*c+a*c4-b*c4;
   let v2 = a4*b-a4*c-a*c4+b*c4;
   let v3 = -(a4*b)-a*b4+a4*c+b4*c;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X697(orbit) {
   /* begin vars */
   let a2=a*a;
   let c2=c*c;
   let b2=b*b;
   let a4=a2*a2;
   let c4=c2*c2;
   let b4=b2*b2;
   /* end vars */
   let v1 = a2/(a*b4-b4*c+a*c4-b*c4);
   let v2 = b2/(a4*b-a4*c-a*c4+b*c4);
   let v3 = c2/(-(a4*b)-a*b4+a4*c+b4*c);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X698(orbit) {
   /* begin vars */
   let a2=a*a;
   let c2=c*c;
   let b2=b*b;
   let a4=a2*a2;
   let c4=c2*c2;
   let b4=b2*b2;
   /* end vars */
   let v1 = a2*b4-b4*c2+a2*c4-b2*c4;
   let v2 = a4*b2-a4*c2-a2*c4+b2*c4;
   let v3 = -(a4*b2)-a2*b4+a4*c2+b4*c2;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X699(orbit) {
   /* begin vars */
   let a2=a*a;
   let c2=c*c;
   let b2=b*b;
   let a4=a2*a2;
   let c4=c2*c2;
   let b4=b2*b2;
   /* end vars */
   let v1 = a2/(a2*b4-b4*c2+a2*c4-b2*c4);
   let v2 = b2/(a4*b2-a4*c2-a2*c4+b2*c4);
   let v3 = c2/(-(a4*b2)-a2*b4+a4*c2+b4*c2);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X700(orbit) {
   /* begin vars */
   let a2=a*a;
   let b2=b*b;
   let c2=c*c;
   let a4=a2*a2;
   let b3=b2*b;
   let c4=c2*c2;
   let c3=c2*c;
   let b4=b2*b2;
   let a3=a2*a;
   /* end vars */
   let v1 = a3*b4-b4*c3+a3*c4-b3*c4;
   let v2 = a4*b3-a4*c3-a3*c4+b3*c4;
   let v3 = -(a4*b3)-a3*b4+a4*c3+b4*c3;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X701(orbit) {
   /* begin vars */
   let a2=a*a;
   let b2=b*b;
   let c2=c*c;
   let a4=a2*a2;
   let b3=b2*b;
   let c4=c2*c2;
   let c3=c2*c;
   let b4=b2*b2;
   let a3=a2*a;
   /* end vars */
   let v1 = a2/(a3*b4-b4*c3+a3*c4-b3*c4);
   let v2 = b2/(a4*b3-a4*c3-a3*c4+b3*c4);
   let v3 = c2/(-(a4*b3)-a3*b4+a4*c3+b4*c3);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X702(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let c4=c2*c2;
   let b4=b2*b2;
   let a4=a2*a2;
   /* end vars */
   let v1 = 2/a4-1/b4-1/c4;
   let v2 = -(1/a4)+2/b4-1/c4;
   let v3 = -(1/a4)-1/b4+2/c4;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X703(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let c4=c2*c2;
   let b4=b2*b2;
   let a4=a2*a2;
   /* end vars */
   let v1 = a2/(a4*b4+a4*c4-2*b4*c4);
   let v2 = b2/(a4*b4-2*a4*c4+b4*c4);
   let v3 = c2/(-2*a4*b4+a4*c4+b4*c4);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X704(orbit) {
   /* begin vars */
   let a2=a*a;
   let c2=c*c;
   let c3=c2*c;
   let b2=b*b;
   let b3=b2*b;
   let a3=a2*a;
   let a4=a2*a2;
   let c5=c2*c3;
   let b5=b2*b3;
   let c4=c2*c2;
   let b4=b2*b2;
   let a5=a2*a3;
   /* end vars */
   let v1 = a5*b4+a5*c4-b5*c4-b4*c5;
   let v2 = a4*b5-a5*c4+b5*c4-a4*c5;
   let v3 = -(a5*b4)-a4*b5+a4*c5+b4*c5;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X705(orbit) {
   /* begin vars */
   let a2=a*a;
   let c2=c*c;
   let c3=c2*c;
   let b2=b*b;
   let b3=b2*b;
   let a3=a2*a;
   let a4=a2*a2;
   let c5=c2*c3;
   let b5=b2*b3;
   let c4=c2*c2;
   let b4=b2*b2;
   let a5=a2*a3;
   /* end vars */
   let v1 = a2/(a5*b4+a5*c4-b5*c4-b4*c5);
   let v2 = b2/(a4*b5-a5*c4+b5*c4-a4*c5);
   let v3 = c2/(-(a5*b4)-a4*b5+a4*c5+b4*c5);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X706(orbit) {
   /* begin vars */
   let a2=a*a;
   let c2=c*c;
   let c4=c2*c2;
   let b2=b*b;
   let b4=b2*b2;
   let a4=a2*a2;
   let c6=c2*c4;
   let b6=b2*b4;
   let a6=a2*a4;
   /* end vars */
   let v1 = a6*b4+a6*c4-b6*c4-b4*c6;
   let v2 = a4*b6-a6*c4+b6*c4-a4*c6;
   let v3 = -(a6*b4)-a4*b6+a4*c6+b4*c6;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X707(orbit) {
   /* begin vars */
   let a2=a*a;
   let c2=c*c;
   let c4=c2*c2;
   let b2=b*b;
   let b4=b2*b2;
   let a4=a2*a2;
   let c6=c2*c4;
   let b6=b2*b4;
   let a6=a2*a4;
   /* end vars */
   let v1 = a2/(a6*b4+a6*c4-b6*c4-b4*c6);
   let v2 = b2/(a4*b6-a6*c4+b6*c4-a4*c6);
   let v3 = c2/(-(a6*b4)-a4*b6+a4*c6+b4*c6);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X708(orbit) {
   /* begin vars */
   let a2=a*a;
   let c2=c*c;
   let c4=c2*c2;
   let c6=c2*c4;
   let b2=b*b;
   let b4=b2*b2;
   let b6=b2*b4;
   let a4=a2*a2;
   let a6=a2*a4;
   let c7=c*c6;
   let b7=b*b6;
   let a7=a*a6;
   /* end vars */
   let v1 = a7*b4+a7*c4-b7*c4-b4*c7;
   let v2 = a4*b7-a7*c4+b7*c4-a4*c7;
   let v3 = -(a7*b4)-a4*b7+a4*c7+b4*c7;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X709(orbit) {
   /* begin vars */
   let a2=a*a;
   let c2=c*c;
   let c4=c2*c2;
   let c6=c2*c4;
   let b2=b*b;
   let b4=b2*b2;
   let b6=b2*b4;
   let a4=a2*a2;
   let a6=a2*a4;
   let c7=c*c6;
   let b7=b*b6;
   let a7=a*a6;
   /* end vars */
   let v1 = a2/(a7*b4+a7*c4-b7*c4-b4*c7);
   let v2 = b2/(a4*b7-a7*c4+b7*c4-a4*c7);
   let v3 = c2/(-(a7*b4)-a4*b7+a4*c7+b4*c7);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X710(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let c4=c2*c2;
   let b4=b2*b2;
   let a4=a2*a2;
   /* end vars */
   let v1 = (a2-b*c)*(a2+b*c)*(a4+b2*c2)*(b4+c4);
   let v2 = (b2-a*c)*(b2+a*c)*(b4+a2*c2)*(a4+c4);
   let v3 = (a4+b4)*(-(a*b)+c2)*(a*b+c2)*(a2*b2+c4);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X711(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let c4=c2*c2;
   let b4=b2*b2;
   let a4=a2*a2;
   /* end vars */
   let v1 = a2/((a2-b*c)*(a2+b*c)*(a4+b2*c2)*(b4+c4));
   let v2 = b2/((b2-a*c)*(b2+a*c)*(b4+a2*c2)*(a4+c4));
   let v3 = c2/((a4+b4)*(-(a*b)+c2)*(a*b+c2)*(a2*b2+c4));
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X712(orbit) {
   /* begin vars */
   let a2=a*a;
   let c2=c*c;
   let b2=b*b;
   let a3=a2*a;
   let c3=c2*c;
   let b3=b2*b;
   /* end vars */
   let v1 = a*b3-b3*c+a*c3-b*c3;
   let v2 = a3*b-a3*c-a*c3+b*c3;
   let v3 = -(a3*b)-a*b3+a3*c+b3*c;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X713(orbit) {
   /* begin vars */
   let a2=a*a;
   let c2=c*c;
   let b2=b*b;
   let a3=a2*a;
   let c3=c2*c;
   let b3=b2*b;
   /* end vars */
   let v1 = a2/(a*b3-b3*c+a*c3-b*c3);
   let v2 = b2/(a3*b-a3*c-a*c3+b*c3);
   let v3 = c2/(-(a3*b)-a*b3+a3*c+b3*c);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X714(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   /* end vars */
   let v1 = (b+c)*(a2*b2-a2*b*c+a2*c2-b2*c2);
   let v2 = (a+c)*(a2*b2-a*b2*c-a2*c2+b2*c2);
   let v3 = (a+b)*(-(a2*b2)+a2*c2-a*b*c2+b2*c2);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X715(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   /* end vars */
   let v1 = a2/((b+c)*(a2*b2-a2*b*c+a2*c2-b2*c2));
   let v2 = b2/((a+c)*(a2*b2-a*b2*c-a2*c2+b2*c2));
   let v3 = c2/((a+b)*(-(a2*b2)+a2*c2-a*b*c2+b2*c2));
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X716(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let c3=c2*c;
   let b3=b2*b;
   let a3=a2*a;
   /* end vars */
   let v1 = 2/a3-1/b3-1/c3;
   let v2 = -(1/a3)+2/b3-1/c3;
   let v3 = -(1/a3)-1/b3+2/c3;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X717(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let c3=c2*c;
   let b3=b2*b;
   let a3=a2*a;
   /* end vars */
   let v1 = a2/(a3*b3+a3*c3-2*b3*c3);
   let v2 = b2/(a3*b3-2*a3*c3+b3*c3);
   let v3 = c2/(-2*a3*b3+a3*c3+b3*c3);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X718(orbit) {
   /* begin vars */
   let c2=c*c;
   let a2=a*a;
   let b2=b*b;
   let c4=c2*c2;
   let a3=a2*a;
   let b4=b2*b2;
   let c3=c2*c;
   let b3=b2*b;
   let a4=a2*a2;
   /* end vars */
   let v1 = (b+c)*(a4*b2-a4*b*c+a4*c2-b3*c3);
   let v2 = (a+c)*(a2*b4-a*b4*c+b4*c2-a3*c3);
   let v3 = (a+b)*(-(a3*b3)+a2*c4-a*b*c4+b2*c4);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X719(orbit) {
   /* begin vars */
   let c2=c*c;
   let a2=a*a;
   let b2=b*b;
   let c4=c2*c2;
   let a3=a2*a;
   let b4=b2*b2;
   let c3=c2*c;
   let b3=b2*b;
   let a4=a2*a2;
   /* end vars */
   let v1 = a2/((b+c)*(a4*b2-a4*b*c+a4*c2-b3*c3));
   let v2 = b2/((a+c)*(a2*b4-a*b4*c+b4*c2-a3*c3));
   let v3 = c2/((a+b)*(-(a3*b3)+a2*c4-a*b*c4+b2*c4));
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X720(orbit) {
   /* begin vars */
   let a2=a*a;
   let c2=c*c;
   let c3=c2*c;
   let b2=b*b;
   let b3=b2*b;
   let a3=a2*a;
   let c5=c2*c3;
   let b5=b2*b3;
   let a5=a2*a3;
   /* end vars */
   let v1 = a5*b3+a5*c3-b5*c3-b3*c5;
   let v2 = a3*b5-a5*c3+b5*c3-a3*c5;
   let v3 = -(a5*b3)-a3*b5+a3*c5+b3*c5;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X721(orbit) {
   /* begin vars */
   let a2=a*a;
   let c2=c*c;
   let c3=c2*c;
   let b2=b*b;
   let b3=b2*b;
   let a3=a2*a;
   let c5=c2*c3;
   let b5=b2*b3;
   let a5=a2*a3;
   /* end vars */
   let v1 = a2/(a5*b3+a5*c3-b5*c3-b3*c5);
   let v2 = b2/(a3*b5-a5*c3+b5*c3-a3*c5);
   let v3 = c2/(-(a5*b3)-a3*b5+a3*c5+b3*c5);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X722(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let c4=c2*c2;
   let b4=b2*b2;
   let a4=a2*a2;
   /* end vars */
   let v1 = (b+c)*(a2-b*c)*(b2-b*c+c2)*(a4+a2*b*c+b2*c2);
   let v2 = (a+c)*(b2-a*c)*(a2-a*c+c2)*(b4+a*b2*c+a2*c2);
   let v3 = (a+b)*(a2-a*b+b2)*(-(a*b)+c2)*(a2*b2+a*b*c2+c4);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X723(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let c4=c2*c2;
   let b4=b2*b2;
   let a4=a2*a2;
   /* end vars */
   let v1 = a2/((b+c)*(a2-b*c)*(b2-b*c+c2)*(a4+a2*b*c+b2*c2));
   let v2 = b2/((a+c)*(b2-a*c)*(a2-a*c+c2)*(b4+a*b2*c+a2*c2));
   let v3 = c2/((a+b)*(a2-a*b+b2)*(-(a*b)+c2)*(a2*b2+a*b*c2+c4));
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X724(orbit) {
   /* begin vars */
   let a2=a*a;
   let c2=c*c;
   let c4=c2*c2;
   let c6=c2*c4;
   let b2=b*b;
   let b4=b2*b2;
   let b6=b2*b4;
   let a4=a2*a2;
   let a6=a2*a4;
   let a3=a2*a;
   let c7=c*c6;
   let b7=b*b6;
   let c3=c2*c;
   let b3=b2*b;
   let a7=a*a6;
   /* end vars */
   let v1 = a7*b3+a7*c3-b7*c3-b3*c7;
   let v2 = a3*b7-a7*c3+b7*c3-a3*c7;
   let v3 = -(a7*b3)-a3*b7+a3*c7+b3*c7;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X725(orbit) {
   /* begin vars */
   let a2=a*a;
   let c2=c*c;
   let c4=c2*c2;
   let c6=c2*c4;
   let b2=b*b;
   let b4=b2*b2;
   let b6=b2*b4;
   let a4=a2*a2;
   let a6=a2*a4;
   let a3=a2*a;
   let c7=c*c6;
   let b7=b*b6;
   let c3=c2*c;
   let b3=b2*b;
   let a7=a*a6;
   /* end vars */
   let v1 = a2/(a7*b3+a7*c3-b7*c3-b3*c7);
   let v2 = b2/(a3*b7-a7*c3+b7*c3-a3*c7);
   let v3 = c2/(-(a7*b3)-a3*b7+a3*c7+b3*c7);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X726(orbit) {
   /* begin vars */
   let a2=a*a;
   let c2=c*c;
   let b2=b*b;
   /* end vars */
   let v1 = a*b2-b2*c+a*c2-b*c2;
   let v2 = a2*b-a2*c-a*c2+b*c2;
   let v3 = -(a2*b)-a*b2+a2*c+b2*c;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X727(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   /* end vars */
   let v1 = a2/(a*b2-b2*c+a*c2-b*c2);
   let v2 = b2/(a2*b-a2*c-a*c2+b*c2);
   let v3 = c2/(-(a2*b)-a*b2+a2*c+b2*c);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X728(orbit) {
   /* begin vars */

   /* end vars */
   let v1 = a*Math.pow(-a+b+c,3);
   let v2 = b*Math.pow(a-b+c,3);
   let v3 = Math.pow(a+b-c,3)*c;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X729(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   /* end vars */
   let v1 = a2/(a2*b2+a2*c2-2*b2*c2);
   let v2 = b2/(a2*b2-2*a2*c2+b2*c2);
   let v3 = c2/(-2*a2*b2+a2*c2+b2*c2);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X730(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let c3=c2*c;
   let b3=b2*b;
   let a3=a2*a;
   /* end vars */
   let v1 = a3*b2+a3*c2-b3*c2-b2*c3;
   let v2 = a2*b3-a3*c2+b3*c2-a2*c3;
   let v3 = -(a3*b2)-a2*b3+a2*c3+b2*c3;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X731(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let c3=c2*c;
   let b3=b2*b;
   let a3=a2*a;
   /* end vars */
   let v1 = a2/(a3*b2+a3*c2-b3*c2-b2*c3);
   let v2 = b2/(a2*b3-a3*c2+b3*c2-a2*c3);
   let v3 = c2/(-(a3*b2)-a2*b3+a2*c3+b2*c3);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X732(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   /* end vars */
   let v1 = (a2-b*c)*(a2+b*c)*(b2+c2);
   let v2 = (b2-a*c)*(b2+a*c)*(a2+c2);
   let v3 = (a2+b2)*(-(a*b)+c2)*(a*b+c2);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X733(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   /* end vars */
   let v1 = a2/((a2-b*c)*(a2+b*c)*(b2+c2));
   let v2 = b2/((b2-a*c)*(b2+a*c)*(a2+c2));
   let v3 = c2/((a2+b2)*(-(a*b)+c2)*(a*b+c2));
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X734(orbit) {
   /* begin vars */
   let c2=c*c;
   let c3=c2*c;
   let b2=b*b;
   let b3=b2*b;
   let a2=a*a;
   let a3=a2*a;
   let c5=c2*c3;
   let b5=b2*b3;
   let a5=a2*a3;
   /* end vars */
   let v1 = a5*b2+a5*c2-b5*c2-b2*c5;
   let v2 = a2*b5-a5*c2+b5*c2-a2*c5;
   let v3 = -(a5*b2)-a2*b5+a2*c5+b2*c5;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X735(orbit) {
   /* begin vars */
   let c2=c*c;
   let c3=c2*c;
   let b2=b*b;
   let b3=b2*b;
   let a2=a*a;
   let a3=a2*a;
   let c5=c2*c3;
   let b5=b2*b3;
   let a5=a2*a3;
   /* end vars */
   let v1 = a2/(a5*b2+a5*c2-b5*c2-b2*c5);
   let v2 = b2/(a2*b5-a5*c2+b5*c2-a2*c5);
   let v3 = c2/(-(a5*b2)-a2*b5+a2*c5+b2*c5);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X736(orbit) {
   /* begin vars */
   let c2=c*c;
   let c4=c2*c2;
   let b2=b*b;
   let b4=b2*b2;
   let a2=a*a;
   let a4=a2*a2;
   let c6=c2*c4;
   let b6=b2*b4;
   let a6=a2*a4;
   /* end vars */
   let v1 = a6*b2+a6*c2-b6*c2-b2*c6;
   let v2 = a2*b6-a6*c2+b6*c2-a2*c6;
   let v3 = -(a6*b2)-a2*b6+a2*c6+b2*c6;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X737(orbit) {
   /* begin vars */
   let c2=c*c;
   let c4=c2*c2;
   let b2=b*b;
   let b4=b2*b2;
   let a2=a*a;
   let a4=a2*a2;
   let c6=c2*c4;
   let b6=b2*b4;
   let a6=a2*a4;
   /* end vars */
   let v1 = a2/(a6*b2+a6*c2-b6*c2-b2*c6);
   let v2 = b2/(a2*b6-a6*c2+b6*c2-a2*c6);
   let v3 = c2/(-(a6*b2)-a2*b6+a2*c6+b2*c6);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X738(orbit) {
   /* begin vars */

   /* end vars */
   let v1 = a/Math.pow(-a+b+c,3);
   let v2 = b/Math.pow(a-b+c,3);
   let v3 = c/Math.pow(a+b-c,3);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X739(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   /* end vars */
   let v1 = a2/(a*b+a*c-2*b*c);
   let v2 = b2/(a*b-2*a*c+b*c);
   let v3 = c2/(-2*a*b+a*c+b*c);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X740(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   /* end vars */
   let v1 = (b+c)*(a2-b*c);
   let v2 = (a+c)*(b2-a*c);
   let v3 = (a+b)*(-(a*b)+c2);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X741(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   /* end vars */
   let v1 = a2/((b+c)*(a2-b*c));
   let v2 = b2/((a+c)*(b2-a*c));
   let v3 = c2/((a+b)*(-(a*b)+c2));
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X742(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let c3=c2*c;
   let b3=b2*b;
   let a3=a2*a;
   /* end vars */
   let v1 = a3*b+a3*c-b3*c-b*c3;
   let v2 = a*b3-a3*c+b3*c-a*c3;
   let v3 = -(a3*b)-a*b3+a*c3+b*c3;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X743(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let c3=c2*c;
   let b3=b2*b;
   let a3=a2*a;
   /* end vars */
   let v1 = a2/(a3*b+a3*c-b3*c-b*c3);
   let v2 = b2/(a*b3-a3*c+b3*c-a*c3);
   let v3 = c2/(-(a3*b)-a*b3+a*c3+b*c3);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X744(orbit) {
   /* begin vars */
   let c2=c*c;
   let a2=a*a;
   let b2=b*b;
   let c4=c2*c2;
   let a3=a2*a;
   let b4=b2*b2;
   let c3=c2*c;
   let b3=b2*b;
   let a4=a2*a2;
   /* end vars */
   let v1 = (b+c)*(a4-b3*c+b2*c2-b*c3);
   let v2 = (a+c)*(b4-a3*c+a2*c2-a*c3);
   let v3 = (a+b)*(-(a3*b)+a2*b2-a*b3+c4);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X745(orbit) {
   /* begin vars */
   let c2=c*c;
   let a2=a*a;
   let b2=b*b;
   let c4=c2*c2;
   let a3=a2*a;
   let b4=b2*b2;
   let c3=c2*c;
   let b3=b2*b;
   let a4=a2*a2;
   /* end vars */
   let v1 = a2/((b+c)*(a4-b3*c+b2*c2-b*c3));
   let v2 = b2/((a+c)*(b4-a3*c+a2*c2-a*c3));
   let v3 = c2/((a+b)*(-(a3*b)+a2*b2-a*b3+c4));
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X746(orbit) {
   /* begin vars */
   let c2=c*c;
   let c3=c2*c;
   let b2=b*b;
   let b3=b2*b;
   let a2=a*a;
   let a3=a2*a;
   let c5=c2*c3;
   let b5=b2*b3;
   let a5=a2*a3;
   /* end vars */
   let v1 = a5*b+a5*c-b5*c-b*c5;
   let v2 = a*b5-a5*c+b5*c-a*c5;
   let v3 = -(a5*b)-a*b5+a*c5+b*c5;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X747(orbit) {
   /* begin vars */
   let c2=c*c;
   let c3=c2*c;
   let b2=b*b;
   let b3=b2*b;
   let a2=a*a;
   let a3=a2*a;
   let c5=c2*c3;
   let b5=b2*b3;
   let a5=a2*a3;
   /* end vars */
   let v1 = a2/(a5*b+a5*c-b5*c-b*c5);
   let v2 = b2/(a*b5-a5*c+b5*c-a*c5);
   let v3 = c2/(-(a5*b)-a*b5+a*c5+b*c5);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X748(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let c3=c2*c;
   let b3=b2*b;
   let a3=a2*a;
   /* end vars */
   let v1 = a3-2*a*b*c;
   let v2 = b3-2*a*b*c;
   let v3 = -2*a*b*c+c3;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X749(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let c3=c2*c;
   let b3=b2*b;
   let a3=a2*a;
   /* end vars */
   let v1 = a2/(a3-2*a*b*c);
   let v2 = b2/(b3-2*a*b*c);
   let v3 = c2/(-2*a*b*c+c3);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X750(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let c3=c2*c;
   let b3=b2*b;
   let a3=a2*a;
   /* end vars */
   let v1 = a3+2*a*b*c;
   let v2 = b3+2*a*b*c;
   let v3 = 2*a*b*c+c3;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X751(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let c3=c2*c;
   let b3=b2*b;
   let a3=a2*a;
   /* end vars */
   let v1 = a2/(a3+2*a*b*c);
   let v2 = b2/(b3+2*a*b*c);
   let v3 = c2/(2*a*b*c+c3);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X752(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let c3=c2*c;
   let b3=b2*b;
   let a3=a2*a;
   /* end vars */
   let v1 = 2*a3-b3-c3;
   let v2 = -a3+2*b3-c3;
   let v3 = -a3-b3+2*c3;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X753(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let c3=c2*c;
   let b3=b2*b;
   let a3=a2*a;
   /* end vars */
   let v1 = a2/(2*a3-b3-c3);
   let v2 = b2/(-a3+2*b3-c3);
   let v3 = c2/(-a3-b3+2*c3);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X754(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let c4=c2*c2;
   let b4=b2*b2;
   let a4=a2*a2;
   /* end vars */
   let v1 = 2*a4-b4-c4;
   let v2 = -a4+2*b4-c4;
   let v3 = -a4-b4+2*c4;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X755(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let c4=c2*c2;
   let b4=b2*b2;
   let a4=a2*a2;
   /* end vars */
   let v1 = a2/(2*a4-b4-c4);
   let v2 = b2/(-a4+2*b4-c4);
   let v3 = c2/(-a4-b4+2*c4);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X756(orbit) {
   /* begin vars */

   /* end vars */
   let v1 = a*(b+c)*(b+c);
   let v2 = b*(a+c)*(a+c);
   let v3 = (a+b)*(a+b)*c;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X757(orbit) {
   /* begin vars */

   /* end vars */
   let v1 = a/(b+c)*(b+c);
   let v2 = b/(a+c)*(a+c);
   let v3 = c/(a+b)*(a+b);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X758(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let c3=c2*c;
   let b3=b2*b;
   let a3=a2*a;
   /* end vars */
   let v1 = a3*(b+c)-a*(b3+c3);
   let v2 = b3*(a+c)-b*(a3+c3);
   let v3 = -((a3+b3)*c)+(a+b)*c3;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X759(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let c3=c2*c;
   let b3=b2*b;
   let a3=a2*a;
   /* end vars */
   let v1 = a2/(a3*(b+c)-a*(b3+c3));
   let v2 = b2/(b3*(a+c)-b*(a3+c3));
   let v3 = c2/(-((a3+b3)*c)+(a+b)*c3);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X760(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let c4=c2*c2;
   let b4=b2*b2;
   let a4=a2*a2;
   /* end vars */
   let v1 = a4*(b+c)-a*(b4+c4);
   let v2 = b4*(a+c)-b*(a4+c4);
   let v3 = -((a4+b4)*c)+(a+b)*c4;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X761(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let c4=c2*c2;
   let b4=b2*b2;
   let a4=a2*a2;
   /* end vars */
   let v1 = a2/(a4*(b+c)-a*(b4+c4));
   let v2 = b2/(b4*(a+c)-b*(a4+c4));
   let v3 = c2/(-((a4+b4)*c)+(a+b)*c4);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X762(orbit) {
   /* begin vars */

   /* end vars */
   let v1 = a*Math.pow(b+c,3);
   let v2 = b*Math.pow(a+c,3);
   let v3 = Math.pow(a+b,3)*c;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X763(orbit) {
   /* begin vars */

   /* end vars */
   let v1 = a/Math.pow(b+c,3);
   let v2 = b/Math.pow(a+c,3);
   let v3 = c/Math.pow(a+b,3);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X764(orbit) {
   /* begin vars */

   /* end vars */
   let v1 = a*Math.pow(b-c,3);
   let v2 = b*Math.pow(-a+c,3);
   let v3 = Math.pow(a-b,3)*c;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X765(orbit) {
   /* begin vars */

   /* end vars */
   let v1 = a/(b-c)*(b-c);
   let v2 = b/(-a+c)*(-a+c);
   let v3 = c/(a-b)*(a-b);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X766(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let c4=c2*c2;
   let b4=b2*b2;
   let a3=a2*a;
   let c3=c2*c;
   let b3=b2*b;
   let a4=a2*a2;
   /* end vars */
   let v1 = a4*(b3+c3)-a3*(b4+c4);
   let v2 = b4*(a3+c3)-b3*(a4+c4);
   let v3 = -((a4+b4)*c3)+(a3+b3)*c4;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X767(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let c4=c2*c2;
   let b4=b2*b2;
   let a3=a2*a;
   let c3=c2*c;
   let b3=b2*b;
   let a4=a2*a2;
   /* end vars */
   let v1 = a2/(a4*(b3+c3)-a3*(b4+c4));
   let v2 = b2/(b4*(a3+c3)-b3*(a4+c4));
   let v3 = c2/(-((a4+b4)*c3)+(a3+b3)*c4);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X768(orbit) {
   /* begin vars */
   let a2=a*a;
   let c2=c*c;
   let b2=b*b;
   let a3=a2*a;
   let c3=c2*c;
   let b3=b2*b;
   /* end vars */
   let v1 = (b-c)*(a*b3+a*b2*c+b3*c+a*b*c2+b2*c2+a*c3+b*c3);
   let v2 = (-a+c)*(a3*b+a3*c+a2*b*c+a2*c2+a*b*c2+a*c3+b*c3);
   let v3 = (a-b)*(a3*b+a2*b2+a*b3+a3*c+a2*b*c+a*b2*c+b3*c);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X769(orbit) {
   /* begin vars */
   let a2=a*a;
   let c2=c*c;
   let b2=b*b;
   let a3=a2*a;
   let c3=c2*c;
   let b3=b2*b;
   /* end vars */
   let v1 = a2/((b-c)*(a*b3+a*b2*c+b3*c+a*b*c2+b2*c2+a*c3+b*c3));
   let v2 = b2/((-a+c)*(a3*b+a3*c+a2*b*c+a2*c2+a*b*c2+a*c3+b*c3));
   let v3 = c2/((a-b)*(a3*b+a2*b2+a*b3+a3*c+a2*b*c+a*b2*c+b3*c));
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X770(orbit) {
   /* begin vars */
   let c2=c*c;
   let a2=a*a;
   let b2=b*b;
   let SA=(b2+c2-a2)/2;
   let SC=(a2+b2-c2)/2;
   let SB=(c2+a2-b2)/2;
   let c4=c2*c2;
   let SA2=SA*SA;
   let b4=b2*b2;
   let SC2=SC*SC;
   let SB2=SB*SB;
   let a4=a2*a2;
   /* end vars */
   let v1 = a*(a-b-c)*(b-c)*(a4*b*c-b*(b-c)*(b-c)*c*(b+c)*(b+c)+4*a2*SB*SC+4*SA*(SB2+SC2));
   let v2 = b*(-a+b-c)*(-a+c)*(a*b4*c-a*c*(-a+c)*(-a+c)*(a+c)*(a+c)+4*b2*SA*SC+4*SB*(SA2+SC2));
   let v3 = (a-b)*c*(-a-b+c)*(-(a*(a-b)*(a-b)*b*(a+b)*(a+b))+a*b*c4+4*c2*SA*SB+4*(SA2+SB2)*SC);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X771(orbit) {
   /* begin vars */
   let c2=c*c;
   let a2=a*a;
   let b2=b*b;
   let SA=(b2+c2-a2)/2;
   let SC=(a2+b2-c2)/2;
   let SB=(c2+a2-b2)/2;
   let c4=c2*c2;
   let SA2=SA*SA;
   let b4=b2*b2;
   let SC2=SC*SC;
   let SB2=SB*SB;
   let a4=a2*a2;
   /* end vars */
   let v1 = a/((a-b-c)*(b-c)*(a4*b*c-b*(b-c)*(b-c)*c*(b+c)*(b+c)+4*a2*SB*SC+4*SA*(SB2+SC2)));
   let v2 = b/((-a+b-c)*(-a+c)*(a*b4*c-a*c*(-a+c)*(-a+c)*(a+c)*(a+c)+4*b2*SA*SC+4*SB*(SA2+SC2)));
   let v3 = c/((a-b)*(-a-b+c)*(-(a*(a-b)*(a-b)*b*(a+b)*(a+b))+a*b*c4+4*c2*SA*SB+4*(SA2+SB2)*SC));
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X772(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let c3=c2*c;
   let b3=b2*b;
   let a3=a2*a;
   /* end vars */
   let v1 = (b-c)*(a3*b3+a3*b2*c+a3*b*c2+a3*c3+b3*c3);
   let v2 = (-a+c)*(a3*b3+a2*b3*c+a*b3*c2+a3*c3+b3*c3);
   let v3 = (a-b)*(a3*b3+a3*c3+a2*b*c3+a*b2*c3+b3*c3);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X773(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let c3=c2*c;
   let b3=b2*b;
   let a3=a2*a;
   /* end vars */
   let v1 = a2/((b-c)*(a3*b3+a3*b2*c+a3*b*c2+a3*c3+b3*c3));
   let v2 = b2/((-a+c)*(a3*b3+a2*b3*c+a*b3*c2+a3*c3+b3*c3));
   let v3 = c2/((a-b)*(a3*b3+a3*c3+a2*b*c3+a*b2*c3+b3*c3));
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X774(orbit) {
   /* begin vars */
   let a2=a*a;
   let c2=c*c;
   let b2=b*b;
   let SA=(b2+c2-a2)/2;
   let SC=(a2+b2-c2)/2;
   let SB=(c2+a2-b2)/2;
   let SA2=SA*SA;
   let SC2=SC*SC;
   let SB2=SB*SB;
   /* end vars */
   let v1 = a*(a2*SB*SC+SA*(SB2+SC2));
   let v2 = b*(b2*SA*SC+SB*(SA2+SC2));
   let v3 = c*(c2*SA*SB+(SA2+SB2)*SC);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X775(orbit) {
   /* begin vars */
   let a2=a*a;
   let c2=c*c;
   let b2=b*b;
   let SA=(b2+c2-a2)/2;
   let SC=(a2+b2-c2)/2;
   let SB=(c2+a2-b2)/2;
   let SA2=SA*SA;
   let SC2=SC*SC;
   let SB2=SB*SB;
   /* end vars */
   let v1 = a/(a2*SB*SC+SA*(SB2+SC2));
   let v2 = b/(b2*SA*SC+SB*(SA2+SC2));
   let v3 = c/(c2*SA*SB+(SA2+SB2)*SC);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X776(orbit) {
   /* begin vars */
   let c2=c*c;
   let c3=c2*c;
   let a2=a*a;
   let b2=b*b;
   let b3=b2*b;
   let a3=a2*a;
   let c5=c2*c3;
   let a4=a2*a2;
   let b5=b2*b3;
   let c4=c2*c2;
   let b4=b2*b2;
   let a5=a2*a3;
   /* end vars */
   let v1 = (b-c)*(a5*b3+a5*b2*c+a5*b*c2+a5*c3-b4*c4);
   let v2 = (-a+c)*(a3*b5+a2*b5*c+a*b5*c2+b5*c3-a4*c4);
   let v3 = (a-b)*(-(a4*b4)+a3*c5+a2*b*c5+a*b2*c5+b3*c5);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X777(orbit) {
   /* begin vars */
   let c2=c*c;
   let c3=c2*c;
   let a2=a*a;
   let b2=b*b;
   let b3=b2*b;
   let a3=a2*a;
   let c5=c2*c3;
   let a4=a2*a2;
   let b5=b2*b3;
   let c4=c2*c2;
   let b4=b2*b2;
   let a5=a2*a3;
   /* end vars */
   let v1 = a2/((b-c)*(a5*b3+a5*b2*c+a5*b*c2+a5*c3-b4*c4));
   let v2 = b2/((-a+c)*(a3*b5+a2*b5*c+a*b5*c2+b5*c3-a4*c4));
   let v3 = c2/((a-b)*(-(a4*b4)+a3*c5+a2*b*c5+a*b2*c5+b3*c5));
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X778(orbit) {
   /* begin vars */
   let c2=c*c;
   let c4=c2*c2;
   let a2=a*a;
   let b2=b*b;
   let b4=b2*b2;
   let a4=a2*a2;
   let c6=c2*c4;
   let b6=b2*b4;
   let a6=a2*a4;
   /* end vars */
   let v1 = (b-c)*(b+c)*(a6*b2+a6*c2-b4*c4);
   let v2 = (-a+c)*(a+c)*(a2*b6+b6*c2-a4*c4);
   let v3 = (a-b)*(a+b)*(-(a4*b4)+a2*c6+b2*c6);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X779(orbit) {
   /* begin vars */
   let c2=c*c;
   let c4=c2*c2;
   let a2=a*a;
   let b2=b*b;
   let b4=b2*b2;
   let a4=a2*a2;
   let c6=c2*c4;
   let b6=b2*b4;
   let a6=a2*a4;
   /* end vars */
   let v1 = a2/((b-c)*(b+c)*(a6*b2+a6*c2-b4*c4));
   let v2 = b2/((-a+c)*(a+c)*(a2*b6+b6*c2-a4*c4));
   let v3 = c2/((a-b)*(a+b)*(-(a4*b4)+a2*c6+b2*c6));
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X780(orbit) {
   /* begin vars */
   let c2=c*c;
   let c4=c2*c2;
   let c6=c2*c4;
   let a2=a*a;
   let a3=a2*a;
   let a4=a2*a2;
   let b2=b*b;
   let b4=b2*b2;
   let b6=b2*b4;
   let c3=c2*c;
   let b3=b2*b;
   let a6=a2*a4;
   let c7=c*c6;
   let a5=a2*a3;
   let b7=b*b6;
   let c5=c2*c3;
   let b5=b2*b3;
   let a7=a*a6;
   /* end vars */
   let v1 = (b-c)*(a7*b3+a7*b2*c+a7*b*c2+a7*c3-b6*c4-b5*c5-b4*c6);
   let v2 = (-a+c)*(a3*b7+a2*b7*c+a*b7*c2+b7*c3-a6*c4-a5*c5-a4*c6);
   let v3 = (a-b)*(-(a6*b4)-a5*b5-a4*b6+a3*c7+a2*b*c7+a*b2*c7+b3*c7);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X781(orbit) {
   /* begin vars */
   let c2=c*c;
   let c4=c2*c2;
   let c6=c2*c4;
   let a2=a*a;
   let a3=a2*a;
   let a4=a2*a2;
   let b2=b*b;
   let b4=b2*b2;
   let b6=b2*b4;
   let c3=c2*c;
   let b3=b2*b;
   let a6=a2*a4;
   let c7=c*c6;
   let a5=a2*a3;
   let b7=b*b6;
   let c5=c2*c3;
   let b5=b2*b3;
   let a7=a*a6;
   /* end vars */
   let v1 = a2/((b-c)*(a7*b3+a7*b2*c+a7*b*c2+a7*c3-b6*c4-b5*c5-b4*c6));
   let v2 = b2/((-a+c)*(a3*b7+a2*b7*c+a*b7*c2+b7*c3-a6*c4-a5*c5-a4*c6));
   let v3 = c2/((a-b)*(-(a6*b4)-a5*b5-a4*b6+a3*c7+a2*b*c7+a*b2*c7+b3*c7));
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X782(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let c4=c2*c2;
   let b4=b2*b2;
   let a4=a2*a2;
   /* end vars */
   let v1 = (a4-b2*c2)*(a4+b2*c2)*(b4-c4);
   let v2 = (b4-a2*c2)*(b4+a2*c2)*(-a4+c4);
   let v3 = (a4-b4)*(-(a2*b2)+c4)*(a2*b2+c4);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X783(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let c4=c2*c2;
   let b4=b2*b2;
   let a4=a2*a2;
   /* end vars */
   let v1 = a2/((a4-b2*c2)*(a4+b2*c2)*(b4-c4));
   let v2 = b2/((b4-a2*c2)*(b4+a2*c2)*(-a4+c4));
   let v3 = c2/((a4-b4)*(-(a2*b2)+c4)*(a2*b2+c4));
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X784(orbit) {
   /* begin vars */
   let a2=a*a;
   let c2=c*c;
   let b2=b*b;
   /* end vars */
   let v1 = (b-c)*(a*b2+a*b*c+b2*c+a*c2+b*c2);
   let v2 = (-a+c)*(a2*b+a2*c+a*b*c+a*c2+b*c2);
   let v3 = (a-b)*(a2*b+a*b2+a2*c+a*b*c+b2*c);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X785(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   /* end vars */
   let v1 = a2/((b-c)*(a*b2+a*b*c+b2*c+a*c2+b*c2));
   let v2 = b2/((-a+c)*(a2*b+a2*c+a*b*c+a*c2+b*c2));
   let v3 = c2/((a-b)*(a2*b+a*b2+a2*c+a*b*c+b2*c));
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X786(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   /* end vars */
   let v1 = (b-c)*(a2*b2+a2*b*c+a2*c2+b2*c2);
   let v2 = (-a+c)*(a2*b2+a*b2*c+a2*c2+b2*c2);
   let v3 = (a-b)*(a2*b2+a2*c2+a*b*c2+b2*c2);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X787(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   /* end vars */
   let v1 = a2/((b-c)*(a2*b2+a2*b*c+a2*c2+b2*c2));
   let v2 = b2/((-a+c)*(a2*b2+a*b2*c+a2*c2+b2*c2));
   let v3 = c2/((a-b)*(a2*b2+a2*c2+a*b*c2+b2*c2));
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X788(orbit) {
   /* begin vars */
   let a2=a*a;
   let c2=c*c;
   let b2=b*b;
   let a3=a2*a;
   let c3=c2*c;
   let b3=b2*b;
   /* end vars */
   let v1 = 1/b3-1/c3;
   let v2 = -(1/a3)+1/c3;
   let v3 = 1/a3-1/b3;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X789(orbit) {
   /* begin vars */
   let a2=a*a;
   let c2=c*c;
   let b2=b*b;
   /* end vars */
   let v1 = 1/(a*(b-c)*(b2+b*c+c2));
   let v2 = 1/(b*(-a+c)*(a2+a*c+c2));
   let v3 = 1/((a-b)*(a2+a*b+b2)*c);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X790(orbit) {
   /* begin vars */
   let c2=c*c;
   let a2=a*a;
   let b2=b*b;
   let c4=c2*c2;
   let a3=a2*a;
   let b4=b2*b2;
   let c3=c2*c;
   let b3=b2*b;
   let a4=a2*a2;
   /* end vars */
   let v1 = (b-c)*(a4*b2+a4*b*c+a4*c2-b3*c3);
   let v2 = (-a+c)*(a2*b4+a*b4*c+b4*c2-a3*c3);
   let v3 = (a-b)*(-(a3*b3)+a2*c4+a*b*c4+b2*c4);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X791(orbit) {
   /* begin vars */
   let c2=c*c;
   let a2=a*a;
   let b2=b*b;
   let c4=c2*c2;
   let a3=a2*a;
   let b4=b2*b2;
   let c3=c2*c;
   let b3=b2*b;
   let a4=a2*a2;
   /* end vars */
   let v1 = a2/((b-c)*(a4*b2+a4*b*c+a4*c2-b3*c3));
   let v2 = b2/((-a+c)*(a2*b4+a*b4*c+b4*c2-a3*c3));
   let v3 = c2/((a-b)*(-(a3*b3)+a2*c4+a*b*c4+b2*c4));
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X792(orbit) {
   /* begin vars */
   let c2=c*c;
   let c3=c2*c;
   let a2=a*a;
   let b2=b*b;
   let b3=b2*b;
   let a3=a2*a;
   let c5=c2*c3;
   let a4=a2*a2;
   let b5=b2*b3;
   let c4=c2*c2;
   let b4=b2*b2;
   let a5=a2*a3;
   /* end vars */
   let v1 = (b-c)*(a5*b2+a5*b*c+a5*c2-b4*c3-b3*c4);
   let v2 = (-a+c)*(a2*b5+a*b5*c+b5*c2-a4*c3-a3*c4);
   let v3 = (a-b)*(-(a4*b3)-a3*b4+a2*c5+a*b*c5+b2*c5);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X793(orbit) {
   /* begin vars */
   let c2=c*c;
   let c3=c2*c;
   let a2=a*a;
   let b2=b*b;
   let b3=b2*b;
   let a3=a2*a;
   let c5=c2*c3;
   let a4=a2*a2;
   let b5=b2*b3;
   let c4=c2*c2;
   let b4=b2*b2;
   let a5=a2*a3;
   /* end vars */
   let v1 = a2/((b-c)*(a5*b2+a5*b*c+a5*c2-b4*c3-b3*c4));
   let v2 = b2/((-a+c)*(a2*b5+a*b5*c+b5*c2-a4*c3-a3*c4));
   let v3 = c2/((a-b)*(-(a4*b3)-a3*b4+a2*c5+a*b*c5+b2*c5));
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X794(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let c4=c2*c2;
   let b4=b2*b2;
   let a4=a2*a2;
   /* end vars */
   let v1 = (b-c)*(a2-b*c)*(b2+b*c+c2)*(a4+a2*b*c+b2*c2);
   let v2 = (-a+c)*(b2-a*c)*(a2+a*c+c2)*(b4+a*b2*c+a2*c2);
   let v3 = (a-b)*(a2+a*b+b2)*(-(a*b)+c2)*(a2*b2+a*b*c2+c4);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X795(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let c4=c2*c2;
   let b4=b2*b2;
   let a4=a2*a2;
   /* end vars */
   let v1 = a2/((b-c)*(a2-b*c)*(b2+b*c+c2)*(a4+a2*b*c+b2*c2));
   let v2 = b2/((-a+c)*(b2-a*c)*(a2+a*c+c2)*(b4+a*b2*c+a2*c2));
   let v3 = c2/((a-b)*(a2+a*b+b2)*(-(a*b)+c2)*(a2*b2+a*b*c2+c4));
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X796(orbit) {
   /* begin vars */
   let c2=c*c;
   let c4=c2*c2;
   let c6=c2*c4;
   let a2=a*a;
   let a3=a2*a;
   let a4=a2*a2;
   let b2=b*b;
   let b4=b2*b2;
   let b6=b2*b4;
   let c3=c2*c;
   let b3=b2*b;
   let a6=a2*a4;
   let c7=c*c6;
   let a5=a2*a3;
   let b7=b*b6;
   let c5=c2*c3;
   let b5=b2*b3;
   let a7=a*a6;
   /* end vars */
   let v1 = (b-c)*(a7*b2+a7*b*c+a7*c2-b6*c3-b5*c4-b4*c5-b3*c6);
   let v2 = (-a+c)*(a2*b7+a*b7*c+b7*c2-a6*c3-a5*c4-a4*c5-a3*c6);
   let v3 = (a-b)*(-(a6*b3)-a5*b4-a4*b5-a3*b6+a2*c7+a*b*c7+b2*c7);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X797(orbit) {
   /* begin vars */
   let c2=c*c;
   let c4=c2*c2;
   let c6=c2*c4;
   let a2=a*a;
   let a3=a2*a;
   let a4=a2*a2;
   let b2=b*b;
   let b4=b2*b2;
   let b6=b2*b4;
   let c3=c2*c;
   let b3=b2*b;
   let a6=a2*a4;
   let c7=c*c6;
   let a5=a2*a3;
   let b7=b*b6;
   let c5=c2*c3;
   let b5=b2*b3;
   let a7=a*a6;
   /* end vars */
   let v1 = a2/((b-c)*(a7*b2+a7*b*c+a7*c2-b6*c3-b5*c4-b4*c5-b3*c6));
   let v2 = b2/((-a+c)*(a2*b7+a*b7*c+b7*c2-a6*c3-a5*c4-a4*c5-a3*c6));
   let v3 = c2/((a-b)*(-(a6*b3)-a5*b4-a4*b5-a3*b6+a2*c7+a*b*c7+b2*c7));
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X798(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let c3=c2*c;
   let b3=b2*b;
   let a3=a2*a;
   /* end vars */
   let v1 = a3*(b2-c2);
   let v2 = b3*(-a2+c2);
   let v3 = (a2-b2)*c3;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X799(orbit) {
   /* begin vars */
   let a2=a*a;
   let c2=c*c;
   let b2=b*b;
   /* end vars */
   let v1 = 1/(a*(b2-c2));
   let v2 = 1/(b*(-a2+c2));
   let v3 = 1/((a2-b2)*c);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X800(orbit) {
   /* begin vars */
   let a2=a*a;
   let c2=c*c;
   let b2=b*b;
   let SA=(b2+c2-a2)/2;
   let SC=(a2+b2-c2)/2;
   let SB=(c2+a2-b2)/2;
   let SA2=SA*SA;
   let SC2=SC*SC;
   let SB2=SB*SB;
   /* end vars */
   let v1 = a2*(a2*SB*SC+SA*(SB2+SC2));
   let v2 = b2*(b2*SA*SC+SB*(SA2+SC2));
   let v3 = c2*(c2*SA*SB+(SA2+SB2)*SC);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X801(orbit) {
   /* begin vars */
   let a2=a*a;
   let c2=c*c;
   let b2=b*b;
   let SA=(b2+c2-a2)/2;
   let SC=(a2+b2-c2)/2;
   let SB=(c2+a2-b2)/2;
   let SA2=SA*SA;
   let SC2=SC*SC;
   let SB2=SB*SB;
   /* end vars */
   let v1 = 1/(a2*SB*SC+SA*(SB2+SC2));
   let v2 = 1/(b2*SA*SC+SB*(SA2+SC2));
   let v3 = 1/(c2*SA*SB+(SA2+SB2)*SC);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X802(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let c3=c2*c;
   let b3=b2*b;
   let a3=a2*a;
   /* end vars */
   let v1 = (b-c)*(a3*b+a3*c-b2*c2);
   let v2 = (-a+c)*(a*b3+b3*c-a2*c2);
   let v3 = (a-b)*(-(a2*b2)+a*c3+b*c3);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X803(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let c3=c2*c;
   let b3=b2*b;
   let a3=a2*a;
   /* end vars */
   let v1 = a2/((b-c)*(a3*b+a3*c-b2*c2));
   let v2 = b2/((-a+c)*(a*b3+b3*c-a2*c2));
   let v3 = c2/((a-b)*(-(a2*b2)+a*c3+b*c3));
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X804(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let c4=c2*c2;
   let b4=b2*b2;
   let a4=a2*a2;
   /* end vars */
   let v1 = (b2-c2)*(a4-b2*c2);
   let v2 = (-a2+c2)*(b4-a2*c2);
   let v3 = (a2-b2)*(-(a2*b2)+c4);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X805(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let c4=c2*c2;
   let b4=b2*b2;
   let a4=a2*a2;
   /* end vars */
   let v1 = a2/((b2-c2)*(a4-b2*c2));
   let v2 = b2/((-a2+c2)*(b4-a2*c2));
   let v3 = c2/((a2-b2)*(-(a2*b2)+c4));
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X806(orbit) {
   /* begin vars */
   let c2=c*c;
   let c3=c2*c;
   let a2=a*a;
   let b2=b*b;
   let b3=b2*b;
   let a3=a2*a;
   let c5=c2*c3;
   let a4=a2*a2;
   let b5=b2*b3;
   let c4=c2*c2;
   let b4=b2*b2;
   let a5=a2*a3;
   /* end vars */
   let v1 = (b-c)*(a5*b+a5*c-b4*c2-b3*c3-b2*c4);
   let v2 = (-a+c)*(a*b5+b5*c-a4*c2-a3*c3-a2*c4);
   let v3 = (a-b)*(-(a4*b2)-a3*b3-a2*b4+a*c5+b*c5);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X807(orbit) {
   /* begin vars */
   let c2=c*c;
   let c3=c2*c;
   let a2=a*a;
   let b2=b*b;
   let b3=b2*b;
   let a3=a2*a;
   let c5=c2*c3;
   let a4=a2*a2;
   let b5=b2*b3;
   let c4=c2*c2;
   let b4=b2*b2;
   let a5=a2*a3;
   /* end vars */
   let v1 = a2/((b-c)*(a5*b+a5*c-b4*c2-b3*c3-b2*c4));
   let v2 = b2/((-a+c)*(a*b5+b5*c-a4*c2-a3*c3-a2*c4));
   let v3 = c2/((a-b)*(-(a4*b2)-a3*b3-a2*b4+a*c5+b*c5));
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X808(orbit) {
   /* begin vars */
   let c2=c*c;
   let c4=c2*c2;
   let a2=a*a;
   let b2=b*b;
   let b4=b2*b2;
   let a4=a2*a2;
   let c6=c2*c4;
   let b6=b2*b4;
   let a6=a2*a4;
   /* end vars */
   let v1 = (b2-c2)*(a6-b4*c2-b2*c4);
   let v2 = (-a2+c2)*(b6-a4*c2-a2*c4);
   let v3 = (a2-b2)*(-(a4*b2)-a2*b4+c6);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X809(orbit) {
   /* begin vars */
   let c2=c*c;
   let c4=c2*c2;
   let a2=a*a;
   let b2=b*b;
   let b4=b2*b2;
   let a4=a2*a2;
   let c6=c2*c4;
   let b6=b2*b4;
   let a6=a2*a4;
   /* end vars */
   let v1 = a2/((b2-c2)*(a6-b4*c2-b2*c4));
   let v2 = b2/((-a2+c2)*(b6-a4*c2-a2*c4));
   let v3 = c2/((a2-b2)*(-(a4*b2)-a2*b4+c6));
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X810(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let SC=(a2+b2-c2)/2;
   let c3=c2*c;
   let SB=(c2+a2-b2)/2;
   let b3=b2*b;
   let SA=(b2+c2-a2)/2;
   let a3=a2*a;
   /* end vars */
   let v1 = a3*(b2-c2)*SA;
   let v2 = b3*(-a2+c2)*SB;
   let v3 = (a2-b2)*c3*SC;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X811(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let SC=(a2+b2-c2)/2;
   let SB=(c2+a2-b2)/2;
   let SA=(b2+c2-a2)/2;
   /* end vars */
   let v1 = 1/(a*(b2-c2)*SA);
   let v2 = 1/(b*(-a2+c2)*SB);
   let v3 = 1/((a2-b2)*c*SC);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X812(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   /* end vars */
   let v1 = (b-c)*(a2-b*c);
   let v2 = (-a+c)*(b2-a*c);
   let v3 = (a-b)*(-(a*b)+c2);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X813(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   /* end vars */
   let v1 = a2/((b-c)*(a2-b*c));
   let v2 = b2/((-a+c)*(b2-a*c));
   let v3 = c2/((a-b)*(-(a*b)+c2));
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X814(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let c3=c2*c;
   let b3=b2*b;
   let a3=a2*a;
   /* end vars */
   let v1 = (b-c)*(a3-b2*c-b*c2);
   let v2 = (-a+c)*(b3-a2*c-a*c2);
   let v3 = (a-b)*(-(a2*b)-a*b2+c3);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X815(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let c3=c2*c;
   let b3=b2*b;
   let a3=a2*a;
   /* end vars */
   let v1 = a2/((b-c)*(a3-b2*c-b*c2));
   let v2 = b2/((-a+c)*(b3-a2*c-a*c2));
   let v3 = c2/((a-b)*(-(a2*b)-a*b2+c3));
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X816(orbit) {
   /* begin vars */
   let c2=c*c;
   let a2=a*a;
   let b2=b*b;
   let c4=c2*c2;
   let a3=a2*a;
   let b4=b2*b2;
   let c3=c2*c;
   let b3=b2*b;
   let a4=a2*a2;
   /* end vars */
   let v1 = (b-c)*(a4-b3*c-b2*c2-b*c3);
   let v2 = (-a+c)*(b4-a3*c-a2*c2-a*c3);
   let v3 = (a-b)*(-(a3*b)-a2*b2-a*b3+c4);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X817(orbit) {
   /* begin vars */
   let c2=c*c;
   let a2=a*a;
   let b2=b*b;
   let c4=c2*c2;
   let a3=a2*a;
   let b4=b2*b2;
   let c3=c2*c;
   let b3=b2*b;
   let a4=a2*a2;
   /* end vars */
   let v1 = a2/((b-c)*(a4-b3*c-b2*c2-b*c3));
   let v2 = b2/((-a+c)*(b4-a3*c-a2*c2-a*c3));
   let v3 = c2/((a-b)*(-(a3*b)-a2*b2-a*b3+c4));
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X818(orbit) {
   /* begin vars */
   let c2=c*c;
   let c3=c2*c;
   let a2=a*a;
   let b2=b*b;
   let b3=b2*b;
   let a3=a2*a;
   let c5=c2*c3;
   let a4=a2*a2;
   let b5=b2*b3;
   let c4=c2*c2;
   let b4=b2*b2;
   let a5=a2*a3;
   /* end vars */
   let v1 = (b-c)*(a5-b4*c-b3*c2-b2*c3-b*c4);
   let v2 = (-a+c)*(b5-a4*c-a3*c2-a2*c3-a*c4);
   let v3 = (a-b)*(-(a4*b)-a3*b2-a2*b3-a*b4+c5);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X819(orbit) {
   /* begin vars */
   let c2=c*c;
   let c3=c2*c;
   let a2=a*a;
   let b2=b*b;
   let b3=b2*b;
   let a3=a2*a;
   let c5=c2*c3;
   let a4=a2*a2;
   let b5=b2*b3;
   let c4=c2*c2;
   let b4=b2*b2;
   let a5=a2*a3;
   /* end vars */
   let v1 = a2/((b-c)*(a5-b4*c-b3*c2-b2*c3-b*c4));
   let v2 = b2/((-a+c)*(b5-a4*c-a3*c2-a2*c3-a*c4));
   let v3 = c2/((a-b)*(-(a4*b)-a3*b2-a2*b3-a*b4+c5));
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X820(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let SC=(a2+b2-c2)/2;
   let SB=(c2+a2-b2)/2;
   let SA=(b2+c2-a2)/2;
   let c3=c2*c;
   let b3=b2*b;
   let SC2=SC*SC;
   let SB2=SB*SB;
   let SA2=SA*SA;
   let a3=a2*a;
   /* end vars */
   let v1 = a3*SA2*(b2*SB2+c2*SC2);
   let v2 = b3*SB2*(a2*SA2+c2*SC2);
   let v3 = c3*(a2*SA2+b2*SB2)*SC2;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X821(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let SC=(a2+b2-c2)/2;
   let SB=(c2+a2-b2)/2;
   let SA=(b2+c2-a2)/2;
   let SC2=SC*SC;
   let SB2=SB*SB;
   let SA2=SA*SA;
   /* end vars */
   let v1 = 1/(a*SA2*(b2*SB2+c2*SC2));
   let v2 = 1/(b*SB2*(a2*SA2+c2*SC2));
   let v3 = 1/(c*(a2*SA2+b2*SB2)*SC2);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X822(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let SC=(a2+b2-c2)/2;
   let SB=(c2+a2-b2)/2;
   let SA=(b2+c2-a2)/2;
   let SC2=SC*SC;
   let c3=c2*c;
   let SB2=SB*SB;
   let b3=b2*b;
   let SA2=SA*SA;
   let a3=a2*a;
   /* end vars */
   let v1 = a3*(b2-c2)*SA2;
   let v2 = b3*(-a2+c2)*SB2;
   let v3 = (a2-b2)*c3*SC2;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X823(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let SC=(a2+b2-c2)/2;
   let SB=(c2+a2-b2)/2;
   let SA=(b2+c2-a2)/2;
   let SC2=SC*SC;
   let SB2=SB*SB;
   let SA2=SA*SA;
   /* end vars */
   let v1 = 1/(a*(b2-c2)*SA2);
   let v2 = 1/(b*(-a2+c2)*SB2);
   let v3 = 1/((a2-b2)*c*SC2);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X824(orbit) {
   /* begin vars */
   let a2=a*a;
   let c2=c*c;
   let b2=b*b;
   let a3=a2*a;
   let c3=c2*c;
   let b3=b2*b;
   /* end vars */
   let v1 = b3-c3;
   let v2 = -a3+c3;
   let v3 = a3-b3;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X825(orbit) {
   /* begin vars */
   let a2=a*a;
   let c2=c*c;
   let b2=b*b;
   let a3=a2*a;
   let c3=c2*c;
   let b3=b2*b;
   /* end vars */
   let v1 = a2/(b3-c3);
   let v2 = b2/(-a3+c3);
   let v3 = c2/(a3-b3);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X826(orbit) {
   /* begin vars */
   let a2=a*a;
   let c2=c*c;
   let b2=b*b;
   let a4=a2*a2;
   let c4=c2*c2;
   let b4=b2*b2;
   /* end vars */
   let v1 = b4-c4;
   let v2 = -a4+c4;
   let v3 = a4-b4;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X827(orbit) {
   /* begin vars */
   let a2=a*a;
   let c2=c*c;
   let b2=b*b;
   let a4=a2*a2;
   let c4=c2*c2;
   let b4=b2*b2;
   /* end vars */
   let v1 = a2/(b4-c4);
   let v2 = b2/(-a4+c4);
   let v3 = c2/(a4-b4);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X828(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let SC=(a2+b2-c2)/2;
   let SB=(c2+a2-b2)/2;
   let SA=(b2+c2-a2)/2;
   let SC2=SC*SC;
   let c3=c2*c;
   let SB2=SB*SB;
   let b3=b2*b;
   let SA2=SA*SA;
   let a3=a2*a;
   /* end vars */
   let v1 = a3*SA2*(b3*SB2+c3*SC2);
   let v2 = b3*SB2*(a3*SA2+c3*SC2);
   let v3 = c3*(a3*SA2+b3*SB2)*SC2;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X829(orbit) {
   /* begin vars */
   let a2=a*a;
   let c2=c*c;
   let b2=b*b;
   let SC=(a2+b2-c2)/2;
   let SB=(c2+a2-b2)/2;
   let SA=(b2+c2-a2)/2;
   let a3=a2*a;
   let SC2=SC*SC;
   let c3=c2*c;
   let SB2=SB*SB;
   let b3=b2*b;
   let SA2=SA*SA;
   /* end vars */
   let v1 = 1/(a*SA2*(b3*SB2+c3*SC2));
   let v2 = 1/(b*SB2*(a3*SA2+c3*SC2));
   let v3 = 1/(c*(a3*SA2+b3*SB2)*SC2);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X830(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   /* end vars */
   let v1 = a*(b-c)*(a2+b2+b*c+c2);
   let v2 = b*(-a+c)*(a2+b2+a*c+c2);
   let v3 = (a-b)*c*(a2+a*b+b2+c2);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X831(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   /* end vars */
   let v1 = a/((b-c)*(a2+b2+b*c+c2));
   let v2 = b/((-a+c)*(a2+b2+a*c+c2));
   let v3 = c/((a-b)*(a2+a*b+b2+c2));
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X832(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let c4=c2*c2;
   let b4=b2*b2;
   let a4=a2*a2;
   /* end vars */
   let v1 = a4*(b-c)+a*(b4-c4);
   let v2 = b4*(-a+c)+b*(-a4+c4);
   let v3 = (a4-b4)*c+(a-b)*c4;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X833(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let c4=c2*c2;
   let b4=b2*b2;
   let a4=a2*a2;
   /* end vars */
   let v1 = a2/(a4*(b-c)+a*(b4-c4));
   let v2 = b2/(b4*(-a+c)+b*(-a4+c4));
   let v3 = c2/((a4-b4)*c+(a-b)*c4);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X834(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let c3=c2*c;
   let b3=b2*b;
   let a3=a2*a;
   /* end vars */
   let v1 = a3*(b2-c2)+a2*(b3-c3);
   let v2 = b3*(-a2+c2)+b2*(-a3+c3);
   let v3 = (a3-b3)*c2+(a2-b2)*c3;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X835(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let c3=c2*c;
   let b3=b2*b;
   let a3=a2*a;
   /* end vars */
   let v1 = a2/(a3*(b2-c2)+a2*(b3-c3));
   let v2 = b2/(b3*(-a2+c2)+b2*(-a3+c3));
   let v3 = c2/((a3-b3)*c2+(a2-b2)*c3);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X836(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let SC=(a2+b2-c2)/2;
   let SB=(c2+a2-b2)/2;
   let SA=(b2+c2-a2)/2;
   let SC2=SC*SC;
   let SB2=SB*SB;
   let SA2=SA*SA;
   /* end vars */
   let v1 = a2*SA2*(b*SB2+c*SC2);
   let v2 = b2*SB2*(a*SA2+c*SC2);
   let v3 = c2*(a*SA2+b*SB2)*SC2;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X837(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let SC=(a2+b2-c2)/2;
   let SB=(c2+a2-b2)/2;
   let SA=(b2+c2-a2)/2;
   let SC2=SC*SC;
   let SB2=SB*SB;
   let SA2=SA*SA;
   /* end vars */
   let v1 = 1/(SA2*(b*SB2+c*SC2));
   let v2 = 1/(SB2*(a*SA2+c*SC2));
   let v3 = 1/((a*SA2+b*SB2)*SC2);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X838(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let c4=c2*c2;
   let b4=b2*b2;
   let a3=a2*a;
   let c3=c2*c;
   let b3=b2*b;
   let a4=a2*a2;
   /* end vars */
   let v1 = a4*(b3-c3)+a3*(b4-c4);
   let v2 = b4*(-a3+c3)+b3*(-a4+c4);
   let v3 = (a4-b4)*c3+(a3-b3)*c4;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X839(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let c4=c2*c2;
   let b4=b2*b2;
   let a3=a2*a;
   let c3=c2*c;
   let b3=b2*b;
   let a4=a2*a2;
   /* end vars */
   let v1 = a2/(a4*(b3-c3)+a3*(b4-c4));
   let v2 = b2/(b4*(-a3+c3)+b3*(-a4+c4));
   let v3 = c2/((a4-b4)*c3+(a3-b3)*c4);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X840(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let c3=c2*c;
   let b3=b2*b;
   let a3=a2*a;
   /* end vars */
   let v1 = a2/(2*a3-2*a2*b+a*b2-b3-2*a2*c+b2*c+a*c2+b*c2-c3);
   let v2 = b2/(-a3+a2*b-2*a*b2+2*b3+a2*c-2*b2*c+a*c2+b*c2-c3);
   let v3 = c2/(-a3+a2*b+a*b2-b3+a2*c+b2*c-2*a*c2-2*b*c2+2*c3);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X841(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let SA=(b2+c2-a2)/2;
   let SC=(a2+b2-c2)/2;
   let SB=(c2+a2-b2)/2;
   let c4=c2*c2;
   let b4=b2*b2;
   let SA2=SA*SA;
   let SC2=SC*SC;
   let SB2=SB*SB;
   let a4=a2*a2;
   /* end vars */
   let v1 = a2/(-(a4*Math.pow(SA,3))+SA*SB*SC*(SB2+10*SB*SC+SC2)+a2*(5*SA2*(SB-SC)*(SB-SC)-4*SB2*SC2));
   let v2 = b2/(-(b4*Math.pow(SB,3))+SA*SB*SC*(SA2+10*SA*SC+SC2)+b2*(5*SB2*(-SA+SC)*(-SA+SC)-4*SA2*SC2));
   let v3 = c2/(SA*SB*(SA2+10*SA*SB+SB2)*SC-c4*Math.pow(SC,3)+c2*(-4*SA2*SB2+5*(SA-SB)*(SA-SB)*SC2));
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X842(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let SC=(a2+b2-c2)/2;
   let SB=(c2+a2-b2)/2;
   let SA=(b2+c2-a2)/2;
   let SC2=SC*SC;
   let SB2=SB*SB;
   let SA2=SA*SA;
   /* end vars */
   let v1 = a2/(-(a2*(SA2+SB*SC))+2*SA*(SB2+SC2));
   let v2 = b2/(-(b2*(SB2+SA*SC))+2*SB*(SA2+SC2));
   let v3 = c2/(2*(SA2+SB2)*SC-c2*(SA*SB+SC2));
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X843(orbit) {
   /* begin vars */
   let a2=a*a;
   let c2=c*c;
   let b2=b*b;
   let SA=(b2+c2-a2)/2;
   let SC=(a2+b2-c2)/2;
   let SB=(c2+a2-b2)/2;
   let SA2=SA*SA;
   let SC2=SC*SC;
   let SB2=SB*SB;
   /* end vars */
   let v1 = a2/(2*(a2-SA)*SA+SB2-4*SB*SC+SC2);
   let v2 = b2/(SA2+2*(b2-SB)*SB-4*SA*SC+SC2);
   let v3 = c2/(SA2-4*SA*SB+SB2+2*(c2-SC)*SC);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X844(orbit) {
   /* begin vars */
   let sc=(a+b-c)/2;
   let s=(a+b+c)/2;
   let Sqrt=Math.sqrt;
   let sb=(c+a-b)/2;
   let sa=(b+c-a)/2;
   /* end vars */
   let v1 = a*(sa*sb*Sqrt((a*b)/(s*sc))+sa*Sqrt((a*c)/(s*sb))*sc-Sqrt((b*c)/(s*sa))*sb*sc);
   let v2 = b*(sa*sb*Sqrt((a*b)/(s*sc))-sa*Sqrt((a*c)/(s*sb))*sc+Sqrt((b*c)/(s*sa))*sb*sc);
   let v3 = c*(-(sa*sb*Sqrt((a*b)/(s*sc)))+sa*Sqrt((a*c)/(s*sb))*sc+Sqrt((b*c)/(s*sa))*sb*sc);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X845(orbit) {
   /* begin vars */
   let sc=(a+b-c)/2;
   let sb=(c+a-b)/2;
   let sa=(b+c-a)/2;
   let Sqrt=Math.sqrt;
   /* end vars */
   let v1 = a*(Sqrt((Math.pow(sa,3)*Math.pow(sb,3))/(a*b))+Sqrt((Math.pow(sa,3)*Math.pow(sc,3))/(a*c))-Sqrt((Math.pow(sb,3)*Math.pow(sc,3))/(b*c)));
   let v2 = b*(Sqrt((Math.pow(sa,3)*Math.pow(sb,3))/(a*b))-Sqrt((Math.pow(sa,3)*Math.pow(sc,3))/(a*c))+Sqrt((Math.pow(sb,3)*Math.pow(sc,3))/(b*c)));
   let v3 = c*(-Sqrt((Math.pow(sa,3)*Math.pow(sb,3))/(a*b))+Sqrt((Math.pow(sa,3)*Math.pow(sc,3))/(a*c))+Sqrt((Math.pow(sb,3)*Math.pow(sc,3))/(b*c)));
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X846(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let SC=(a2+b2-c2)/2;
   let SB=(c2+a2-b2)/2;
   let SA=(b2+c2-a2)/2;
   /* end vars */
   let v1 = a*(a*b+a*c+b*c+2*SA);
   let v2 = b*(a*b+a*c+b*c+2*SB);
   let v3 = c*(a*b+a*c+b*c+2*SC);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X847(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let SC=(a2+b2-c2)/2;
   let SB=(c2+a2-b2)/2;
   let SA=(b2+c2-a2)/2;
   let SC2=SC*SC;
   let SB2=SB*SB;
   let SA2=SA*SA;
   let area=triAreaHeron(a,b,c);
   /* end vars */
   let v1 = (b2*c2)/(SA*(-4*area*area+SA2));
   let v2 = (a2*c2)/(SB*(-4*area*area+SB2));
   let v3 = (a2*b2)/(SC*(-4*area*area+SC2));
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X848(orbit) {
   /* begin vars */
   let cosC=lawOfCosines(c,a,b);
   let sinC=getSin(cosC);
   let cosB=lawOfCosines(b,a,c);
   let sinB=getSin(cosB);
   let cosA=lawOfCosines(a,b,c);
   let sinA=getSin(cosA);
   let cotCp=getCotPrime(c,a,b);
   let cotC=cosC/sinC;
   let cotBp=getCotPrime(b,c,a);
   let cotB=cosB/sinB;
   let cotAp=getCotPrime(a,b,c);
   let cotA=cosA/sinA;
   /* end vars */
   let v1 = 1/(cotA-cotAp);
   let v2 = 1/(cotB-cotBp);
   let v3 = 1/(cotC-cotCp);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X849(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let c3=c2*c;
   let b3=b2*b;
   let a3=a2*a;
   /* end vars */
   let v1 = a3/(b+c)*(b+c);
   let v2 = b3/(a+c)*(a+c);
   let v3 = c3/(a+b)*(a+b);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X850(orbit) {
   /* begin vars */
   let a2=a*a;
   let c2=c*c;
   let b2=b*b;
   /* end vars */
   let v1 = (b2-c2)/a2;
   let v2 = (-a2+c2)/b2;
   let v3 = (a2-b2)/c2;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X851(orbit) {
   /* begin vars */
   let c2=c*c;
   let a2=a*a;
   let b2=b*b;
   let c4=c2*c2;
   let a3=a2*a;
   let b4=b2*b2;
   let c3=c2*c;
   let b3=b2*b;
   let a4=a2*a2;
   /* end vars */
   let v1 = a*(b+c)*(a4-a2*b2+a2*b*c-b3*c-a2*c2+2*b2*c2-b*c3);
   let v2 = b*(a+c)*(-(a2*b2)+b4-a3*c+a*b2*c+2*a2*c2-b2*c2-a*c3);
   let v3 = (a+b)*c*(-(a3*b)+2*a2*b2-a*b3-a2*c2+a*b*c2-b2*c2+c4);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X852(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let SC=(a2+b2-c2)/2;
   let SB=(c2+a2-b2)/2;
   let SA=(b2+c2-a2)/2;
   let SC2=SC*SC;
   let SB2=SB*SB;
   let SA2=SA*SA;
   /* end vars */
   let v1 = a2*SA*(SA2*SB2+SA2*SC2-2*SB2*SC2);
   let v2 = b2*SB*(SA2*SB2-2*SA2*SC2+SB2*SC2);
   let v3 = c2*SC*(-2*SA2*SB2+SA2*SC2+SB2*SC2);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X853(orbit) {
   /* begin vars */
   let c2=c*c;
   let c3=c2*c;
   let b2=b*b;
   let b3=b2*b;
   let a2=a*a;
   let a3=a2*a;
   let c5=c2*c3;
   let c4=c2*c2;
   let b5=b2*b3;
   let b4=b2*b2;
   let a4=a2*a2;
   let a5=a2*a3;
   /* end vars */
   let v1 = a2*(a-b-c)*(a5*b2-a4*b3-a3*b4+a2*b5+a4*b2*c-a2*b4*c+a5*c2+a4*b*c2-2*b5*c2-a4*c3+2*b4*c3-a3*c4-a2*b*c4+2*b3*c4+a2*c5-2*b2*c5);
   let v2 = b2*(-a+b-c)*(a5*b2-a4*b3-a3*b4+a2*b5-a4*b2*c+a2*b4*c-2*a5*c2+a*b4*c2+b5*c2+2*a4*c3-b4*c3+2*a3*c4-a*b2*c4-b3*c4-2*a2*c5+b2*c5);
   let v3 = (-a-b+c)*c2*(-2*a5*b2+2*a4*b3+2*a3*b4-2*a2*b5+a5*c2-a4*b*c2-a*b4*c2+b5*c2-a4*c3-b4*c3-a3*c4+a2*b*c4+a*b2*c4-b3*c4+a2*c5+b2*c5);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X854(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let SC=(a2+b2-c2)/2;
   let SB=(c2+a2-b2)/2;
   let SA=(b2+c2-a2)/2;
   /* end vars */
   let v1 = -((a2*c2*(-(a2*SA)+b2*SB))/((a+b-c)*(-a+b+c)))+(a2*b2*(a2*SA-c2*SC))/((a-b+c)*(-a+b+c));
   let v2 = (b2*c2*(-(a2*SA)+b2*SB))/((a+b-c)*(a-b+c))-(a2*b2*(-(b2*SB)+c2*SC))/((a-b+c)*(-a+b+c));
   let v3 = -((b2*c2*(a2*SA-c2*SC))/((a+b-c)*(a-b+c)))+(a2*c2*(-(b2*SB)+c2*SC))/((a+b-c)*(-a+b+c));
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X855(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let SC=(a2+b2-c2)/2;
   let SB=(c2+a2-b2)/2;
   let SA=(b2+c2-a2)/2;
   /* end vars */
   let v1 = -((a*c*(-(a2*SA)+b2*SB))/((a+b-c)*(-a+b+c)))+(a*b*(a2*SA-c2*SC))/((a-b+c)*(-a+b+c));
   let v2 = (b*c*(-(a2*SA)+b2*SB))/((a+b-c)*(a-b+c))-(a*b*(-(b2*SB)+c2*SC))/((a-b+c)*(-a+b+c));
   let v3 = -((b*c*(a2*SA-c2*SC))/((a+b-c)*(a-b+c)))+(a*c*(-(b2*SB)+c2*SC))/((a+b-c)*(-a+b+c));
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X856(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let SC=(a2+b2-c2)/2;
   let SB=(c2+a2-b2)/2;
   let SA=(b2+c2-a2)/2;
   /* end vars */
   let v1 = a*SA*(-(b*c*(b+c)*SB*SC)+a2*SA*(b*SB+c*SC));
   let v2 = b*SB*(-(a*c*(a+c)*SA*SC)+b2*SB*(a*SA+c*SC));
   let v3 = c*SC*(-(a*b*(a+b)*SA*SB)+c2*(a*SA+b*SB)*SC);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X857(orbit) {
   /* begin vars */
   let a2=a*a;
   let c2=c*c;
   let b2=b*b;
   let a3=a2*a;
   let SC=(a2+b2-c2)/2;
   let c3=c2*c;
   let SB=(c2+a2-b2)/2;
   let b3=b2*b;
   let SA=(b2+c2-a2)/2;
   /* end vars */
   let v1 = a2*(b+c)*SA-b3*SB-c3*SC;
   let v2 = -(a3*SA)+b2*(a+c)*SB-c3*SC;
   let v3 = -(a3*SA)-b3*SB+(a+b)*c2*SC;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X858(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let c4=c2*c2;
   let b4=b2*b2;
   let a4=a2*a2;
   /* end vars */
   let v1 = c2*(a4-a2*b2+b4-c4)+b2*(a4-b4-a2*c2+c4);
   let v2 = c2*(a4-a2*b2+b4-c4)+a2*(-a4+b4-b2*c2+c4);
   let v3 = b2*(a4-b4-a2*c2+c4)+a2*(-a4+b4-b2*c2+c4);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X859(orbit) {
   /* begin vars */
   let a2=a*a;
   let c2=c*c;
   let b2=b*b;
   let SA=(b2+c2-a2)/2;
   let SC=(a2+b2-c2)/2;
   let SB=(c2+a2-b2)/2;
   /* end vars */
   let v1 = a2*(a+b)*(a+c)*(a*b*c-b*SB-c*SC);
   let v2 = (a+b)*b2*(b+c)*(a*b*c-a*SA-c*SC);
   let v3 = (a+c)*(b+c)*c2*(a*b*c-a*SA-b*SB);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X860(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let SC=(a2+b2-c2)/2;
   let SB=(c2+a2-b2)/2;
   let SA=(b2+c2-a2)/2;
   /* end vars */
   let v1 = (b+c)*(b*c-2*SA)*SB*SC;
   let v2 = (a+c)*SA*(a*c-2*SB)*SC;
   let v3 = (a+b)*SA*SB*(a*b-2*SC);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X861(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let SC=(a2+b2-c2)/2;
   let SB=(c2+a2-b2)/2;
   let SA=(b2+c2-a2)/2;
   /* end vars */
   let v1 = a*(-a+b+c)*((a+b-c)*c*(-(a2*SA)+b2*SB)-b*(a-b+c)*(a2*SA-c2*SC));
   let v2 = b*(a-b+c)*(-((a+b-c)*c*(-(a2*SA)+b2*SB))+a*(-a+b+c)*(-(b2*SB)+c2*SC));
   let v3 = (a+b-c)*c*(b*(a-b+c)*(a2*SA-c2*SC)-a*(-a+b+c)*(-(b2*SB)+c2*SC));
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X862(orbit) {
   /* begin vars */
   let a2=a*a;
   let c2=c*c;
   let b2=b*b;
   let SA=(b2+c2-a2)/2;
   let SC=(a2+b2-c2)/2;
   let SB=(c2+a2-b2)/2;
   /* end vars */
   let v1 = a*(b+c)*(-a2+b*c)*SB*SC;
   let v2 = b*(a+c)*(-b2+a*c)*SA*SC;
   let v3 = (a+b)*c*(a*b-c2)*SA*SB;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X863(orbit) {
   /* begin vars */
   let a2=a*a;
   let c2=c*c;
   let b2=b*b;
   let a3=a2*a;
   let SC=(a2+b2-c2)/2;
   let SB=(c2+a2-b2)/2;
   let c3=c2*c;
   let b3=b2*b;
   let SA=(b2+c2-a2)/2;
   /* end vars */
   let v1 = a2*(-(SA/b3)-SA/c3)+SB/b+SC/c;
   let v2 = SA/a+b2*(-(SB/a3)-SB/c3)+SC/c;
   let v3 = SA/a+SB/b+c2*(-(SC/a3)-SC/b3);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X864(orbit) {
   /* begin vars */
   let a2=a*a;
   let c2=c*c;
   let b2=b*b;
   let a4=a2*a2;
   let SC=(a2+b2-c2)/2;
   let SB=(c2+a2-b2)/2;
   let c4=c2*c2;
   let b4=b2*b2;
   let SA=(b2+c2-a2)/2;
   /* end vars */
   let v1 = a2*(-(SA/b4)-SA/c4)+SB/b2+SC/c2;
   let v2 = SA/a2+b2*(-(SB/a4)-SB/c4)+SC/c2;
   let v3 = SA/a2+SB/b2+c2*(-(SC/a4)-SC/b4);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X865(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let SC=(a2+b2-c2)/2;
   let SB=(c2+a2-b2)/2;
   let SA=(b2+c2-a2)/2;
   let SC2=SC*SC;
   let SB2=SB*SB;
   let SA2=SA*SA;
   /* end vars */
   let v1 = a2*(SB-SC)*(SB-SC)*(a2*SB*SC-SA*(SA2+SB*SC));
   let v2 = b2*(-SA+SC)*(-SA+SC)*(b2*SA*SC-SB*(SB2+SA*SC));
   let v3 = c2*(SA-SB)*(SA-SB)*(c2*SA*SB-SC*(SA*SB+SC2));
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X866(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let SC=(a2+b2-c2)/2;
   let c3=c2*c;
   let SB=(c2+a2-b2)/2;
   let b3=b2*b;
   let SA=(b2+c2-a2)/2;
   let a3=a2*a;
   /* end vars */
   let v1 = a*(b-c)*(a3*(b-c)*SA-b3*c*SB+b*c3*SC+a*b*c*(b*SB-c*SC));
   let v2 = b*(-a+c)*(a3*c*SA+b3*(-a+c)*SB-a*c3*SC+a*b*c*(-(a*SA)+c*SC));
   let v3 = (a-b)*c*(-(a3*b*SA)+a*b3*SB+a*b*c*(a*SA-b*SB)+(a-b)*c3*SC);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X867(orbit) {
   /* begin vars */
   let a2=a*a;
   let c2=c*c;
   let b2=b*b;
   let a3=a2*a;
   let SC=(a2+b2-c2)/2;
   let c3=c2*c;
   let SB=(c2+a2-b2)/2;
   let b3=b2*b;
   let SA=(b2+c2-a2)/2;
   /* end vars */
   let v1 = (b-c)*(a2*(-b+c)*SA+b3*SB-c3*SC+a*(-(b2*SB)+c2*SC));
   let v2 = (-a+c)*(-(a3*SA)+b2*(a-c)*SB+c3*SC+b*(a2*SA-c2*SC));
   let v3 = (a-b)*(a3*SA-b3*SB+c*(-(a2*SA)+b2*SB)+(-a+b)*c2*SC);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X868(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let SC=(a2+b2-c2)/2;
   let SB=(c2+a2-b2)/2;
   let SA=(b2+c2-a2)/2;
   let SC2=SC*SC;
   let SB2=SB*SB;
   let SA2=SA*SA;
   /* end vars */
   let v1 = (SB-SC)*(SB-SC)*(-SA2+SB*SC);
   let v2 = (-SA+SC)*(-SA+SC)*(-SB2+SA*SC);
   let v3 = (SA-SB)*(SA-SB)*(SA*SB-SC2);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X869(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let c3=c2*c;
   let b3=b2*b;
   let a3=a2*a;
   /* end vars */
   let v1 = a3*(b2+b*c+c2);
   let v2 = b3*(a2+a*c+c2);
   let v3 = (a2+a*b+b2)*c3;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X870(orbit) {
   /* begin vars */
   let a2=a*a;
   let c2=c*c;
   let b2=b*b;
   /* end vars */
   let v1 = 1/(a*(b2+b*c+c2));
   let v2 = 1/(b*(a2+a*c+c2));
   let v3 = 1/((a2+a*b+b2)*c);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X871(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let c3=c2*c;
   let b3=b2*b;
   let a3=a2*a;
   /* end vars */
   let v1 = 1/(a3*(b2+b*c+c2));
   let v2 = 1/(b3*(a2+a*c+c2));
   let v3 = 1/((a2+a*b+b2)*c3);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X872(orbit) {
   /* begin vars */

   /* end vars */
   let v1 = Math.pow(a,3)*(b+c)*(b+c);
   let v2 = Math.pow(b,3)*(a+c)*(a+c);
   let v3 = (a+b)*(a+b)*Math.pow(c,3);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X873(orbit) {
   /* begin vars */

   /* end vars */
   let v1 = 1/(a*(b+c)*(b+c));
   let v2 = 1/(b*(a+c)*(a+c));
   let v3 = 1/((a+b)*(a+b)*c);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X874(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   /* end vars */
   let v1 = (a2-b*c)/(a*b-a*c);
   let v2 = (b2-a*c)/(-(a*b)+b*c);
   let v3 = (-(a*b)+c2)/(a*c-b*c);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X875(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   /* end vars */
   let v1 = (a2*(a*b-a*c))/(a2-b*c);
   let v2 = (b2*(-(a*b)+b*c))/(b2-a*c);
   let v3 = ((a*c-b*c)*c2)/(-(a*b)+c2);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X876(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   /* end vars */
   let v1 = (a*b-a*c)/(a2-b*c);
   let v2 = (-(a*b)+b*c)/(b2-a*c);
   let v3 = (a*c-b*c)/(-(a*b)+c2);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X877(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let SC=(a2+b2-c2)/2;
   let SB=(c2+a2-b2)/2;
   let SA=(b2+c2-a2)/2;
   let SC2=SC*SC;
   let SB2=SB*SB;
   let SA2=SA*SA;
   /* end vars */
   let v1 = (SA-SB)*SB*(SA-SC)*SC*(SA2-SB*SC);
   let v2 = SA*(-SA+SB)*(SB-SC)*SC*(SB2-SA*SC);
   let v3 = SA*SB*(-SA+SC)*(-SB+SC)*(-(SA*SB)+SC2);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X878(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let SC=(a2+b2-c2)/2;
   let SB=(c2+a2-b2)/2;
   let SA=(b2+c2-a2)/2;
   let SC2=SC*SC;
   let SB2=SB*SB;
   let SA2=SA*SA;
   /* end vars */
   let v1 = (a2*SA*(SB-SC))/(SA2-SB*SC);
   let v2 = (b2*SB*(-SA+SC))/(SB2-SA*SC);
   let v3 = (c2*(SA-SB)*SC)/(-(SA*SB)+SC2);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X879(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let SC=(a2+b2-c2)/2;
   let SB=(c2+a2-b2)/2;
   let SA=(b2+c2-a2)/2;
   let SC2=SC*SC;
   let SB2=SB*SB;
   let SA2=SA*SA;
   /* end vars */
   let v1 = (SA*(SB-SC))/(SA2-SB*SC);
   let v2 = (SB*(-SA+SC))/(SB2-SA*SC);
   let v3 = ((SA-SB)*SC)/(-(SA*SB)+SC2);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X880(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let c4=c2*c2;
   let b4=b2*b2;
   let a4=a2*a2;
   /* end vars */
   let v1 = (a4-b2*c2)/(a2*(b2-c2));
   let v2 = (b4-a2*c2)/(b2*(-a2+c2));
   let v3 = (-(a2*b2)+c4)/((a2-b2)*c2);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X881(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let c4=c2*c2;
   let b4=b2*b2;
   let a4=a2*a2;
   /* end vars */
   let v1 = (a4*(b2-c2))/(a4-b2*c2);
   let v2 = (b4*(-a2+c2))/(b4-a2*c2);
   let v3 = ((a2-b2)*c4)/(-(a2*b2)+c4);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X882(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let c4=c2*c2;
   let b4=b2*b2;
   let a4=a2*a2;
   /* end vars */
   let v1 = (a2*(b2-c2))/(a4-b2*c2);
   let v2 = (b2*(-a2+c2))/(b4-a2*c2);
   let v3 = ((a2-b2)*c2)/(-(a2*b2)+c4);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X883(orbit) {
   /* begin vars */
   let a2=a*a;
   let c2=c*c;
   let b2=b*b;
   /* end vars */
   let v1 = (b2-a*(b+c)+c2)/((a-b-c)*(b-c));
   let v2 = (a2-b*(a+c)+c2)/((-a+b-c)*(-a+c));
   let v3 = (a2+b2-(a+b)*c)/((a-b)*(-a-b+c));
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X884(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   /* end vars */
   let v1 = (a2*(a-b-c)*(b-c))/(b2-a*(b+c)+c2);
   let v2 = (b2*(-a+b-c)*(-a+c))/(a2-b*(a+c)+c2);
   let v3 = ((a-b)*(-a-b+c)*c2)/(a2+b2-(a+b)*c);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X885(orbit) {
   /* begin vars */
   let a2=a*a;
   let c2=c*c;
   let b2=b*b;
   /* end vars */
   let v1 = ((a-b-c)*(b-c))/(b2-a*(b+c)+c2);
   let v2 = ((-a+b-c)*(-a+c))/(a2-b*(a+c)+c2);
   let v3 = ((a-b)*(-a-b+c))/(a2+b2-(a+b)*c);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X886(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   /* end vars */
   let v1 = 1/(a2*(b2-c2)*(-2*b2*c2+a2*(b2+c2)));
   let v2 = 1/(b2*(-a2+c2)*(-2*a2*c2+b2*(a2+c2)));
   let v3 = 1/((a2-b2)*c2*(-2*a2*b2+(a2+b2)*c2));
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X887(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let c4=c2*c2;
   let b4=b2*b2;
   let a4=a2*a2;
   /* end vars */
   let v1 = a4*(b2-c2)*(-2*b2*c2+a2*(b2+c2));
   let v2 = b4*(-a2+c2)*(-2*a2*c2+b2*(a2+c2));
   let v3 = (a2-b2)*(-2*a2*b2+(a2+b2)*c2)*c4;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X888(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   /* end vars */
   let v1 = a2*(b2-c2)*(-2*b2*c2+a2*(b2+c2));
   let v2 = b2*(-a2+c2)*(-2*a2*c2+b2*(a2+c2));
   let v3 = (a2-b2)*c2*(-2*a2*b2+(a2+b2)*c2);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X889(orbit) {
   /* begin vars */

   /* end vars */
   let v1 = 1/(a*(b-c)*(-2*b*c+a*(b+c)));
   let v2 = 1/(b*(-a+c)*(-2*a*c+b*(a+c)));
   let v3 = 1/((a-b)*c*(-2*a*b+(a+b)*c));
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X890(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let c3=c2*c;
   let b3=b2*b;
   let a3=a2*a;
   /* end vars */
   let v1 = a3*(b-c)*(-2*b*c+a*(b+c));
   let v2 = b3*(-a+c)*(-2*a*c+b*(a+c));
   let v3 = (a-b)*(-2*a*b+(a+b)*c)*c3;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X891(orbit) {
   /* begin vars */

   /* end vars */
   let v1 = a*(b-c)*(-2*b*c+a*(b+c));
   let v2 = b*(-a+c)*(-2*a*c+b*(a+c));
   let v3 = (a-b)*c*(-2*a*b+(a+b)*c);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X892(orbit) {
   /* begin vars */
   let a2=a*a;
   let c2=c*c;
   let b2=b*b;
   /* end vars */
   let v1 = 1/((b2-c2)*(-2*a2+b2+c2));
   let v2 = 1/((-a2+c2)*(a2-2*b2+c2));
   let v3 = 1/((a2-b2)*(a2+b2-2*c2));
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X893(orbit) {
   /* begin vars */
   let c2=c*c;
   let a2=a*a;
   let b2=b*b;
   /* end vars */
   let v1 = 1/(b2*(a2+b*c)*c2);
   let v2 = 1/(a2*(b2+a*c)*c2);
   let v3 = 1/(a2*b2*(a*b+c2));
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X894(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   /* end vars */
   let v1 = a2+b*c;
   let v2 = b2+a*c;
   let v3 = a*b+c2;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X895(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let SC=(a2+b2-c2)/2;
   let SB=(c2+a2-b2)/2;
   let SA=(b2+c2-a2)/2;
   /* end vars */
   let v1 = a2/((a2-2*SA)*SB*SC);
   let v2 = b2/(SA*(b2-2*SB)*SC);
   let v3 = c2/(SA*SB*(c2-2*SC));
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X896(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   /* end vars */
   let v1 = a*(2*a2-b2-c2);
   let v2 = b*(-a2+2*b2-c2);
   let v3 = c*(-a2-b2+2*c2);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X897(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   /* end vars */
   let v1 = a/(2*a2-b2-c2);
   let v2 = b/(-a2+2*b2-c2);
   let v3 = c/(-a2-b2+2*c2);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X898(orbit) {
   /* begin vars */

   /* end vars */
   let v1 = a/((b-c)*(-2*b*c+a*(b+c)));
   let v2 = b/((-a+c)*(-2*a*c+b*(a+c)));
   let v3 = c/((a-b)*(-2*a*b+(a+b)*c));
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X899(orbit) {
   /* begin vars */

   /* end vars */
   let v1 = a*(-2/a+1/b+1/c);
   let v2 = b*(1/a-2/b+1/c);
   let v3 = (1/a+1/b-2/c)*c;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X900(orbit) {
   /* begin vars */

   /* end vars */
   let v1 = (b-c)*(-2*a+b+c);
   let v2 = (-a+c)*(a-2*b+c);
   let v3 = (a-b)*(a+b-2*c);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X901(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   /* end vars */
   let v1 = a2/((b-c)*(-2*a+b+c));
   let v2 = b2/((-a+c)*(a-2*b+c));
   let v3 = c2/((a-b)*(a+b-2*c));
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X902(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   /* end vars */
   let v1 = a2*(-2*a+b+c);
   let v2 = b2*(a-2*b+c);
   let v3 = (a+b-2*c)*c2;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X903(orbit) {
   /* begin vars */

   /* end vars */
   let v1 = 1/(-2*a+b+c);
   let v2 = 1/(a-2*b+c);
   let v3 = 1/(a+b-2*c);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X904(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let c3=c2*c;
   let b3=b2*b;
   let a3=a2*a;
   /* end vars */
   let v1 = a3/(a2+b*c);
   let v2 = b3/(b2+a*c);
   let v3 = c3/(a*b+c2);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X905(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let SC=(a2+b2-c2)/2;
   let SB=(c2+a2-b2)/2;
   let SA=(b2+c2-a2)/2;
   /* end vars */
   let v1 = a*(b-c)*SA;
   let v2 = b*(-a+c)*SB;
   let v3 = (a-b)*c*SC;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X906(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let SC=(a2+b2-c2)/2;
   let c3=c2*c;
   let SB=(c2+a2-b2)/2;
   let b3=b2*b;
   let SA=(b2+c2-a2)/2;
   let a3=a2*a;
   /* end vars */
   let v1 = (2*a3*SA)/(b-c);
   let v2 = (2*b3*SB)/(-a+c);
   let v3 = (2*c3*SC)/(a-b);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X907(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   /* end vars */
   let v1 = a2/((b2-c2)*(3*a2+b2+c2));
   let v2 = b2/((-a2+c2)*(a2+3*b2+c2));
   let v3 = c2/((a2-b2)*(a2+b2+3*c2));
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X908(orbit) {
   /* begin vars */
   let a2=a*a;
   let c2=c*c;
   let b2=b*b;
   let SA=(b2+c2-a2)/2;
   let SC=(a2+b2-c2)/2;
   let SB=(c2+a2-b2)/2;
   /* end vars */
   let v1 = a*b*c-b*SB-c*SC;
   let v2 = a*b*c-a*SA-c*SC;
   let v3 = a*b*c-a*SA-b*SB;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X909(orbit) {
   /* begin vars */
   let a2=a*a;
   let c2=c*c;
   let b2=b*b;
   let SA=(b2+c2-a2)/2;
   let SC=(a2+b2-c2)/2;
   let SB=(c2+a2-b2)/2;
   /* end vars */
   let v1 = a2/(a*b*c-b*SB-c*SC);
   let v2 = b2/(a*b*c-a*SA-c*SC);
   let v3 = c2/(a*b*c-a*SA-b*SB);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X910(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let c3=c2*c;
   let b3=b2*b;
   let a3=a2*a;
   /* end vars */
   let v1 = a*(2*a3-a2*b-b3-a2*c+b2*c+b*c2-c3);
   let v2 = b*(-a3-a*b2+2*b3+a2*c-b2*c+a*c2-c3);
   let v3 = c*(-a3+a2*b+a*b2-b3-a*c2-b*c2+2*c3);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X911(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let c3=c2*c;
   let b3=b2*b;
   let a3=a2*a;
   /* end vars */
   let v1 = a3/(2*a3-a2*(b+c)-(b-c)*(b-c)*(b+c));
   let v2 = b3/(2*b3-b2*(a+c)-(-a+c)*(-a+c)*(a+c));
   let v3 = c3/(-((a-b)*(a-b)*(a+b))-(a+b)*c2+2*c3);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X912(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let SC=(a2+b2-c2)/2;
   let SB=(c2+a2-b2)/2;
   let SA=(b2+c2-a2)/2;
   /* end vars */
   let v1 = a*SA*(a*((a+b-c)*(a-b+c)*(b+c)-2*a*SA)-4*SB*SC);
   let v2 = b*SB*(b*((a+b-c)*(a+c)*(-a+b+c)-2*b*SB)-4*SA*SC);
   let v3 = c*SC*(-4*SA*SB+c*((a+b)*(a-b+c)*(-a+b+c)-2*c*SC));
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X913(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let SC=(a2+b2-c2)/2;
   let SB=(c2+a2-b2)/2;
   let SA=(b2+c2-a2)/2;
   /* end vars */
   let v1 = a2/(2*a*SA*((a+b-c)*(a-b+c)*(b+c)-2*a*SA)-8*SA*SB*SC);
   let v2 = b2/(2*b*SB*((a+b-c)*(a+c)*(-a+b+c)-2*b*SB)-8*SA*SB*SC);
   let v3 = c2/(-8*SA*SB*SC+2*c*SC*((a+b)*(a-b+c)*(-a+b+c)-2*c*SC));
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X914(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let SC=(a2+b2-c2)/2;
   let SB=(c2+a2-b2)/2;
   let SA=(b2+c2-a2)/2;
   /* end vars */
   let v1 = 2*a*SA*((a+b-c)*(a-b+c)*(b+c)-2*a*SA)-8*SA*SB*SC;
   let v2 = 2*b*SB*((a+b-c)*(a+c)*(-a+b+c)-2*b*SB)-8*SA*SB*SC;
   let v3 = -8*SA*SB*SC+2*c*SC*((a+b)*(a-b+c)*(-a+b+c)-2*c*SC);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X915(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let SC=(a2+b2-c2)/2;
   let SB=(c2+a2-b2)/2;
   let SA=(b2+c2-a2)/2;
   /* end vars */
   let v1 = a/(SA*(a*((a+b-c)*(a-b+c)*(b+c)-2*a*SA)-4*SB*SC));
   let v2 = b/(SB*(b*((a+b-c)*(a+c)*(-a+b+c)-2*b*SB)-4*SA*SC));
   let v3 = c/(SC*(-4*SA*SB+c*((a+b)*(a-b+c)*(-a+b+c)-2*c*SC)));
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X916(orbit) {
   /* begin vars */
   let c2=c*c;
   let c3=c2*c;
   let b2=b*b;
   let b3=b2*b;
   let a2=a*a;
   let a3=a2*a;
   let c5=c2*c3;
   let b5=b2*b3;
   let SC=(a2+b2-c2)/2;
   let SB=(c2+a2-b2)/2;
   let a5=a2*a3;
   let SA=(b2+c2-a2)/2;
   /* end vars */
   let v1 = a2*SA*(a5+2*a3*SA-a*(SB-SC)*(SB-SC)-2*(b3*SB+c3*SC));
   let v2 = b2*SB*(b5+2*b3*SB-b*(-SA+SC)*(-SA+SC)-2*(a3*SA+c3*SC));
   let v3 = c2*SC*(c5-c*(SA-SB)*(SA-SB)-2*(a3*SA+b3*SB)+2*c3*SC);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X917(orbit) {
   /* begin vars */
   let c2=c*c;
   let c3=c2*c;
   let b2=b*b;
   let b3=b2*b;
   let a2=a*a;
   let a3=a2*a;
   let c5=c2*c3;
   let b5=b2*b3;
   let SC=(a2+b2-c2)/2;
   let SB=(c2+a2-b2)/2;
   let a5=a2*a3;
   let SA=(b2+c2-a2)/2;
   /* end vars */
   let v1 = 1/(SA*(a5+2*a3*SA-a*(SB-SC)*(SB-SC)-2*(b3*SB+c3*SC)));
   let v2 = 1/(SB*(b5+2*b3*SB-b*(-SA+SC)*(-SA+SC)-2*(a3*SA+c3*SC)));
   let v3 = 1/(SC*(c5-c*(SA-SB)*(SA-SB)-2*(a3*SA+b3*SB)+2*c3*SC));
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X918(orbit) {
   /* begin vars */
   let a2=a*a;
   let c2=c*c;
   let b2=b*b;
   /* end vars */
   let v1 = (b-c)*(b2-a*(b+c)+c2);
   let v2 = (-a+c)*(a2-b*(a+c)+c2);
   let v3 = (a-b)*(a2+b2-(a+b)*c);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X919(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   /* end vars */
   let v1 = a2/((b-c)*(b2-a*(b+c)+c2));
   let v2 = b2/((-a+c)*(a2-b*(a+c)+c2));
   let v3 = c2/((a-b)*(a2+b2-(a+b)*c));
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X920(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let SC=(a2+b2-c2)/2;
   let SB=(c2+a2-b2)/2;
   let SA=(b2+c2-a2)/2;
   let SC2=SC*SC;
   let SB2=SB*SB;
   let SA2=SA*SA;
   /* end vars */
   let v1 = a*(-(a2*SA2)+b2*SB2+c2*SC2);
   let v2 = b*(a2*SA2-b2*SB2+c2*SC2);
   let v3 = c*(a2*SA2+b2*SB2-c2*SC2);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X921(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let SC=(a2+b2-c2)/2;
   let SB=(c2+a2-b2)/2;
   let SA=(b2+c2-a2)/2;
   let SC2=SC*SC;
   let SB2=SB*SB;
   let SA2=SA*SA;
   /* end vars */
   let v1 = a/(-(a2*SA2)+b2*SB2+c2*SC2);
   let v2 = b/(a2*SA2-b2*SB2+c2*SC2);
   let v3 = c/(a2*SA2+b2*SB2-c2*SC2);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X922(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let c3=c2*c;
   let b3=b2*b;
   let a3=a2*a;
   /* end vars */
   let v1 = a3*(2*a2-b2-c2);
   let v2 = b3*(-a2+2*b2-c2);
   let v3 = (-a2-b2+2*c2)*c3;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X923(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let c3=c2*c;
   let b3=b2*b;
   let a3=a2*a;
   /* end vars */
   let v1 = a3/(2*a2-b2-c2);
   let v2 = b3/(-a2+2*b2-c2);
   let v3 = c3/(-a2-b2+2*c2);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X924(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let SC=(a2+b2-c2)/2;
   let SB=(c2+a2-b2)/2;
   let SA=(b2+c2-a2)/2;
   let SC2=SC*SC;
   let SB2=SB*SB;
   let SA2=SA*SA;
   let area=triAreaHeron(a,b,c);
   /* end vars */
   let v1 = a2*(-4*area*area+SA2)*(SB-SC);
   let v2 = b2*(-4*area*area+SB2)*(-SA+SC);
   let v3 = c2*(SA-SB)*(-4*area*area+SC2);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X925(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let SC=(a2+b2-c2)/2;
   let SB=(c2+a2-b2)/2;
   let SA=(b2+c2-a2)/2;
   let SC2=SC*SC;
   let SB2=SB*SB;
   let SA2=SA*SA;
   let area=triAreaHeron(a,b,c);
   /* end vars */
   let v1 = 1/((-4*area*area+SA2)*(SB-SC));
   let v2 = 1/((-4*area*area+SB2)*(-SA+SC));
   let v3 = 1/((SA-SB)*(-4*area*area+SC2));
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X926(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   /* end vars */
   let v1 = a2*(a-b-c)*(b-c)*(-b2+a*(b+c)-c2);
   let v2 = b2*(-a+b-c)*(-a+c)*(-a2+b*(a+c)-c2);
   let v3 = (a-b)*(-a-b+c)*(-a2-b2+(a+b)*c)*c2;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X927(orbit) {
   /* begin vars */
   let a2=a*a;
   let c2=c*c;
   let b2=b*b;
   /* end vars */
   let v1 = 1/((a-b-c)*(b-c)*(-b2+a*(b+c)-c2));
   let v2 = 1/((-a+b-c)*(-a+c)*(-a2+b*(a+c)-c2));
   let v3 = 1/((a-b)*(-a-b+c)*(-a2-b2+(a+b)*c));
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X928(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let SC=(a2+b2-c2)/2;
   let SB=(c2+a2-b2)/2;
   let SA=(b2+c2-a2)/2;
   let SC2=SC*SC;
   let c3=c2*c;
   let SB2=SB*SB;
   let b3=b2*b;
   let SA2=SA*SA;
   let a3=a2*a;
   /* end vars */
   let v1 = a2*(b-c)*(-2*a2*b*c+a3*(b+c)-a*(b-c)*(b-c)*(b+c)+2*(SA2-SB*SC));
   let v2 = b2*(-a+c)*(-2*a*b2*c+b3*(a+c)-b*(-a+c)*(-a+c)*(a+c)+2*(SB2-SA*SC));
   let v3 = (a-b)*c2*(-((a-b)*(a-b)*(a+b)*c)-2*a*b*c2+(a+b)*c3+2*(-(SA*SB)+SC2));
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X929(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let SC=(a2+b2-c2)/2;
   let SB=(c2+a2-b2)/2;
   let SA=(b2+c2-a2)/2;
   let SC2=SC*SC;
   let c3=c2*c;
   let SB2=SB*SB;
   let b3=b2*b;
   let SA2=SA*SA;
   let a3=a2*a;
   /* end vars */
   let v1 = 1/((b-c)*(-2*a2*b*c+a3*(b+c)-a*(b-c)*(b-c)*(b+c)+2*(SA2-SB*SC)));
   let v2 = 1/((-a+c)*(-2*a*b2*c+b3*(a+c)-b*(-a+c)*(-a+c)*(a+c)+2*(SB2-SA*SC)));
   let v3 = 1/((a-b)*(-((a-b)*(a-b)*(a+b)*c)-2*a*b*c2+(a+b)*c3+2*(-(SA*SB)+SC2)));
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X930(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let SC=(a2+b2-c2)/2;
   let SB=(c2+a2-b2)/2;
   let SA=(b2+c2-a2)/2;
   let SC2=SC*SC;
   let SB2=SB*SB;
   let SA2=SA*SA;
   let area=triAreaHeron(a,b,c);
   /* end vars */
   let v1 = 1/((b2-c2)*(12*area*area-SA2));
   let v2 = 1/((-a2+c2)*(12*area*area-SB2));
   let v3 = 1/((a2-b2)*(12*area*area-SC2));
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X931(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   /* end vars */
   let v1 = a/((a2+a*b+a*c+2*b*c)*(b2-c2));
   let v2 = b/((a*b+b2+2*a*c+b*c)*(-a2+c2));
   let v3 = c/((a2-b2)*(2*a*b+a*c+b*c+c2));
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X932(orbit) {
   /* begin vars */

   /* end vars */
   let v1 = (a*(a-b)*(a-c))/(a*b+a*c-b*c);
   let v2 = (b*(-a+b)*(b-c))/(a*b-a*c+b*c);
   let v3 = (c*(-a+c)*(-b+c))/(-(a*b)+a*c+b*c);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X933(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let SC=(a2+b2-c2)/2;
   let SB=(c2+a2-b2)/2;
   let SA=(b2+c2-a2)/2;
   /* end vars */
   let v1 = a2/(SA*(SB-SC)*(a2*SA+2*SB*SC));
   let v2 = b2/(SB*(-SA+SC)*(b2*SB+2*SA*SC));
   let v3 = c2/((SA-SB)*SC*(2*SA*SB+c2*SC));
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X934(orbit) {
   /* begin vars */

   /* end vars */
   let v1 = a*(a-b)*(a+b-c)*(a+b-c)*(-a+c)*(a-b+c)*(a-b+c);
   let v2 = (a-b)*b*(b-c)*(a+b-c)*(a+b-c)*(-a+b+c)*(-a+b+c);
   let v3 = (b-c)*c*(-a+c)*(a-b+c)*(a-b+c)*(-a+b+c)*(-a+b+c);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X935(orbit) {
   /* begin vars */
   let a2=a*a;
   let c2=c*c;
   let b2=b*b;
   let SA=(b2+c2-a2)/2;
   let SC=(a2+b2-c2)/2;
   let SB=(c2+a2-b2)/2;
   /* end vars */
   let v1 = (SB*SC)/((SB-SC)*(SA*(a2+SA)-3*SB*SC));
   let v2 = (SA*SC)/((-SA+SC)*(SB*(b2+SB)-3*SA*SC));
   let v3 = (SA*SB)/((SA-SB)*(-3*SA*SB+SC*(c2+SC)));
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X936(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let c3=c2*c;
   let b3=b2*b;
   let a3=a2*a;
   /* end vars */
   let v1 = a*(a3-a2*b-a*b2+b3-a2*c+2*a*b*c+3*b2*c-a*c2+3*b*c2+c3);
   let v2 = b*(a3-a2*b-a*b2+b3+3*a2*c+2*a*b*c-b2*c+3*a*c2-b*c2+c3);
   let v3 = c*(a3+3*a2*b+3*a*b2+b3-a2*c+2*a*b*c-b2*c-a*c2-b*c2+c3);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X937(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let c3=c2*c;
   let b3=b2*b;
   let a3=a2*a;
   /* end vars */
   let v1 = a/(a3-a2*b-a*b2+b3-a2*c+2*a*b*c+3*b2*c-a*c2+3*b*c2+c3);
   let v2 = b/(a3-a2*b-a*b2+b3+3*a2*c+2*a*b*c-b2*c+3*a*c2-b*c2+c3);
   let v3 = c/(a3+3*a2*b+3*a*b2+b3-a2*c+2*a*b*c-b2*c-a*c2-b*c2+c3);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X938(orbit) {
   /* begin vars */
   let b2=b*b;
   let a2=a*a;
   let c2=c*c;
   let b3=b2*b;
   let b4=b2*b2;
   let a3=a2*a;
   let a4=a2*a2;
   let c4=c2*c2;
   let c3=c2*c;
   /* end vars */
   let v1 = Math.pow(a-b,3)*(a+b)-2*a*(a+b)*(a+b)*c-2*(a-b)*b*c2+2*a*c3-c4;
   let v2 = -a4+2*a3*b-2*a2*(b-c)*c+Math.pow(b-c,3)*(b+c)-2*a*b*(b+c)*(b+c);
   let v3 = -b4+2*b3*c-2*a*b2*(-a+c)+Math.pow(-a+c,3)*(a+c)-2*b*c*(a+c)*(a+c);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X939(orbit) {
   /* begin vars */
   let b2=b*b;
   let a2=a*a;
   let c2=c*c;
   let b3=b2*b;
   let b4=b2*b2;
   let a3=a2*a;
   let a4=a2*a2;
   let c4=c2*c2;
   let c3=c2*c;
   /* end vars */
   let v1 = a2/(Math.pow(a-b,3)*(a+b)-2*a*(a+b)*(a+b)*c-2*(a-b)*b*c2+2*a*c3-c4);
   let v2 = b2/(-a4+2*a3*b-2*a2*(b-c)*c+Math.pow(b-c,3)*(b+c)-2*a*b*(b+c)*(b+c));
   let v3 = c2/(-b4+2*b3*c-2*a*b2*(-a+c)+Math.pow(-a+c,3)*(a+c)-2*b*c*(a+c)*(a+c));
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X940(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   /* end vars */
   let v1 = a*(a2+a*b+a*c+2*b*c);
   let v2 = b*(a*b+b2+2*a*c+b*c);
   let v3 = c*(2*a*b+a*c+b*c+c2);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X941(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   /* end vars */
   let v1 = a/(a2+a*b+a*c+2*b*c);
   let v2 = b/(a*b+b2+2*a*c+b*c);
   let v3 = c/(2*a*b+a*c+b*c+c2);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X942(orbit) {
   /* begin vars */
   let a2=a*a;
   let c2=c*c;
   let b2=b*b;
   let SA=(b2+c2-a2)/2;
   let SC=(a2+b2-c2)/2;
   let SB=(c2+a2-b2)/2;
   /* end vars */
   let v1 = b*(a+c)*SB+(a+b)*c*SC;
   let v2 = a*(b+c)*SA+(a+b)*c*SC;
   let v3 = a*(b+c)*SA+b*(a+c)*SB;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X943(orbit) {
   /* begin vars */
   let a2=a*a;
   let c2=c*c;
   let b2=b*b;
   let SA=(b2+c2-a2)/2;
   let SC=(a2+b2-c2)/2;
   let SB=(c2+a2-b2)/2;
   /* end vars */
   let v1 = a2/(b*(a+c)*SB+(a+b)*c*SC);
   let v2 = b2/(a*(b+c)*SA+(a+b)*c*SC);
   let v3 = c2/(a*(b+c)*SA+b*(a+c)*SB);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X944(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let c3=c2*c;
   let b3=b2*b;
   let SC=(a2+b2-c2)/2;
   let SB=(c2+a2-b2)/2;
   let SA=(b2+c2-a2)/2;
   let a3=a2*a;
   /* end vars */
   let v1 = a3*(b+c)-a*(b-c)*(b-c)*(b+c)+a2*(-2*b*c+2*SA)-2*SB*SC;
   let v2 = b3*(a+c)-b*(-a+c)*(-a+c)*(a+c)+b2*(-2*a*c+2*SB)-2*SA*SC;
   let v3 = -((a-b)*(a-b)*(a+b)*c)+(a+b)*c3-2*SA*SB+c2*(-2*a*b+2*SC);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X945(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let c3=c2*c;
   let b3=b2*b;
   let SC=(a2+b2-c2)/2;
   let SB=(c2+a2-b2)/2;
   let SA=(b2+c2-a2)/2;
   let a3=a2*a;
   /* end vars */
   let v1 = a2/(a3*(b+c)-a*(b-c)*(b-c)*(b+c)+a2*(-2*b*c+2*SA)-2*SB*SC);
   let v2 = b2/(b3*(a+c)-b*(-a+c)*(-a+c)*(a+c)+b2*(-2*a*c+2*SB)-2*SA*SC);
   let v3 = c2/(-((a-b)*(a-b)*(a+b)*c)+(a+b)*c3-2*SA*SB+c2*(-2*a*b+2*SC));
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X946(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let c3=c2*c;
   let b3=b2*b;
   let SC=(a2+b2-c2)/2;
   let SB=(c2+a2-b2)/2;
   let SA=(b2+c2-a2)/2;
   let a3=a2*a;
   /* end vars */
   let v1 = a3*(b+c)-a*(b-c)*(b-c)*(b+c)-2*a2*(b*c-SA)+4*SB*SC;
   let v2 = b3*(a+c)-b*(-a+c)*(-a+c)*(a+c)-2*b2*(a*c-SB)+4*SA*SC;
   let v3 = -((a-b)*(a-b)*(a+b)*c)+(a+b)*c3+4*SA*SB-2*c2*(a*b-SC);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X947(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let c3=c2*c;
   let b3=b2*b;
   let SC=(a2+b2-c2)/2;
   let SB=(c2+a2-b2)/2;
   let SA=(b2+c2-a2)/2;
   let a3=a2*a;
   /* end vars */
   let v1 = a2/(a3*(b+c)-a*(b-c)*(b-c)*(b+c)-2*a2*(b*c-SA)+4*SB*SC);
   let v2 = b2/(b3*(a+c)-b*(-a+c)*(-a+c)*(a+c)-2*b2*(a*c-SB)+4*SA*SC);
   let v3 = c2/(-((a-b)*(a-b)*(a+b)*c)+(a+b)*c3+4*SA*SB-2*c2*(a*b-SC));
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X948(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let c3=c2*c;
   let b3=b2*b;
   let a3=a2*a;
   /* end vars */
   let v1 = (a+b-c)*(a-b+c)*(a3-a2*(b+c)-(b-c)*(b-c)*(b+c)+a*(b+c)*(b+c));
   let v2 = (a+b-c)*(-a+b+c)*(b3-b2*(a+c)-(-a+c)*(-a+c)*(a+c)+b*(a+c)*(a+c));
   let v3 = (a-b+c)*(-a+b+c)*(-((a-b)*(a-b)*(a+b))+(a+b)*(a+b)*c-(a+b)*c2+c3);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X949(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let c3=c2*c;
   let b3=b2*b;
   let a3=a2*a;
   /* end vars */
   let v1 = a2/((a+b-c)*(a-b+c)*(a3-a2*(b+c)-(b-c)*(b-c)*(b+c)+a*(b+c)*(b+c)));
   let v2 = b2/((a+b-c)*(-a+b+c)*(b3-b2*(a+c)-(-a+c)*(-a+c)*(a+c)+b*(a+c)*(a+c)));
   let v3 = c2/((a-b+c)*(-a+b+c)*(-((a-b)*(a-b)*(a+b))+(a+b)*(a+b)*c-(a+b)*c2+c3));
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X950(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let c3=c2*c;
   let b3=b2*b;
   let a3=a2*a;
   /* end vars */
   let v1 = (a-b-c)*(2*a3+a2*b+b3+a2*c-b2*c-b*c2+c3);
   let v2 = (-a+b-c)*(a3+a*b2+2*b3-a2*c+b2*c-a*c2+c3);
   let v3 = (-a-b+c)*(a3-a2*b-a*b2+b3+a*c2+b*c2+2*c3);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X951(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let c3=c2*c;
   let b3=b2*b;
   let a3=a2*a;
   /* end vars */
   let v1 = a2/((a-b-c)*(2*a3+a2*b+b3+a2*c-b2*c-b*c2+c3));
   let v2 = b2/((-a+b-c)*(a3+a*b2+2*b3-a2*c+b2*c-a*c2+c3));
   let v3 = c2/((-a-b+c)*(a3-a2*b-a*b2+b3+a*c2+b*c2+2*c3));
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X952(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let c3=c2*c;
   let b3=b2*b;
   let SC=(a2+b2-c2)/2;
   let SB=(c2+a2-b2)/2;
   let SA=(b2+c2-a2)/2;
   let a3=a2*a;
   /* end vars */
   let v1 = a3*(b+c)-a*(b-c)*(b-c)*(b+c)+a2*(-2*b*c+SA)-2*SB*SC;
   let v2 = b3*(a+c)-b*(-a+c)*(-a+c)*(a+c)+b2*(-2*a*c+SB)-2*SA*SC;
   let v3 = -((a-b)*(a-b)*(a+b)*c)+(a+b)*c3-2*SA*SB+c2*(-2*a*b+SC);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X953(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let c3=c2*c;
   let b3=b2*b;
   let SC=(a2+b2-c2)/2;
   let SB=(c2+a2-b2)/2;
   let SA=(b2+c2-a2)/2;
   let a3=a2*a;
   /* end vars */
   let v1 = a2/(a3*(b+c)-a*(b-c)*(b-c)*(b+c)+a2*(-2*b*c+SA)-2*SB*SC);
   let v2 = b2/(b3*(a+c)-b*(-a+c)*(-a+c)*(a+c)+b2*(-2*a*c+SB)-2*SA*SC);
   let v3 = c2/(-((a-b)*(a-b)*(a+b)*c)+(a+b)*c3-2*SA*SB+c2*(-2*a*b+SC));
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X954(orbit) {
   /* begin vars */
   let a2=a*a;
   let c2=c*c;
   let b2=b*b;
   let c3=c2*c;
   let b3=b2*b;
   let a3=a2*a;
   let SA=(b2+c2-a2)/2;
   let c5=c2*c3;
   let b5=b2*b3;
   let SC=(a2+b2-c2)/2;
   let SB=(c2+a2-b2)/2;
   let a5=a2*a3;
   /* end vars */
   let v1 = a*(a5-a*(SB-SC)*(SB-SC)+2*(b5+c5-c*(SA-SB)*(SA-SB)+2*b3*SB-b*(SA-SC)*(SA-SC)+2*c3*SC));
   let v2 = b*(b5-b*(-SA+SC)*(-SA+SC)+2*(a5+c5+2*a3*SA-c*(-SA+SB)*(-SA+SB)-a*(SB-SC)*(SB-SC)+2*c3*SC));
   let v3 = c*(c5-c*(SA-SB)*(SA-SB)+2*(a5+b5+2*a3*SA+2*b3*SB-b*(-SA+SC)*(-SA+SC)-a*(-SB+SC)*(-SB+SC)));
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X955(orbit) {
   /* begin vars */
   let a2=a*a;
   let c2=c*c;
   let b2=b*b;
   let c3=c2*c;
   let b3=b2*b;
   let a3=a2*a;
   let SA=(b2+c2-a2)/2;
   let c5=c2*c3;
   let b5=b2*b3;
   let SC=(a2+b2-c2)/2;
   let SB=(c2+a2-b2)/2;
   let a5=a2*a3;
   /* end vars */
   let v1 = a/(a5-a*(SB-SC)*(SB-SC)+2*(b5+c5-c*(SA-SB)*(SA-SB)+2*b3*SB-b*(SA-SC)*(SA-SC)+2*c3*SC));
   let v2 = b/(b5-b*(-SA+SC)*(-SA+SC)+2*(a5+c5+2*a3*SA-c*(-SA+SB)*(-SA+SB)-a*(SB-SC)*(SB-SC)+2*c3*SC));
   let v3 = c/(c5-c*(SA-SB)*(SA-SB)+2*(a5+b5+2*a3*SA+2*b3*SB-b*(-SA+SC)*(-SA+SC)-a*(-SB+SC)*(-SB+SC)));
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X956(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let SC=(a2+b2-c2)/2;
   let SB=(c2+a2-b2)/2;
   let SA=(b2+c2-a2)/2;
   /* end vars */
   let v1 = a*b*c*(-a+b+c)+a2*SA;
   let v2 = a*b*c*(a-b+c)+b2*SB;
   let v3 = a*b*(a+b-c)*c+c2*SC;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X957(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let SC=(a2+b2-c2)/2;
   let SB=(c2+a2-b2)/2;
   let SA=(b2+c2-a2)/2;
   /* end vars */
   let v1 = a2/(a*b*c*(-a+b+c)+a2*SA);
   let v2 = b2/(a*b*c*(a-b+c)+b2*SB);
   let v3 = c2/(a*b*(a+b-c)*c+c2*SC);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X958(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   /* end vars */
   let v1 = a*(a-b-c)*(a2+a*b+a*c+2*b*c);
   let v2 = b*(-a+b-c)*(a*b+b2+2*a*c+b*c);
   let v3 = c*(-a-b+c)*(2*a*b+a*c+b*c+c2);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X959(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   /* end vars */
   let v1 = a/((a-b-c)*(a2+a*b+a*c+2*b*c));
   let v2 = b/((-a+b-c)*(a*b+b2+2*a*c+b*c));
   let v3 = c/((-a-b+c)*(2*a*b+a*c+b*c+c2));
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X960(orbit) {
   /* begin vars */
   let a2=a*a;
   let c2=c*c;
   let b2=b*b;
   /* end vars */
   let v1 = a*(-a+b+c)*(a*b+b2+a*c+c2);
   let v2 = b*(a-b+c)*(a2+a*b+b*c+c2);
   let v3 = (a+b-c)*c*(a2+b2+a*c+b*c);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X961(orbit) {
   /* begin vars */
   let a2=a*a;
   let c2=c*c;
   let b2=b*b;
   /* end vars */
   let v1 = a/((-a+b+c)*(a*b+b2+a*c+c2));
   let v2 = b/((a-b+c)*(a2+a*b+b*c+c2));
   let v3 = c/((a+b-c)*(a2+b2+a*c+b*c));
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X962(orbit) {
   /* begin vars */
   let b2=b*b;
   let a2=a*a;
   let c2=c*c;
   let b3=b2*b;
   let b4=b2*b2;
   let a3=a2*a;
   let a4=a2*a2;
   let c4=c2*c2;
   let c3=c2*c;
   /* end vars */
   let v1 = (a-b)*Math.pow(a+b,3)+2*a*(a-b)*(a-b)*c+2*b*(a+b)*c2-2*a*c3-c4;
   let v2 = -a4-2*a3*b+2*a*b*(b-c)*(b-c)+2*a2*c*(b+c)+(b-c)*Math.pow(b+c,3);
   let v3 = -b4-2*b3*c+2*b*c*(-a+c)*(-a+c)+2*a*b2*(a+c)+(-a+c)*Math.pow(a+c,3);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X963(orbit) {
   /* begin vars */
   let b2=b*b;
   let a2=a*a;
   let c2=c*c;
   let b3=b2*b;
   let b4=b2*b2;
   let a3=a2*a;
   let a4=a2*a2;
   let c4=c2*c2;
   let c3=c2*c;
   /* end vars */
   let v1 = a2/((a-b)*Math.pow(a+b,3)+2*a*(a-b)*(a-b)*c+2*b*(a+b)*c2-2*a*c3-c4);
   let v2 = b2/(-a4-2*a3*b+2*a*b*(b-c)*(b-c)+2*a2*c*(b+c)+(b-c)*Math.pow(b+c,3));
   let v3 = c2/(-b4-2*b3*c+2*b*c*(-a+c)*(-a+c)+2*a*b2*(a+c)+(-a+c)*Math.pow(a+c,3));
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X964(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let c4=c2*c2;
   let c3=c2*c;
   let b3=b2*b;
   let b4=b2*b2;
   let a3=a2*a;
   let a4=a2*a2;
   /* end vars */
   let v1 = a4+a3*(b+c)+a2*(b+c)*(b+c)+b*c*(b+c)*(b+c)+a*(b+c)*(b2+b*c+c2);
   let v2 = b4+b3*(a+c)+b2*(a+c)*(a+c)+a*c*(a+c)*(a+c)+b*(a+c)*(a2+a*c+c2);
   let v3 = a*b*(a+b)*(a+b)+(a+b)*(a2+a*b+b2)*c+(a+b)*(a+b)*c2+(a+b)*c3+c4;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X965(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let c4=c2*c2;
   let c3=c2*c;
   let b3=b2*b;
   let b4=b2*b2;
   let a3=a2*a;
   let a4=a2*a2;
   /* end vars */
   let v1 = a*(a4-a3*(b+c)+2*b*c*(b+c)*(b+c)+a*Math.pow(b+c,3)-a2*(b2+c2));
   let v2 = b*(b4-b3*(a+c)+2*a*c*(a+c)*(a+c)+b*Math.pow(a+c,3)-b2*(a2+c2));
   let v3 = c*(2*a*b*(a+b)*(a+b)+Math.pow(a+b,3)*c-(a2+b2)*c2-(a+b)*c3+c4);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X966(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let SC=(a2+b2-c2)/2;
   let SB=(c2+a2-b2)/2;
   let SA=(b2+c2-a2)/2;
   /* end vars */
   let v1 = b*c+a*(b+c)+SA;
   let v2 = a*c+b*(a+c)+SB;
   let v3 = a*b+(a+b)*c+SC;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X967(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let SC=(a2+b2-c2)/2;
   let SB=(c2+a2-b2)/2;
   let SA=(b2+c2-a2)/2;
   /* end vars */
   let v1 = a2/(b*c+a*(b+c)+SA);
   let v2 = b2/(a*c+b*(a+c)+SB);
   let v3 = c2/(a*b+(a+b)*c+SC);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X968(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   /* end vars */
   let v1 = a*(a2-2*a*(b+c)-(b+c)*(b+c));
   let v2 = b*(b2-2*b*(a+c)-(a+c)*(a+c));
   let v3 = c*(-(a+b)*(a+b)-2*(a+b)*c+c2);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X969(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   /* end vars */
   let v1 = a/(a2-2*a*(b+c)-(b+c)*(b+c));
   let v2 = b/(b2-2*b*(a+c)-(a+c)*(a+c));
   let v3 = c/(-(a+b)*(a+b)-2*(a+b)*c+c2);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X970(orbit) {
   /* begin vars */
   let a2=a*a;
   let c2=c*c;
   let b2=b*b;
   let a4=a2*a2;
   let c3=c2*c;
   let b3=b2*b;
   let c4=c2*c2;
   let b4=b2*b2;
   let a3=a2*a;
   /* end vars */
   let v1 = a2*(a3*(b+c)*(b+c)+a2*(b+c)*(b2+c2)-(b+c)*(b4+c4)-a*(b4+2*b3*c+2*b*c3+c4));
   let v2 = b2*(b3*(a+c)*(a+c)+b2*(a+c)*(a2+c2)-(a+c)*(a4+c4)-b*(a4+2*a3*c+2*a*c3+c4));
   let v3 = c2*(-((a+b)*(a4+b4))-(a4+2*a3*b+2*a*b3+b4)*c+(a+b)*(a2+b2)*c2+(a+b)*(a+b)*c3);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X971(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let c4=c2*c2;
   let c3=c2*c;
   let b3=b2*b;
   let b4=b2*b2;
   let a3=a2*a;
   let a4=a2*a2;
   /* end vars */
   let v1 = a*(a4*(b+c)-(b-c)*(b-c)*Math.pow(b+c,3)-2*a3*(b2-b*c+c2)+2*a*(b-c)*(b-c)*(b2+b*c+c2));
   let v2 = b*(b4*(a+c)-(-a+c)*(-a+c)*Math.pow(a+c,3)-2*b3*(a2-a*c+c2)+2*b*(-a+c)*(-a+c)*(a2+a*c+c2));
   let v3 = c*(-((a-b)*(a-b)*Math.pow(a+b,3))+2*(a-b)*(a-b)*(a2+a*b+b2)*c-2*(a2-a*b+b2)*c3+(a+b)*c4);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X972(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let c4=c2*c2;
   let c3=c2*c;
   let b3=b2*b;
   let b4=b2*b2;
   let a3=a2*a;
   let a4=a2*a2;
   /* end vars */
   let v1 = a/(a4*(b+c)-(b-c)*(b-c)*Math.pow(b+c,3)-2*a3*(b2-b*c+c2)+2*a*(b-c)*(b-c)*(b2+b*c+c2));
   let v2 = b/(b4*(a+c)-(-a+c)*(-a+c)*Math.pow(a+c,3)-2*b3*(a2-a*c+c2)+2*b*(-a+c)*(-a+c)*(a2+a*c+c2));
   let v3 = c/(-((a-b)*(a-b)*Math.pow(a+b,3))+2*(a-b)*(a-b)*(a2+a*b+b2)*c-2*(a2-a*b+b2)*c3+(a+b)*c4);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X973(orbit) {
   /* begin vars */
   let b2=b*b;
   let a2=a*a;
   let c2=c*c;
   let SC=(a2+b2-c2)/2;
   let SB=(c2+a2-b2)/2;
   let SA=(b2+c2-a2)/2;
   let b4=b2*b2;
   let a4=a2*a2;
   let SC2=SC*SC;
   let c4=c2*c2;
   let SB2=SB*SB;
   let SA2=SA*SA;
   /* end vars */
   let v1 = a2*(a2*SA+2*SB*SC)*(SA*SB*(SA2+7*SA*SB+2*SB2)*SC+c4*Math.pow(SC,3)+c2*(SA2*SB2+(SA2+6*SA*SB+SB2)*SC2));
   let v2 = b2*(b2*SB+2*SA*SC)*(a4*Math.pow(SA,3)+SA*SB*SC*(SB2+7*SB*SC+2*SC2)+a2*(SB2*SC2+SA2*(SB2+6*SB*SC+SC2)));
   let v3 = c2*(2*SA*SB+c2*SC)*(b4*Math.pow(SB,3)+SA*SB*SC*(2*SA2+7*SA*SC+SC2)+b2*(SA2*SC2+SB2*(SA2+6*SA*SC+SC2)));
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X974(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let SC=(a2+b2-c2)/2;
   let SB=(c2+a2-b2)/2;
   let SA=(b2+c2-a2)/2;
   let SC2=SC*SC;
   let b4=b2*b2;
   let SB2=SB*SB;
   let a4=a2*a2;
   let SA2=SA*SA;
   let c4=c2*c2;
   /* end vars */
   let v1 = a2*SA*(c4*Math.pow(SC,4)+Math.pow(SB,3)*SC*(-2*SA2+2*SA*SB-3*SA*SC+SB*SC)+c2*SA*(SA*Math.pow(SB,3)+(SA-3*SB)*Math.pow(SC,3)));
   let v2 = b2*SB*(a4*Math.pow(SA,4)+SA*Math.pow(SC,3)*(-3*SA*SB-2*SB2+SA*SC+2*SB*SC)+a2*SB*(Math.pow(SA,3)*(SB-3*SC)+SB*Math.pow(SC,3)));
   let v3 = c2*SC*(b4*Math.pow(SB,4)+b2*SC*(Math.pow(SA,3)*SC+Math.pow(SB,3)*(-3*SA+SC))+Math.pow(SA,3)*SB*(SA*SB+2*SA*SC-3*SB*SC-2*SC2));
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X975(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let c3=c2*c;
   let b3=b2*b;
   let a3=a2*a;
   /* end vars */
   let v1 = a*(a3+a2*(b+c)+Math.pow(b+c,3)+a*(b2+4*b*c+c2));
   let v2 = b*(b3+b2*(a+c)+Math.pow(a+c,3)+b*(a2+4*a*c+c2));
   let v3 = c*(Math.pow(a+b,3)+(a2+4*a*b+b2)*c+(a+b)*c2+c3);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X976(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let c3=c2*c;
   let b3=b2*b;
   let a3=a2*a;
   /* end vars */
   let v1 = a*(a3+(b+c)*(b2+c2));
   let v2 = b*(b3+(a+c)*(a2+c2));
   let v3 = c*((a+b)*(a2+b2)+c3);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X977(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let c3=c2*c;
   let b3=b2*b;
   let a3=a2*a;
   /* end vars */
   let v1 = a/(a3+(b+c)*(b2+c2));
   let v2 = b/(b3+(a+c)*(a2+c2));
   let v3 = c/((a+b)*(a2+b2)+c3);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X978(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   /* end vars */
   let v1 = a*(a2*(b+c)-b*c*(b+c)+a*(b2-b*c+c2));
   let v2 = b*(b2*(a+c)-a*c*(a+c)+b*(a2-a*c+c2));
   let v3 = c*(-(a*b*(a+b))+(a2-a*b+b2)*c+(a+b)*c2);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X979(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   /* end vars */
   let v1 = a/(a2*(b+c)-b*c*(b+c)+a*(b2-b*c+c2));
   let v2 = b/(b2*(a+c)-a*c*(a+c)+b*(a2-a*c+c2));
   let v3 = c/(-(a*b*(a+b))+(a2-a*b+b2)*c+(a+b)*c2);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X980(orbit) {
   /* begin vars */
   let a2=a*a;
   let c2=c*c;
   let b2=b*b;
   /* end vars */
   let v1 = a*((a*b+a*c+b*c)*(b2+c2)+a2*(b2+b*c+c2));
   let v2 = b*((a*b+a*c+b*c)*(a2+c2)+b2*(a2+a*c+c2));
   let v3 = c*((a2+b2)*(a*b+a*c+b*c)+(a2+a*b+b2)*c2);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X981(orbit) {
   /* begin vars */
   let a2=a*a;
   let c2=c*c;
   let b2=b*b;
   /* end vars */
   let v1 = a/((a*b+a*c+b*c)*(b2+c2)+a2*(b2+b*c+c2));
   let v2 = b/((a*b+a*c+b*c)*(a2+c2)+b2*(a2+a*c+c2));
   let v3 = c/((a2+b2)*(a*b+a*c+b*c)+(a2+a*b+b2)*c2);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X982(orbit) {
   /* begin vars */
   let a2=a*a;
   let c2=c*c;
   let b2=b*b;
   /* end vars */
   let v1 = a*(b2-b*c+c2);
   let v2 = b*(a2-a*c+c2);
   let v3 = (a2-a*b+b2)*c;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X983(orbit) {
   /* begin vars */
   let a2=a*a;
   let c2=c*c;
   let b2=b*b;
   /* end vars */
   let v1 = a/(b2-b*c+c2);
   let v2 = b/(a2-a*c+c2);
   let v3 = c/(a2-a*b+b2);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X984(orbit) {
   /* begin vars */
   let a2=a*a;
   let c2=c*c;
   let b2=b*b;
   /* end vars */
   let v1 = a*(b2+b*c+c2);
   let v2 = b*(a2+a*c+c2);
   let v3 = (a2+a*b+b2)*c;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X985(orbit) {
   /* begin vars */
   let a2=a*a;
   let c2=c*c;
   let b2=b*b;
   /* end vars */
   let v1 = a/(b2+b*c+c2);
   let v2 = b/(a2+a*c+c2);
   let v3 = c/(a2+a*b+b2);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X986(orbit) {
   /* begin vars */
   let a2=a*a;
   let c2=c*c;
   let b2=b*b;
   let a3=a2*a;
   let c3=c2*c;
   let b3=b2*b;
   /* end vars */
   let v1 = a*(b3+a*(b2+b*c+c2)+c3);
   let v2 = b*(a3+b*(a2+a*c+c2)+c3);
   let v3 = c*(a3+b3+(a2+a*b+b2)*c);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X987(orbit) {
   /* begin vars */
   let a2=a*a;
   let c2=c*c;
   let b2=b*b;
   let a3=a2*a;
   let c3=c2*c;
   let b3=b2*b;
   /* end vars */
   let v1 = a/(b3+a*(b2+b*c+c2)+c3);
   let v2 = b/(a3+b*(a2+a*c+c2)+c3);
   let v3 = c/(a3+b3+(a2+a*b+b2)*c);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X988(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let c3=c2*c;
   let b3=b2*b;
   let a3=a2*a;
   /* end vars */
   let v1 = a*(a3-a2*(b+c)-(3*a+b+c)*(b2+c2));
   let v2 = b*(b3-b2*(a+c)-(a+3*b+c)*(a2+c2));
   let v3 = c*(-((a2+b2)*(a+b+3*c))-(a+b)*c2+c3);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X989(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let c3=c2*c;
   let b3=b2*b;
   let a3=a2*a;
   /* end vars */
   let v1 = a/(a3-a2*(b+c)-(3*a+b+c)*(b2+c2));
   let v2 = b/(b3-b2*(a+c)-(a+3*b+c)*(a2+c2));
   let v3 = c/(-((a2+b2)*(a+b+3*c))-(a+b)*c2+c3);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X990(orbit) {
   /* begin vars */
   let c2=c*c;
   let c3=c2*c;
   let b2=b*b;
   let b3=b2*b;
   let a2=a*a;
   let a3=a2*a;
   let c5=c2*c3;
   let c4=c2*c2;
   let b4=b2*b2;
   let b5=b2*b3;
   let a4=a2*a2;
   let a5=a2*a3;
   /* end vars */
   let v1 = a*(a5-2*a3*b*c-a4*(b+c)+(b-c)*(b-c)*Math.pow(b+c,3)-a*(b-c)*(b-c)*(b2+c2));
   let v2 = b*(b5-2*a*b3*c-b4*(a+c)+(-a+c)*(-a+c)*Math.pow(a+c,3)-b*(-a+c)*(-a+c)*(a2+c2));
   let v3 = c*((a-b)*(a-b)*Math.pow(a+b,3)-(a-b)*(a-b)*(a2+b2)*c-2*a*b*c3-(a+b)*c4+c5);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X991(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let c3=c2*c;
   let b3=b2*b;
   let a3=a2*a;
   /* end vars */
   let v1 = a2*(a3*(b+c)-a*(b+c)*(b2+c2)-a2*(b2-b*c+c2)+(b-c)*(b3-c3));
   let v2 = b2*(b3*(a+c)-b*(a+c)*(a2+c2)-b2*(a2-a*c+c2)+(-a+c)*(-a3+c3));
   let v3 = c2*((a-b)*(a3-b3)-(a+b)*(a2+b2)*c-(a2-a*b+b2)*c2+(a+b)*c3);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X992(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let c3=c2*c;
   let b3=b2*b;
   let a3=a2*a;
   /* end vars */
   let v1 = a*(a3*(b+c)-a*b*c*(b+c)-b*c*(b+c)*(b+c)+a2*(b2+c2));
   let v2 = b*(b3*(a+c)-a*b*c*(a+c)-a*c*(a+c)*(a+c)+b2*(a2+c2));
   let v3 = c*(-(a*b*(a+b)*(a+b))-a*b*(a+b)*c+(a2+b2)*c2+(a+b)*c3);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X993(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let c3=c2*c;
   let b3=b2*b;
   let a3=a2*a;
   /* end vars */
   let v1 = a*(a3-b*c*(b+c)-a*(b2+c2));
   let v2 = b*(b3-a*c*(a+c)-b*(a2+c2));
   let v3 = c*(-(a*b*(a+b))-(a2+b2)*c+c3);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X994(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let c3=c2*c;
   let b3=b2*b;
   let a3=a2*a;
   /* end vars */
   let v1 = a/(a3-b*c*(b+c)-a*(b2+c2));
   let v2 = b/(b3-a*c*(a+c)-b*(a2+c2));
   let v3 = c/(-(a*b*(a+b))-(a2+b2)*c+c3);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X995(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   /* end vars */
   let v1 = a2*(a*b+b2+a*c-b*c+c2);
   let v2 = b2*(a2+a*b-a*c+b*c+c2);
   let v3 = (a2-a*b+b2+a*c+b*c)*c2;
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X996(orbit) {
   /* begin vars */
   let a2=a*a;
   let c2=c*c;
   let b2=b*b;
   /* end vars */
   let v1 = 1/(a*b+b2+a*c-b*c+c2);
   let v2 = 1/(a2+a*b-a*c+b*c+c2);
   let v3 = 1/(a2-a*b+b2+a*c+b*c);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X997(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let c3=c2*c;
   let b3=b2*b;
   let a3=a2*a;
   /* end vars */
   let v1 = a*(a3-a*(b-c)*(b-c)-a2*(b+c)+(b+c)*(b2+c2));
   let v2 = b*(b3-b*(-a+c)*(-a+c)-b2*(a+c)+(a+c)*(a2+c2));
   let v3 = c*((a+b)*(a2+b2)-(a-b)*(a-b)*c-(a+b)*c2+c3);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X998(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let c3=c2*c;
   let b3=b2*b;
   let a3=a2*a;
   /* end vars */
   let v1 = a/(a3-a*(b-c)*(b-c)-a2*(b+c)+(b+c)*(b2+c2));
   let v2 = b/(b3-b*(-a+c)*(-a+c)-b2*(a+c)+(a+c)*(a2+c2));
   let v3 = c/((a+b)*(a2+b2)-(a-b)*(a-b)*c-(a+b)*c2+c3);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}

function barys_X999(orbit) {
   /* begin vars */
   let c2=c*c;
   let b2=b*b;
   let a2=a*a;
   let SC=(a2+b2-c2)/2;
   let SB=(c2+a2-b2)/2;
   let SA=(b2+c2-a2)/2;
   /* end vars */
   let v1 = a2*(2*b*c-SA);
   let v2 = b2*(2*a*c-SB);
   let v3 = c2*(2*a*b-SC);
   let barys = [v1,v2,v3];
   return barys_to_cartesian(orbit, barys);
}