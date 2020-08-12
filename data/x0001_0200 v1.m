(* ::Package:: *)

{
{"X(1)", Hold@{1, 1, 1}, "INCENTER"},
{"X(2)", Hold@{b c, c a, a b}, "CENTROID"},
{"X(3)", Hold@{cosA, cosB, cosC}, "CIRCUMCENTER"},
{"X(4)", Hold@{secA, secB, secC}, "ORTHOCENTER"},
{"X(5)", Hold@{cosB cosC+sinB sinC, cosC cosA+sinC sinA, cosA cosB+sinA sinB}, "NINE-POINT CENTER"},
{"X(6)", Hold@{a, b, c}, "SYMMEDIAN / LEMOINE / GREBE POINT"},
{"X(7)", Hold@{b c/(b+c-a), c a/(c+a-b), a b/(a+b-c)}, "GERGONNE POINT"},
{"X(8)", Hold@{(b+c-a)/a, (c+a-b)/b, (a+b-c)/c},"NAGEL POINT"},
{"X(9)", Hold@{b+c-a, c+a-b, a+b-c}, "MITTENPUNKT"},
{"X(10)", Hold@{b c (b+c), c a (c+a), a b (a+b)}, "SPIEKER CENTER"},
{"X(11)", Hold@{1-cosB cosC-sinB sinC, 1-cosC cosA-sinC sinA, 1-cosA cosB-sinA sinB}, "FEUERBACH POINT"},
{"X(12)", Hold@{1+cosB cosC+sinB sinC, 1+cosC cosA+sinC sinA, 1+cosA cosB+sinA sinB}, "{X(1),X(5)}-HARMONIC CONJUGATE OF X(11)"},
{"X(13)", Hold@{cscApPi3, cscBpPi3, cscCpPi3}, "1st ISOGONIC CENTER (FERMAT / TORRICELLI POINT)"},
{"X(14)", Hold@{cscAmPi3, cscBmPi3, cscCmPi3}, "2nd ISOGONIC CENTER"},
{"X(15)", Hold@{sinApPi3, sinBpPi3, sinCpPi3}, "1st ISODYNAMIC POINT"},
{"X(16)", Hold@{sinAmPi3, sinBmPi3, sinCmPi3}, "2nd ISODYNAMIC POINT"},
{"X(17)", Hold@{cscApPi6, cscBpPi6, cscCpPi6}, "1st NAPOLEON POINT"},
{"X(18)", Hold@{cscAmPi6, cscBmPi6, cscCmPi6}, "2nd NAPOLEON POINT"},
{"X(19)", Hold@{1/(b2+c2-a2), 1/(a2+c2-b2), 1/(a2+b2-c2)}, "CLAWSON POINT"}, 
{"X(20)", Hold@{cosA-cosB cosC, cosB-cosC cosA, cosC-cosA cosB}, "DE LONGCHAMPS POINT"}, 
{"X(21)", Hold@{(b+c-a)/(b+c), (a+c-b)/(a+c), (a+b-c)/(a+b)}, "SCHIFFLER POINT"},
{"X(22)", Hold@{a (b4+c4-a4), b (c4+a4-b4), c (a4+b4-c4)}, "EXETER POINT"},
{"X(23)", Hold@{a (b4+c4-a4-b2 c2), b (a4+c4-b4-a2 c2), c (b4+a4-c4-b2 a2)}, "FAR-OUT POINT"},
{"X(24)", Hold@{secA cos2A, secB cos2B, secC cos2C}, "PERSPECTOR OF ABC AND ORTHIC-OF-ORTHIC TRIANGLE"},
{"X(25)", Hold@{a/(b2+c2-a2), b/(c2+a2-b2), c/(a2+b2-c2)}, "HOMOTHETIC CENTER OF ORTHIC AND TANGENTIAL TRIANGLES"},
{"X(26)", Hold@{a (b2 cos2B+c2 cos2C-a2 cos2A), b (a2 cos2A+c2 cos2C-b2 cos2B), c (a2 cos2A+b2 cos2B-c2 cos2C)}, "CIRCUMCENTER OF THE TANGENTIAL TRIANGLE"},
{"X(27)", Hold@{secA/(b+c), secB/(c+a), secC/(a+b)}, "CEVAPOINT OF ORTHOCENTER AND CLAWSON CENTER"}, 
{"X(28)", Hold@{tanA/(b+c), tanB/(c+a), tanC/(a+b)}, "CEVAPOINT OF X(19) AND X(25)"}, 
{"X(29)", Hold@{secA/(cosB+cosC), secB/(cosC+cosA), secC/(cosA+cosB)}, "CEVAPOINT OF INCENTER AND ORTHOCENTER"},
{"X(30)", Hold@{b+c-a, c+a-b, a+b-c}(* WRONG: copy mitten to avoid infinity *), "EULER INFINITY POINT"},
{"X(31)", Hold@{a2, b2, c2}, "2nd POWER POINT"},
{"X(32)", Hold@{a3, b3, c3}, "3rd POWER POINT"},
{"X(33)", Hold@{1+secA, 1+secB, 1+secC}, "PERSPECTOR OF THE ORTHIC AND INTANGENTS TRIANGLES"}, 
{"X(34)", Hold@{1-secA, 1-secB, 1-secC}, "X(4)-BETH CONJUGATE OF X(4)"}, 
{"X(35)", Hold@{1+2 cosA, 1+2 cosB, 1+2 cosC}, "{X(1),X(3)}-HARMONIC CONJUGATE OF X(36)"}, 
{"X(36)", Hold@{1-2 cosA, 1-2 cosB, 1-2 cosC}, "INVERSE-IN-CIRCUMCIRCLE OF INCENTER"}, 
{"X(37)", Hold@{b+c, c+a, a+b}, "CROSSPOINT OF INCENTER AND CENTROID"}, 
{"X(38)", Hold@{b2+c2, c2+a2,a2+b2}, "CROSSPOINT OF X(1) AND X(75)"}, 
{"X(39)", Hold@{a (b2+c2), b (c2+a2), c (a2+b2)}, "BROCARD MIDPOINT"},
{"X(40)", Hold@{cosB+cosC-cosA-1, cosA+cosC-cosB-1, cosA+cosB-cosC-1}, "BEVAN POINT"},
{"X(41)", Hold@{a2 (b+c-a), b2 (c+a-b), c2 (a+b-c)}, "X(6)-CEVA CONJUGATE OF X(31)"}, 
{"X(42)", Hold@{a (b+c), b (c+a), c (a+b)}, "CROSSPOINT OF INCENTER AND SYMMEDIAN POINT"},
{"X(43)", Hold@{a b+a c-b c, b c+b a-c a, c a+c b-a b}, "X(6)-CEVA CONJUGATE OF X(1)"}, 
{"X(44)", Hold@{b+c-2 a, c+a-2 b, a+b-2 c}, "X(6)-LINE CONJUGATE OF X(1)"}, 
{"X(45)", Hold@{2 b+2 c-a, 2 c+2 a-b, 2 a+2 b-c}, "X(9)-BETH CONJUGATE OF X(1)"}, 
{"X(46)", Hold@{cosB+cosC-cosA, cosA+cosC-cosB, cosA+cosB-cosC}, "X(4)-CEVA CONJUGATE OF X(1)"}, 
{"X(47)", Hold@{cos2A, cos2B, cos2C}, "X(110)-BETH CONJUGATE OF X(34)"},
{"X(48)", Hold@{tanB+tanC, tanA+tanC, tanA+tanB}, "CROSSPOINT OF X(1) AND X(63)"}, 
{"X(49)", Hold@{cos3A, cos3B, cos3C}, "CENTER OF SINE-TRIPLE-ANGLE CIRCLE"},
{"X(50)", Hold@{sin3A, sin3B, sin3C}, "X(74)-CEVA CONJUGATE OF X(184)"},
{"X(51)", Hold@{a2 (cosB cosC+sinB sinC), b2 (cosC cosA+sinC sinA), c2 (cosA cosB+sinA sinB)}, "CENTROID OF ORTHIC TRIANGLE"}, 
{"X(52)", Hold@{cos2A (cosB cosC+sinB sinC), cos2B (cosC cosA+sinC sinA), cos2C (cosA cosB+sinA sinB)}, "ORTHOCENTER OF ORTHIC TRIANGLE"}, 
{"X(53)", Hold@{tanA (cosB cosC+sinB sinC), tanB (cosC cosA+sinC sinA), tanC (cosA cosB+sinA sinB)}, "SYMMEDIAN POINT OF ORTHIC TRIANGLE"}, 
{"X(54)", Hold@{1/(cosB cosC+sinB sinC), 1/(cosC cosA+sinC sinA), 1/(cosA cosB+sinA sinB)}, "KOSNITA POINT"}, 
{"X(55)", Hold@{a (b+c-a), b (c+a-b), c (a+b-c)}, "INSIMILICENTER(CIRCUMCIRCLE, INCIRCLE)"},
{"X(56)", Hold@{a/(b+c-a), b/(c+a-b), c/(a+b-c)}, "EXSIMILICENTER(CIRCUMCIRCLE, INCIRCLE)"}, 
{"X(57)", Hold@{1/(b+c-a), 1/(c+a-b), 1/(a+b-c)}, "ISOGONAL CONJUGATE OF X(9)"}, 
{"X(58)", Hold@{a/(b+c), b/(c+a), c/(a+b)}, "ISOGONAL CONJUGATE OF X(10)"},
{"X(59)", Hold@{1/(1-(cosB cosC+sinB sinC)), 1/(1-(cosC cosA+sinC sinA)), 1/(1-(cosA cosB+sinA sinB))}, "ISOGONAL CONJUGATE OF X(11)"},
{"X(60)", Hold@{1/(1+(cosB cosC+sinB sinC)), 1/(1+(cosC cosA+sinC sinA)), 1/(1+(cosA cosB+sinA sinB))}, "ISOGONAL CONJUGATE OF X(12)"},
{"X(61)", Hold@{(sinA cPi6+sPi6 cosA), (sinB cPi6+sPi6 cosB), (sinC cPi6+sPi6 cosC)}, "ISOGONAL CONJUGATE OF X(17)"}, 
{"X(62)", Hold@{(sinA cPi6-sPi6 cosA), (sinB cPi6-sPi6 cosB), (sinC cPi6-sPi6 cosC)}, "ISOGONAL CONJUGATE OF X(18)"},
{"X(63)", Hold@{cotA, cotB, cotC}, "ISOGONAL CONJUGATE OF X(19)"},
{"X(64)", Hold@{1/(cosA-cosB cosC), 1/(cosB-cosC cosA), 1/(cosC-cosA cosB)}, "ISOGONAL CONJUGATE OF X(20)"}, 
{"X(65)", Hold@{cosB+cosC, cosC+cosA, cosA+cosB}, "ORTHOCENTER OF THE INTOUCH TRIANGLE"},
{"X(66)", Hold@{b c/(b4+c4-a4), c a/(c4+a4-b4), a b/(a4+b4-c4)}, "ISOGONAL CONJUGATE OF X(22)"}, 
{"X(67)", Hold@{b c/(b4+c4-a4-b2 c2), c a/(c4+a4-b4-c2 a2), a b/(a4+b4-c4-a2 b2)}, "ISOGONAL CONJUGATE OF X(23)"}, 
{"X(68)", Hold@{cosA sec2A, cosB sec2B, cosC sec2C}, "PRASOLOV POINT"},
{"X(69)", Hold@{(cosA)/a2, (cosB)/b2, (cosC)/c2}, "SYMMEDIAN POINT OF THE ANTICOMPLEMENTARY TRIANGLE"},
{"X(70)", Hold@{b c/(b2 cos2B+c2 cos2C-a2 cos2A), c a/(c2 cos2C+a2 cos2A-b2 cos2B), a b/(a2 cos2A+b2 cos2B-c2 cos2C)}, "ISOGONAL CONJUGATE OF X(26)"},
{"X(71)", Hold@{(b+c) cosA, (c+a) cosB, (a+b) cosC}, "ISOGONAL CONJUGATE OF X(27)"}, 
{"X(72)", Hold@{(b+c) cotA, (c+a) cotB, (a+b) cotC}, "ISOGONAL CONJUGATE OF X(28)"}, 
{"X(73)", Hold@{secB+secC, secC+secA, secA+secB}, "CROSSPOINT OF INCENTER AND CIRCUMCENTER"},
{"X(74)", Hold@{1/(cosA-2 cosB cosC), 1/(cosB-2 cosC cosA), 1/(cosC-2 cosA cosB)}, "ISOGONAL CONJUGATE OF EULER INFINITY POINT"}, 
{"X(75)", Hold@{1/a2, 1/b2, 1/c2}, "ISOTOMIC CONJUGATE OF INCENTER"},
{"X(76)", Hold@{1/a3, 1/b3, 1/c3}, "3rd BROCARD POINT"},
{"X(77)", Hold@{1/(1+secA), 1/(1+secB), 1/(1+secC)}, "ISOGONAL CONJUGATE OF X(33)"},
{"X(78)", Hold@{1/(1-secA), 1/(1-secB), 1/(1-secC)}, "ISOGONAL CONJUGATE OF X(34)"},
{"X(79)", Hold@{1/(1+2 cosA), 1/(1+2 cosB), 1/(1+2 cosC)}, "ISOGONAL CONJUGATE OF X(35)"},
{"X(80)", Hold@{1/(1-2 cosA), 1/(1-2 cosB), 1/(1-2 cosC)}, "REFLECTION OF INCENTER IN FEUERBACH POINT"}, 
{"X(81)", Hold@{1/(b+c), 1/(c+a), 1/(a+b)}, "CEVAPOINT OF INCENTER AND SYMMEDIAN POINT"},
{"X(82)", Hold@{1/(b2+c2), 1/(c2+a2), 1/(a2+b2)}, "ISOGONAL CONJUGATE OF X(38)"},
{"X(83)", Hold@{b c/(b2+c2), a c/(c2+a2), a b/(a2+b2)}, "CEVAPOINT OF CENTROID AND SYMMEDIAN POINT"},
{"X(84)", Hold@{1/(cosB+cosC-cosA-1), 1/(cosA+cosC-cosB-1), 1/(cosA+cosB-cosC-1)}, "ISOGONAL CONJUGATE OF X(40)"},
{"X(85)", Hold@{b2 c2/(b+c-a), c2 a2/(c+a-b), a2 b2/(a+b-c)}, "ISOTOMIC CONJUGATE OF X(9)"},
{"X(86)", Hold@{(b c)/(b+c), (c a)/(c+a), (a b)/(a+b)}, "CEVAPOINT OF INCENTER AND CENTROID"},
{"X(87)", Hold@{1/(a b+a c-b c), 1/(b c+b a-c a),  1/(c a+c b-a b)}, "X(2)-CROSS CONJUGATE OF X(1)"},
  (* drawing billiard *)
{"X(88)", Hold@{1/(b+c-2 a), 1/(c+a-2 b), 1/(a+b-2 c)}, "ISOGONAL CONJUGATE OF X(44)"}, 
{"X(89)", Hold@{1/(2 b+2 c-a), 1/(2 c+2 a-b), 1/(2 a+2 b-c)}, "ISOGONAL CONJUGATE OF X(45)"},
{"X(90)", Hold@{1/(cosB+cosC-cosA), 1/(cosC+cosA-cosB), 1/(cosA+cosB-cosC)}, "X(3)-CROSS CONJUGATE OF X(1)"},
{"X(91)", Hold@{sec2A, sec2B, sec2C}, "ISOGONAL CONJUGATE OF X(47)"},
{"X(92)", Hold@{csc2A, csc2B, csc2C}, "CEVAPOINT OF INCENTER AND CLAWSON POINT"},
{"X(93)", Hold@{sec3A, sec3B, sec3C}, "ISOGONAL CONJUGATE OF X(49)"},
{"X(94)", Hold@{csc3A, csc3B, csc3C}, "ISOGONAL CONJUGATE OF X(50)"},
{"X(95)", Hold@{b2 c2 /(cosB cosC+sinB sinC), c2 a2 /(cosC cosA+sinC sinA), a2 b2 /(cosA cosB+sinA sinB)}, "CEVAPOINT OF CENTROID AND CIRCUMCENTER"},
{"X(96)", Hold@{sec2A /(cosB cosC+sinB sinC), sec2B /(cosC cosA+sinC sinA), sec2C /(cosA cosB+sinA sinB)}, "ISOGONAL CONJUGATE OF X(52)"},
{"X(97)", Hold@{cotA/(cosB cosC+sinB sinC),cotB/(cosC cosA+sinC sinA),cotC/(cosA cosB+sinA sinB)}, "ISOGONAL CONJUGATE OF X(53)"}, 
{"X(98)", Hold@{b c/(b4+c4-a2 b2-a2 c2), c a/(c4+a4-b2 c2-b2 a2), a b/(a4+b4-c2 a2-c2 b2)}, "TARRY POINT"},
{"X(99)", Hold@{(b c)/(b2-c2), (c a)/(c2-a2), (a b)/(a2-b2)}, "STEINER POINT"}, 
{"X(100)", Hold@{1/(b-c), 1/(c-a), 1/(a-b)}, "ANTICOMPLEMENT OF FEUERBACH POINT"},
{ "X(101)", Hold@{a/(b-c),b/(c-a),c/(a-b)},"\[CapitalLambda](INCENTER,ORTHOCENTER)\[CapitalPsi](INCENTER,SYMMEDIAN POINT)" },
{ "X(102)", Hold@{a/(2 a4-(b+c)a3-(b-c)^2 a2 (b-c)^2 (b+c)a-(b2-c2)^2),b/(2 b4-(c+a)b3-(c-a)^2 b2 (c-a)^2 (c+a)b-(c2-a2)^2),c/(2 c4-(a+b)c3-(a-b)^2 c2 (a-b)^2 (a+b)c-(a2-b2)^2)},"\[CapitalLambda](INCENTER,ORTHOCENTER)" },
{ "X(103)", Hold@{a/(2 a3-a2(b+c)-(b-c)^2 (b+c)),b/(2 b3-b2(c+a)-(c-a)^2 (c+a)),c/(2 c3-c2(a+b)-(a-b)^2 (a+b))},"ANTIPODE OF X(101)" },
{ "X(104)", Hold@{1/(b3+c3-(a2+b c)(b+c)+2 a b c),1/(c3+a3-(b2+c a)(c+a)+2 b c a),1/(a3+b3-(c2+a b)(a+b)+2 c a b)},"ANTIPODE OF X(100)" },
{ "X(105)", Hold@{1/(b2+c2-a(b+c)),1/(c2+a2-b(c+a)),1/(a2+b2-c(a+b))},"\[CapitalLambda](INCENTER,SYMMEDIAN POINT)" },
{ "X(106)", Hold@{a/(2 a-b-c),b/(2 b-c-a),c/(2 c-a-b)},"\[CapitalLambda](INCENTER,CENTROID)" },
{ "X(107)", Hold@{b c/((b2-c2)(b2+c2-a2)^2),c a/((c2-a2)(c2+a2-b2)^2),a b/((a2-b2)(a2+b2-c2)^2)},"\[CapitalPsi](SYMMEDIAN POINT,ORTHOCENTER)" },
{ "X(108)", Hold@{1/((b-c)(b+c-a)(b2+c2-a2)),1/((c-a)(c+a-b)(c2+a2-b2)),1/((a-b)(a+b-c)(a2+b2-c2))},"\[CapitalPsi](CIRCUMCENTER,INCENTER)" },
{ "X(109)", Hold@{a/((b-c)(b+c-a)),b/((c-a)(c+a-b)),c/((a-b)(a+b-c))},"\[CapitalPsi](INCENTER,CIRCUMCENTER)" },
{ "X(110)", Hold@{a/(b2-c2),b/(c2-a2),c/(a2-b2)},"FOCUS OF KIEPERT PARABOLA" },
{ "X(111)", Hold@{a/(2 a2-b2-c2),b/(2 b2-c2-a2),c/(2 c2-a2-b2)},"PARRY POINT" },
{ "X(112)", Hold@{a/((b2-c2)(b2+c2-a2)),b/((c2-a2)(c2+a2-b2)),c/((a2-b2)(a2+b2-c2))},"\[CapitalPsi](ORTHOCENTER,SYMMEDIAN POINT)" },
{"X(113)", Hold@{(sinB sinC)(sinC/(cosC-2 cosA cosB)+sinB/(cosB-2 cosA cosC)),(sinC sinA)(sinA/(cosA-2 cosB cosC)+sinC/(cosC-2 cosB cosA)),(sinA sinB)(sinB/(cosB-2 cosC cosA)+sinA/(cosA-2 cosC cosB))},"JERABEK ANTIPODE" },
{ "X(114)", Hold@{(b4+c4-a2 b2-a2 c2)(2 a4+b4+c4-a2 b2-a2 c2-2 b2 c2)/a,(c4+a4-b2 c2-b2 a2)(2 b4+c4+a4-b2 c2-b2 a2-2 c2 a2)/b,(a4+b4-c2 a2-c2 b2)(2 c4+a4+b4-c2 a2-c2 b2-2 a2 b2)/c},"KIEPERT ANTIPODE" },
{ "X(115)", Hold@{b c (b2-c2)^2,c a (c2-a2)^2,a b (a2-b2)^2},"CENTER OF KIEPERT HYPERBOLA" },
{ "X(116)", Hold@{b c((b-c)^2(b2+b c+c2-a b-a c)),c a((c-a)^2(c2+c a+a2-b c-b a)),a b((a-b)^2(a2+a b+b2-c a-c b))},"MIDPOINT OF X(4) AND X(103)" },
{ "X(117)", Hold@{(b2 c)/(c(secB-secC)+a(secB-secA))+(c2 b)/(b(secC-secB)+a(secC-secA)),(c2 a)/(a(secC-secA)+b(secC-secB))+(a2 c)/(c(secA-secC)+b(secA-secB)),(a2 b)/(b(secA-secB)+c(secA-secC))+(b2 a)/(a(secB-secA)+c(secB-secC))},"MIDPOINT OF X(4) AND X(109)" },
{ "X(118)", Hold@{(b3 c)/((b-c)/tanA+(b-a)/tanC)+(c3 b)/((c-b)/tanA+(c-a)/tanB),(c3 a)/((c-a)/tanB+(c-b)/tanA)+(a3 c)/((a-c)/tanB+(a-b)/tanC),(a3 b)/((a-b)/tanC+(a-c)/tanB)+(b3 a)/((b-a)/tanC+(b-c)/tanA)},"MIDPOINT OF X(4) AND X(101)" },
{ "X(119)", Hold@{cscA(cosB+cosC-1)(sin2B+sin2C+2(cosA-1)(sinB+sinC)),cscB(cosC+cosA-1)(sin2C+sin2A+2(cosB-1)(sinC+sinA)),cscC(cosA+cosB-1)(sin2A+sin2B+2(cosC-1)(sinA+sinB))},"X(119) = FEUERBACH ANTIPODE" },
{ "X(120)", Hold@{b c(2 a b c-(b+c)(a2+(b-c)^2))(b2+c2-a b-a c),c a(2 b c a-(c+a)(b2+(c-a)^2))(c2+a2-b c-b a),a b(2 c a b-(a+b)(c2+(a-b)^2))(a2+b2-c a-c b)},"X(105)-OF-MEDIAL-TRIANGLE" },
{ "X(121)", Hold@{b c(b+c-2 a)(b3+c3+a(b2+c2)-2 b c(b+c)),c a(c+a-2 b)(c3+a3+b(c2+a2)-2 c a(c+a)),a b(a+b-2 c)(a3+b3+c(a2+b2)-2 a b(a+b))},"X(106)-OF-MEDIAL-TRIANGLE" },
{ "X(122)", Hold@{(b2-c2)^2 (cosA-cosB cosC)/(tanA*tanA),(c2-a2)^2 (cosB-cosC cosA)/(tanB*tanB),(a2-b2)^2 (cosC-cosA cosB)/(tanC*tanC)},"X(107)-OF-MEDIAL-TRIANGLE" },
{ "X(123)", Hold@{ cscA(secB-secC)(secA(sin2B-sin2C)+sinC tanC-sinB tanB), cscB(secC-secA)(secB(sin2C-sin2A)+sinA tanA-sinC tanC), cscC(secA-secB)(secC(sin2A-sin2B)+sinB tanB-sinA tanA)},"X(108)-OF-MEDIAL-TRIANGLE" },
{ "X(124)", Hold@{b c(b+c-a)(b-c)^2 ((b+c)(b2+c2-a2-b c)+a b c),c a(c+a-b)(c-a)^2 ((c+a)(c2+a2-b2-c a)+b c a),a b(a+b-c)(a-b)^2 ((a+b)(a2+b2-c2-a b)+c a b)},"X(109)-OF-MEDIAL-TRIANGLE" },
{ "X(125)", Hold@{b c(b2+c2-a2)(b2-c2)^2,c a(c2+a2-b2)(c2-a2)^2,a b(a2+b2-c2)(a2-b2)^2},"CENTER OF JERABEK HYPERBOLA" },
{ "X(126)", Hold@{b c(2 a2-b2-c2)(b4+c4+a2(b2+c2)-4 b2 c2),c a(2 b2-c2-a2)(c4+a4+b2(c2+a2)-4 c2 a2),a b(2 c2-a2-b2)(a4+b4+c2(a2+b2)-4 a2 b2)},"X(111)-OF-MEDIAL-TRIANGLE" },
{ "X(127)", Hold@{b c(sin2B-sin2C)((b2-c2)sin2A-b2 sin2B+c2 sin2C),c a(sin2C-sin2A)((c2-a2)sin2B-c2 sin2C+a2 sin2A),a b(sin2A-sin2B)((a2-b2)sin2C-a2 sin2A+b2 sin2B)},"X(112)-OF-MEDIAL-TRIANGLE" },
{ "X(128)", Hold@{secA(cos2B+cos2C)(1+2 cos2A)(cos2A+2 cos2B cos2C),secB(cos2C+cos2A)(1+2 cos2B)(cos2B+2 cos2C cos2A),secC(cos2A+cos2B)(1+2 cos2C)(cos2C+2 cos2A cos2B)},"X(74)-OF-ORTHIC-TRIANGLE" },
{ "X(129)", Hold@{secA sin2A(sin2B+sin2C)(sin2B^4+ sin2C^4-(sin2A sin2B)^2-(sin2A sin2C)^2)(sin2A^4+sin2A^2(sin2B sin2C-sin2B^2-sin2C^2)+(sin2B sin2C)(sin2B-sin2C)^2),secB sin2B(sin2C+sin2A)(sin2C^4+ sin2A^4-(sin2B sin2C)^2-(sin2B sin2A)^2)(sin2B^4+sin2B^2(sin2C sin2A-sin2C^2-sin2A^2)+(sin2C sin2A)(sin2C-sin2A)^2),secC sin2C(sin2A+sin2B)(sin2A^4+ sin2B^4-(sin2C sin2A)^2-(sin2C sin2B)^2)(sin2C^4+sin2C^2(sin2A sin2B-sin2A^2-sin2B^2)+(sin2A sin2B)(sin2A-sin2B)^2)},"X(98)-OF-ORTHIC-TRIANGLE" },
{ "X(130)", Hold@{sinA(sin2B+sin2C)(sin2B-sin2C)^2 (sin2A^2+sin2B sin2C),sinB(sin2C+sin2A)(sin2C-sin2A)^2 (sin2B^2+sin2C sin2A),sinC(sin2A+sin2B)(sin2A-sin2B)^2 (sin2C^2+sin2A sin2B)},"X(99)-OF-ORTHIC-TRIANGLE" },
{ "X(131)", Hold@{secA(2 sumT2-sumS2(sec2B+sec2C))(sumT2-S sec2A),secB(2 sumT2-sumS2(sec2C+sec2A))(sumT2-S sec2B),secC(2 sumT2-sumS2(sec2A+sec2B))(sumT2-S sec2C)},"INTERSECTION OF LINES X(3)X(125) AND X(4)X(135)" },
{ "X(132)", Hold@{secA (sin2A^2+(sin2B-sin2C)^2+sin2A(sin2A-sin2B-sin2C)) (sin2B^2+sin2C^2-sin2A sin2B-sin2A sin2C),secB (sin2B^2+(sin2C-sin2A)^2+sin2B(sin2B-sin2C-sin2A)) (sin2C^2+sin2A^2-sin2B sin2C-sin2B sin2A),secC (sin2C^2+(sin2A-sin2B)^2+sin2C(sin2C-sin2A-sin2B)) (sin2A^2+sin2B^2-sin2C sin2A-sin2C sin2B)},"X(2)X(107)\:2229X(4)X(32)" },
{ "X(133)", Hold@{secA((sin2B-sin2C)^2+sin2A sin2B+sin2A sin2C-2 sin2B sin2C)(2 sin2A-sin2B-sin2C),secB((sin2C-sin2A)^2+sin2B sin2C+sin2B sin2A-2 sin2C sin2A)(2 sin2B-sin2C-sin2A),secC((sin2A-sin2B)^2+sin2C sin2A+sin2C sin2B-2 sin2A sin2B)(2 sin2C-sin2A-sin2B)},"INTERSECTION OF LINES X(4)X(74) AND X(5)X(122)" },
{ "X(134)", Hold@{secA sin2A(sin2B^2-sin2C^2)(sin2B^2+sin2C^2-sin2A^2)^2 (sin2C(sin2A^2-sin2B^2)(sin2A^2+sin2B^2-sin2C^2)^2 - sin2B(sin2A^2-sin2C^2)(sin2A^2+ sin2C^2-sin2B^2)^2),secB sin2B(sin2C^2-sin2A^2)(sin2C^2+sin2A^2-sin2B^2)^2 (sin2A(sin2B^2-sin2C^2)(sin2B^2+sin2C^2-sin2A^2)^2 - sin2C(sin2B^2-sin2A^2)(sin2B^2+ sin2A^2-sin2C^2)^2),secC sin2C(sin2A^2-sin2B^2)(sin2A^2+sin2B^2-sin2C^2)^2 (sin2B(sin2C^2-sin2A^2)(sin2C^2+sin2A^2-sin2B^2)^2 - sin2A(sin2C^2-sin2B^2)(sin2C^2+ sin2B^2-sin2A^2)^2)},"X(107)-OF-ORTHIC-TRIANGLE" },
{ "X(135)", Hold@{secA(sin2B/(sec2C-sec2A)+sin2C/(sec2A-sec2B)),secB(sin2C/(sec2A-sec2B)+sin2A/(sec2B-sec2C)),secC(sin2A/(sec2B-sec2C)+sin2B/(sec2C-sec2A))},"INTERSECTION OF LINE X(4)X(131) AND X(25)X(114)" },
{ "X(136)", Hold@{((b2-c2)^2 (a4+b4+c4-2 a2 b2-2 a2 c2)/(a2-b2-c2))/a,((c2-a2)^2 (b4+c4+a4-2 b2 c2-2 b2 a2)/(b2-c2-a2))/b,((a2-b2)^2 (c4+a4+b4-2 c2 a2-2 c2 b2)/(c2-a2-b2))/c},"INTERSECTION OF LINE X(4)X(110) AND X(25)X(132)" },
{ "X(137)", Hold@{secA(sin2B+sin2C)(sin2B-sin2C)^2 (sin2A^2-sin2B^2-sin2C^2-sin2B sin2C),secB(sin2C+sin2A)(sin2C-sin2A)^2 (sin2B^2-sin2C^2-sin2A^2-sin2C sin2A),secC(sin2A+sin2B)(sin2A-sin2B)^2 (sin2C^2-sin2A^2-sin2B^2-sin2A sin2B)},"X(110)-OF-ORTHIC-TRIANGLE" },
{ "X(138)", Hold@{secA(sin2B/(2 sin2B^2-sin2C^2-sin2A^2)+sin2C/(2 sin2C^2-sin2A^2-sin2B^2)),secB(sin2C/(2 sin2C^2-sin2A^2-sin2B^2)+sin2A/(2 sin2A^2-sin2B^2-sin2C^2)),secC(sin2A/(2 sin2A^2-sin2B^2-sin2C^2)+sin2B/(2 sin2B^2-sin2C^2-sin2A^2))},"X(111)-OF-ORTHIC-TRIANGLE" },
{ "X(139)", Hold@{secA(sin2B+sin2C)(sin2B-sin2C)^2(sin2B^2+sin2C^2- sin2A^2)(sin2B^4+sin2C^4-sin2A^4+(sin2B sin2C)(sin2B^2+sin2C^2-sin2A^2)),secB(sin2C+sin2A)(sin2C-sin2A)^2(sin2C^2+sin2A^2- sin2B^2)(sin2C^4+sin2A^4-sin2B^4+(sin2C sin2A)(sin2C^2+sin2A^2-sin2B^2)),secC(sin2A+sin2B)(sin2A-sin2B)^2(sin2A^2+sin2B^2- sin2C^2)(sin2A^4+sin2B^4-sin2C^4+(sin2A sin2B)(sin2A^2+sin2B^2-sin2C^2))},"X(112)-OF-ORTHIC-TRIANGLE" },
{ "X(140)", Hold@{2 a3-3 a(b2+c2)+(b2-c2)^2/a,2 b3-3 b(c2+a2)+(c2-a2)^2/b,2 c3-3 c(a2+b2)+(a2-b2)^2/c},"MIDPOINT OF X(3) AND X(5)" },
{ "X(141)", Hold@{ b c(b2+c2), c a(c2+a2), a b(a2+b2)},"COMPLEMENT OF SYMMEDIAN POINT" },
{ "X(142)", Hold@{ b+c-((b-c)^2)/a, c+a-((c-a)^2)/b, a+b-((a-b)^2)/c},"COMPLEMENT OF X(9)" },
{ "X(143)", Hold@{ a (a2 b2+a2 c2+2 b2 c2-b4-c4)(a4+b4+c4-2 a2 b2-2 a2 c2-b2 c2), b (b2 c2+b2 a2+2 c2 a2-c4-a4)(b4+c4+a4-2 b2 c2-2 b2 a2-c2 a2), c (c2 a2+c2 b2+2 a2 b2-a4-b4)(c4+a4+b4-2 c2 a2-2 c2 b2-a2 b2)},"NINE-POINT CENTER OF ORTHIC TRIANGLE" },
{ "X(144)", Hold@{3 a-2(b+c)-((b-c)^2)/a,3 b-2(c+a)-((c-a)^2)/b,3 c-2(a+b)-((a-b)^2)/c},"ANTICOMPLEMENT OF X(7)" },
{ "X(145)", Hold@{ b c(3 a- b-c), c a(3 b- c-a), a b(3 c- a-b)},"ANTICOMPLEMENT OF NAGEL POINT" },
{ "X(146)", Hold@{a^9+a^7(b2+c2)-a^5(8 b4-9 b2 c2+8 c4)+2 a3 (b2+c2)(4 b4-7 b2 c2+4 c4)-a (b2-c2)^2 (b4+9 b2 c2+c4)-(b2-c2)^4 (b2+c2)/a,b^9+b^7(c2+a2)-b^5(8 c4-9 c2 a2+8 a4)+2 b3 (c2+a2)(4 c4-7 c2 a2+4 a4)-b (c2-a2)^2 (c4+9 c2 a2+a4)-(c2-a2)^4 (c2+a2)/b,c^9+c^7(a2+b2)-c^5(8 a4-9 a2 b2+8 b4)+2 c3 (a2+b2)(4 a4-7 a2 b2+4 b4)-c (a2-b2)^2 (a4+9 a2 b2+b4)-(a2-b2)^4 (a2+b2)/c},"REFLECTION OF X(20) IN X(110)" },
{ "X(147)", Hold@{b c(a8+(b2+c2)a6-(2 b4+3 b2 c2+2 c4)a4+(b6+b4 c2+b2 c4+c6)a2-b8+b6 c2+b2 c6-c8),c a(b8+(c2+a2)b6-(2 c4+3 c2 a2+2 a4)b4+(c6+c4 a2+c2 a4+a6)b2-c8+c6 a2+c2 a6-a8),a b(c8+(a2+b2)c6-(2 a4+3 a2 b2+2 b4)c4+(a6+a4 b2+a2 b4+b6)c2-a8+a6 b2+a2 b6-b8)},"TARRY POINT OF ANTICOMPLEMENTARY TRIANGLE" },
{ "X(148)", Hold@{b c(a4-(b2-c2)^2+b2 c2-a2 b2-a2 c2),c a(b4-(c2-a2)^2+c2 a2-b2 c2-b2 a2),a b(c4-(a2-b2)^2+a2 b2-c2 a2-c2 b2)},"STEINER POINT OF ANTICOMPLEMENTARY TRIANGLE" },
{ "X(149)", Hold@{b c(b3+c3-a3+(a2-b c)(b+c)+a(b c-b2-c2)),c a(c3+a3-b3+(b2-c a)(c+a)+b(c a-c2-a2)),a b(a3+b3-c3+(c2-a b)(a+b)+c(a b-a2-b2))},"REFLECTION OF X(20) IN X(104)" },
{ "X(150)", Hold@{b c(b4+c4-a4+a(b c2+c b2-b3-c3)-b c(a2+b2+c2)+(b+c)a3),c a(c4+a4-b4+b(c a2+a c2-c3-a3)-c a(b2+c2+a2)+(c+a)b3),a b(a4+b4-c4+c(a b2+b a2-a3-b3)-a b(c2+a2+b2)+(a+b)c3)},"REFLECTION OF X(20) IN X(103)" },
{ "X(151)", Hold@{(1/a)(a+b+c)^2 (a^10-a^9(b+c)+12 a^5 b (b-c)^2 c(b+c)+a (b-c)^6 (b+c)^3+2 a^7(b+c)(b2-4 b c+c2)-(b2-c2)^4 (b2-b c+c2)+a^8 (b2+3 b c+c2)-2 a3 (b-c)^4 (b+c)(b2+4 b c+c2)+4 a4 (b2-c2)^2 (2 b2-3 b c+2 c2)+2 a6(-4 b4+b3 c+8 b2 c2+b c3-4 c4)-a2 (b2-c2)^2 (b4-6 b3 c+14 b2 c2-6 b c3+c4)),(1/b)(b+c+a)^2 (b^10-b^9(c+a)+12 b^5 c (c-a)^2 a(c+a)+b (c-a)^6 (c+a)^3+2 b^7(c+a)(c2-4 c a+a2)-(c2-a2)^4 (c2-c a+a2)+b^8 (c2+3 c a+a2)-2 b3 (c-a)^4 (c+a)(c2+4 c a+a2)+4 b4 (c2-a2)^2 (2 c2-3 c a+2 a2)+2 b6(-4 c4+c3 a+8 c2 a2+c a3-4 a4)-b2 (c2-a2)^2 (c4-6 c3 a+14 c2 a2-6 c a3+a4)),(1/c)(c+a+b)^2 (c^10-c^9(a+b)+12 c^5 a (a-b)^2 b(a+b)+c (a-b)^6 (a+b)^3+2 c^7(a+b)(a2-4 a b+b2)-(a2-b2)^4 (a2-a b+b2)+c^8 (a2+3 a b+b2)-2 c3 (a-b)^4 (a+b)(a2+4 a b+b2)+4 c4 (a2-b2)^2 (2 a2-3 a b+2 b2)+2 c6(-4 a4+a3 b+8 a2 b2+a b3-4 b4)-c2 (a2-b2)^2 (a4-6 a3 b+14 a2 b2-6 a b3+b4))},"REFLECTION OF X(20) IN X(109)" },
{ "X(152)", Hold@{(1/a)(a8-a^7 (b+c)+a (b-c)^4 (b+c)^3-a4 b c(b2-6 b c+c2)-(b-c)^4 (b+c)^2 (b2+b c+c2)+a^6 (2 b2+b c+2 c2)-a5(b+c)(5 b2-6 b c+5 c2)+a3 (b-c)^2(b+c)(5 b2+6 b c+5 c2)-a2 (b-c)^2 (2 b4+5 b3 c+10 b2 c2+5 b c3+2 c4)),(1/b)(b8-b^7 (c+a)+b (c-a)^4 (c+a)^3-b4 c a(c2-6 c a+a2)-(c-a)^4 (c+a)^2 (c2+c a+a2)+b^6 (2 c2+c a+2 a2)-b5(c+a)(5 c2-6 c a+5 a2)+b3 (c-a)^2(c+a)(5 c2+6 c a+5 a2)-b2 (c-a)^2 (2 c4+5 c3 a+10 c2 a2+5 c a3+2 a4)),(1/c)(c8-c^7 (a+b)+c (a-b)^4 (a+b)^3-c4 a b(a2-6 a b+b2)-(a-b)^4 (a+b)^2 (a2+a b+b2)+c^6 (2 a2+a b+2 b2)-c5(a+b)(5 a2-6 a b+5 b2)+c3 (a-b)^2(a+b)(5 a2+6 a b+5 b2)-c2 (a-b)^2 (2 a4+5 a3 b+10 a2 b2+5 a b3+2 b4))},"REFLECTION OF X(20) IN X(101)" },
{ "X(153)", Hold@{(1/a)(a^7-a6(b+c)-(b-c)^4(b+c)^3-a5(b2-7 b c+c2)+a4(b+c) (b2-6 b c+c2)+a (b2-c2)^2(b2-5 b c+c2)+a2 (b-c)^2(b+c)(b2+6 b c+c2)-a3(b4+2 b3 c-10 b2 c2+2 b c3+c4)),(1/b)(b^7-b6(c+a)-(c-a)^4(c+a)^3-b5(c2-7 c a+a2)+b4(c+a) (c2-6 c a+a2)+b (c2-a2)^2(c2-5 c a+a2)+b2 (c-a)^2(c+a)(c2+6 c a+a2)-b3(c4+2 c3 a-10 c2 a2+2 c a3+a4)),(1/c)(c^7-c6(a+b)-(a-b)^4(a+b)^3-c5(a2-7 a b+b2)+c4(a+b) (a2-6 a b+b2)+c (a2-b2)^2(a2-5 a b+b2)+c2 (a-b)^2(a+b)(a2+6 a b+b2)-c3(a4+2 a3 b-10 a2 b2+2 a b3+b4))},"REFLECTION OF X(20) IN X(100)" },
{ "X(154)", Hold@{a(tanB+tanC-tanA),b(tanC+tanA-tanB),c(tanA+tanB-tanC)},"X(3)-CEVA CONJUGATE OF X(6)" },
{ "X(155)", Hold@{ cosA(cosB^2+cosC^2-cosA^2), cosB(cosC^2+cosA^2-cosB^2), cosC(cosA^2+cosB^2-cosC^2)},"EIGENCENTER OF ORTHIC TRIANGLE" },
{ "X(156)", Hold@{a(a2 (a2-b2)^3+(-3 a6+2 a4 b2+b6)c2+(3 a4-2 b4)c4+(-a2+b2)c6),b(b2 (b2-c2)^3+(-3 b6+2 b4 c2+c6)a2+(3 b4-2 c4)a4+(-b2+c2)a6),c(c2 (c2-a2)^3+(-3 c6+2 c4 a2+a6)b2+(3 c4-2 a4)b4+(-c2+a2)b6)},"X(5)-OF-TANGENTIAL-TRIANGLE" },
{ "X(157)", Hold@{a(b3 cosB+c3 cosC-a3 cosA),b(c3 cosC+a3 cosA-b3 cosB),c(a3 cosA+b3 cosB-c3 cosC)},"X(6)-OF-TANGENTIAL-TRIANGLE" },
{ "X(158)", Hold@{secA^2,secB^2,secC^2},"X(19)-CROSS CONJUGATE OF X(92)" },
{ "X(159)", Hold@{a((a2+b2+c2)sin2A+(c2-b2-a2)sin2B+(b2-c2-a2)sin2C),b((b2+c2+a2)sin2B+(a2-c2-b2)sin2C+(c2-a2-b2)sin2A),c((c2+a2+b2)sin2C+(b2-a2-c2)sin2A+(a2-b2-c2)sin2B)},"X(9)-OF-TANGENTIAL-TRIANGLE" },
{ "X(160)", Hold@{a((b2+c2)sin2A+(c2-a2)sin2B+(b2-a2)sin2C),b((c2+a2)sin2B+(a2-b2)sin2C+(c2-b2)sin2A),c((a2+b2)sin2C+(b2-c2)sin2A+(a2-c2)sin2B)},"X(37)-OF-TANGENTIAL-TRIANGLE" },
{ "X(161)", Hold@{a((a2+b2+c2)sin2A^2+(c2-b2-a2)sin2B^2+(b2-c2-a2)sin2C^2),b((b2+c2+a2)sin2B^2+(a2-c2-b2)sin2C^2+(c2-a2-b2)sin2A^2),c((c2+a2+b2)sin2C^2+(b2-a2-c2)sin2A^2+(a2-b2-c2)sin2B^2)},"X(63)-OF-TANGENTIAL-TRIANGLE" },
{ "X(162)", Hold@{1/((b2-c2)(b2+c2-a2)),1/((c2-a2)(c2+a2-b2)),1/((a2-b2)(a2+b2-c2))},"CEVAPOINT OF X(108) AND X(109)" },
{ "X(163)", Hold@{a2/(b2-c2),b2/(c2-a2),c2/(a2-b2)},"TRILINEAR PRODUCT X(6)*X(110)" },
{ "X(164)", Hold@{sinHalfB+sinHalfC-sinHalfA,sinHalfC+sinHalfA-sinHalfB,sinHalfA+sinHalfB-sinHalfC},"INCENTER OF EXCENTRAL TRIANGLE" },
{ "X(165)", Hold@{3 a2-2 a(b+c)-(b-c)^2,3 b2-2 b(c+a)-(c-a)^2,3 c2-2 c(a+b)-(a-b)^2},"CENTROID OF THE EXCENTRAL TRIANGLE" },
{ "X(166)", Hold@{tanHalfA/(cosHalfB+cosHalfC-cosHalfA)-tanHalfB/(cosHalfC+cosHalfA-cosHalfB)-tanHalfC/(cosHalfA+cosHalfB-cosHalfC),tanHalfB/(cosHalfC+cosHalfA-cosHalfB)-tanHalfC/(cosHalfA+cosHalfB-cosHalfC)-tanHalfA/(cosHalfB+cosHalfC-cosHalfA),tanHalfC/(cosHalfA+cosHalfB-cosHalfC)-tanHalfA/(cosHalfB+cosHalfC-cosHalfA)-tanHalfB/(cosHalfC+cosHalfA-cosHalfB)},"GERGONNE POINT OF EXCENTRAL TRIANGLE" },
{ "X(167)", Hold@{sinHalfB (cosHalfC+cosHalfA-cosHalfB)/cosHalfB+sinHalfC (cosHalfA+cosHalfB-cosHalfC)/cosHalfC-sinHalfA (cosHalfB+cosHalfC-cosHalfA)/cosHalfA,sinHalfC (cosHalfA+cosHalfB-cosHalfC)/cosHalfC+sinHalfA (cosHalfB+cosHalfC-cosHalfA)/cosHalfA-sinHalfB (cosHalfC+cosHalfA-cosHalfB)/cosHalfB,sinHalfA (cosHalfB+cosHalfC-cosHalfA)/cosHalfA+sinHalfB (cosHalfC+cosHalfA-cosHalfB)/cosHalfB-sinHalfC (cosHalfA+cosHalfB-cosHalfC)/cosHalfC},"NAGEL POINT OF EXCENTRAL TRIANGLE" },
{ "X(168)", Hold@{b/(1-sinHalfB)+c/(1-sinHalfC)-a/(1-sinHalfA),c/(1-sinHalfC)+a/(1-sinHalfA)-b/(1-sinHalfB),a/(1-sinHalfA)+b/(1-sinHalfB)-c/(1-sinHalfC)},"MITTENPUNKT OF EXCENTRAL TRIANGLE" },
{ "X(169)", Hold@{a3-a2(b+c)+a(b2+c2)-(b-c)^2(b+c),b3-b2(c+a)+b(c2+a2)-(c-a)^2(c+a),c3-c2(a+b)+c(a2+b2)-(a-b)^2(a+b)},"X(85)-CEVA CONJUGATE OF X(1)" },
{ "X(170)", Hold@{-tanHalfA/(cosHalfA^2)+tanHalfB/(cosHalfB^2)+tanHalfC/(cosHalfC^2),-tanHalfB/(cosHalfB^2)+tanHalfC/(cosHalfC^2)+tanHalfA/(cosHalfA^2),-tanHalfC/(cosHalfC^2)+tanHalfA/(cosHalfA^2)+tanHalfB/(cosHalfB^2)},"X(9)-ALEPH CONJUGATE OF X(9)" },
{ "X(171)", Hold@{a2+b c,b2+c a,c2+a b},"{X(2),X(31)}-HARMONIC CONJUGATE OF X(238)" },
{ "X(172)", Hold@{a3+a b c,b3+b c a,c3+c a b},"TRILINEAR PRODUCT X(6)*X(171)" },
{ "X(173)", Hold@{cosHalfB+cosHalfC-cosHalfA,cosHalfC+cosHalfA-cosHalfB,cosHalfA+cosHalfB-cosHalfC},"CONGRUENT ISOSCELIZERS POINT" },
{ "X(174)", Hold@{1/cosHalfA,1/cosHalfB,1/cosHalfC},"YFF CENTER OF CONGRUENCE" },
{ "X(175)", Hold@{(cosHalfB cosHalfC)/cosHalfA-1,(cosHalfC cosHalfA)/cosHalfB-1,(cosHalfA cosHalfB)/cosHalfC-1},"ISOPERIMETRIC POINT" },
{ "X(176)", Hold@{(cosHalfB cosHalfC)/cosHalfA+1,(cosHalfC cosHalfA)/cosHalfB+1,(cosHalfA cosHalfB)/cosHalfC+1},"EQUAL DETOUR POINT" },
{ "X(177)", Hold@{(cosHalfB+cosHalfC)/cosHalfA,(cosHalfC+cosHalfA)/cosHalfB,(cosHalfA+cosHalfB)/cosHalfC},"1st MID-ARC POINT" },
{ "X(178)", Hold@{(cosHalfB+cosHalfC)/a,(cosHalfC+cosHalfA)/b,(cosHalfA+cosHalfB)/c},"2nd MID-ARC POINT" },
{ "X(179)", Hold@{1/(cosQuarterA^4),1/(cosQuarterB^4),1/(cosQuarterC^4)},"1st AJIMA-MALFATTI POINT" },
{ "X(180)", Hold@{1/(1 + 2(cosQuarterC cosQuarterA/cosQuarterB)^2)+1/(1 + 2(cosQuarterB cosQuarterA/cosQuarterC)^2)-1/(1 + 2(cosQuarterC cosQuarterB/cosQuarterA)^2),1/(1 + 2(cosQuarterA cosQuarterB/cosQuarterC)^2)+1/(1 + 2(cosQuarterC cosQuarterB/cosQuarterA)^2)-1/(1 + 2(cosQuarterA cosQuarterC/cosQuarterB)^2),1/(1 + 2(cosQuarterB cosQuarterC/cosQuarterA)^2)+1/(1 + 2(cosQuarterA cosQuarterC/cosQuarterB)^2)-1/(1 + 2(cosQuarterB cosQuarterA/cosQuarterC)^2)},"2nd AJIMA-MALFATTI POINT" },
{ "X(181)", Hold@{a (b+c)^2 /(b+c-a),b (c+a)^2 /(c+a-b),c (a+b)^2 /(a+b-c)},"APOLLONIUS POINT" },
{ "X(182)", Hold@{a(a4-a2 b2-a2 c2-2 b2 c2),b(b4-b2 c2-b2 a2-2 c2 a2),c(c4-c2 a2-c2 b2-2 a2 b2)},"MIDPOINT OF BROCARD DIAMETER" },
{ "X(183)", Hold@{1/tanA+tanOmega,1/tanB+tanOmega,1/tanC+tanOmega},"TRILINEAR PRODUCT X(75)X(182)" },
{ "X(184)", Hold@{a2 cosA,b2 cosB,c2 cosC},"INVERSE OF X(125) IN THE BROCARD CIRCLE" },
{ "X(185)", Hold@{cosA(cosB^2+cosC^2),cosB(cosC^2+cosA^2),cosC(cosA^2+cosB^2)},"NAGEL POINT OF THE ORTHIC TRIANGLE" },
{ "X(186)", Hold@{4 cosA-1/cosA,4 cosB-1/cosB,4 cosC-1/cosC},"INVERSE-IN-CIRCUMCIRCLE OF X(4)" },
{ "X(187)", Hold@{a(2 a2-b2-c2),b(2 b2-c2-a2),c(2 c2-a2-b2)},"INVERSE-IN-CIRCUMCIRCLE OF X(6) (SCHOUTE CENTER)" },
{ "X(188)", Hold@{1/sinHalfA,1/sinHalfB,1/sinHalfC},"2nd MID-ARC POINT OF ANTICOMPLEMENTARY TRIANGLE" },
{ "X(189)", Hold@{(b c)/(cosB+cosC-cosA-1),(c a)/(cosC+cosA-cosB-1),(a b)/(cosA+cosB-cosC-1)},"CYCLOCEVIAN CONJUGATE OF X(8)" },
{ "X(190)", Hold@{b c/(b-c),c a/(c-a),a b/(a-b)},"YFF PARABOLIC POINT" },
{ "X(191)", Hold@{(b+c-a)(b c+c a+a b)+b3+c3-a3,(c+a-b)(c a+a b+b c)+c3+a3-b3,(a+b-c)(a b+b c+c a)+a3+b3-c3},"X(10)-CEVA CONJUGATE OF X(1)" },
{ "X(192)", Hold@{b c(c a+a b-b c),c a(a b+b c-c a),a b(b c+c a-a b)},"X(1)-CEVA CONJUGATE OF X(2)" },
{ "X(193)", Hold@{3 a-(b2+c2)/a,3 b-(c2+a2)/b,3 c-(a2+b2)/c},"X(4)-CEVA CONJUGATE OF X(2)" },
{ "X(194)", Hold@{b c(a2 b2+a2 c2-b2 c2),c a(b2 c2+b2 a2-c2 a2),a b(c2 a2+c2 b2-a2 b2)},"X(6)-CEVA CONJUGATE OF X(2)" },
{ "X(195)", Hold@{a(a8+b8+c8-4 a6(b2+c2)+a4(6 b4+6 c4+5 b2 c2)-a2(4 b6+4 c6-b4 c2-b2 c4)-2 b2 c2(b4+c4-b2 c2)),b(b8+c8+a8-4 b6(c2+a2)+b4(6 c4+6 a4+5 c2 a2)-b2(4 c6+4 a6-c4 a2-c2 a4)-2 c2 a2(c4+a4-c2 a2)),c(c8+a8+b8-4 c6(a2+b2)+c4(6 a4+6 b4+5 a2 b2)-c2(4 a6+4 b6-a4 b2-a2 b4)-2 a2 b2(a4+b4-a2 b2))},"X(5)-CEVA CONJUGATE OF X(3)" },
{ "X(196)", Hold@{(cosB+cosC-cosA-1)tanHalfA/cosA,(cosC+cosA-cosB-1)tanHalfB/cosB,(cosA+cosB-cosC-1)tanHalfC/cosC},"X(7)-CEVA CONJUGATE OF X(4)" },
{ "X(197)", Hold@{a(-a2 tanHalfA+b2 tanHalfB+c2 tanHalfC),b(-b2 tanHalfB+c2 tanHalfC+a2 tanHalfA),c(-c2 tanHalfC+a2 tanHalfA+b2 tanHalfB)},"X(8)-CEVA CONJUGATE OF X(6)" },
{ "X(198)", Hold@{a(a3+a2(b+c)-a (b+c)^2-(b-c)^2 (b+c)),b(b3+b2(c+a)-b (c+a)^2-(c-a)^2 (c+a)),c(c3+c2(a+b)-c (a+b)^2-(a-b)^2 (a+b))},"X(9)-CEVA CONJUGATE OF X(6)" },
{ "X(199)", Hold@{a(b4+c4-a4+(b2+c2-a2)(b c+c a+a b)),b(c4+a4-b4+(c2+a2-b2)(c a+a b+b c)),c(a4+b4-c4+(a2+b2-c2)(a b+b c+c a))},"X(10)-CEVA CONJUGATE OF X(6)" },
{ "X(200)", Hold@{(b+c-a)^2,(c+a-b)^2,(a+b-c)^2},"X(8)-CEVA CONJUGATE OF X(9)" },
{ "X(376)", Hold@{2 cosA-cosB cosC,2 cosB-cosC cosA,2 cosC-cosA cosB},"CENTROID OF THE ANTIPEDAL TRIANGLE OF X(2)" }
}

