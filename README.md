# math-trilins-to-js

R code to convert trilinear file produced with mathematica into js

Usage: at linux terminal type:

$ Rscript p5js.R 'fname_math_in.cform' [barys]

Examples

1) Rscript p5js.R 'data/x0001_0200 v2b.cform'
2) Rscript p5js.R 'data_moses/x0001_1000 v1.cform' barys

Output: 'js/fname_math_in.js' (should be renamed trilins_final.js' when moving it to the app)
