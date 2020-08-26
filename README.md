# math-trilins-to-js

R code to convert trilinear file produced with mathematica into js

Usage: at linux terminal type:

$ Rscript p5js.R 'fname_math_in.cform'

e.g. Rscript psjs.R 'x0001_0200 v2b.cform'

1) creates intermediate file 'fname_math_in.csv'
2) using the csv creates 'js/fname_math_in.js'