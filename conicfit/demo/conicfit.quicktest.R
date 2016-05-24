library(conicfit)

# demo for ellipse fitting
# created by Prof. Chernov 2014

XY <- matrix(c(1,7,2,6,5,8,7,7,9,5,3,7,6,2,8,4),8,2,byrow=TRUE)
#  XY is the coordinates of the data points.
# This is a benchmark example from the journal paper
#       W. Gander, G. H. Golub, and R. Strebel,
#       "Least squares fitting of circles and ellipses"
#       BIT, volume 34, (1994), pages 558-578
  
ParGini <- matrix(c(0,0,2,1,0),ncol=1)

# Parameters of the initial ellipse are chosen arbitrarily:
#   center (0,0),  axes 2 and 1,  angle of tilt 0

LambdaIni <- 0.1

#  initial value for the control parameter "lambda" set to 0.1

tmp <- fit.ellipseLMG(XY,ParGini,LambdaIni)
print(tmp)
# ParG =
#     2.6996
#     3.8160
#     6.5187
#     3.0319
#     0.3596
# RSS = 1.3733
# iters = 15
# TF = 0

#  geometric ellipse fit returns parameters of the best fitting ellipse
#    correct answers:  
#  center at (2.6996,3.180), axes (6.5187,3.0319), angle of tilt 0.3596

ParAini <- GtoA(ParGini)
#    A B C D E  F
# 0.25 0 1 0 0 -1

#  algebraic parameters of the initial ellipse

tmp <- fit.conicLMA(XY,c(ParAini),LambdaIni)
ParA <- tmp$ParA
RSS <- tmp$RSS
iters <- tmp$iters
print(tmp)
# ParA =
#     0.0551
#    -0.0908
#     0.1588
#     0.0489
#    -0.9669
#     0.1620
# RSS = 1.3733
# iters = 18
# code = 1

#  geometric conic fit using algebraic parameters 
#  returns algebraic parameters of the best fitting conic

tmp <- AtoG(ParA)
print(tmp)
# ParG =
#     2.6996
#     3.8160
#     6.5187
#     3.0319
#     0.3596
# code = 1

#  convert the algebraic to geometric parameters
RSS
iters
# 1.373306
# 18