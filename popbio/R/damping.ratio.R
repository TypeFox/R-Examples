damping.ratio<-function(A)
{
   ev <- eigen(A)
   ## Use second largest magnitude in case of ties
   ## Also, add round for imprimitive matrices.   For example,
   ## without rounding, Mod(ev$values) of this matrix
   ## A<-matrix(c(0,0,2,.3,0,0,0,.6,0), nrow=3,byrow=TRUE) is   
   ## 0.7113786608980130 0.7113786608980126 and damping.ratio would be 1.
   ## with rounding, only one modulus .711 and damping.ration is NA
   dr<-rle(round(Mod(ev$values), 5  ))$values
   dr<-dr[1]/dr[2]
   dr
}
