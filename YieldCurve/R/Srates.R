`Srates`<- function( Coeff, maturity, whichRate="Forward" )
  {
    if(ncol(Coeff)==1) Coeff<-matrix(as.vector(Coeff),1,nrow(Coeff))
    Curve <- xts(matrix( 0, nrow(Coeff), length(maturity) ), order.by=time(Coeff))
    colnames(Curve) <- make.names(maturity)
    Coeff <- as.matrix(Coeff)
    
    switch(whichRate,
      Forward =
      {
        CurveForward <- Curve
        for(i in 1:nrow(Coeff))
          {
            CurveForward[i,] <- Coeff[i,1]+
              Coeff[i,2] * .beta1Forward( maturity, Coeff[i,5] ) +
              Coeff[i,3] * .beta2Forward( maturity, Coeff[i,5] ) +
              Coeff[i,4] * .beta2Forward( maturity, Coeff[i,6] )
          }
          FinalCurve<-CurveForward
      },
      Spot =
      {
         for(i in 1:nrow(Coeff))
          {
            Curve[i,] <- Coeff[i,1] +
              Coeff[i,2] * .beta1Spot( maturity, Coeff[i,5] ) +
              Coeff[i,3] * .beta2Spot( maturity, Coeff[i,5] ) +
              Coeff[i,4] * .beta2Spot( maturity, Coeff[i,6] )
          }
          FinalCurve <- Curve
      })
    reclass( FinalCurve, Coeff )
  }
