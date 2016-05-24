print.GPCE.quad <- function(x, ...) {
  if (!is.null(x$Args)){
    print(paste("Problem definition:"),...)
    if (!is.null(x$Args$InputDim)) {print(paste("Number of random input is",x$Args$InputDim),...)}
    if (!is.null(x$Args$InputDistrib)) {print(paste("Input distribution(s) :",x$Args$InputDistrib),...)}
    if (!is.null(x$Args$ExpPoly)) {print(paste("Basis of expansion :",x$Args$ExpPoly),...)}
    if (!is.null(x$Args$QuadType)) {
      print(paste("Quadrature construction : ",x$Args$QuadType),...)
      if (x$Args$QuadType=='FULL'){
        for (ii in 1:x$Args$InputDim){{print(paste("Gauss-",x$Args$QuadPoly[ii]," quadrature in dimension",ii),...)}}}
      if (x$Args$QuadType=='SPARSE'){
        for (ii in 1:x$Args$InputDim){{print(paste(x$Args$QuadPoly[ii]," sparse-quadrature in dimension",ii),...)}}}
    }
    if (!is.null(x$Design$QuadSize)) {print(paste("Number of quadrature sample is",x$Design$QuadSize),...)}
  }
  if (!is.null(x$Moments)) {
    print("PCE Results:")
    if (!is.null(x$Args$p)) {print(paste("Expansion up to order p =",x$Args$p,"with",getM(x$Args$InputDim,x$Args$p),"terms."),...)}    
    if (!is.null(x$Moments$PCEMean)) {print(paste("PCE mean is",x$Moments$PCEMean),...)}
    if (!is.null(x$Moments$PCEVar))  {print(paste("PCE variance is",x$Moments$PCEVar),...)}  
    if (!is.null(x$Moments$PCESkew)) {print(paste("PCE skewness is",x$Moments$PCESkew),...)}  
    if (!is.null(x$Moments$PCEKurt)) {print(paste("PCE kurtosis is",x$Moments$PCEKurt),...)}  
  } else {stop('There are no PCE moments in the results.')}
      
  if (!is.null(x$Sensitivity)) {
    for (ii in 1:length(x$Sensitivity$Values)){
      if (nchar(x$Sensitivity$Names[ii]) == 1){
        print(paste("Sobol indices for variable",x$Sensitivity$Names[ii],"is", x$Sensitivity$Values[ii]),...)
      }else{
        print(paste("Sobol indices for interaction of variables",x$Sensitivity$Names[ii],"is", x$Sensitivity$Values[ii]),...)        
      }
    }
  }
}
