
bceps = 0.0001

logLbc = function(bc, y, x, 
    logy=log(y), 
    xvx = crossprod(x),     
    xvxinv = solve(xvx),
    sumlogy = sum(logy,na.rm=TRUE)
){
  
  if(abs(bc)< bceps){
    ybc = logy
  } else if (abs(bc-1)>bceps){
    ybc = (exp(bc*logy)-1)/bc
  } else {
    ybc = y
  }
  twologj = 2*(sumlogy)*(bc-1)
  
  
  betaHat = xvxinv %*% crossprod(x, ybc)
  ssq = ybc - x %*% betaHat
  
  length(y)*log(sum(ssq^2,na.rm=TRUE)) - twologj
  
}

optBoxCox = function(y,x,boxcoxSeq){
  y = as.matrix(y)
  if(ncol(y)>1) warning('only y should have one column when using boxcox')
  logy = log(y)
  xvx = crossprod(x)

  boxcox = optimize(
    logLbc, interval=c(-1.5, 2.5),
    y=y, x=x,
    logy = logy, xvx = xvx, 
    xvxinv=solve(xvx),
    sumlogy = sum(logy,na.rm=TRUE)
)$min

Nboxcox = ceiling(boxcoxSeq['len']-1)/2
Sboxcox = sort(unique(round((boxcox + boxcoxSeq)/(100*bceps)))*100*bceps)
Sboxcox
}


forBoxCox = function(y,x,seqBoxCox){

Sboxcox = optBoxCox(y,x,seqBoxCox)

theOnes = abs(Sboxcox-1) < bceps
theZeros = abs(Sboxcox) < bceps

Ymat = matrix(NA, nrow(y), length(Sboxcox))
colnames(Ymat) = as.character(Sboxcox)

logy = log(y)
for(D in Sboxcox)
  Ymat[,as.character(D)] = (exp(D*logy)-1)/D
colnames(Ymat) = paste("bc",colnames(Ymat),sep='')

Ymat[,theZeros] = logy
Ymat[,theOnes] = y

sumlogy = sum(logy,na.rm=TRUE)

jacobian = 2*(sumlogy)*(
      Sboxcox -1
      )


return(list(
        y = Ymat,
        boxcox=Sboxcox,
        jacobian=jacobian
        ))
}