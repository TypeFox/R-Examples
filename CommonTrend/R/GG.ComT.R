GG.ComT <- function (data, rank, k)
{ 

  begin= k+1
  max = nrow(data)
  width = ncol(data)
  length=(max-k)

  cotest = ca.jo(data, type = "trace", ecdet = "none", K=k,    spec="transitory")
  equation = cajorls(cotest, r = rank , reg.number = NULL)

  alpha=matrix(NA,width,rank)
  alpha [,1:rank] =t(equation$rlm$coefficients[1:rank,1:width])#### otherwise, when rank is one, the data structure will change 
  oalpha <- qr.Q(qr(alpha),complete=TRUE) [ ,(rank+1):width] #### (width-rank) is the number of common trends

  beta= equation$beta
  obeta <- qr.Q(qr(beta),complete=TRUE)[,(rank+1):width]
  
  loading=obeta %*% ginv (t(oalpha)%*%obeta) #  Kasa's method and GG's method are different here.
  
  Pure.T= matrix(NA,(width-rank),(max-k) )
  Pure.T[1:(width-rank),]= t(oalpha) %*% t (data[begin:max,]) #  Kasa's method and GG's method are different here.
  Com.T= loading%*%Pure.T
  
  station.p = alpha %*% ginv ( t(beta)%*% alpha) %*% t(beta) %*% t (data[begin:max,])
  
  result=list (method="Gonzalo and Grange(1995)",length = length, lag.chosen=k,
               beta = beta, othog.beta = obeta, alpha=alpha,othog.alpha=oalpha,
               common.trend=Com.T,  pure.trend=Pure.T, 
               loading.vector= loading,  stationary=station.p, data.used=t(data[(k+1):max,]) )
  
  class (result) <- "ComT"
  return(result)
}
