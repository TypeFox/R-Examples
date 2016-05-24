#Calculates likelihood
#@numb Number of covariates on bb
#@Rfun -- what R function to use to calculate
#@See documentation for description of other parameters

like <- function(param, y, x, n, Zb, Zw, numb, erho, esigma, ebeta,
                 ealphab, ealphaw, Rfun){

                                        #Transform parameters
  Bb0 <- param[1]
  Bw0 <- param[2]
  sb0 <- param[3]
  sw0 <- param[4]
  rho0 <- param[5]
  Bb0v <- param[6:(5+numb)]
  Bw0v <- param[(numb+6):length(param)]
  sb=exp(sb0)
  sw=exp(sw0)
  Zb <- as.matrix(Zb)
  Zw <- as.matrix(Zw)
  bb=Bb0*(.25+sb^2) + .5 +
    as.matrix(apply(Zb,2, function(x) x - mean(x)))%*%as.matrix	(Bb0v)
  bw=Bw0*(.25+sw^2) + .5 +
    as.matrix(apply(Zw,2, function(x) x - mean(x)))%*%as.matrix(Bw0v)
  rho=(exp(2*rho0)-1)/(exp(2*rho0) +1)
  sigb2 <- sb^2
  sigw2 <- sw^2
  sigbw = rho*sb*sw

#print(c(mean(Bb0),mean(Bw0),sb0,sw0,rho0,mean(Bb0v),mean(Bw0v)))
#print(c(mean(bb),mean(bw),sb,sw,rho))

                                        #Create Demographic Categories
  homoindx <- ifelse(x==0, 1, 0)
  homoindx <- ifelse(x==1, 2, homoindx)
  enumtol=.0001
  cT0 <- y<enumtol & homoindx==0
  cT1 <- y>(1-enumtol) & homoindx==0
  ok <- ifelse(homoindx==0 & cT0==0 & cT1==0,T, F)

#Compute likelihood for different categories

#0<T<1, 0<X<1
  omx <- 1-x
  mu = bb*x + bw*omx
  epsilon = y - mu
  s2 = sigb2*(x^2) + sigw2*(omx^2) + 2*sigbw*x*omx
  omega = sigb2*x + sigbw*omx
  ebb = bb + (omega/s2)*epsilon
  vbb = sigb2 - (omega^2)/s2
  vbb = ifelse(vbb<1*10^-32, .0001, vbb)
  bounds <- bounds1(x, y, n)
  s <- ifelse(vbb>=0 & vbb!=Inf & !is.na(vbb),sqrt(vbb),NaN)
  res <- NULL
  b.s = (bounds[ok,][,2]-ebb[ok])/s[ok]
  as = (bounds[ok,][,1]-ebb[ok])/s[ok]
  res[ok] <- log(pnorm(as, lower.tail=F) - pnorm(b.s, lower.tail=F))
#res[ok] <- ifelse(abs(res[ok])==Inf, log(1*10^-15),res[ok])
#res[ok] <- ifelse(abs(res[ok])==Inf, NaN,res[ok])
#res[ok] <- log(pnorm(bounds[ok,2], mean=ebb[ok], sd=s[ok]) -
##pnorm(bounds[ok,1],#mean=ebb[ok], sd=s[ok]))s
#res[ok] <- ifelse(res[ok]==NA,-999,res[ok])
#print(summary(res))
  R <- NULL
  bs <- as.matrix(cbind(bb, bw))
  R[ok] <- .createR(ok,Rfun,bb,bw,sb,sw,rho,x)

#print(summary(R))
	
#lliki  <- -.5*(log(s2[ok])+epsilon[ok]^2/s2[ok]) + res[ok] - R[ok]
  llik.het <- -.5*sum((log(s2[ok])+(epsilon[ok]^2)/(s2[ok])))
  llik.het <- llik.het + sum(res[ok]) - sum(R[ok])

#Homogenous precincts

#X=0
  wh <- homoindx==1
  llik.wh=0
  if (sum(wh)>0){
    epsilon= y[wh]-bw[wh]
    llik.whi = -.5*(log(sigw2)+(epsilon^2)/(sigw2))
    llik.wh=-.5*sum((log(sigw2)+(epsilon^2)/(sigw2)))
    bnds=cbind(rep(0,sum(wh)),rep(1,sum(wh)))
    Ebb = bb[wh]+rho*(sb/sw)*epsilon
    vbb = sigb2*(1-rho^2)
    vbb = ifelse(vbb<1*10^-32, .0001, vbb)
    s <- ifelse(vbb>=0 & vbb!=Inf & !is.na(vbb),sqrt(vbb),NaN)
    b.s = (bnds[,2]-Ebb)/s
    as = (bnds[,1]-Ebb)/s
    res <- log(pnorm(as, lower.tail=F) - pnorm(b.s, lower.tail=F))
    #res <- log(pnorm(bnds[,2], mean=Ebb, sd=s) - pnorm(bnds[,1],
    #mean=Ebb, sd=s))
    R[wh] <- .createR(wh,Rfun,bb,bw,sb,sw,rho,x)
    llik.wh = llik.wh + sum(res)-sum(R[wh])
  }

#X=1
  bl <- homoindx==2
  llik.bl=0
  if (sum(bl)>0){
    epsilon = y[bl] - bb[bl]
    llik.bl = -.5*sum((log(sigb2)+(epsilon^2)/(sigb2)))
    bnds = cbind(rep(0, sum(bl)),rep(1,sum(bl)))
    Ebb=bw[bl] + rho*(sw/sb)*epsilon
    vbb=sigw2*(1-rho^2)
	vbb = ifelse(vbb<1*10^-32, .0001, vbb)
    s <- ifelse(vbb>=0 & vbb!=Inf & !is.na(vbb),sqrt(vbb),NaN)
    b.s = (bnds[,2]-Ebb)/s
    as = (bnds[,1]-Ebb)/s
    res <- log(pnorm(as, lower.tail=F) - pnorm(b.s, lower.tail=F))
    #res <- log(pnorm(bnds[,2], mean=Ebb, sd=s) - pnorm(bnds[,1],
    #mean=Ebb, sd=s)) #res[ok] <- ifelse(abs(res[ok])==Inf,
    #NaN,res[ok])
    R[bl] <- .createR(bl,Rfun,bb,bw,sb,sw,rho,x)
    llik.bl = llik.bl + sum(res)-sum(R[bl])
  }

#T=0, 0<X<1
  llik.cT0=0
  if(sum(cT0)>0){
    bb.cT0 = bs[cT0,][1]
    bw.cT0 = bs[cT0,][2]
    sigma = matrix(c(sigb2,sigbw,sigbw,sigw2),nrow=2)
    if(sum(cT0)==1){
      first = log(dmvnorm(c(0,0),mean=bs[cT0,], sigma=sigma))
      second <- .createR(cT0,Rfun,bb,bw,sb,sw,rho,x)
      llik.cT0=sum(first)-sum(second)
    }
    else{
      first = apply(bs[cT0,], 1,
        function (x) log(dmvnorm(c(0,0),mean=as.vector(x),
                                 sigma=sigma)))
      second <- NULL
      second <- .createR(cT0,Rfun,bb,bw,sb,sw,rho,x)
      llik.cT0=sum(first)-sum(second)
    }
  }

#T=1, 0<X<1
  llik.cT1=0
  if(sum(cT1)>0){
    bb.cT1 = bs[cT1,][1]
    bw.cT1 = bs[cT1,][2]
    sigma=matrix(c(sigb2,sigbw,sigbw,sigw2),nrow=2)
    if(sum(cT1)==1){
      first = log(dmvnorm(c(1,1),mean=bs[cT1,], sigma=sigma))
    #qi <- pmvnorm(lower=c(-bb.cT1/sb,-bw.cT1/sw),
     #upper=c(-bb.cT1/sb + 1/sb, -bw.cT1/sw + 1/sw),
      #mean=c(0,0), #corr=matrix(c(1,rho,rho,1), nrow=2))
    #qi <- ifelse(qi<0 | qi==0, 1*10^-322,qi)
    #second = ifelse(qi<0|qi==0, -999,log(qi))
      second <- .createR(cT1,Rfun,bb,bw,sb,sw,rho,x)
      llik.cT1=sum(first)-sum(second)
    }
    if(sum(cT1)>1){
      first = apply(as.matrix(bs[cT1,]), 1,
        function (x) log(dmvnorm(c(1,1),mean=as.vector(x),
                                 sigma=sigma)))
      second <- NULL
      second <- .createR(cT1,Rfun,bb,bw,sb,sw,rho,x)
#for (i in 1:length(bb.cT1)){
#qi <- pmvnorm(lower=c(-bb.cT1[i]/sb,-bw.cT1[i]/sw), upper=c(-bb.cT1[i]/sb + 1/sb, -bw.cT1[i]/sw + 1/sw), mean=c#(0,0), corr=matrix(c(1,rho,rho,1), nrow=2))
#qi <- ifelse(qi<0 | qi==0, 1*10^-322,qi)
#second[i] = ifelse(qi<0|qi==0, -999,log(qi))
#}
      llik.cT1=sum(first)-sum(second)
    }
  }
  llik=llik.het + llik.bl + llik.wh + llik.cT0 + llik.cT1


#Priors
  prior=0
  lpdfnorm = log(dnorm(rho0,0,sd=erho))
  if (esigma>0) prior = prior-(1/(2*esigma^2))*(sigb2+sigw2);
  if(erho>0) prior = prior +lpdfnorm;
  if(ebeta>0 & (mean(bb)<0)) prior = prior -.5*((mean(bb)^2)/ebeta);
  if(ebeta>0 & mean(bb)>1) prior = prior -.5*((mean(bb)-1)^2/ebeta);
  if(ebeta>0 & mean(bw)<0) prior = prior -.5*((mean(bw)^2)/ebeta);
  if(ebeta>0 & mean(bw)>1) prior = prior -.5*((mean(bw)-1)^2/ebeta);
  if(sum(is.na(ealphab))==0) prior=prior +
    sum(dmvnorm(Bb0v, ealphab[,1], sigma=diag(ealphab[,2]^2), log=T));
  if(sum(is.na(ealphaw))==0) prior=prior +
    sum(dmvnorm(Bw0v, ealphaw[,1], sigma=diag(ealphaw[,2]), log=T));
  llik = llik + prior
                                        #print(-llik)
  if(is.na(llik)|abs(llik)==Inf) llik = NaN
  return(-llik)
}


#Alternative res
#s = sqrt(vbb)
#b.s = (bounds[,2]-ebb)/s
#as = (bounds[,1]-ebb)/s
#log(pnorm(b.s) - pnorm(as))

#for (i in 1:75) {
#res[i] <- log(pmvnorm(lower=c(as[i]), upper=c(b.s[i]), mean=c(0), #sigma=as.matrix(1)))
#}

#for (i in 1:75) {
#res[i] <- log(pmvnorm(lower=c(bounds[,1][i]), upper=c(bounds[,2][i]), #mean=ebb[i], sigma=as.matrix(vbb[i])))
#}
#Alternative R

#R <- ifelse(rho==0,
		#log(pmvnorm(lower=c(0), upper=c(1), mean=bb, #sigma=as.matrix(sigb2))) + log(pmvnorm(lower=c(0), upper=c(1), #mean=bw, sigma=as.matrix(sigw2))), log(pmvnorm(lower=c(0,0), #upper=c(1,1), mean=c(bb, bw), sigma=matrix(c(sigb2,sigbw,sigbw,sigw2), nrow=2))))
		
