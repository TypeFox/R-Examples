##required library(s): numDeriv                                       
## Functions to fit the Negative binomial mixture model with an AR(1) dependence structure
## current options for distribution of random effects (RE): "G" (gamma),"N" (lognormal), "NoN" (non-parametric)
##===========MLE=====================
fitParaAR1 <- function(
                       formula,
                       ## an object of class "formula"
                       ## (or one that can be coerced to that class):
                       ## a symbolic description of the model to be fitted.
                       data,
                       ## a data frame, list or environment (or object coercible
                       ## by as.data.frame to a data frame)
                       ## containing the variables in the model.
                       ID,
                       ## a vector of length n*ni containing patient IDs of observations in data
                       Vcode,
                       ## scan number need to be integers with increment of one, i.e., -1, 0, 1, 2,..
                       p.ini=NULL,     ## initial values for the parameters
                       ## c(log(a), log(th), lgt(dt), b0, b1, ...)
                       IPRT=FALSE,    ## FALSE control: T = print iterations
                       RE="G",    # dist'n for RE G= gamma, N = lognormal
                       i.tol=1.e-75, # tolerance; for integration (ar1.lk)
                       o.tol=1.e-3
                       ## tolerance: for optim 
                       ## enter the number of scans when there is no missing data 
                       ) 
{
  ## Adjust ID and Vcode so that they have the same length as the cleaned data
  ftd <- formulaToDat(formula=formula,data=data,ID=ID,Vcode=Vcode)
  dat <- ftd$dat
  Vcode <- ftd$Vcode
  ## dat = (ID, Y, x1, x2, ...) numeric matrix
  DT <- getDT(dat) ## list

  DT$dif <- c(0,diff(Vcode))   
  DT$dif[DT$ind] <- 0
  ## DT$ind contains the locations of the 1st repeated measures for each patient
  ## DT$diff  = 0 if its the first repeated measure and = 1 ifelse
  DT$dif <- DT$dif[1:DT$totN]    # scan lag

  if (is.null(p.ini))
    {
      p.ini <- rep(0, 4+DT$cn)
      p.ini[4] <- mean(DT$y)   
    }else{
      if (length(p.ini) != 4 + DT$cn) stop("The length of p.ini does not agree to 4 + # covariates")
      }
  
  if (IPRT)
    cat("\n\n estimates: log(a),log(theta),lgt(d),b0,b1,... and the negative of the log-likelihood")

  

  ## R function version
  tt <- optim(p.ini,  ## c(log(a), log(th), lgt(dt), b0, b1, ...)
              ar1.lk, ##ar1.lk likelihood function for gamma/log-normal RE model 
              hessian=TRUE,  
              control=list(reltol=o.tol),
              DT=DT,#input data (output from getDT; lag)
              IPRT=IPRT,
              i.tol=i.tol,
              RE=RE ## Notice that the model option is passed to ar1.lk
              )
  
  
  nlk <- tt$value #neg likelihood
  if (rcond(tt$hessian) < 1E-6){
    vcm <- matrix(NA,nrow=length(p.ini),ncol=length(p.ini))
  }else{
    vcm <- solve(tt$hessian)
  }
  if (is.matrix(vcm)) colnames(vcm) <- rownames(vcm) <- c("log_a", "log_th", "logit_dt","(Intercept)", DT$xnames)
  p.est <- cbind(tt$p, sqrt(diag(vcm)))
  row.names(p.est) <- c("log_a", "log_th", "logit_dt", "(Intercept)", DT$xnames)
  re <- list(opt=tt, nlk=nlk, V=vcm, est=p.est, RE=RE,##idat=data.frame(dat),
             Vcode=Vcode,AR=TRUE,formula=formula)
  class(re) <- "LinearMixedEffectNBFreq"
  return(re)
}



ar1.lk <- function(para,
                   ## c(log(a), log(th), lgt(dt), b0, b1, ...)
                   DT, #input data (output from getDT; lag)
                   IPRT=TRUE,      #print control
                   i.tol=1.e-75,  # tolerance for integration 
                   sig=FALSE,       # if T, the compute full likelihood
                   ##FixN,     # FixN: if uij = uj
                   RE = "G"   #dist'n of the RE; G= gamma, N = lognormal
                   ) 
{
  if(IPRT) cat("\n",para," ")
  ainv <- exp(-para[1])    ## ainv=1/a  
  th1 <- exp(para[2])      ## scalar of gamma /var(G_i) of lognormal
  if (RE=="G") shp <- 1/th1           ## gamma-shape

  ## When G_i ~ LN, var(G_i)=theta
  if (RE=="N") 
    { s.ln <- log(th1+1) ## sigma^2 = log(th1+1) of log-normal
      u.ln <- -s.ln/2    ## mu = -log(th1+1)/2 of log-normal
      s.ln <- sqrt(s.ln) # sigma of log-normal
    }

  dt <- ilgt(para[3])    #inverse logit = delta
  
  tem <- rep(0, DT$totN) #tem = zero vector of length (s*sn) if there is no covariate 
  ## If the number of covariates is greater than 1
  if (DT$cn>0) {
    b <- para[5:(DT$cn+4)] ## = [b1 b2 ...]
    tem <- DT$x%*%b ## tem = b1*x1 + b2*x2 + ... (DT$x is (n*sn) by # covariates)
  }
                                        #r[i,j] = u[i,j]/alpha =exp(- log(alpha)+b0 + b1*x1 + ... )
  th2 <- exp(tem+ ## b1*x1+b2*x2+...
             para[4]- ## b0
             para[1] ## log a
             ) ## th2 is a vector of length (n*sn) 
  
  Dl <- dt^DT$dif ## Dl = 1 for the first repeated measure of each patient and  Dl = d ifelse

  szm <- c(0, th2[-DT$totN]) #r[i, j-1]


  
  U <- Dl*szm  #d*r[i,j-1] 
  V <- szm-U   #(1-d)*r[i, j-1]
  sz2 <- th2-U #r[i, j] - d*r[i, j-1]
  ## tem contains location of the 1st measurements for each patient
  ## Function assumes that DTa contains each patient's measures in a block 
  tem <- DT$ind[1:DT$np]
  sz2[tem] <- th2[tem] #set sz2 = r[i, j] for the first scan
  
  if (any(sz2<=0)) return(1.e15)

  nllk <- 0               #total likelihood
  lki <- rep(0, DT$np)   #likelihood for each individual (when sig=T)

  llk0 <- 0               #likelihood for an observation with all 0 counts
  for (i in 1:DT$np) ## i in 1, ..., # patients
    { ll=DT$ind[i]:(DT$ind[i+1]-1)
      ## ll is a vector, containing the locations of the repeated measures
      ## corresponding to i^th obs.

      if (RE=="G") tem <- integrate(ar1.intg, lower=0, upper=Inf, abs.tol=i.tol,
            a_inv=ainv, sh=shp, sc=th1, 
            y=DT$y[ll], u=U[ll], v=V[ll], s2=sz2[ll])
      
      if (RE=="N") tem <- integrate(ar1.ln.intg, lower=0, upper=Inf, abs.tol=i.tol, 
            a_inv=ainv, mu=u.ln, sig=s.ln, #check
            y=DT$y[ll], u=U[ll], v=V[ll], s2=sz2[ll])
                                        #print(c(DT$ys[i], us[i], tem$v))
      lki[i]=log(tem$value)
      ##if (DT$ni[i]==FixN&DT$ys[i]==0) llk0=lki[i]
      ##}
      nllk=nllk-lki[i]
    }
  if (IPRT) cat(" nllk=", nllk)

  if (sig) return(lki) 
  else return(nllk)
}

##integradient oc to prob(Y|g)*f(g); G~gamma
ar1.intg <- function(x=2,       #G=x
                     a_inv=0.5, #1/alpha 
                     sh=0.5, sc=2, # gamma parameters
                     y=1:3,        #count 
                     u=c(0,1,1),   # d*r[i,j-1] 
                     v=c(0,2,2),   # (1-d)*r[i,j-1]
                     s2=c(3,2,2)   # r[i,j]- d*r[i, j-1]
                     )
  ##y,u,v, s2 need to be of the same length
  
{ 
  Pr=a_inv/(x+a_inv)
  tem=NULL
  for ( pr in  Pr)
    { #print(pr)
      tem <- c(tem,ar1.fun(y, u, v, s2, pr))
    }
  res <- tem*dgamma(x, shape=sh, scale=sc)
  return(res)
}

                                        #example: integrate with gamma density
                                        #integrate(ar1.intg, lower=0, upper=Inf, a_inv=0.5, sh=1/5, sc=5, y=1:3)

##integradient oc to prob(Y|g)*f(g); G~lognormal
ar1.ln.intg=function(x=2, a_inv=0.5, mu=0.5, sig=2, y=1:3, u=c(0,1,1), v=c(0,2,2), s2=c(3,2,2))
{ #see ar1.intg; mu, sig : log normal parameters

  Pr=a_inv/(x+a_inv)
  tem=NULL
  for ( pr in  Pr)
    { #print(pr)
      tem=c(tem,ar1.fun(y, u, v, s2, pr))
    }                                      
  res=tem*dlnorm(x, meanlog=mu, sdlog=sig) 
  return(res)
}

                                        # missing DTa dealt by the approximation
##ar1.fun oc to prob(Y=y|g)
ar1.fun <- function(
                    y=c(0,1,3), #counts
                    U=c(0,1,1), # d*r[i,j-1]; r[i,j-1]=u[i,j]/alpha 
                    V=c(0,2,2), # (1-d)*r[i,j-1]
                    sz2=c(3,2,2), # r[i,j]- d*r[i, j-1]), 
                    pr=0.5      # pr = prob =(1/a)/(g+1/a)
                    )
{
  ## Compute: P(Y_i1=y_i1|G=g) P(Y_i2=y_i2|y_i1,G=g)...P(Y_ini=y_ini|y_i(ni-1),G=g)
  ## where P(Y_i1=y_i1|g_i) = NB(r_i1,pi)
  ##   and P(Y_ij=y_ij|y_i(j-1),g_i) = sum_{k=0}^{min{y_ij,y_ij-1} } P(Z=k|y_i(j-1),g_i)P(eps=y_ij - k |y_i(j-1),g_i)
  n=length(y) ## the number of repeated measure
  
  ## the first repeated measure P(Y_i1=y_i1|g_i)
  pb=dnbinom(y[1], size = sz2[1], prob=pr)
  if (n==1) return(pb)
  for ( i in 2:n)
    { ## P(Y_ij=y_ij|y_i(j-1),g_i)
      k = 0:min(y[i], y[i-1])
      ## betabinomial: three parameters N,u,v
      ## dbb(x, N, u, v) = beta(x+u, N-x+v)/beta(u,v)*choose(N,x)   Wiki-notation is u=alpha, v=beta
      ## x takes values 0,1,..., N 
      pp=dbb(x=k, N=y[i-1], al=U[i], bet=V[i])
      pp=pp*dnbinom(y[i]-k, size = sz2[i], prob=pr)
      p1=sum(pp)
      pb=pb*p1 
    }
  return(pb)
}
##============= end MLE ============

##==========  Index functions  ============
                                        # pr(q(Ynew) >= q(ynew) | Ypre=ypre)
                                        # input: data from one patient

jCP.ar1 <- function(tpar, ## log(a),log(theta),log(delta),b0,...
                                        # parameters (obj is an output from fitSemiAR1 or fitParaAR1) 
                    ypre,## vector of length # pre, containing CEL
                    ynew,  ## vector of length # new, containing CEL
                    y2m=NULL, 
                    XM, # matrix of covariates
                    stp,  
                    RE="G", #options for RE dist'n: G=Gamma, N=lognormal, NoN=nonparametric
                    LG=FALSE,    # if T, return logit(P)
                    MC=FALSE, N.MC=40000, #Monte carlo integration 
                    qfun="sum", #q function
                    oth=NULL,    #if RE=="NoN", oth=obj$gi
                    i.tol=1E-75
                    )
{  
  a <- exp(tpar[1])
  th <- exp(tpar[2])
  dt <- ilgt(tpar[3])
  
  sn <- length(ypre)+length(ynew) ## the number of repeated measures 
  
  u0 <- exp(tpar[4]) ## beta0
  u <- rep(u0, sn)
  ## with covariates
  if (length(tpar)>4) u <- u*exp(XM%*%tpar[-(1:4)]) 
  
  if (MC)
    tem <- MCCP.ar1(ypre=ypre, ynew=ynew, stp=stp, u=u,
                    th=th, a=a, dt=dt, RE=RE, N.MC=N.MC, oth=oth, qfun=qfun)
  else 
    tem <- CP1.ar1(ypre=ypre, ynew=ynew, y2m=y2m, stp=stp, u=u, th=th, a=a, dt=dt, RE=RE, oth=oth, qfun=qfun,i.tol=i.tol)
  if (LG) tem <- lgt(tem)
  return(tem)
}

##subroutines for jCP.ar1
                                        #CP1.ar1 (copied Psum.jun25.R)
                                        #parameters on the original scales

                                        #CP1.ar1(ynew=c(1,0,1), mod="NoN", oth=obj$gi, qfun="max")
CP1.ar1 <- function(ypre,ynew,y2m=NULL, stp, 
                    u,th, a, dt, #parameters on the original scales
                    RE="G",  #G=gamma, NoN=nonparam
                    oth, ## NULL for parametric model
                    qfun,i.tol=1E-75)
  {
    if (all(is.na(ynew))) return (NA)
    
    Qf=match.fun(qfun)
    newQ=Qf(ynew, na.rm=T)
    if (newQ==0) return(1)

    ain=1/a
    if (RE=="G") shp=1/th
    if (RE=="NoN") 
      { tem=sort(round(oth,6))
        oth1=unique(tem)
        gp=table(tem)/length(tem)
                                        #Pr=ain/(ain+oth1)
        Pr=1/(oth1*a+1)
      }
    
    n0=length(ypre)
    n1=length(ynew)
    l1=n0+n1  

    y=c(ypre, ynew) 
    

    sz=u/a
    DT=dt^stp
    
                                        #parameters for ypre
    szm=c(0, sz) 
    szm=szm[-length(szm)]
    u=DT*szm
    v=szm-u
    sz2=sz-u                           #possible combinations for ynew under q()
                                        #y2m=getY2.1(newQ, l1-l0, qfun)
    if (is.null(y2m)) y2m=getY2.1(newQ, l1-n0, qfun) 
    ## l1-n0 = n1 # the number of new-scans
    ## newQ = q(Y_new) 
    if (RE=="G")
      { #Pr(Ypre=Ypre)
        tem <- integrate(ar1.intg, lower=0, upper=Inf,  a_inv=ain, sh=shp, sc=th, 
                         y=ypre, u=u[1:n0], v=v[1:n0], s2=sz2[1:n0],abs.tol = i.tol)
        bot <- tem$v
        ##cat("\n den=",bot)
        tem <- integrate(ar1.2.mmg, lower=0, upper=Inf, rel.tol=0.1, 
                         a_inv=ain, sh=shp, sc=th, 
                         y1=ypre, y2m=y2m, u=u, v=v, sz2=sz2,abs.tol = i.tol)
        top <- tem$v
                                        #print(top)
      }
    if (RE=="N") #added May 18, 2012
      {
        ## YK: June 6 {
        s.ln <- log(th+1)    
        u.ln <- -s.ln/2     
        s.ln <- sqrt(s.ln)
        ## }
        tem <- integrate(ar1.ln.intg, lower=0, upper=Inf,  
                         a_inv=ain, mu=u.ln, sig=s.ln, 
                         y=ypre, u=u[1:n0], v=v[1:n0], s2=sz2[1:n0],abs.tol = i.tol) 
        bot <- tem$v
        tem <- integrate(ar1.ln2.mmg, lower=0, upper=Inf, abs.tol = i.tol, 
                         a_inv=ain, mu=u.ln, sig=s.ln, 
                         y1=ypre, y2m=y2m, u=u, v=v, sz2=sz2)
        top <- tem$v
      }

    if (RE=="NoN")
      { p1.non=p2.non=NULL
        for ( pr in Pr)
          { #p1=ar1.fun(y1, u[1:l0], v[1:l0], sz2[1:l0], pr)
            p1=ar1.fun(ypre, u[1:n0], v[1:n0], sz2[1:n0], pr)

            p1.non=c(p1.non, p1)
            
            p2=0
            for (j in 1:nrow(y2m))
              { #yy=c(y1[l0], y2m[j,]) 
                                        #p2=p2+ar22.fun(yy, u[l0:l1], v[l0:l1], sz2[l0:l1], pr) 
                yy=c(ypre[n0], y2m[j,]) 
                p2=p2+ar22.fun(yy, u[n0:l1], v[n0:l1], sz2[n0:l1], pr) 
              }
            p2.non=c(p2.non, p2*p1) 
          }
        bot=sum(p1.non*gp)
        top=sum(p2.non*gp)
      }

                                        #print(c(bot, top))
    res=1-top/bot
    return(res)
  }


                                        # Pr( Ypre=y1, Ynew[1:k-1] = y2m[i,1:(k-1)], Ynew[k]<=y2m[i,k] | G=x) 
                                        # with G~gamma for all row of y2m
ar1.2.mmg=function(x=0.5, a_inv=2, sh=1/3, sc=3, y1=c(0,1), 
  y2m=getY2.1(), #a matrix, each row is a set of value for Ynew
  u=c(0,1,1,1,1), v=c(0,2,2,2,2), sz2=c(3,2,2,2,2))
{  #print(x)
                                        #print(length(x))
  
  m=nrow(y2m)
  Pr=a_inv/(a_inv+x)
  l0=length(y1)
  l1=ncol(y2m)+l0     

  tem=NULL
  for ( pr in Pr)
    { 
      p1=0
      for (j in 1:m)
        {  
          y=c(y1[l0], y2m[j,])
          ## ar22.fun computes Pr(Yfol1 = yfol1,...,YfolM<=yfolM | Ypre2=ypre2 )
          p1=p1+ar22.fun(y=y, U=u[l0:l1], V=v[l0:l1], sz2=sz2[l0:l1], pr=pr) 
        }
      ## ar1.fun computes Pr(Ypre1 = ypre1 Ypre2=ypre2 )
      p1=p1*ar1.fun(y=y1, U=u[1:l0], V=v[1:l0], sz2=sz2[1:l0], pr=pr) 
      tem=c(tem, p1)
    }
  tem=tem*dgamma(x, shape=sh, scale=sc)
  return(tem)
}

ar1.ln2.mmg=function(x=0.5, a_inv=2, mu=1/3, sig=3, y1=c(0,1), 
  y2m=getY2.1(), #a matrix, each row is a set of value for Ynew
  u=c(0,1,1,1,1), v=c(0,2,2,2,2), sz2=c(3,2,2,2,2))
{  
  m=nrow(y2m)
  l0=length(y1)
  l1=l0+ncol(y2m)
  Pr=a_inv/(a_inv+x)
  ## x = RE. As this function is input of "integrate" function,
  ## this function should be evaluated with vector x
  
  tem=NULL
  for ( pr in Pr)
    { #Pr(Ypre=ypre)
      
      p1=0
      for (j in 1:m) ## run through all the possible of Yfol that return Yfol+ = yfol+
        {   
          y=c(y1[l0], y2m[j,])
          ## ar22.fun computes Pr(Yfol1 = yfol1,...,YfolM<=yfolM | Ypre2=ypre2 )
          p1 = p1 + ar22.fun(y=y, U=u[l0:l1], V=v[l0:l1], sz2=sz2[l0:l1], pr=pr) 
        }
      ## ar1.fun computes Pr(Ypre1 = ypre1 Ypre2=ypre2 )
      p1=p1*ar1.fun(y=y1, U=u[1:l0], V=v[1:l0], sz2=sz2[1:l0], pr=pr)
      tem=c(tem, p1)
                                        #print(pr)
    }
  tem=tem*dlnorm(x, meanlog=mu, sdlog=sig)
  return(tem)
}

## Pr(Y2=y2, Y3=y3, ..., Yn<= yn|Y1=y1)
ar22.fun=function(y=c(0,1,3), U=c(0,1,1), V=c(0,2,2), sz2=c(3,2,2), pr=0.5)
{  n=length(y)   

   pb=1
   for ( i in 2:n)
     { k = 0:min(y[i], y[i-1])
                                        #betabinomial dbb(x, N, u, v)
                                        #source("/home/yinshan/Mydesk/Mylib/Rlib/myfun/beta-binomial.R")
       pp=dbb(x=k, N=y[i-1], al=U[i], bet=V[i])
       if (i<n) { pp=pp*dnbinom(y[i]-k, size = sz2[i], prob=pr)}
       else  { pp=pp*pnbinom(y[i]-k, size = sz2[i], prob=pr) }
       p1=sum(pp)
       pb=pb*p1 
     }
   return(pb)
 }


                                        #MCCP.ar1: monto carlo for conditional probability
                                        #MCCP.ar1=function(ypre, ynew, u, th, a, dt, RE="G", N.MC=1000, oth=gi, qfun=sum)
MCCP.ar1 <- function(ypre, ynew, stp, u, th, a, dt, RE="G", N.MC=1000, oth, qfun="sum") 
{ if (all(is.na(ynew))) return (NA)

  Qfun=match.fun(qfun)
  newQ=Qfun(ynew, na.rm=TRUE)
  if (newQ==0) return(1)

  if (newQ==1) 
    {  res=CP1.ar1(ypre=ypre, ynew=ynew, stp=stp, u=u, th=th, a=a, dt=dt, RE=RE, oth=oth,qfun=qfun)
       return(res)
     }

  ain=1/a
  if (RE=="G") shp=1/th

  if (RE=="N") 
    { 
      s.ln=log(th+1)   
      u.ln=-s.ln/2     
      s.ln=sqrt(s.ln)  
    }

  if (RE=="NoN") 
    { tem=sort(round(oth,6))
      oth1=unique(tem)
      gp=table(tem)/length(tem)
      Pr=1/(1+a*oth1)
    }
  
  n0=length(ypre)
  n1=length(ynew)
  l1=n0+n1  

  sz=u/a

  dt=dt^stp 
  szm=c(0, sz)
  szm=szm[-length(szm)]
  u0=dt*szm 
  v0=szm-u0
                                        #sz2=sz0[ms0]-u0  #YZ remove (May 24, 2012)
  sz2=sz-u0   
  Opre=list(y=ypre, u=u0[1:n0], v=v0[1:n0], s2=sz2[1:n0])
  Onew=list(y0=ypre[n0], u=u0[-(1:n0)], v=v0[-(1:n0)], s2=sz2[-(1:n0)], n1=n1)
  
  
  y=c(ypre, ynew) 
  if (RE=="G")
    { #Pr(Ypre=Ypre)
      tem=integrate(ar1.intg, lower=0, upper=Inf,  
        a_inv=ain, sh=shp, sc=th, 
                                        #y=ypre[ms0], u=u0, v=v0, s2=sz2) #YZ old
        y=ypre, u=u0[1:n0], v=v0[1:n0], s2=sz2[1:n0]) #YZ new (May 24, 2012)
      bot=tem$v

      tem=integrate(ar1.mmg, lower=0, upper=Inf, rel.tol=1.e-2, 
        a_inv=ain, sh=shp, sc=th, 
        O1=Opre, O2=Onew, tot=newQ, N=N.MC, qfun=qfun) #YZ new (May 24, 2012)
      
      top=tem$v
    }
  
  if (RE=="N")
    { #Pr(Ypre=Ypre)
      tem=integrate(ar1.ln.intg, lower=0, upper=Inf,  
        a_inv=ain, mu=u.ln, sig=s.ln, 
        y=ypre, u=u0[1:n0], v=v0[1:n0], s2=sz2[1:n0])
      bot=tem$v

                                        #Pr(Ypre=ypre, Ynew+ < newSum) 
      tem=integrate(ar1.mm.ln, lower=0, upper=Inf, rel.tol=1.e-2, 
        a_inv=ain, mu=u.ln, sig=s.ln, 
        O1=Opre, O2=Onew, tot=newQ, N=N.MC, qfun=qfun) 
      
      top=tem$v
    }
  

  if (RE=="NoN")
    { p1.non=p2.non=NULL
      for ( pr in Pr)
        { 
          p1=ar1.fun(ypre, u0, v0, sz2, pr)  
          p2=p1*mmS.fun(Obj=Onew, pr=pr, nSum=newQ, N=N.MC, qfun=qfun)

          p1.non=c(p1.non, p1)
          p2.non=c(p2.non, p2)
        }
      bot=sum(p1.non*gp)
      top=sum(p2.non*gp)
    }

                                        #print(c(bot, top))
  res=1-top/bot
  return(res)
}

##Pr(q(Ynew) < nSum |y0; pr) using MC

mmS.fun=function(Obj=list(y0=0, u=c(1,1,1), v=c(2,2,2), s2=c(2,2,2), n1=3),
  pr=0.5, nSum=2, N=1000, qfun="sum")
{ Qfun=match.fun(qfun)

                                        #k=length(sz1)
                                        #u=sz1[-k]*dt
                                        #v=sz1[-k]-u
                                        #s2=sz1[-1]-u

  u=Obj$u
  v=Obj$v 
  s2=Obj$s2

  y0=Obj$y0
  YY=NULL
  for (i in 1:Obj$n1)
    {  
      z1=rbb(N, y0, u[i], v[i])
      z2=rnbinom(N, size=s2[i], prob=pr)
      y1=z1+z2
      YY=cbind(YY, y1)
      y0=y1
    }
  if (Obj$n1>1) 
    { tot=apply(YY, 1, Qfun) 
    }
  else tot=YY

  return(mean(tot<nSum))
}


ar1.mmg=function(x, a_inv, sh, sc, O1, O2, tot, N, qfun=sum)  
{  #print(x)
                                        #print(length(x))
  Pr=a_inv/(a_inv+x)
  tem=NULL
  for ( pr in Pr)
    { #Pr(Ypre=ypre)
      p1=ar1.fun(y=O1$y, U=O1$u, V=O1$v, sz2=O1$s2, pr=pr)   
      p2=mmS.fun(Obj=O2, pr=pr, nSum=tot, N=N, qfun=qfun) 
                                        #print(c(pr, p1, p2))
      tem=c(tem,p1*p2)
    }
  tem=tem*dgamma(x, shape=sh, scale=sc)
  return(tem)
}


ar1.mm.ln <- function(x, a_inv, sig, mu, O1, O2, tot, N, qfun=sum)
{  Pr=a_inv/(a_inv+x)
   tem=NULL
   for ( pr in Pr)
     {
       p1=ar1.fun(y=O1$y, U=O1$u, V=O1$v, sz2=O1$s2, pr=pr)   
       p2=mmS.fun(Obj=O2, pr=pr, nSum=tot, N=N, qfun=qfun)
       tem=c(tem,p1*p2)
     }
   tem=tem*dlnorm(x, meanlog=mu, sdlog=sig)
   return(tem)
 }










CP.ar1.se <- function(
                      tpar,
                      ypre,
                      ynew,
                      y2m=NULL,
                      XM, 
                      ## see jCP.ar1
                      stp, 
                      RE="G",   # dist'n of REs G=Gamma, N =lognormal
                      V,  
                      MC=FALSE,       # if true use MC
                      qfun="sum",i.tol=1E-75)
{  if (all(is.na(ynew))) return(c(NA, NA))
   Qf=match.fun(qfun)
   newQ=Qf(ynew, na.rm=T)
   if (newQ==0) return(c(1,0))
   
   p <- jCP.ar1(
                tpar=tpar,## log(a),log(theta),log(delta),b0,...
                ypre=ypre,## vector of length # pre, containing CEL
                ynew=ynew,## vector of length # new, containing CEL
                y2m=y2m, stp=stp,
                XM=XM, RE=RE, MC=MC, qfun=qfun,LG=FALSE,i.tol=i.tol)
   if (!is.null(V)){
     if (MC) mth="simple" else mth="Richardson"
     
     jac <- jacobian(func=jCP.ar1, x=tpar, method=mth,
                     method.args=list(eps=0.01, d=0.01, r=2), ypre=ypre, ynew=ynew,
                     y2m=y2m, XM=XM, stp=stp, RE=RE, LG=TRUE, MC=MC, qfun=qfun,i.tol=i.tol)  
     s2=jac%*%V%*%t(jac)
     s=sqrt(s2) #SE of logit(Phat)
   }else s <- NA
   return(c(p,s)) ## (Phat, SE(logit(Phat)))
 }





fitParaIND <-
  function(formula,     ## an object of class "formula"
           ## (or one that can be coerced to that class): a symbolic description of the model to be fitted.
           data,        ## a data frame, list or environment (or object coercible by as.data.frame to a data frame)
           ## containing the variables in the model.
           ID,          ## a vector of length n*ni containing patient IDs of observations in data
           p.ini=NULL,    # initial value c(log(a), log(th), b0, b1, ...)
           IPRT=FALSE,   # printing control: T = print iterations
           RE="G", #RE options: G = gamma; N = lognormal
           i.tol=1e-75,
           o.tol=1.e-3,  # tolerance: for optim
           COV=TRUE     ## Covariance matrix == TRUE
           )
{ ##########################
  ## Tips for computation ##
##########################
  ## Note that fitParaIND fails to compute the negative log-likelihood values
  ## when the total sum of response counts of a patient is **very large**
  ## This is because the log-likelihood is:
  ## 
  ##         integrate dnbinom(sum_j=1^n_i y_ij; sum_j=1^n_i r_ij,p_i)*f(g) dg
  ## \propto integrate p_i^{sum_j=1^n_i r_ij} (1-p_i)^{sum_j=1^n_i y_ij} *f(g) dg
  ##
  ## since 0 <1-p_i < 1 if sum_j=1^n_i y_ij is **very** large, then (1-p_i)^{sum_j=1^n_i y_ij} = 0.
  ## and the integrated value become zero. e.g., 0.5^1000 = 0 in R
  ## Since log of this quantity is taken,
  ## the evaluated value become log(0) = - Inf
  ##
  cmd <- match.call() 
  ftd <- formulaToDat(formula=formula,data=data,ID=ID)  # dat = (ID, Y, x1, x2, ...) numeric matrix
  dat <- ftd$dat
  DT <- getDT(dat)

  if (is.null(p.ini))
    {
      p.ini=rep(0, 3+DT$cn) 
      p.ini[3]=mean(DT$y)   
    }
  
  if (IPRT) cat("Print log_a log_th log_u0", DT$xnames, " and negative of the log-likelihood at each iteration")
  
  if (RE=="G") ## Gi ~ Gamma(scale=theta,shape=1/theta)
    tt <- optim(p.ini, lk.fun,control=list(reltol=o.tol), hessian=TRUE, dat=DT, Iprt=IPRT,tol=i.tol)
  else if (RE=="N")## Gi ~ logN(mu=1,var(G_i)=theta)
    tt <- optim(p.ini, lk.ln.fun,control=list(reltol=o.tol), hessian=TRUE, dat=DT, Iprt=IPRT,tol=i.tol)
  
  nlk <- tt$value+sum(lgamma(DT$y+1))  #full -loglikelihood
  ## if (IPRT) print(tt$hessian)
  
  if (COV){
    if (rcond(tt$hessian) > 1E-6){
      vcm <- solve(tt$hessian)   # the covariance matrix  = inverse Hessian
    }else vcm <- "hessian is not invertible" 
  }else vcm <- matrix(NA,length(tt$p),length(tt$p))
  p.est <- cbind(tt$p, sqrt(diag(vcm)))
  row.names(p.est) <- rownames(vcm)<-colnames(vcm)<- c("log_a", "log_th", "(Intercept)", DT$xnames)
  colnames(p.est) <- c("estimate","SE")
  
  re <- list(call=cmd, p.ini=p.ini,opt=tt,formula=formula,
             nlk=nlk, V=vcm, est=p.est, RE=RE, ##idat=data.frame(dat),
             AR=FALSE)
  class(re) <- "LinearMixedEffectNBFreq"
  return(re)
}








## Likelihood function for gamma random effect model
lk.fun <-
  function(para, ## parameters (log(a), log(th), b0, b1, ...) 
           dat,  ## dat=list(id=ID, y=ydat, x=xdat, cn=ncov, 
           ##                ys=Ysum, np=N, totN=totN, ni=Ni, ind=IND)
           ## output of getDT()
           Iprt=FALSE, #printing control 
           tol=1.e-75, #control parameter for integration 
           sig=FALSE # if T, compute the full likelihood 
                                        #(i.e., with the constant)
           )
{ if (Iprt) cat(para)
  a <- exp(para[1])
  th1 <- exp(para[2]) #scale
  ainv <- 1/a  
  shp <- 1/th1

  u0 <- exp(para[3])
  th2 <- rep(u0/a, dat$totN)
  if (dat$cn>0) {
    ## if there are covariates
    b <- para[4:(dat$cn+3)]
    tem <- exp(dat$x%*%b)
    th2 <- tem*th2 ## th2 = r_{ij}=exp(X^T beta)/a: Y_{ij}|G_i=gi ~NB(r_{ij},p_i)
  }
  ## -loglikelihood
  ## constant terms of - loglikelihood
  nllk <- sum( - lgamma(dat$y+th2) + lgamma(th2))
  us <- tapply(th2, dat$id, sum)
  lk <- rep(0, dat$np)
  for (i in 1:dat$np)
    { tem=integrate(int.fun, lower=0, upper=Inf, a_inv=ainv, abs.tol=tol,
        ## dat$ys[i] = sum(y_{ij},j=1,...,ni)
        sh=shp, sc=th1, ysum=dat$ys[i],  usum=us[i])
                                        #print(c(dat$ys[i], us[i], tem$v))
      nllk = nllk - log(tem$v)
                                        #sig=T full likelihood
      if (sig) 
        {  ll=dat$ind[i]:(dat$ind[i+1]-1)
           lk[i]=log(tem$v)+sum(lgamma(dat$y[ll]+th2[ll]) - lgamma(th2[ll])-lgamma(dat$y[ll]+1))
         }
    }
  if (Iprt) cat(" nllk=", nllk, "\n")
  if (sig) return(lk) #return full likelihood (vector of length n)
  else return(nllk)   #return log likelihood without constant
}

int.fun <-
  function(x=2, # value of G 
           a_inv=0.5,  # a_inv=1/a
           sh=0.5, sc=2, #shape and scale; sh=1/th; sc=th to have E(G_i)=scale*shape=1
           ysum=2, #sum(y.ij,j=1,..,ni); ysum is a number
           usum=3  #sum(u.ij/a,j=1,...,ni); usum is a number
           )
  ## note p=1-p
{   p=x/(x+a_inv)  
    tem=p^ysum*(1-p)^usum*dgamma(x, shape=sh, scale=sc)
    return(tem)
  }


## ========= log-normal random effects =========
lk.ln.fun <- function(para, dat, Iprt=FALSE, tol=1e-75)
{ if (Iprt) cat("\n",para)
  
  a=exp(para[1])
  th1=exp(para[2]) #scale
  
  s.ln=log(th1+1) ##s^2
  u.ln=-s.ln/2    ##-s^2/2, mean for dlnorm()
  s.ln=sqrt(s.ln) ##s, sd for dlnorm()

  ainv=1/a  

  u0=exp(para[3])
  th2=rep(u0/a, dat$totN)
  if (dat$cn>0) {
    b=para[4:(dat$cn+3)]
    tem=exp(dat$x%*%b)
    th2=tem*th2
  }

  nllk=sum( - lgamma(dat$y+th2) + lgamma(th2))
  us=tapply(th2, dat$id, sum)
  
  for (i in 1:dat$np) 
    { tem=integrate(int1.ln, lower=0, upper=Inf, a_inv=ainv,
        sh=u.ln, sc=s.ln, ysum=dat$ys[i],  usum=us[i], abs.tol=tol)
                                        #print(c(dat$ys[i], us[i], tem$v))
      nllk = nllk - log(tem$v)
    }
  if (Iprt) cat(" nllk=", nllk)
  return(nllk)
}






getY2.1=function(tot=2, ## q(Y_new)
  n=3, ## # new scans
  qfun="sum")
{
  if (tot==0){ return(matrix(rep(0,n),nrow=1 ))}
  if (n==1) 
    { ym=tot-1
      return(matrix(ym, ncol=n))
    }

  N=tot^(n-1)
  ii=0:(N-1)
  
  ym=NULL
  for (j in 1:(n-1))
    { t2=ii%%tot
      ii=ii%/%tot
      ym=cbind(ym, t2)
    }
  if (qfun=="sum") 
    {
      ysum=apply(ym, 1,sum)
      yn=tot-ysum-1
      ym=cbind(ym, yn)
      ym=ym[yn>=0,]

    }else{
      
      yn=rep(tot-1, nrow(ym)) 
      ym=cbind(ym, yn)
    }
  
  return(matrix(ym, ncol=n))
}




































Psum1 <-
  function(Y1, Y2,     #Y1=sum(y.pre), Y2=sum(y.new)
           u1, u2, # u1=mean(Y1), u2=mean(Y2)
           a, Th,           #Var(Gi), under the gamma model, Th=scale
           RE="G",       #"G" = gamma, "N" = log-normal, "GN" = mix of gamma and normal#"U" = log-uniform #"NoN" = non-parametric
           othr=NULL,       # if dist= "GN",
           sn1,i.tol
           ## othr=list(u.n=3, s.n=0.5, p.mx=0.05, sh.mx=NA)
           ## if dist="NoN", 
           ## othr = ghat frequency table of gi (olmeNB$gtb)
           ## for other dist options, othr = NULL 
           )
  { if (Y2==0) {return(1)}
    uVa=c(u1, u2)/a
    if (RE=="NoN")
      { 
        tem <- Psum.non(Y1=Y1, Y2=Y2, u1=u1, u2=u2, a=a, gi=othr)
        return(tem)
      }
    else if (RE=="G") #gamma
      { 
        tem1 <- integrate(cum.fun, lower=0, upper=Inf, a_inv=1/a, sh=1/Th,
                          sc=Th, y1=Y1, y2=Y2, u1=uVa[1], u2=uVa[2], abs.tol=i.tol,sn1=sn1)
        if (sn1 > 0) tem2 <- integrate(int.fun, lower=0, upper=Inf, a_inv=1/a, sh=1/Th, sc=Th, ysum=Y1, usum=uVa[1], abs.tol=i.tol)
        else{
          tem2 <- list();
          tem2$v <- 1
        }
      }  
    else if (RE=="N")
      {  t1 <- sqrt(log(Th+1)) #sd (sc)
         t2 <- -t1^2/2 #mean (sh)
         tem1 <- integrate(cum1.ln, lower=0, upper=Inf, a_inv=1/a, sh=t2, sc=t1, y1=Y1, y2=Y2, u1=uVa[1], u2=uVa[2], abs.tol=i.tol,sn1=sn1)
         if (sn1 > 0) tem2 <- integrate(int1.ln, lower=0, upper=Inf, a_inv=1/a, sh=t2, sc=t1, ysum=Y1, usum=uVa[1],abs.tol=i.tol)
         else{
           tem2 <- list();
           tem2$v <- 1
         }
       }
    else {
      stop("RE must be G, N or NoN")
    }

    val=tem1$v/tem2$v
    return(1-val)
  }



int1.ln <-
  function(x=2, 
           a_inv=0.5, #see int.fun
           sh=0.5, #mean for dlnorm()
           sc=2,   #sd for dlnorm()
           ysum=2, usum=3 #see int.fun
           )
{   p=x/(x+a_inv)  
    tem=p^ysum*(1-p)^usum*dlnorm(x, meanlog=sh, sdlog=sc)
    return(tem)
  }

Pmax1 <-
  function(Y1=0,                #sum(ypre)
           Y2=1,                #max(ynew)
           u1=3,                #Ex(sum(Ypre))
           u2=c(1.5, 1.5, 1.5), #Ex(Ynew); vector 
           a=0.5, 
           Th=3,                #var(G), no use if RE = "NoN" 
           RE="G",            # "G"=gamma, "N" = lognormal, "NoN = nonparametric
           othr=NULL,            # othr=gi (a vector) if  RE = "NoN", not used otherwise
           sn1,i.tol
           )
{ if (Y2==0) {return(1)}
  uVa1 <- u1/a
  uVa2 <- u2/a
  
  if (RE=="NoN")
    { tem <- Pmax.non(Y1=Y1, Y2=Y2, u1=u1, u2=u2, a=a, gi=othr)
      return(tem)
    }
  else if (RE=="G") 
    {
      tem1 <- integrate(max.fun, lower=0, upper=Inf, a_inv=1/a, sh=1/Th, sc=Th, y1=Y1, maxY2=Y2, u1=uVa1, u2=uVa2, abs.tol=i.tol,sn1=sn1)
      if (sn1 > 0) tem2 <- integrate(int.fun, lower=0, upper=Inf, a_inv=1/a, sh=1/Th, sc=Th, ysum=Y1, usum=uVa1, abs.tol=i.tol)
      else{
        tem2 <- list();tem2$v <- 1
      }
    }  
  else if (RE=="N")
    {  t1=sqrt(log(Th+1)) #sd (sc)
       t2=-t1^2/2 #mean (sh)
       tem1 <- integrate(max.ln, lower=0, upper=Inf, a_inv=1/a, sh=t2, sc=t1, y1=Y1, maxY2=Y2, u1=uVa1, u2=uVa2, abs.tol=i.tol,sn1=sn1)
       if (sn1 > 0) tem2 <- integrate(int1.ln, lower=0, upper=Inf, a_inv=1/a, sh=t2, sc=t1, ysum=Y1, usum=uVa1,abs.tol=i.tol)
       else{
         tem2 <- list();tem2$v <- 1
       }
     }
  else stop("RE must be G, N or NoN")
  val <- tem1$v/tem2$v
  return(1-val)
}





## ================ numerator of CPI ===========================

cum.fun <-
  function(
           ## Pr(Y_i,new+ >= y_i,new+, Y_i,pre+ = y_i,pre+, G_i=g) 
           ## = Pr(Y_i,new+ >= y_i,new+|Y_i,pre+ = y_i,pre+, G_i=g) Pr(Y_i,pre+ = y_i,pre+|G_i=g) Pr(G_i=g)
           ## = Pr(Y_i,new+ >= y_i,new+|G_i=g) Pr(Y_i,pre+ = y_i,pre+|G_i=g) Pr(G_i=g) ## conditional independence
           ## = (1-pnbinom(Y_i,new+;size=u2,prob=p))*dnbinom(Y_ipre+;size=u1,prob=p)*dgamma(gi;shape=sh,scale=1/sh)
           ## where u1 = sum_{j in old scan} exp(beta^T * Xij) and u2 = sum_{j in new scan} exp(beta^T * Xij)  
           x=2,          # value of the random effect
           a_inv=0.5,    # 1/a
           sh=0.5, sc=2, # shape and scale of the Gamma RE'n
           y1=2, y2=2,   # Y1=sum(y.pre), Y2=sum(y.new)
           u1=3, u2=3,    # u1 = r1 = mu1/a; u2= r2 =mu2/a
           sn1
           )
{
  p <- x/(x+a_inv)  ## in this manuscript p = 1-p
  if (sn1>0) tem <- p^y1*(1-p)^u1 else tem <- 1
  tem <- tem*dgamma(x, shape=sh, scale=sc)
  tem <- tem*pnbinom(y2-1, prob=1-p, size=u2)
  return(tem)
}



max.fun <-
  function(x=1,      #value of the random effect
           y1=0,                #sum(ypre); i.e. ypre+
           maxY2=3,             #max(Ynew)
           u1=1.5,              #Ex(Ypre+)
           u2=c(1.5,1.5,1.5),   #Ex(Ynew), vector
           a_inv=0.5,           #1/a
           sh=0.5, sc=2,sn1         #shame and scale of the Gamma dist'n
           )
{   
  P <- x/(x+a_inv) 
  ## pr(Y1+ = y1+ |g) 
  if (sn1 > 0) tem <- P^y1*(1-P)^u1 else tem <- 1
  tem <- tem*dgamma(x, shape=sh, scale=sc)
  for (i in 1:length(u2))
    { 
      tem1 <- pnbinom(maxY2-1, prob=1-P, size=u2[i])
      tem <- tem*tem1
    }
  return(tem) ## change from tem*v Jun 13
}



cum1.ln <-
  function( x=2, #value of the random effect
           a_inv=0.5,             # 1/a
           sh=0.5, sc=2,          # lnorm paramters, sc=s=sqrt(log(Th+1)), sh = u = -s^2/2 
           y1=2, y2=2, u1=3, u2=3,sn1 #see cum.fun
           )
{   p=x/(x+a_inv)  
    if (sn1 > 0) tem <- p^y1*(1-p)^u1 else tem <- 1
    tem <- tem*dlnorm(x, meanlog=sh, sdlog=sc)
    tem <- tem*pnbinom(y2-1, prob=1-p, size=u2)
    return(tem)
  }

max.ln <-
  function(x=1, y1=0, maxY2=3, u1=1.5, u2=c(1.5,1.5,1.5), a_inv=0.5, 
                                        #see max.fun  
           sh=0.5, sc=2,sn1 #sh=mean, sc=sd of the lognormal dist'n
           )
{   
  P <- x/(x+a_inv) 
                                        #pr(Y1+= y1+ |g) 
  if (sn1 > 0) tem <- P^y1*(1-P)^u1 else tem <- 1
  tem <- tem*dlnorm(x, meanlog=sh, sdlog=sc)
  for (i in 1:length(u2)){
    
    tem1 <- pnbinom(maxY2-1, prob=1-P, size=u2[i])
    tem <- tem*tem1
  }
  return(tem)
}

## ===========================================

CP.se <-
  function(tpar,
           ## return
           ## 1) the point estimate of the conditional probability of observing
           ## the response counts as large as the observed ones given the previous counts
           ## 2) its asymptotic standard error
           
           ## estimates of log(alpha, sigma.G, betas); 
           ## e.g., tpar=olmeNB$est[,1] where olmeNB = output of the mle.*.fun function 
           Y1,          ## Y1 = y_i,old+ the sum of the response counts in pres
           Y2,          # Y2 = q(y_new)
           sn1, sn2,  # sn1 and sn2: number of scans in the pre and new sets.
           XM=NULL,       # XM : matrix of covariates
           RE="G",      # distribution of the random effects (G = gamma, N = log-normal)
           V,    # variance-covariance matrix of the parameter estimates; olmeNB$V
           qfun="sum",i.tol=1E-75)     # q() = "sum" or "max" 
{ if (Y2==0) return(c(1,0))
  
  ## the point estimate Phat
  p <- jCP(tpar=tpar, Y1=Y1, Y2=Y2, sn1=sn1, sn2=sn2, XM=XM, RE=RE, qfun=qfun, LG=FALSE,i.tol=i.tol)
  
  ## SE of logit(Phat)
  jac <- jacobian(func=jCP, x=tpar, Y1=Y1, Y2=Y2, sn1=sn1, sn2=sn2, XM=XM, RE=RE, LG=TRUE, qfun=qfun,i.tol=i.tol)
  s2 <- jac%*%V%*%t(jac)
  s <- sqrt(s2) 
  return(c(p,s)) ## c(Phat, SE(logit(Phat)))
}


jCP <-
  function(tpar,
           ## estimates of log(alpha, sigma.G) betas; 
           ## e.g., tpar=olmeNB$est[,1] where olmeNB = output of the mle.*.fun function 
           Y1,
           ## Y1 = y_i,old+ the sum of the response counts in pres
           Y2,
           ## Y2 = q(y_new)
           sn1, sn2,
           ## sn1 and sn2: number of scans in the pre and new sets.
           XM=NULL,
           ## XM : matrix of covariates
           RE="G", 
           ## same as CP.se for description
           LG=FALSE,     #indicatior for logit transformation
           oth=NULL, # see Psum1 and Pmax1
           qfun="sum",i.tol=1E-75) 
{
  ## Return a point estimate
  a <- exp(tpar[1])
  th <- exp(tpar[2])

  sn <- sn1+sn2 ## ni 
  u0 <- exp(tpar[3])
  u <- rep(u0, sn)
  if (length(tpar)>3) u <- u*exp(XM%*%tpar[-(1:3)])
  
  u1 <- sum(u[1:sn1])
  if (qfun=="sum") 
    {
      u2 <- sum(u[-(1:sn1)])
      temp <- Psum1(Y1=Y1, Y2=Y2, u1=u1, u2=u2, a=a, Th=th, RE=RE, othr=oth,sn1=sn1,i.tol=i.tol)
    }else { ## max
      u2 <- u[-(1:sn1)]
      temp <- Pmax1(Y1=Y1, Y2=Y2, u1=u1, u2=u2, a=a, Th=th, RE=RE, othr=oth,sn1=sn1,i.tol=i.tol)
    }
  if (LG) temp <- lgt(temp)
  return(temp)
}


tp.fun <- function(i1, ## idat$ID
                   i2, ## idat$Vcode
                   y ## idat$CEL
                   )
{
  ## transpose the counts into a matrix so that each patient
  ## is a row and each visit is a column.
  ## It also puts an NA at any  visit with no data.
  x1 = table(i1,i2)
  x2 = match(x1,1) # 0->NA
  x = matrix(x2, nrow(x1), ncol(x1))
  dimnames(x) = dimnames(x1)
  i1 = as.numeric(factor(i1))
  i2 = as.numeric(factor(i2))
  for (i in 1:length(i1))
    {
      x[i1[i],i2[i]]=y[i]
      ## x1[i1[i],i2[i]]=prism.na$NASS.D[i]
    }
  ## return i1 by i2 contingency table
  return(x)
}

