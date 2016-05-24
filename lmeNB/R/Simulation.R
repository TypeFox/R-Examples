## Jun 22, 2012
## combining ind.sdt and ar1.sdt

## rNBME.R is copied fro ind.sdt
## main changes: 
## (1) Add option 'NoN' to gdist. Under this option, othrp = gi.
## (2) Add argument 'd', a value between (0,1) or NULL
## if 'is.null(d)',  generate data from the independent model
## if '!is.null(d)', generate data from the AR(1) model
## rNBME.R calls 'rnb.ar1' to generate data from a AR(1) model

## rNBME.R requires 'nr.fun' and 'rbb'

rNBME.R <-
  function(gdist="G",  #dist'n for G, G=gamma, N=lognormal, U=uniform, GN=mix of G and normal, NoN for nonparametric
           n=200,  #number of patients,
           sn=5,   #number of scans per person         
           th=exp(1.3), ## scale parameter of gamma
           ## u[i] = exp(beta) = E(Y_{ij}) = mu_{ij}
           u1=rep(1.5,5),   #PL 
           u2=rep(1.5,5),   #Rx
           a=exp(-0.5), ## Y_ij | Gi=gi ~ NB(r_ij,p_i) where p_i =1/(gi*a + 1 )
           d=NULL, #delta in the AR(1) model
           othrp=list(u.n=3, s.n=0.5, p.mx=0.05, sh.mx=NA))  
                                        #parameters, for the "GN" option; 'gi's for the "NoN" option
{  #save info for gi
  Gpara=list(th=th, dist=gdist)
  
                                        #generate gi 
  if (gdist=="G")  
    { g=rgamma(n, scale=th, shape=1/th) ## E(G)= scale*shape = 1 by assumption
                                        #qg=qgamma(c((1:9)/10, (91:99)/100), shape=1/th, scale=th)
    }
  if (gdist=="N") {
                                        # g=rlnorm()   # rlnorm(u,s) u=-ss/2
                                        # EX(G) =1; Var(G)= exp(ss)-1 =th
    ss=log(th+1)    
    mu=-ss/2
    g=rlnorm(n, mu, sqrt(ss))
    Gpara=list(dist=gdist, mu=mu, s=sqrt(ss))
                                        #qg=qlnorm(c((1:9)/10, (91:99)/100), mu, sqrt(ss))
  }
  if (gdist=="U")
    {  #mean(g)=1, #var(g)=th g~U(-b,a)
      ab=nr.fun(th) #find a,b
      g=exp(runif(n, -ab[2], ab[1]))
      Gpara=list(dist=gdist, ab)
                                        #qg=exp(qunif(c((1:9)/10, (91:99)/100), -ab[2], ab[1]))
    }
  if (gdist=="GN") #mixture of G and normal
                                        # th1*p +(1-p)*th2=th; p, th, th1 specified
                                        # th2=(th-th1*p) /(1-p)
                                        # vth=sample(c(th1,th2), n, replace=T, prob=c(p, 1-p))
                                        # g=dgamma(n*p, shape=1/vth, scale=vth)}
    { othrp$sh.mx=getSH.gn(th, othrp)       #see up.index.R
      
      
      if (othrp$sh.mx  < 0)
        stop("invalid combination parameters in othrp! no positive shape parameter to satisfy E(G_i) = 1!!")

      g=rgamma(n, scale=th, shape=othrp$sh.mx)
      nm=rbinom(1, n, othrp$p.mx) 
      nm=sample(1:n, nm) 
      g.norm=rnorm(nm, othrp$u.n, othrp$s.n)
      g.norm[g.norm<0]=0
      g[nm]=g.norm
      Gpara=list(dist=gdist, othrp)
    }
  if (gdist=="NoN") g=sample(othrp, n, replace=T)

  if (is.null(d)) #independent
    {  g=rep(g, rep(sn,n))
       uu=c(rep(u1,n/2), rep(u2,n/2))
       
                                        #pb=1/(a*g+1) = size/(size+mu)
                                        #sz=uu/a where uu = mu_{ij} = E(Y_{ij})
       y=rnbinom(n*sn, mu=uu*g, size=uu/a)
     }
  else  #AR(1)
    { 
      pr=1/(g*a+1) ## p_i = 1/(g_i*alpha + 1)
                                        #check r.ij-d*r.i(j-1) > 0
      if ( any(u1[2:sn] <= d*u1[1:(sn-1)]) |
          any(u2[2:sn] <= d*u2[1:(sn-1)]) )
        {  stop("u[j] - d*u[(j-1)], j=2, ..., sn, MUST BE POSITIVE")
         }

      y=rnb.ar1(n=n/2, sn=sn, pr=pr[1:(n/2)], sz=u1/a, d=d)
      y=c(y, rnb.ar1(n=n/2, sn=sn, pr=pr[(n/2+1):n], sz=u2/a, d=d))      
    }

  id=rep(1:n, rep(sn, n))## c(1,..,1,2,...,2,3,...,n,...,n)
  vn=rep(1:sn, n) ## c(1,..,sn,1,...,sn,1,...,sn)
  gp=rep(1:2, rep(sn*n/2,2)) ## c(1,.......,1,2,.......,2)
  
  return(list(id=id, vn=vn, gp=gp, y=y, g=g, Ginfo=Gpara ##, qg=qg
              ))
}

rnb.ar1=function(n, sn, pr, sz, d)
{ # n = # of subject, sn = # of visits of each subject
                                        #pr is of length n, u is of length sn
                                        #d is a scalar of positive value

  ## generate n response count at the first time point from
  ## the negative binomial dist'n
  ## Y_i(-1)|G_i=g_i ~ NB(r_i(-1),p_i) = NB(sz,pr)
  
  g1=d*sz[1:(sn-1)]    ## alpha of beta-binomial 
  g2=sz[2:sn]-g1 ## beta of beta-binomial

  y1=rnbinom(n, size=sz[1], prob=pr)
  y=y1
  ## generate n response count at the time points > 1 
  ## Three steps:
  ## (1) Z(y_i(j-1),d) ~ beta-binom(size=y_i(j-1),
  ##                           alpha=d*r_i(j-1),beta=(1-d)*r_i(j-1))
  ## (2) e_ij ~ NB(r_ij-d*r_ij,p_i)
  ## (3) Y_ij = Z(y_i(j-1),d) + e_ij
  for (i in 1:(sn-1))
    { z1=rbb(n, y1, g1[i], g2[i])
      z2=rnbinom(n, size=g2[i], prob=pr)
      y2=z1+z2
      y=c(y, y2)
      y1=y2
    }
  tem=rep(1:n, sn)
  y=y[order(tem)] #order y according subject
  return(y)
}

getSH.gn <-
  function(th=exp(1.3), #scale
           othr=list(u.n=3, s.n=0.5, p.mx=0.05) # normal parameters: u.n = mean, s.n=SD, p.mx=proportion
           ) 
{ sha=(1-othr$p.mx*othr$u.n)/(1-othr$p.mx)/th 
  return(sha)
}


nr.fun <-
  function(th)
{ a=log(2)+log(th+1) 
  b=exp(a)
                                        #need solve f1=0 and f2=0 (see nr.pdf)
  f1=exp(a)-exp(-b)-(a+b)
  f2=exp(2*a)-exp(-2*b)-2*(a+b)*(th+1)
  
  ch=max(abs(f1), abs(f2))
  while (ch>0.0001)
    {
                                        #print(c(a,b, f1,f2) )
      df11=exp(a)-1
      df12=exp(-b)-1

      df21=2*exp(2*a)-2*th-2
      df22=2*exp(-2*b)-2*th-2

      m=matrix(c(df11,df21, df12, df22),2,2)
      v=-c(f1,f2)

      tem=solve(m,v)

      a=a+tem[1]
      b=b+tem[2]
      f1=exp(a)-exp(-b)-(a+b)
      f2=exp(2*a)-exp(-2*b)-2*(a+b)*(th+1)
      
      ch=max(abs(f1), abs(f2))
    }

                                        #print(c(f1,f2))
                                        #print(c(a,b))
  return(c(a,b))
}

