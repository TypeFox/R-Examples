BiVarEst <- function(Z, kopule="gumbel", rodiny=c("weibull","weibull"))
{
  gammafit <- function(time, event) {
    gammaLik <- function(x) {
      ie <- which(event==0)
      it <- which(event==1)
      if(length(it)>0) dg <- dgamma(time[it], shape=x[1], scale=x[2], log=TRUE) else {
        print("count of noncensored times is 0")
      }
      if(length(ie)>0) sg <- pgamma(time[ie],shape=x[1],scale=x[2],lower.tail=FALSE,log.p=TRUE) else sg = 0
#      print(paste("shape =",x[1],"  scale =",x[2],"  gammaLik =",-(sum(dg)+sum(sg))))
      return(-(sum(dg) + sum(sg)))
    }
    no = length(time)
    s = log(sum(time)/no) - sum(log(time))/no
    sx = (3 - s + sqrt((s-3)*(s-3) + 24*s))/(12*s)
    xx = c(sx, max(time)/2)
#    print("x0:"); print(xx)
    return(nmkb(xx, gammaLik, lower=c(0.01,0.1), upper=c(sx*10, max(time))))
  }

  LiBiCop <- function(x)
  {
#  print(x)
  Lik = 0
  d11 <- which(Z[,3]==1 & Z[,4]==1)
  d00 <- which(Z[,3]==0 & Z[,4]==0)
  d01 <- which(Z[,3]==0 & Z[,4]==1)
  d10 <- which(Z[,3]==1 & Z[,4]==0)
  if(length(d11)>0){
    if (rodiny[1] == "weibull") {
      du = dweibull(Z[,1][d11],shape=x[1],scale=x[2])
      u = pweibull(Z[,1][d11],shape=x[1],scale=x[2])
    } else {if (rodiny[1] == "lnorm") {
      du = dlnorm(Z[,1][d11],meanlog=x[1],sdlog=x[2])
      u = plnorm(Z[,1][d11],meanlog=x[1],sdlog=x[2])
    } else {if(rodiny[1]=="norm") {
      du = dnorm(Z[,1][d11],mean=x[1],sd=x[2])
      u = pnorm(Z[,1][d11],mean=x[1],sd=x[2])
    } else {if(rodiny[1]=="gamma"){
      du = dgamma(Z[,1][d11],shape=x[1],scale=x[2])
      u = pgamma(Z[,1][d11],shape=x[1],scale=x[2])
    } else {print("distribution = ?")}}}}
    # v
    if (rodiny[2] == "weibull") {
      dv = dweibull(Z[,2][d11],shape=x[3],scale=x[4])
      v = pweibull(Z[,2][d11],shape=x[3],scale=x[4])
    } else {if (rodiny[2] == "lnorm") {
      dv = dlnorm(Z[,2][d11],meanlog=x[3],sdlog=x[4])
      v = plnorm(Z[,2][d11],meanlog=x[3],sdlog=x[4])
    } else {
      if (rodiny[2]=="norm") {
        dv = dnorm(Z[,2][d11],mean=x[3],sd=x[4])
        v = pnorm(Z[,2][d11],mean=x[3],sd=x[4])
      } else {if(rodiny[2]=="gamma"){
        dv = dgamma(Z[,2][d11],shape=x[3],scale=x[4])
        v = pgamma(Z[,2][d11],shape=x[3],scale=x[4])
      } else {print("distribution = ?")}}}}
    a = x[length(x)]
    if (kopule=="gumbel") {
      # have to be a >= 1
      if (a < 1) a <- 1.01
      gu = -log(u)
      gv = -log(v)
      gua = gu^a
      gva = gv^a
      guv = gua+gva
      guv2 = guv^(2/a-2)
      guv1 = guv^(1/a-2)
      exa = exp(-(guv)^(1/a))
      cip = gua*du*exa*gva*dv
      cid = guv2 - (1-a)*guv1 # watch out for cid < or = 0
      Lik = Lik + sum(log((cip*cid)/(u*v*gu*gv)))
    } else {
      if (kopule=="clayton") {
          Lik = Lik + sum(log((1+a) * (u*v)^(-a-1) * du * dv * (1/v^a + 1/u^a -1)^(-1/a-2)))
      } else {
        if (kopule=="frank") {
 #         print(paste("a =",a))
          if (abs(a) < 0.01) a <- a + sign(a)*0.01
          aux = exp(-a*u)-1
          avy = exp(-a*v)-1
          axy = exp(-a*(u+v))
          ea = exp(-a)-1
          me = aux*avy/ea+1
          Lik = Lik + sum(log(a*du*dv*axy*(aux*avy/(ea*ea*me*me)-1/(ea*me))))
#          print(paste("d11:",Lik))
        } else {print("kopule = ?")}}}
  } # d11 > 0
  if(length(d00)>0){
    if (rodiny[1] == "weibull") {
      u = pweibull(Z[,1][d00],shape=x[1],scale=x[2],lower.tail=FALSE)
    } else {if (rodiny[1] == "lnorm") {
      u = plnorm(Z[,1][d00],meanlog=x[1],sdlog=x[2],lower.tail=FALSE)
    } else {if(rodiny[1]=="norm") {
      u = pnorm(Z[,1][d00],mean=x[1],sd=x[2],lower.tail=FALSE)
    } else {if(rodiny[1]=="gamma"){
      u = pgamma(Z[,1][d00],shape=x[1],scale=x[2],lower.tail=FALSE)
    } else {print("distribution = ?")}}}}
    # v
    if (rodiny[2] == "weibull") {
      v = pweibull(Z[,2][d00],shape=x[3],scale=x[4],lower.tail=FALSE)
    } else {if (rodiny[2] == "lnorm") {
      v = plnorm(Z[,2][d00],meanlog=x[3],sdlog=x[4],lower.tail=FALSE)
    } else {
      if (rodiny[2]=="norm") {
        v = pnorm(Z[,2][d00],mean=x[3],sd=x[4],lower.tail=FALSE)
      } else {if(rodiny[2]=="gamma"){
        v = pgamma(Z[,2][d00],shape=x[3],scale=x[4],lower.tail=FALSE)
      } else {print("distribution = ?")}}}}
    a = x[length(x)]
    if (kopule=="gumbel") {
      gua = (-log(u))^a
      gva = (-log(v))^a
      Lik = Lik - sum((gva+gua)^(1/a))
    } else {
      if (kopule=="clayton") {
        for (iu in 1:length(u)) {
          Lik = Lik + sum(-log((1/v[iu])^a + (1/u[iu])^a -1)/a)
        }
      } else {
        if (kopule=="frank") {
            aux = exp(-a*u)-1
            avy = exp(-a*v)-1
            ea = exp(-a)-1
            Lik = Lik + sum(log(-log(((aux*avy)/ea)+1)/a))
  #          print(paste("d00:",Lik))
        } else {print("kopule = ?")}}}
  } # d00 > 0
  if(length(d01)>0){
    if (rodiny[1] == "weibull") {
      du = dweibull(Z[,1][d01],shape=x[1],scale=x[2])
      u = pweibull(Z[,1][d01],shape=x[1],scale=x[2],lower.tail=FALSE)
    } else {if (rodiny[1] == "lnorm") {
      du = dlnorm(Z[,1][d01],meanlog=x[1],sdlog=x[2])
      u = plnorm(Z[,1][d01],meanlog=x[1],sdlog=x[2],lower.tail=FALSE)
    } else {if(rodiny[1]=="norm") {
      du = dnorm(Z[,1][d01],mean=x[1],sd=x[2])
      u = pnorm(Z[,1][d01],mean=x[1],sd=x[2],lower.tail=FALSE)
    } else {if(rodiny[1]=="gamma"){
      du = dgamma(Z[,1][d01],shape=x[1],scale=x[2])
      u = pgamma(Z[,1][d01],shape=x[1],scale=x[2],lower.tail=FALSE)
    } else {print("distribution = ?")}}}}
    # v
    if (rodiny[2] == "weibull") {
      dv = dweibull(Z[,2][d01],shape=x[3],scale=x[4])
      v = pweibull(Z[,2][d01],shape=x[3],scale=x[4],lower.tail=FALSE)
    } else {if (rodiny[2] == "lnorm") {
      dv = dlnorm(Z[,2][d01],meanlog=x[3],sdlog=x[4])
      v = plnorm(Z[,2][d01],meanlog=x[3],sdlog=x[4],lower.tail=FALSE)
    } else {
      if (rodiny[2]=="norm") {
        dv = dnorm(Z[,2][d01],mean=x[3],sd=x[4])
        v = pnorm(Z[,2][d01],mean=x[3],sd=x[4],lower.tail=FALSE)
      } else {if(rodiny[2]=="gamma"){
        dv = dgamma(Z[,2][d01],shape=x[3],scale=x[4])
        v = pgamma(Z[,2][d01],shape=x[3],scale=x[4],lower.tail=FALSE)
      } else {print("distribution = ?")}}}}
    a = x[length(x)]
    if (kopule=="gumbel") {
        gu = -log(u)
        gv = -log(v)
        gua = gu^a
        gva = gv^a
        guv = gua+gva
        guv1 = guv^(1/a-1)
        exa = exp(-(guv)^(1/a))
        Lik = Lik + sum(log((exa*guv1*(-gva*dv))/(v*log(v))))
    } else {
      if (kopule=="clayton") {
          Lik = Lik + sum(log((v^(-a-1) * dv * (1/v^a + 1/u^a -1)^(-1/a-1))))
      } else {
        if (kopule=="frank") {
            aux = exp(-a*u)-1
            avy = exp(-a*v)-1
            ea = exp(-a)-1
            me = aux*avy/ea+1
            Lik = Lik + sum(log((aux*(dv*exp(-a*v)))/(ea*me)))
  #          print(paste("d01:",Lik))
        } else {print("kopule = ?")}}}
  } # d01 > 0
  if(length(d10)>0){
    if (rodiny[1] == "weibull") {
      du = dweibull(Z[,1][d10],shape=x[1],scale=x[2])
      u = pweibull(Z[,1][d10],shape=x[1],scale=x[2],lower.tail=FALSE)
    } else {if (rodiny[1] == "lnorm") {
      du = dlnorm(Z[,1][d10],meanlog=x[1],sdlog=x[2])
      u = plnorm(Z[,1][d10],meanlog=x[1],sdlog=x[2],lower.tail=FALSE)
    } else {if(rodiny[1]=="norm") {
      du = dnorm(Z[,1][d10],mean=x[1],sd=x[2])
      u = pnorm(Z[,1][d10],mean=x[1],sd=x[2],lower.tail=FALSE)
    } else {if(rodiny[1]=="gamma"){
      du = dgamma(Z[,1][d10],shape=x[1],scale=x[2])
      u = pgamma(Z[,1][d10],shape=x[1],scale=x[2],lower.tail=FALSE)
    } else {print("distribution = ?")}}}}
    # v
    if (rodiny[2] == "weibull") {
      dv = dweibull(Z[,2][d10],shape=x[3],scale=x[4])
      v = pweibull(Z[,2][d10],shape=x[3],scale=x[4],lower.tail=FALSE)
    } else {if (rodiny[2] == "lnorm") {
      dv = dlnorm(Z[,2][d10],meanlog=x[3],sdlog=x[4])
      v = plnorm(Z[,2][d10],meanlog=x[3],sdlog=x[4],lower.tail=FALSE)
    } else {
      if (rodiny[2]=="norm") {
        dv = dnorm(Z[,2][d10],mean=x[3],sd=x[4])
        v = pnorm(Z[,2][d10],mean=x[3],sd=x[4],lower.tail=FALSE)
      } else {if(rodiny[2]=="gamma"){
        dv = dgamma(Z[,2][d10],shape=x[3],scale=x[4])
        v = pgamma(Z[,2][d10],shape=x[3],scale=x[4],lower.tail=FALSE)
      } else {print("distribution = ?")}}}}
    a = x[length(x)]
    if (kopule=="gumbel") {
        gu = -log(u)
        gv = -log(v)
        gua = gu^a
        gva = gv^a
        guv = gua+gva
        guv1 = guv^(1/a-1)
        exa = exp(-(guv)^(1/a))
        Lik = Lik +  sum(log((exa*guv1*(-gua*du))/(u*log(u))))
    } else {
      if (kopule=="clayton") {
          Lik = Lik + sum(log((u^(-a-1) * du * (1/v^a + 1/u^a -1)^(-1/a-1))))
      } else {
        if (kopule=="frank") {
            aux = exp(-a*u)-1
            avy = exp(-a*v)-1
            ea = exp(-a)-1
            me = aux*avy/ea+1
 #           print(paste("d10:",Lik))
            Lik = Lik + sum(log((avy*(du*exp(-a*u)))/(ea*me)))
        } else {print("kopule = ?")}}}
  } # d10 > 0
#  print(paste("Lik=",Lik))
  return(-Lik)
 }
  # starting value
  fam = c(NA,NA)
  for (i in 1:2) {
    if(rodiny[i]=="weibull") fam[i] = "weibull" else {
      if(rodiny[i]=="norm") fam[i] = "gaussian" else {
        if(rodiny[i]=="lnorm") fam[i] = "lognormal" else fam[i] = "gamma"}}
  }
  if (fam[1]!="gamma") wx = survreg(Surv(Z[,1],Z[,3])~1,dist=fam[1]) else wx = gammafit(Z[,1],Z[,3])
  if (fam[2]!="gamma") wy = survreg(Surv(Z[,2],Z[,4])~1,dist=fam[2]) else wy = gammafit(Z[,2],Z[,4])
  # prs
  if(fam[1]=="weibull") {ax0 = exp(wx$coefficients); bx0 = 1/wx$scale} else {
if(fam[1]=="gamma") {ax0=wx$par[2]; bx0 = wx$par[1]} else {
  if(fam[1]=="gaussian") {bx0 = wx$coefficients; ax0 = wx$scale} else {bx0 = wx$coefficients; ax0 = wx$scale}}}
  if(fam[2]=="weibull") {ay0 = exp(wy$coefficients); by0 = 1/wy$scale} else {
if(fam[2]=="gamma") {ay0=wy$par[2]; by0 = wy$par[1]} else {
if(fam[2]=="gaussian") {by0 = wy$coefficients; ay0 = wy$scale} else {by0 = wy$coefficients; ay0 = wy$scale}}}
  af = 2
  prs = c(bx0, ax0, by0, ay0, af)
  Lw = prs/5
  Upp = prs*10
  if(fam[1]=="gaussian") {Lw[1] = min(prs[1] - 3*prs[2], Z[, 1]); Upp[1] = max(-prs[1] + 3*prs[2], Z[, 1])}
  if(fam[2]=="gaussian") {Lw[3] = min(prs[3] - 3*prs[4], Z[, 2]); Upp[3] = max(-prs[3] + 3*prs[4], Z[, 2])}
  if(fam[1]=="lognormal" & prs[1] < 0) {Lw[1] = prs[1] - 3*prs[2]; Upp[1] = -prs[1] + 3*prs[2]}
  if(fam[2]=="lognormal" & prs[3] < 0) {Lw[3] = prs[3] - 3*prs[4]; Upp[3] = -prs[3] + 3*prs[4]}
  if(kopule=="frank") {Lw[5] = -af*10; Upp[5]=af*10}
 # print(Lw); print(prs); print(Upp)
#  Z = Nxyd
 return(nmkb(par=prs,fn=LiBiCop,lower=Lw,upper=Upp))
}
