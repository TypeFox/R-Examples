################################################################################
######################### Modelos T Mistos Censurados ##########################
################################################################################


################################################################################
#################### Momentos da Distribuição t Truncada #######################
################################################################################

#GB = GenzBretz(maxpts = 1e4, abseps = 1e-9, releps = 0)

cdfNI<-function(x,mu,sigma2,nu,type="Normal"){
resp<-matrix(0,length(x),1)


if(type=="Normal"){
resp<-pnorm(x,mu,sqrt(sigma2))
                  }

if(type=="T"){
z=(x-mu)/sqrt(sigma2)
resp=pt(z,df=nu)
             }
return(resp)
}


################################################################################
#################### Momentos da Distribuição t Truncada #######################
################################################################################

TT.moment = function(a,b,R,nu)
{

  require(mvtnorm)

  GB = GenzBretz(maxpts = 5e4, abseps = 1e-9, releps = 0)

  p = length(a)

  if(p==1){

  if(a== -Inf)  a <- -1e12
  if(b== Inf)   b <- 1e12

    G1<- 0.5*(gamma((nu-1)/2)*nu^(nu/2))/((cdfNI(b,0,1,nu,"T")-cdfNI(a,0,1,nu,"T"))*gamma(nu/2)*gamma(1/2))
    EX<- ((G1*((nu+a^2)^(-(nu-1)/2)-(nu+b^2)^(-(nu-1)/2))))
    EXX<- nu/(nu-2)+(G1*(a*(nu+a^2)^(-(nu-1)/2)-b*(nu+b^2)^(-(nu-1)/2)))

          }

  else{

  a = ifelse(a==-Inf,rep(-1e12,p),a)

  b = ifelse(b== Inf,rep( 1e12,p),b)

  al0 = pmvt(lower = a, upper = b, sigma = R, df = nu, algorithm = GB)[1]

  ### pdf & cdf

  la1 = (nu-2)/nu; la2 = (nu-4)/nu

  da = (nu-1)/(nu+a^2); db = (nu-1)/(nu+b^2)

  f1a = sqrt(la1)*dt(sqrt(la1)*a,df=nu-2)

  f1b = sqrt(la1)*dt(sqrt(la1)*b,df=nu-2)

  f2 = matrix(NA, p, p)

  G1a = G1b = rep(NA, p)

  G2 = matrix(NA, p, p)



  H = matrix(0,p,p)

  for(r in 1:(p-1))

  {

    temp = R[-r,r]

    S1 = R[-r,-r] - temp %*% t(R[r,-r])

    mua = temp * a[r]; low = a[-r]-mua; upp = b[-r]-mua

    G1a[r] = ifelse(p==2,pt(upp/sqrt(S1/da[r]),df=nu-1)-pt(low/sqrt(S1/da[r]),df=nu-1)

        ,pmvt(lower = low, upper = upp, sigma = S1/da[r], df = nu-1, algorithm = GB)[1])

    mub = temp * b[r]; low = a[-r]-mub; upp = b[-r]-mub

    G1b[r] = ifelse(p==2,pt(upp/sqrt(S1/db[r]),df=nu-1)-pt(low/sqrt(S1/db[r]),df=nu-1)

        ,pmvt(lower = low, upper = upp, sigma = S1/db[r], df = nu-1, algorithm = GB)[1])

    for(s in (r+1):p)

    {

      rs = c(r,s)

      pdf.aa = dmvt(c(a[r],a[s]),sigma=R[rs,rs]/la2,df=nu-4, log =F)

      pdf.ab = dmvt(c(a[r],b[s]),sigma=R[rs,rs]/la2,df=nu-4, log =F)

      pdf.ba = dmvt(c(b[r],a[s]),sigma=R[rs,rs]/la2,df=nu-4, log =F)

      pdf.bb = dmvt(c(b[r],b[s]),sigma=R[rs,rs]/la2,df=nu-4, log =F)

      if(p==2){cdf.aa=cdf.ab=cdf.ba=cdf.bb=1}

      if(p>2)

      {

        tmp = R[-rs,rs]%*%solve(R[rs,rs])

        mu.aa = c(tmp%*%c(a[r],a[s]))

        mu.ab = c(tmp%*%c(a[r],b[s]))

        mu.ba = c(tmp%*%c(b[r],a[s]))

        mu.bb = c(tmp%*%c(b[r],b[s]))

        daa = (nu-2)/(nu+(a[r]^2-2*R[r,s]*a[r]*a[s]+a[s]^2)/(1-R[r,s]^2))

        dab = (nu-2)/(nu+(a[r]^2-2*R[r,s]*a[r]*b[s]+b[s]^2)/(1-R[r,s]^2))

        dba = (nu-2)/(nu+(b[r]^2-2*R[r,s]*b[r]*a[s]+a[s]^2)/(1-R[r,s]^2))

        dbb = (nu-2)/(nu+(b[r]^2-2*R[r,s]*b[r]*b[s]+b[s]^2)/(1-R[r,s]^2))

        R21 = R[-rs,-rs] - R[-rs,rs]%*%solve(R[rs,rs]) %*% R[rs,-rs]

        cdf.aa = ifelse(p==3,pt((b[-rs]-mu.aa)/sqrt(R21/daa),df=nu-2)-pt((a[-rs]-mu.aa)/sqrt(R21/daa),df=nu-2)

            ,pmvt(lower = a[-rs]-mu.aa, upper = b[-rs]-mu.aa, sigma = R21/daa, df=nu-2, algorithm = GB)[1])

        cdf.ab = ifelse(p==3,pt((b[-rs]-mu.ab)/sqrt(R21/dab),df=nu-2)-pt((a[-rs]-mu.ab)/sqrt(R21/dab),df=nu-2)

            ,pmvt(lower = a[-rs]-mu.ab, upper = b[-rs]-mu.ab, sigma = R21/dab, df=nu-2, algorithm = GB)[1])

        cdf.ba = ifelse(p==3,pt((b[-rs]-mu.ba)/sqrt(R21/dba),df=nu-2)-pt((a[-rs]-mu.ba)/sqrt(R21/dba),df=nu-2)

            ,pmvt(lower = a[-rs]-mu.ba, upper = b[-rs]-mu.ba, sigma = R21/dba, df=nu-2, algorithm = GB)[1])

        cdf.bb = ifelse(p==3,pt((b[-rs]-mu.bb)/sqrt(R21/dbb),df=nu-2)-pt((a[-rs]-mu.bb)/sqrt(R21/dbb),df=nu-2)

            ,pmvt(lower = a[-rs]-mu.bb, upper = b[-rs]-mu.bb, sigma = R21/dbb, df=nu-2, algorithm = GB)[1])

      }

      H[r,s] = H[s,r] = pdf.aa*cdf.aa - pdf.ab*cdf.ab - pdf.ba*cdf.ba + pdf.bb*cdf.bb

    }

  }

  ##last part do loop
  r <- p
  temp = R[-r,r]

  S1 = R[-r,-r] - temp %*% t(R[r,-r])

  mua = temp * a[r]; low = a[-r]-mua; upp = b[-r]-mua

  G1a[r] = ifelse(p==2,pt(upp/sqrt(S1/da[r]),df=nu-1)-pt(low/sqrt(S1/da[r]),df=nu-1)

        ,pmvt(lower = low, upper = upp, sigma = S1/da[r], df = nu-1, algorithm = GB)[1])

  mub = temp * b[r]; low = a[-r]-mub; upp = b[-r]-mub

  G1b[r] = ifelse(p==2,pt(upp/sqrt(S1/db[r]),df=nu-1)-pt(low/sqrt(S1/db[r]),df=nu-1)

        ,pmvt(lower = low, upper = upp, sigma = S1/db[r], df = nu-1, algorithm = GB)[1])

  qa = f1a*G1a; qb = f1b*G1b

  EX = c(R %*% (qa-qb)) / al0 / la1

  H = H / la2

  D = matrix(0,p,p)

  diag(D) = a * qa - b * qb - diag(R%*%H)

  al1 = pmvt(lower = a, upper = b, sigma = R/la1, df=nu-2, algorithm = GB)[1]

  EXX = (al1 * R + R %*% (H + D) %*% R) / al0 / la1

  }

return(list(EX=EX,EXX=EXX))

}


##Calculate the first to moments when mu not 0 and Sigma not R

Mtmvt <- function(mu,Sigma,nu,lower,upper){

p=length(lower)

if(p==1){

  if(lower== -Inf)  lower <- -1e12
  if(upper== Inf)  upper <- 1e12

  a1<-(lower-mu)/sqrt(Sigma)
  b1<-(upper-mu)/sqrt(Sigma)
  M <- TT.moment(a1, b1, 1, nu)
  Ey<- mu+sqrt(Sigma)*M$EX
  Eyy<- mu^2+Sigma*M$EXX+2*mu*sqrt(Sigma)*M$EX
  Vary<- Eyy - Ey^2
}

else{
          Lambda <- diag(1/sqrt(diag(Sigma)))

   if(length(which(upper == Inf)) != 0)  upper[which(upper == Inf)] <- 1e12
   b <- as.vector(diag(1/sqrt(diag(Sigma))) %*% (upper - mu))

   if(length(which(lower == -Inf)) != 0) lower[which(lower == -Inf)] <- -1e12
   a <- as.vector(diag(1/sqrt(diag(Sigma))) %*% (lower - mu))

   R <- Lambda %*% Sigma %*% Lambda

   M <- TT.moment(a, b, R, nu)

   Ey <- mu + solve(Lambda) %*% M$EX

   Eyy <- mu %*% t(mu) + solve(Lambda) %*% M$EX %*% t(mu) + mu %*% t(M$EX) %*% solve(Lambda) + solve(Lambda) %*% M$EXX %*% solve(Lambda)
   Vary<- Eyy- Ey%*%t(Ey)
}

   return(list(Ey=Ey,Eyy=Eyy,Vary=Vary))
}

