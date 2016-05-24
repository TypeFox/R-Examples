# R commands: slice sampling for the truncated multivariate t (TMVT) distribution
TT.GS = function(n, mu=rep(0,nrow(S)), S=diag(length(mu)), nu=2, lower=rep(-Inf, length(mu)), upper=rep(Inf, length(mu)))

{
  
#  require(mvtnorm)
  
  p=length(mu)

  ## Verify error at parameters specification
  if(length(lower) != length(upper)) stop("The lengths of the 'lower' and 'upper' truncated values must be equal!")
  for(i in 1:length(lower))
  {
    if(upper[i] <= lower[i]) stop("The lower limit is larger than upper limit for truncation!")
  }
  if(length(lower) != nrow(S) | length(upper) != nrow(S) ) stop("The dimension of scale-covariance matrix must be equal to the length of the 'lower/upper' truncated value!")
  if(nu <=0) stop("The degree of freedom must be larger than zero!")
  if(det(S)<=0) stop("The S matrix must be inversible!")
  
  TT.GS.sp = function(n,R,nu,lower,upper)
    
  {
        
    #initial value
    x=qt(runif(rep(1,p),pt(lower,df=nu),pt(upper,df=nu)),df=nu)
    
    if(n<1) return(t(x))    
    
    X = matrix(NA, n, p)
    
    R.inv = solve(R)
    
    for(i in 1:n)
      
    {
      
      delta = sum(colSums(x*R.inv)*x)
      
      y = runif(1,0,exp(-.5*(nu+p)*log(1+delta/nu)))
      
      kap = nu*(y^(-2/(nu+p))-1)
      
      for(j in 1:p)
        
      {
        
        ss = x[-j]%*%R.inv[-j,-j]%*%x[-j]
        
        mj = - sum(R.inv[-j,j]*x[-j]) / R.inv[j,j]
        
        tj = sqrt(mj^2 + (kap - ss) / R.inv[j,j])
        
        xij = runif(1,max(lower[j],mj-tj),min(upper[j],mj+tj))
        
        X[i,j] = xij
        
        x[j] = xij
        
      }
      
    }
    
    return(X)
    
  }
  
  s = sqrt(diag(S))
  
  R = S/outer(s,s,"*")
  
  Z = TT.GS.sp(n,R,nu,lower=(lower-mu)/s,upper=(upper-mu)/s)
  
  X = t(mu + t(Z) * s)
  
  return(X)
  
}

# R commands: calculation of the first two moments of the TMVT distribution
TT.moment = function(R=diag(length(lower)), nu=5, lower=rep(-Inf, nrow(R)), upper=rep(Inf, nrow(R)))
    
{
  
    ## Verify error at parameters specification
    if(length(lower) != length(upper)) stop("The lengths of the 'lower' and 'upper' truncated values must be equal!")
    for(i in 1:length(lower))
    {
      if(upper[i] <= lower[i]) stop("The lower limit is larger than upper limit for truncation!")
    }
    if(length(lower) != nrow(R)) stop("The dimension of scale-covariance matrix must be equal to the length of the 'lower/upper' truncated value!")
    if(nu <= 2) stop("The first moment exists only when the degree of freedom is larger than 2!")
    if(det(R)<=0) stop("The R matrix must be inversible!")
    
    nu=ifelse(nu<3,3,round(nu))
    
#    require(mvtnorm)
  
  GB = GenzBretz(maxpts = 5e4, abseps = 1e-9, releps = 0)
  a = lower; b = upper
  p = length(a)
  
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
  
  for(r in 1:p)
    
  {
    
    temp = R[-r,r]
    
    S1 = R[-r,-r] - temp %*% t(R[r,-r])
    
    mua = temp * a[r]; low = a[-r]-mua; upp = b[-r]-mua
    
    G1a[r] = ifelse(p==2,pt(upp/sqrt(S1/da[r]),df=nu-1)-pt(low/sqrt(S1/da[r]),df=nu-1)
                    
                    ,pmvt(lower = low, upper = upp, sigma = S1/da[r], df = nu-1, algorithm = GB)[1])
    
    mub = temp * b[r]; low = a[-r]-mub; upp = b[-r]-mub
    
    G1b[r] = ifelse(p==2,pt(upp/sqrt(S1/db[r]),df=nu-1)-pt(low/sqrt(S1/db[r]),df=nu-1)
                    
                    ,pmvt(lower = low, upper = upp, sigma = S1/db[r], df = nu-1, algorithm = GB)[1])
    
  }
  
  qa = f1a*G1a; qb = f1b*G1b
  
  EX = c(R %*% (qa-qb)) / al0 / la1
  
  if(nu>4){
  
  H = matrix(0,p,p)
  
  for(r in 1:(p-1))
    
  {
    
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
  
  H = H / la2
  
  D = matrix(0,p,p)
  
  diag(D) = a * qa - b * qb - diag(R%*%H)
  
  al1 = pmvt(lower = a, upper = b, sigma = R/la1, df=nu-2, algorithm = GB)[1]
  
  EXX = (al1 * R + R %*% (H + D) %*% R) / al0 / la1
  } else {
    EXX=matrix(NA,p,p)
    cat('Warning message:','\n')
    cat('The theoretical second moment exists only when the degrees of freedom is larger than 4!')
  }
  
  return(list(EX=EX,EXX=EXX))
  
}
