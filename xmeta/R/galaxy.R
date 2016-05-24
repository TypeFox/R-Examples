galaxy=function(data, rhow, type, method, k, L, estimator, side, maxiter){    
  galaxy.cl=function(y1, s1, y2, s2, L, estimator, side, maxiter){
    v1=s1^2
    v2=s2^2
    
    y1.org=y1
    y2.org=y2
    v1.org=v1
    v2.org=v2
    
    rho.sign=sign(cor(y1, y2))
    if(rho.sign==0){rho.sign=1}
    if (rho.sign<0){y2=-y2}
    rho.w=0
    v12=cbind(v1, rho.w*sqrt(v1*v2), v2)
    
    LC=y1+L*y2
    LC.var=v1+L^2*v2+2*rho.w*sqrt(v1*L^2*v2)
    
    LCi=LC
    LCVi=LC.var
    y1i=y1
    y2i=y2
    v1i=v1
    v2i=v2
    v12i=v12
    
    if (missing(side)){side=NULL}
    estimator = match.arg(estimator, c("L0", "R0", "Q0"))
    
    if (is.null(side)){
      bbb=rma(LC, LC.var, mods=sqrt(LC.var))$b[2]
      if (bbb < 0) {
        side = "right"
      }
      else {
        side = "left"
      }
    }
    else {side=match.arg(side, c("left", "right"))}
    
    if (side=="right"){
      LCi=-1*LCi
      y1i=-1*y1i
      y2i=-1*y2i
    }
    idix = sort(LCi, index.return = TRUE)$ix
    LCi = LCi[idix]
    LCVi = LCVi[idix]
    y1i=y1i[idix]
    y2i=y2i[idix]
    v1i=v1i[idix]
    v2i=v2i[idix]
    v12i=v12i[idix,]
    
    k = length(LCi)
    k0.sav = -1
    k0 = 0
    iter = 0
    
    junk1=mvmeta(y1i, v1i, method = "reml")
    junk2=mvmeta(y2i, v2i, method = "reml")
    res_galaxy = c(junk1$coefficients,  junk2$coefficients)   
    cov_galaxy=c(vcov(junk1), vcov(junk2))
    
    while (abs(k0 - k0.sav) > 0) {
      k0.sav = k0
      iter = iter + 1
      if (iter > maxiter) 
        stop("Galaxy algorithm did not converge.")
      LCi.t = LCi[1:(k - k0)]
      LCVi.t = LCVi[1:(k - k0)]
      
      y1.t=y1i[1:(k-k0)]
      v1.t=v1i[1:(k-k0)]
      y2.t=y2i[1:(k-k0)]
      v2.t=v2i[1:(k-k0)]
      v12.t=v12i[1:(k-k0),]
      
      #res = mvmeta(cbind(y1.t, y2.t), S=v12.t, method = "reml")
      res1 = mvmeta(y1.t, S=v1.t, method = "reml")
      res2 = mvmeta(y2.t, S=v2.t, method = "reml")
      
      b = c(res1$coefficients, res2$coefficients)
      
      LCi.c = LCi - sum(b)
      y1i.c=y1i-b[1]
      y2i.c=y2i-b[2]
      LCi.c.r = rank(abs(LCi.c), ties.method = "first")
      LCi.c.r.s = sign(LCi.c) * LCi.c.r
      if (estimator == "R0") {
        k0 = (k - max(-1 * LCi.c.r.s[LCi.c.r.s < 0])) - 1
      }
      if (estimator == "L0") {
        Sr <- sum(LCi.c.r.s[LCi.c.r.s > 0])
        k0 <- (4 * Sr - k * (k + 1))/(2 * k - 1)
        varSr <- 1/24 * (k * (k + 1) * (2 * k + 1) + 10 * 
                           k0^3 + 27 * k0^2 + 17 * k0 - 18 * k * k0^2 - 
                           18 * k * k0 + 6 * k^2 * k0)
      }
      if (estimator == "Q0") {
        Sr <- sum(LCi.c.r.s[LCi.c.r.s > 0])
        k0 <- k - 1/2 - sqrt(2 * k^2 - 4 * Sr + 1/4)
        varSr <- 1/24 * (k * (k + 1) * (2 * k + 1) + 10 * 
                           k0^3 + 27 * k0^2 + 17 * k0 - 18 * k * k0^2 - 
                           18 * k * k0 + 6 * k^2 * k0)
      }
      k0 = max(0, k0)
      k0 = round(k0)
    }
    LCi.c = LCi.c - sum(b)
    y1i.c = y1i.c-b[1]
    y2i.c = y2i.c-b[2]
    
    x.rma=rma(LCi, LCVi)
    if(k0==0){
      LCi.fill = x.rma$yi.f
      LCVi.fill = x.rma$vi.f
      
      y1.fill=y1
      v1.fill=v1
      y2.fill=y2
      v2.fill=v2
      v12.fill=v12
      
      junk1=mvmeta(y1.org, v1.org, method = "reml")
      junk2=mvmeta(y2.org, v2.org, method = "reml")
      res_galaxy=c(junk1$coefficients, junk2$coefficients)
      cov_galaxy=c(vcov(junk1), vcov(junk2))
    }
    
    if(k0>0){
      LCi.fill = c(x.rma$yi.f, -1 * LCi.c[(k - k0 + 1):k])
      LCVi.fill = c(x.rma$vi.f, LCVi[(k - k0 + 1):k])
      
      y1.fill=c(y1i, -1 * y1i.c[(k - k0 + 1):k])
      v1.fill=c(v1i, v1i[(k - k0 + 1):k])
      y2.fill=c(y2i, -1 * y2i.c[(k - k0 + 1):k])
      v2.fill=c(v2i, v2i[(k - k0 + 1):k])
      
      junk1=mvmeta(y1.fill, v1.fill, method = "reml")
      junk2=mvmeta(y2.fill, v2.fill, method = "reml")
      res_galaxy=c(junk1$coefficients, junk2$coefficients)
      if (side=="right"){res_galaxy=-res_galaxy}
      if (rho.sign<0){res_galaxy[2]=-res_galaxy[2]} 
      cov_galaxy=c(vcov(junk1), vcov(junk2))
    }
    return(list(res_galaxy=res_galaxy, cov_galaxy=cov_galaxy, k0=k0, side=side, method=method, type=type))
  }
  
  if (missing(data)){data=NULL}
  if (is.null(data)){
    stop("The dataset must be specified.")
  }
  
  if (missing(k)){k=NULL}
  if (is.null(k)){
    stop("The number of outcomes must be specified.")
  }  
  
  if (missing(type)){type=NULL}
  if (is.null(type)){
    stop('The type of outcome must be specified. Please input "continuous" or "binary". ')
  }      
  
  if (type=="binary"){
    stop("The method for meta-analysis with binary outcome is currently
         under development. ")
  }
  
  if (missing(method)){method=NULL}
  if (is.null(method)){
    stop("A method must be specified. ")
  }
  if (missing(rhow)){rhow=NULL}
  if ((is.null(rhow))&(method=="nn.reml")){
    stop('The within-study correlations are required for the input method "nn.reml".')
  }
  if (k>2){
    stop('The method for MMA with more than 2 outcomes are currently under development. ')
  }
    
  if (type=="continuous"){
    y1=data$y1
    s1=data$s1
    y2=data$y2
    s2=data$s2
    
    if(method=="galaxy.cl"){
      res=galaxy.cl(y1, s1, y2, s2, L, estimator, side, maxiter)}
  
    if ((method%in%c("galaxy.cl"))!=1){
      stop("The input method is not available for continuous data")
    }
}
class(res) = c("galaxy")
return(res)
}




