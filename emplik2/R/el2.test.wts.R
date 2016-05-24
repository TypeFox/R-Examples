el2.test.wts <- function (u,v,wu,wv,mu0,nu0,indicmat,mean) {

#If mean is not a scalar then stop
    if (length(mean) != 1) 
        stop("mean must be a scalar")

#Calculate scalars to be used in calculations
    sumwu <- sum(wu)
    sumwv <- sum(wv)
    nu <- length(u)
    nv <- length(v)
  
#Calculate matrix and vectors to be used in calculations
    indic4mu <- nu0 %*% t(indicmat) 
    indic4nu <- mu0 %*% indicmat  
#Calculate the delta for lam search
    du <- 0.02 * sumwu/abs(sum(indic4mu))
    dv <- 0.02 * sumwv/abs(sum(indic4nu))
    dd <- min(du,dv)

#Define lamfun, where lamfun(true lam) = 0
    lamfun <- function(lam, wu, wv, sumwu, sumwv, indic4mu,
      indic4nu, indicmat) {
    mu <- wu/(sumwu+lam*indic4mu)
    nu <- wv/(sumwv+lam*indic4nu)
    return(mu %*% indicmat %*% t(nu))
      }

#Find upper and lower bounds on lam, for uniroot
    if (lamfun(0, wu, wv, sumwu, sumwv, indic4mu,
      indic4nu,indicmat) == 0)
      lam0 <- 0 else 

    {if (lamfun(0, wu, wv, sumwu, sumwv, indic4mu,
      indic4nu,indicmat) > 0)
       {lo <- 0
       up <- dd 
       while (lamfun(up, wu, wv, sumwu, sumwv, indic4mu,
         indic4nu,indicmat) > 0)
         {up <- up + dd} } else
       
           {up <- 0
           lo <- -dd  
           while (lamfun(lo, wu, wv, sumwu, sumwv, indic4mu,
             indic4nu,indicmat) < 0)
             {lo <- lo - dd}}}

#Find lam using uniroot
 lam <- uniroot(lamfun, lower = lo, upper = up,
    tol = 1e-09, wu=wu, wv=wv, sumwu=sumwu, sumwv=sumwv,
    indic4mu=indic4mu, indic4nu=indic4nu, indicmat=indicmat)$root

#Calculate updated mu1,nu1 using the lagrangian lam
    mu1 <- wu/(sumwu + lam * nu0 %*% t(indicmat))
    nu1 <- wv/(sumwv + lam * mu0 %*% indicmat)

#list the original data & weights plus the p_i, lam0, and mean
    list(u=u, wu=wu, jumpu=mu1, v=v, wv=wv, jumpv=nu1, lam=lam, mean)

}