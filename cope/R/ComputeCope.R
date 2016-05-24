#' Compute CoPE sets.
#' 
#' Computes CoPE sets for the data Y using the algorithm from Sommerfeld, Sain
#' and Schwartzman (2015).
#'
#' @param Z A list with components "x", "y" and "z". Here, x and y are the 
#'          coordinates of the grid on which the data is observed and z is an 
#'          array with dimensions c(length(x),length(y),n), containing the
#'          data. n is the number of observations.
#' @param level The level of interest.
#' @param X The design matrix of the linear model. If NULL, it is set to 
#'          matrix(rep(1,dim(Y)[3]),ncol=1) corresponding to i.i.d. data.
#' @param tau The scaling factor ||v||_2*pi_n^(1/2). The default 1/sqrt(n) 
#'            corresponds to i.i.d. data.
#' @param alpha The significance level. Inclusion for CoPE sets holds with 
#'              probability 1-alpha.
#' @param N The number of bootstrap realizations to generate for determining
#'          the threshold.
#' @param mu The true parameter function. Use the default NULL if unknown. 
#' @param mask Pixels outside the mask (i.e. where mask is ==NA) are ignored.
#'   
#' @return An object of class cope. A list containing the following
#' \itemize{
#'  \item{x, y: }{The grid coordinates from the input.}
#'  \item{mu, level, tau, X, N, alpha, mask: }
#'  {The corresponding values from the input.}
#'  \item{mu_hat, norm_est: }{The estimatot for mu and its normalized version.}
#'  \item{a_MB, a_MB_true, a_Tay, a_Tay_true: }{Thresholds for the CoPE sets 
#'  determined using the multiplier bootstrap or Taylor's method and the 
#'  estimated or the true contour, respectively.}
#'  \item{incl_MB, incl_MB_true, incl_Tay, incl_Tay_true: }{Booleans indicating 
#'  whether inclusion of the excursion set in the upper CoPE set and inclusion
#'  of the lower CoPE set in the excursion set holds, when CoPE sets are
#'  determined by a_MB, a_MB_true, a_Tay or a_Tay_true, respectively. Only 
#'  available if mu is given.}
#' }
#' 
#' @export
#' 
#' @references M. Sommerfeld, S. Sain and A. Schwartzman. Confidence regions for 
#'             excursion sets in asymptotically Gaussian
#'             random fields, with an application to climate. Preprint, 2015. 
#'        
#' @examples
#' # An example using the ToyNoise and ToySignal of this package.
#' n = 30
#' Data = ToyNoise1(n = n)
#' Data$z = Data$z + rep(ToySignal()$z, n)
#' CopeSet = ComputeCope(Data,level=4/3, mu=ToySignal()$z)
#' PlotCope(CopeSet)
ComputeCope = function(Z,level,
                X=NULL,
                tau=1/sqrt(dim(Z$z)[3]),
                alpha=0.1,
                N=1000,
                mu=NULL,
                mask=NULL){
  x = Z$x
  y = Z$y
  Y = Z$z
  n = dim(Y)[3]
  #Compute the estimator beta_hat.
  
  
  
  if(is.null(X)){
    beta_hat = array(0, c(length(x),length(y), 1))
    beta_hat[,,1] = apply(Y,1:2,mean)
    R = Y - rep(beta_hat[,,1],n)
  } else{
    beta_hat = array(0, c(length(x),length(y), ncol(X)))
    M = MASS::ginv(t(X) %*% X) %*% t(X)
    for(i in 1:length(x))
      for(j in 1:length(y))
        beta_hat[i,j,] = M %*% Y[i,j,]
    
    #Compute the residuals.
    R = array(0,c(length(x),length(y),n))
    for(i in 1:length(x))
      for(j in 1:length(y))
        R[i,j,] = Y[i,j,] - X %*% beta_hat[i,j,]
  }
  
  #Ignore values not on mask.
  if(!is.null(mask)) beta_hat[,,1][is.na(mask)] = level-10
  
  
  #Compute empirical variance and normalize.
  sigma_hat = apply(R,1:2,sd)
  R_tilde = R / rep(sigma_hat,n)
  
  #Compute a using Multiplier Bootstrap.
  a_MB = quantile(MBContour(x=x,y=y,R=R_tilde,N=N,f=beta_hat[,,1],level=level),
                  probs=1-alpha,type=8)
  #Compute a using Taylors method.
  Tay_fun = TaylorContour(x=x,y=y,f=beta_hat[,,1],level=level,R=R_tilde)
  a_Tay = 0
  while(2*Tay_fun(a_Tay)>alpha) a_Tay = a_Tay + 0.01  #The two is for the absolute value.
  
  #Compute threshold a with the true boundary if available.
  if(is.null(mu)) {a_MB_true = NA; a_Tay_true = NA} else{
    #Compute a using Multiplier Bootstrap.
    a_MB_true = quantile(MBContour(x=x,y=y,R=R_tilde,N=N,f=mu,level=level),
                         probs=1-alpha,type=8)
    #Compute a using Taylors method.
    Tay_fun = TaylorContour(x=x,y=y,f=mu,level=level,R=R_tilde) #The two is for the absolute value.
    a_Tay_true = 0
    while(2*Tay_fun(a_Tay_true)>alpha) a_Tay_true = a_Tay_true + 0.01
  }
  
  #Determine normalized estimators.
  norm_est =  (beta_hat[,,1] - level) / tau / sigma_hat
  
  #Determine whether inclusion holds.
  if(is.null(mu)){incl_MB = NA; incl_Tay = NA; incl_MB_true = NA; incl_Tay_true = NA} else{
    incl_MB = SubsetContour(x,y,norm_est,a_MB,mu,level) & SubsetContour(x,y,mu,level,norm_est,-a_MB)
    incl_Tay = SubsetContour(x,y,norm_est,a_Tay,mu,level) & SubsetContour(x,y,mu,level,norm_est,-a_Tay)
    incl_MB_true = SubsetContour(x,y,norm_est,a_MB_true,mu,level) & SubsetContour(x,y,mu,level,norm_est,-a_MB_true)
    incl_Tay_true = SubsetContour(x,y,norm_est,a_Tay_true,mu,level) & SubsetContour(x,y,mu,level,norm_est,-a_Tay_true)
  }
  
  result = structure(
    list(x=x,y=y,mu=mu,level=level,mu_hat=beta_hat[,,1],tau=tau,X=X,alpha=alpha,N=N,norm_est=norm_est, 
         a_MB=a_MB,a_Tay=a_Tay,a_MB_true=a_MB_true,a_Tay_true=a_Tay_true,
         incl_MB=incl_MB, incl_Tay=incl_Tay,incl_MB_true=incl_MB_true,incl_Tay_true=incl_Tay_true, mask=mask), class = "cope")
}

#' Plots CoPE sets. 
#'
#' @param cope An object of class cope to be plotted.
#' @param plot.taylor  Boolean indicating whether the CoPE sets with the threshold 
#'                      obtained by Taylor's method should be plotted. Default is
#'                      FALSE. 
#' @param use.true.function  Boolean indicating whether the threshold obtained 
#'                            from the true function should be used. Default is 
#'                            FALSE.  
#' @param map If TRUE plot the cope set on a map of the world. The coordinates 
#'            in this case are interpreted as longitude and latitude.  
#' @param ... Additional graphical parameters passed to fields::image.plot.
#' @export
PlotCope = function(cope,plot.taylor=FALSE, use.true.function = FALSE, map=FALSE, ...){
  
  if(map){
    plot.function = ImageMap
  } else{
    plot.function = fields::image.plot
  }
  
  if(use.true.function & is.null(cope$mu)){
    print("True function is not available. Estimated boundary will be used.")
  } 
  
  x = cope$x
  y = cope$y
  
  if(use.true.function & !is.null(cope$mu)){
    if(plot.taylor) a = cope$a_Tay_true else a = cope$a_MB_true
    plot.function(x,y,cope$mu_hat,horizontal=FALSE, mask=cope$mask, ...)
    DrawContour(x,y,cope$mu,level=cope$level,col="purple")
    DrawContour(x,y,cope$norm_est,level=a,col="darkred")
    DrawContour(x,y,cope$norm_est,level=-a,col="darkgreen")
  } else{
    if(plot.taylor) a = cope$a_Tay else a = cope$a_MB
    plot.function(x,y,cope$mu_hat,horizontal=FALSE, mask=cope$mask, ...)
    if(!is.null(cope$mu)){
      DrawContour(x,y,cope$mu,level=cope$level,col="purple")
    } else{
      DrawContour(x,y,cope$mu_hat,level=cope$level,col="purple")
    }
    DrawContour(x,y,cope$norm_est,level=a,col="darkred")
    DrawContour(x,y,cope$norm_est,level=-a,col="darkgreen")
  }
}