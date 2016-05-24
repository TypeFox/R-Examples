#'@title Asymmetry curve based on depths
#'
#'@description Produces an asymmetry curve estimated from given data.
#'
#'  @param x  The data as a matrix or data frame. If it is a matrix or data frame, then each row is viewed as one multivariate observation. 
#'  @param y  Additional matrix of multivariate data.
#'  @param alpha  An ordered vector containing indices of central regins used for asymmetry curve calculation.
#'  @param method Character string which determines the depth function used. The method can be "Projection" (the default), "Mahalanobis", "Euclidean", "Tukey" or 'LP'.  For details see \code{\link{depth}.}
#' @param movingmedian  Logical. For default FALSE only one depth median is used to compute asymmetry norm. If TRUE - for every central area, a new depth median will be used - this approach needs much more time.
#' @param name Name of set X - used in plot legend
#' @param name_y Name of set Y - used in plot legend
#' @param ... Any additional parameters for function depth
#'
#' @details 
#'  
#'  For sample depth function  \eqn{ D({x},{{{Z}}^{n}}) } ,  \eqn{ {x}\in {{{R}}^{d}} } ,  \eqn{ d\ge 2 } ,  \eqn{ {Z}^{n}=\{{{{z}}_{1}},...,{{{z}}_{n}}\}\subset {{{R}}^{d}} } ,   \eqn{ {{D}_{\alpha }}({{{Z}}^{n}}) }  denoting  \eqn{ \alpha- }  central region, we can define {the asymmetry curve}  \eqn{ AC(\alpha )=\left( \alpha ,\left\| {{c}^{-1}}(\{{\bar{z}}-med|{{D}_{\alpha }}({{{Z}}^{n}})\}) \right\| \right)\subset {{{R}}^{2}} } , for   \eqn{ \alpha \in [0,1]  }  being nonparametric scale and asymmetry functional correspondingly, where  \eqn{ c- } denotes constant,  \eqn{ {\bar{z}}- }  denotes mean vector, denotes multivariate median induced by depth function and  \eqn{ vol- }  denotes a volume.
#'  
#'  Asymmetrycurve takes uses function convhulln from package geometry for computing a volume of convex hull containing central region.
#'  
#'  
#'  @author Daniel Kosiorowski, Mateusz Bocian, Anna Wegrzynkiewicz and Zygmunt Zawadzki from Cracow University of Economics.
#'  
#'  
#'  @references 
#'  
#' Serfling  R. J.  Multivariate Symmetry and Asymmetry, \emph{Encyclopedia of Statistical Science}, S Kotz, C.B. Read, N. Balakrishnan, B. Vidakovic (eds), 2nd, ed., John Wiley.
#'
#'Liu, R.Y., Parelius, J.M. and Singh, K. (1999), Multivariate analysis by data depth: Descriptive statistics, graphics and inference (with discussion), \emph{Ann. Statist.}, \bold{27}, 783--858.
#'
#'Chaudhuri, P. (1996), On a Geometric Notion of Quantiles for Multivariate Data, \emph{Journal of the American Statistical Association}, 862--872.
#'
#'Dyckerhoff, R. (2004), Data Depths Satisfying the Projection Property, \emph{Allgemeines Statistisches Archiv.},  \bold{88}, 163--190.
#'
#'  
#'  @seealso \code{\link{scaleCurve}}, \code{\link{depth}}
#'  
#'  @export
#'  
#'  @examples
#'  
#' #EXAMPLE 1
#' require(sn)
#' xi = c(0,0)
#' alpha <- c(2,-5)
#' Omega <- diag(2)*5
#' 
#' n = 500
#' X = mvrnorm(n, xi, Omega)  # normal distribution
#' Y = rmst(n, xi, Omega, alpha, nu=1)
#' asymmetryCurve(X,Y,name = "NORM",name_y = "S_T(2,-5,10)")
#'  
#'  
#' #EXAMPLE 2
#' data(under5.mort)
#' data(inf.mort)
#' data(maesles.imm)
#' data1990=cbind(under5.mort[,1],inf.mort[,1],maesles.imm[,1])
#' data2011=cbind(under5.mort[,22],inf.mort[,22],maesles.imm[,22])
#' as1990=asymmetryCurve(data1990,name='scale curve 1990')
#' as2011=asymmetryCurve(data2011,name='scale curve 2011')
#' figure=getPlot(as1990 %+% as2011)+ggtitle('Scale curves')
#' figure
#'  
#'  @keywords
#'  multivariate
#'  nonparametric
#'  robust
#'  depth function
#'  asymmetry
#'
asymmetryCurve<-function(x, y = NULL, alpha = seq(0,1,0.01), method = "Projection",
	movingmedian = FALSE, name = "X", name_y = "Y",...)
{
  x = na.omit(x)
  
  if(nrow(x)< NCOL(x)*10) stop("Too small sample!")
	if(!is.matrix(x)) stop("X must be a matrix!")
	if(!is.null(y)) if(!is.matrix(y)) stop("Y must be a matrix!")
  x = na.omit(x)
	
	depth_est <- depth(x,x,method=method, name=name) 
  
  if(!movingmedian) x_est = .asCurve(x, depth_est, alpha, method, ...) #function in as_curve.R
	else x_est = .asCurveMM(x,alpha,method,...)  #function in as_curve_mm.R
	
	asc = new("AsymmetryCurve", x_est[,2], depth = depth_est, alpha = x_est[,1])
  
  if(!is.null(y)){
    asc = asc %+% asymmetryCurve(y, y =NULL, alpha, method,
                                         movingmedian, name = name_y, name_y = "Y",...)
  }
	return(asc)
}

.asCurve<-function(X, depth_est = NULL,alpha = NULL,method = "Projection",...)
{
  dim_X <- dim(X)[2] 
  
  if(is.null(depth_est)) depth_est = depth(X,X,method=method,...) 
  
  
  median =  depthMedian(X,method=method,...)
  
  if((length(median)/dim(X)[2])!=1)
  {
    median = colMeans(median)
  }
  
  k <-length(alpha)
  vol = 1:k 
  alpha_est = 1:k 
  means = matrix(nrow = k, ncol = dim_X) 
  
  
  for(i in 1:k)
  {
    tmp_X <- X[depth_est >= alpha[i],]
    np <- dim(as.matrix(tmp_X))[1] 
    
    if ((np > ((2*(dim_X+2))^2))&(np > 100))
    {
      vol[i] <- convhulln(tmp_X,options = "FA")$vol
      alpha_est[i] <- alpha[i]
      means[i,] <- colMeans(tmp_X)
      
    }
    else
    {
      vol[i]=NA 
      alpha_est[i] <- NA  
    }
  }
  
  
  #vol <- vol[!is.na(vol)]
  #alpha_est <- alpha_est[!is.na(alpha_est)]
  #means <- matrix(means[!is.na(means)],ncol=dim_X)
  
  k = length(vol)
  kmedian = matrix(rep(median,k),byrow=TRUE,ncol=dim_X) # taka ma?a optymalizacja
  
  n=(means - kmedian)
  nn = 2*sqrt(rowSums(n^2))/(vol^(1/dim_X)) 
  
  matrix(c(rev(1-alpha_est),rev(nn)),ncol=2)
  
}


####################################
####################################
.asCurveMM<-function(X, depth_est = NULL,alpha = NULL,method,...)
{
  
  dim_X <- dim(X)[2]
  
  if(is.null(depth_est)) depth_est = depth(X,X,method=method,...) 
  if(is.null(alpha)) alpha = seq(0,1,0.01)
  
  
  k <-length(alpha)
  vol = 1:k 
  alpha_est = 1:k
  means = matrix(nrow = k, ncol = dim_X) 
  medians = matrix(nrow = k, ncol = dim_X)
  
  for(i in 1:k)
  {
    tmp_X <- X[depth_est >= alpha[i],]
    np <- dim(as.matrix(tmp_X))[1] 
    
    if ((np > ((2*(dim_X+2))^2))&(np > 100))
    { 
      vol[i] <- convhulln(tmp_X,options = "FA")$vol
      alpha_est[i] <- alpha[i]
      means[i,] <- colMeans(tmp_X)
      
      tmp_depth_est = depth(tmp_X,tmp_X,method,...)
      
      med = tmp_X[tmp_depth_est==max(tmp_depth_est),]
      
      if((length(med)/dim(X)[2])!=1)
      {
        
        
        
        med = colMeans(med)
        
      }
      medians[i,] = med
      
      
    }
    else
    {
      vol[i]=NA 
      alpha_est[i] <- NA  
    }
  }
  
  
  vol <- vol[!is.na(vol)]
  alpha_est <- alpha_est[!is.na(alpha_est)]
  means <- matrix(means[!is.na(means)],ncol=dim_X)
  medians <- matrix(medians[!is.na(medians)],ncol=dim_X)
  #koniec usuwani NA
  
  k = length(vol)
  kvol = matrix(rep(vol,dim_X),ncol=dim_X) 
  n=(means - medians)#/(kvol^(1/dim_X))
  nn = 2*sqrt(rowSums(n^2))/(vol^(1/dim_X)) 
  matrix(c(rev(1-alpha_est),rev(nn)),ncol=2)
  
}



