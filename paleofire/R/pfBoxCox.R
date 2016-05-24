#' Box-Cox transformation of Charcoal series
#' 
#' Box-Cox transformation of charcoal series, the maximum likelihood estimation
#' of lambda is derived from the boxcox.R function in the Venables and Ripley
#' MASS library included in R 2.6.1
#' 
#' 
#' @param serie A vector of charcoal values.
#' @param alpha Numeric, the "shift" parameter, default=0.01.
#' @param type Character, the Box-Cox transformation formulation, can be either
#' "BoxCox1964" (default) for the original Box & Cox (1964) formulation, or
#' "JohnDraper" for the John & Draper (1980) modulus transformation.
#' @return
#' 
#' \item{X}{Vector of transformed charcoal values}
#' @author P. Bartlein
#' @seealso \code{\link{pfTransform}}
#' @references Venables, W. N., Ripley, B. D., & Venables, W. N. (1994). Modern
#' applied statistics with S-PLUS (Vol. 250). New York: Springer-verlag. \cr
#' \cr Box, G.E.P. & Cox, D. R.(1964) An analysis of transformations, Journal
#' of the Royal Statistical Society, Series B, 26, 211-252. \cr \cr John, J. A.
#' & Draper N. R. (1980) Analternative family of transformations, Applied
#' Statistics, 29, 190-197.
#' @examples
#' 
#' # Select a site
#' ID=pfSiteSel(site_name=="Pas-de-Fond")
#' 
#' # Extract data
#' A=pfExtract(ID)
#' 
#' B=pfBoxCox(A[,4],0.1)
#' plot(B,type="l")
#' 
#' 
pfBoxCox=function(serie,alpha=0.01,type="BoxCox1964")
{
  types=c("BoxCox1964", "JohnDraper")
  warntype=type[(type %in% types)==FALSE]
  if(length(warntype)!=0){stop(paste(warntype, "is not a valid type for pfBoxCox", sep=" "))}
  
  # initial minimax rescaling of data, and addition of "shift" parameter (alpha)
  # value used in DP1
  if (alpha=="alternative"){
    alpha <- 0.5*min(serie[serie != 0])  # alternative alpha: 0.5 time the smallest nonzero value of quant
  }
  quant2=serie+alpha
  
  # maximum likelihood estimation of lambda 
  # derived from the boxcox.R function in the Venables and Ripley MASS library included in R 2.6.1
  
  npts <- 201 # number of estimates of lambda
  y <- quant2
  n <- length(y)
  logy <- log(y)
  ydot <- exp(mean(logy))
  lasave <- matrix(1:npts)
  liksave <- matrix(1:npts)
  for (i in 1:npts) {
    la <- -2.0+(i-1)*(4/(npts-1))
    if (la != 0.0) yt <- (y^la-1)/la else yt <- logy*(1+(la*logy)/2*(1+(la*logy)/3*(1+(la*logy)/4)))
    zt <- yt/ydot^(la-1)
    loglik <- -n/2*log(sum((zt - mean(zt))^2 ))
    lasave[i] <- la
    liksave[i] <- loglik
  }
  #print(cbind(lasave,liksave))
  cbind(1,liksave[which.max(liksave)],lasave[which.max(liksave)])
  laopt=as.character(lasave[which.max(liksave)])
  lafit <- lasave[which.max(liksave)]
  
  # Box-Cox transformation of data
  
  # Box-Cox original
  if (type=="BoxCox1964"){
    if (lafit == 0.0) tquant = log(quant2) else tquant = (quant2^lafit - 1)/lafit
  }
  # Modulus transformation
  if (type=="JohnDraper"){
    s=sign(quant2)
    s[s==0]=1
    
    if (lafit == 0.0)  tquant = s*(log(abs(quant2)+1))
    else tquant =s*( ((abs(quant2)+1)^lafit-1) /lafit )
  }
  
  return(tquant)
}
