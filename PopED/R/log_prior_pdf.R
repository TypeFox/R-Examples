#' Compute the natural log of the PDF for the parameters in an E-family design
#' 
#' @inheritParams ed_laplace_ofv
#' @inheritParams evaluate.fim
#' @inheritParams create.poped.database
#' @inheritParams Doptim
#' @param alpha  A parameter vector.
#' @param return_hessian Should the hessian be returned?
#' 
#' @example tests/testthat/examples_fcn_doc/warfarin_optimize.R
#' @example tests/testthat/examples_fcn_doc/examples_log_prior_pdf.R
#' @export
#' @keywords internal
log_prior_pdf <- function(alpha, bpopdescr, ddescr,return_gradient=F,return_hessian=F){
  #returns the logarithm of the probability density for alpha,given the prior
  #and if required the gradient
  #priordescr=matrix(c(bpopdescr[bpopdescr[,1]!=0,], ddescr[ddescr[,1]!=0,]),nrow=1,byrow=T)
  priordescr <- rbind(bpopdescr,ddescr)
  priordescr <- priordescr[priordescr[,1]!=0,,drop=FALSE]
  if(any(priordescr[,1]>4 | priordescr[,1]==3)){
    stop(sprintf('Specified prior distribution not supported'))
  }
  mu=priordescr[,2]	#means
  sigma=priordescr[,3] #variances
  n=priordescr[,1]==1
  u=priordescr[,1]==2
  ln=priordescr[,1]==4
  
  #normal
  #temp=(-log(sqrt(2*pi*sigma))-(alpha-mu)^2/(2*sigma))*n
  temp_n <- log(dnorm(alpha[n],mean=mu[n],sd=sqrt(sigma[n]))) 
  #p=sum(temp[n])
  p=sum(temp_n)
  
  #uniform
  #temp=log((alpha>=mu-sigma/2 & alpha<=mu+sigma/2)*(1/sigma))
  #p=p+sum(temp[u])
  temp_u=log((alpha[u]>=mu[u]-sigma[u]/2 & alpha[u]<=mu[u]+sigma[u]/2)*(1/sigma[u]))
  p=p+sum(temp_u)
  
  #log normal
  #temp=log(1/alpha*own_normpdf(-log(mu)+log(alpha),0,sqrt(sigma)))
  #p=p+sum(temp[ln])
  #temp_ln=log(1/alpha[ln]*own_normpdf(-log(mu[ln])+log(alpha[ln]),0,sqrt(sigma[ln])))
  temp_ln <- log(dlnorm(alpha[ln],log(mu[ln]^2/sqrt(sigma[ln]+mu[ln]^2)),sqrt(log(sigma[ln]/mu[ln]^2+1))))
  p=p+sum(temp_ln)
  
  grad <- NULL
  if(return_gradient){
    grad=zeros(length(alpha),1)
    #normal
    grad[n]=-(alpha[n]-mu[n])/sigma[n]
    #uniform
    grad[u]=0
    #log normal
    temp=-(sigma+log(alpha)-log(mu))/(alpha*sigma)
    grad[ln]=temp[ln]
  }
 
  hess <- NULL
  if(return_hessian){
    diaghess=zeros(length(alpha), 1)
    diaghess[n]=-1/sigma[n]
    diaghess[u]=0
    temp=(-1+sigma+log(alpha)-log(mu))/(alpha^2*sigma)
    diaghess[ln]=temp[ln]
    hess=diag_matlab(diaghess)
  }
  
  #return(list( p= p,grad=grad,hess=hess)) 
  ret_args <- "p"
  if(!is.null(grad) || !is.null(hess)) ret_args <- "list(p=p"
  if(!is.null(grad)) ret_args <- paste(ret_args,",grad=grad",sep="")
  if(!is.null(hess)) ret_args <- paste(ret_args,",hess=hess",sep="")
  if(!is.null(grad) || !is.null(hess)) ret_args <- paste(ret_args,")",sep="")
        
  return(eval(parse(text=ret_args))) 
  
}

# own_normpdf <- function(x,mu,sigma){
#   y = exp(-0.5 * ((x - mu)/sigma)^2) / (sqrt(2*pi) * sigma)
#   return( y ) 
# }
