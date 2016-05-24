#' Estimation of dichotomous logistic Rasch model (Rasch, 1960)
#'
#' This function estimates the dichotomous Rasch model by Rasch
#' (1960).
#'
#' Parameters are estimated by CML.
#'
#' @param data Data matrix or data frame; rows represent observations
#' (persons), columns represent the items.
#' @param desmat Design matrix; if missing, the design matrix for a dichotomous Rasch model will be created automatically.
#' @param start starting values for parameter estimation. If missing, a vector of 0 is used as starting values.
#' @param control list with control parameters for the estimation process e.g. the convergence criterion. For details please see the help pages to the R built-in function \code{optim}
#'
#' @return \item{data}{data matrix according to the input} \item{design}{design
#' matrix either according to the input or according to the automatically
#' generated matrix} \item{logLikelihood}{conditional log-likelihood}
#' \item{estpar}{estimated basic item parameters}
#' \item{estpar_se}{estimated standard errors for basic item
#' parameters} \item{itempar}{estimated item parameters}
#' \item{itempar_se}{estimated standard errors for item parameters}
#' \item{hessian}{Hessian matrix} \item{convergence}{convergence of solution
#' (see help files in \code{\link{optim}})} \item{fun_calls}{number of function
#' calls (see help files in \code{\link{optim}})}
#' @author Christine Hohensinn
#' @references Fischer, G. H. (1974). Einfuehrung in die Theorie psychologischer Tests
#' [Introduction to test theory]. Bern: Huber.
#'
#' Rasch, G. (1960). Probabalistic models for some intelligence and attainment tests. Danmarks paedagogiske institut.
#' @keywords continuous rating scale model
#'
#' @rdname drm
#'
#' @export
#'
#' @examples
#'
#' #estimate Rasch model parameters
#' data(reason)
#' res_drm <- DRM(reason.test[,1:11])
#'
#' summary(res_drm)
#'
#'
DRM <-
  function(data, desmat, start, control){

call <- match.call()

data <- as.matrix(data)
rv <- rowSums(data)
si <- colSums(data)
dat.r <- data[!(rv %in% c(0,ncol(data))),]

rv.r <- rowSums(dat.r)
si.r <- colSums(dat.r)
fr.r <- tabulate(rv.r, nbins=(ncol(data)-1))

if(missing(start)){
  para <- rep(0,(ncol(data)-1))
} else if(start == "improved"){
  para <-(nrow(dat.r)-si.r[1:(ncol(dat.r)-1)])/si.r[1:(ncol(dat.r)-1)]-(sum(log((nrow(dat.r)-si.r)/si.r))/length(si.r))
} else if(start != "improved" & is.numeric(start)==FALSE){
  stop("start values have to be numeric!")
} else {
  para <- start
}

if(missing(desmat)){
  desmat <- rbind(diag(1,nrow=length(para)),-1)
}

cmlLoglik <- function(sigs, si, roh, desmat){
  sigs2 <- desmat %*%sigs
  term1 <- as.vector(si%*%sigs2)
  epsi  <- exp(sigs2)
  fgmat <- matrix(0,ncol=length(epsi), nrow=length(epsi))
  fgmat[1,]<- cumsum(epsi)
  gam <- gamfunk(epsi, fgmat)[1:(ncol(data)-1)]
  term2 <- as.vector(roh%*%log(gam))
  term1-term2
}

if(missing(control)){
  control <- list(fnscale=-1, maxit=1000)}
res <- optim(para, fn=cmlLoglik, si=si.r, roh=fr.r,desmat=desmat, control=control, hessian=TRUE, method="L-BFGS-B")

# normalize item parameters
itempar <- as.vector(desmat %*% res$par)

#SE
estpar_se <- sqrt(diag(solve(res$hessian*(-1))))

#itempar_se <- sqrt(diag(desmat %*% solve(res$hessian*(-1)) %*% t(desmat)))
itempar_se <- c(estpar_se,NA)

if(!is.null(colnames(data))){
  names(itempar) <- paste("beta", colnames(data))
  names(itempar_se) <- paste("SE", colnames(data))
} else {
  names(itempar) <- paste("beta item", 1:length(itempar))
  names(itempar_se) <- paste("SE item", 1:length(itempar))
}


res_all <- list(data=data, design=desmat,logLikelihood=res$value, estpar=res$par, estpar_se=estpar_se, itempar=itempar*(-1), itempar_se=itempar_se, hessian=res$hessian*(-1), convergence=res$convergence, fun_calls=res$counts, call=call)
class(res_all) <- "DRM"
res_all
}


