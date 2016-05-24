#' Mark-Recapture Distance Sampling (MRDS) Removal - FI
#'
#' Mark-Recapture Distance Sampling (MRDS) Analysis of Removal Observer
#' Configuration with Full Independence
#'
#' The mark-recapture data derived from an removal observer distance sampling
#' survey can only derive conditional detection functions (p_j(y)) for both
#' observers (j=1) because technically it assumes that detection probability
#' does not vary by occasion (observer in this case).  It is a conditional
#' detection function because detection probability for observer 1 is
#' conditional on the observations seen by either of the observers. Thus,
#' p_1(y) is estimated by p_1|2(y).
#'
#' If detections by the observers are
#' independent (full independence) then p_1(y)=p_1|2(y) and for the union, full
#' independence means that p(y)=p_1(y) + p_2(y) - p_1(y)*p_2(y) for each
#' distance y.  In fitting the detection functions the likelihood from Laake
#' and Borchers (2004) are used. That analysis does not require the usual
#' distance sampling assumption that perpendicular distances are uniformly
#' distributed based on line placement that is random relative to animal
#' distribution.  However, that assumption is used in computing predicted
#' detection probability which is averaged based on a uniform distribution (see
#' eq 6.11 of Laake and Borchers 2004).
#'
#' For a complete description of each of the calling arguments, see
#' \code{\link{ddf}}.  The argument \code{model} in this function is the same
#' as \code{mrmodel} in \code{ddf}.  The argument \code{dataname} is the name
#' of the dataframe specified by the argument \code{data} in \code{ddf}. The
#' arguments \code{control},\code{meta.data},and \code{method} are defined the
#' same as in \code{ddf}.
#'
#' @export
#' @method ddf rem.fi
#' @param model mark-recapture model specification
#' @param data analysis dataframe
#' @param meta.data list containing settings controlling data structure
#' @param control list containing settings controlling model fitting
#' @param call original function call used to call \code{ddf}
#' @param method analysis method; only needed if this function called from
#'   \code{ddf.io}
#' @return result: an rem.fi model object
#' @author Jeff Laake
#' @seealso \code{\link{ddf.io}},\code{\link{rem.glm}}
#' @references Laake, J.L. and D.L. Borchers. 2004. Methods for incomplete
#'   detection at distance zero. In: Advanced Distance Sampling, eds. S.T.
#'   Buckland, D.R.Anderson, K.P. Burnham, J.L. Laake, D.L. Borchers, and L.
#'   Thomas. Oxford University Press.
#' @keywords Statistical Models
ddf.rem.fi<-function(model,data,meta.data=list(),control=list(),call="",method){
  #  NOTE: gams are only partially implemented

  # The following are dummy glm and gam functions that are defined here to
  # provide the list of arguments for use in the real glm/gam functions. These
  # dummy functions are removed after they are used so the real ones can be
  # used in the model fitting.

  glm <- function(formula,link="logit"){
    if(class(formula)!="formula"){
      if(class(try(as.formula(formula)))!="formula"){
        stop("Invalid formula")
      }
    }else{
      formula <- paste(as.character(formula),collapse="") 
    }
    if(class(link)=="function"){
      link <- substitute(link)
    }
    link <- match.arg(link,c("logit"))
    return(list(fct="glm", formula=formula, link=substitute(link)))
  }


  gam <- function(formula,link="logit"){
    if(class(formula)!="formula"){
      if(class(try(as.formula(formula)))!="formula"){
        stop("Invalid formula")
      }
    }else{
      formula <- paste(as.character(formula),collapse="")
    }
    if(class(link)=="function"){
      link <- substitute(link)
    }
    link <- match.arg(link,c("logit"))
    return(list(fct="gam", formula=formula, link=substitute(link)))
  }

  # Test to make sure that observer not used in mrmodel
  if(length(grep("observer",model))!=0){
    stop("observer cannot be included in models for removal configuration\n")
  }

  # Save current user options and then set design contrasts to treatment style
  save.options <- options()
  options(contrasts=c("contr.treatment","contr.poly"))

  # Set up meta data values
  meta.data <- assign.default.values(meta.data, left=0, width=NA, binned=FALSE,
                                     int.range=NA, point=FALSE)

  # Set up control values
  control <- assign.default.values(control, showit=0,
                                   estimate=TRUE, refit=TRUE, nrefits=25,
                                   initial=NA, lowerbounds=NA, upperbounds=NA,
                                   mono.points=20)

  # Assign model values; this uses temporarily defined functions glm and gam
  modpaste <- paste(model)
  modelvalues <- try(eval(parse(text=modpaste[2:length(modpaste)])))
  if(class(modelvalues)=="try-error"){
    stop("Invalid model specification: ",model)
  }
  rm(glm,gam)

  # Process data if needed
   if(is.data.frame(data)){
     data.list <- process.data(data,meta.data)
     meta.data <- data.list$meta.data
     xmat <- data.list$xmat
   }else{
     xmat <- data
   }

  # Setup default breaks
  if(meta.data$binned){
    meta.data$breaks <- c(max(0,
                           min(as.numeric(levels(as.factor(xmat$distbegin))))),
                          as.numeric(levels(as.factor(xmat$distend))))
  }

  # Create result list with some arguments
  result <- list(call=call, data=data, model=model, meta.data=meta.data,
                 control=control, method="rem.fi")
  class(result) <- c("rem.fi","ddf")

  # Create formula and model frame (if not GAM);
  # xmat2 contains data for observer 2 and xmat for observer 1.
  # This is only done to allow a field for personnel which are rotated in
  # the observer roles of 1/2.
  xmat$offsetvalue <- rep(0,dim(xmat)[1])
  model.formula <- paste("detected",modelvalues$formula)
  xmat2 <- xmat[xmat$observer==2,]
  xmat1 <- xmat[xmat$observer==1,]
  p.formula <- as.formula(model.formula)
  npar <- ncol(model.matrix(p.formula,xmat1))

  #  fit=optim(par=rep(0,npar),lnl.removal,x1=xmat1,x2=xmat2,
  #            models=list(p.formula=p.formula),hessian=TRUE,
  #            control=list(maxit=5000))
  #  fit$hessian=hessian(lnl.removal,x=fit$par,method="Richardson",x1=xmat1,
  #                      x2=xmat2,models=list(p.formula=p.formula))
  GAM <- FALSE
  #if(modelvalues$fct=="gam"){
  #  GAM=TRUE
  #}else
  xmat1 <- create.model.frame(xmat1,as.formula(model.formula),meta.data)
  xmat2 <- create.model.frame(xmat2,as.formula(model.formula),meta.data)
  model.formula <- as.formula(paste(model.formula,"+offset(offsetvalue)"))

  # Fit the conditional detection functions using io.glm
  suppressWarnings(result$mr <- rem.glm(xmat1,model.formula,GAM,
                                        datavec2=xmat2,iterlimit=1))
  # if(GAM)result$mr$data=xmat1

  #  Now use optimx with starting values perturbed by 5%
  fit <- optimx(1.05*result$mr$coefficients, lnl.removal, method="nlminb",
                hessian=TRUE, x1=xmat1, x2=xmat2,
                models=list(p.formula=p.formula))
  topfit.par <- coef(fit, order="value")[1, ]
  details <- attr(fit,"details")[1,]
  fit <- as.list(summary(fit, order="value")[1, ])
  fit$par <- topfit.par
  fit$message <- ""
  names(fit)[names(fit)=="convcode"] <- "conv"
  fit$hessian <- details$nhatend
  result$mr$mr$coefficients <- fit$par
  result$hessian <- fit$hessian

  # Compute the L_omega portion of the likelihood value, AIC and hessian
  cond.det <- predict(result)
  p1 <- cond.det$p1
  p2 <- cond.det$p2
  result$par <- coef(result$mr)
  npar <- length(result$par)
  p.c.omega <- p1^result$mr$data$detected[result$mr$data$observer==1]*
            ((1-p1)*p2)^(1-result$mr$data$detected[result$mr$data$observer==1])
  result$lnl <- sum(log(p.c.omega)) - sum(log(cond.det$fitted))

  #if(GAM){
  #  result$hessian <- result$mr$Vp
  #}else{
  #  result$hessian <- solve(summary(result$mr)$cov.unscaled) 
  #}

  # Compute fitted values
  result$fitted <- predict(result,newdata=xmat,integrate=TRUE)$fitted

  # If this is not TI mode, compute the lnlU2 portion of the likelihood
  # using the conditional detection functions.
  if(method=="rem.fi"){
    distances <- xmat$distance[xmat$observer==1]
    if(!meta.data$binned){
      if(meta.data$point){
        result$lnl<- result$lnl +
                     sum(log(predict(result,newdat=xmat,integrate=FALSE)$fitted*
                     2*distances/meta.data$width^2)) - sum(log(result$fitted))
      }else{
        result$lnl<- result$lnl +
                     sum(log(predict(result,newdat=xmat,
                                     integrate=FALSE)$fitted/meta.data$width))-
                      sum(log(result$fitted))
      }
    }else{
      for(i in 1:(nrow(xmat)/2)){
        int.val <- predict(result,newdata=xmat[(2*(i-1)+1):(2*i),],
                           int.range=as.vector(as.matrix(
                              xmat[(2*(i-1)+1),c("distbegin","distend")])),
                           integrate=TRUE)$fitted
        result$lnl <- result$lnl + log(int.val)
      }
      result$lnl <- result$lnl- sum(log(result$fitted))
    }
  }

  # Compute AIC and Nhat in covered region
  result$Nhat <- NCovered(result)
  result$criterion <- -2*result$lnl + 2*npar
  result$par <- coef(result$mr)

  # Restore user options
  options(save.options)

  # Return result
   return(result)
}

lnl.removal <- function(par,x1,x2,models){
  # compute probabilities
  p.list <- p.removal.mr(par,x1,x2,models)
  # if any are 0 set to a small value
  p11 <- p.list$p11
  p01 <- p.list$p01
  p01[p01==0] <- 1e-6
  # compute negative log-likelihood value and return it
  lnl <- sum((1-x1$detected)*log(p01))+ sum(x1$detected*log(p11))-
         sum(log(p.list$pdot))
  return(-lnl)
}

p.removal.mr <- function(par,x1,x2,models){
  x <- rbind(x1,x2)
  dmrows <- nrow(x1)
  xmat <- model.matrix(models$p.formula,x)
  xmat1 <- xmat[1:dmrows,,drop=FALSE]
  xmat2 <- xmat[(dmrows+1):(2*dmrows),,drop=FALSE]
  p01 <- rem.p01(xmat1,xmat2,beta=par)
  p11 <- rem.p11(xmat1,xmat2,beta=par)
  pdot <- rem.pdot(xmat1,xmat2,beta=par)
  return(list(p11=as.vector(p11),p01=as.vector(p01),pdot=as.vector(pdot)))
}

rem.p01 <- function(xmat1,xmat2,beta){
  p1 <- plogis(xmat1%*%beta)
  p2 <- plogis(xmat2%*%beta)
  return(p2*(1-p1))
}

rem.p11 <- function(xmat1,xmat2=NULL,beta){
  return(plogis(xmat1%*%beta))
}

rem.pdot <- function(xmat1,xmat2,beta){
  p1 <- plogis(xmat1%*%beta)
  p2 <- plogis(xmat2%*%beta)
  return(p1+p2-p1*p2)
}
