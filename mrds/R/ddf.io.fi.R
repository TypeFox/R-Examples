#' Mark-Recapture Distance Sampling (MRDS) IO - FI
#'
#' Mark-Recapture Analysis of Independent Observer Configuration with Full
#' Independence
#'
#' The mark-recapture data derived from an independent observer distance
#' sampling survey can be used to derive conditional detection functions
#' (p_j(y)) for both observers (j=1,2).  They are conditional detection
#' functions because detection probability for observer j is based on seeing or
#' not seeing observations made by observer 3-j. Thus, p_1(y) is estimated by
#' p_1|2(y).
#'
#' If detections by the observers are independent (full
#' independence) then p_1(y)=p_1|2(y),p_2(y)=p_2|1(y) and for the union, full
#' independence means that p(y)=p_1(y) + p_2(y) - p_1(y)*p_2(y) for each
#' distance y.  In fitting the detection functions the likelihood given by eq
#' 6.8 and 6.16 in Laake and Borchers (2004) is used. That analysis does not
#' require the usual distance sampling assumption that perpendicular distances
#' are uniformly distributed based on line placement that is random relative to
#' animal distribution.  However, that assumption is used in computing
#' predicted detection probability which is averaged based on a uniform
#' distribution (see eq 6.11 of Laake and Borchers 2004).
#'
#' For a complete description of each of the calling arguments, see
#' \code{\link{ddf}}.  The argument \code{model} in this function is the same
#' as \code{mrmodel} in \code{ddf}.  The argument \code{dataname} is the name
#' of the dataframe specified by the argument \code{data} in \code{ddf}. The
#' arguments \code{control},\code{meta.data},and \code{method} are defined the
#' same as in \code{ddf}.
#'
#' @export
#' @method ddf io.fi
#' @param model mark-recapture model specification
#' @param data analysis dataframe
#' @param meta.data list containing settings controlling data structure
#' @param control list containing settings controlling model fitting
#' @param call original function call used to call \code{ddf}
#' @param method analysis method; only needed if this function called from
#'   \code{ddf.io}
#' @return result: an io.fi model object
#' @author Jeff Laake
#' @seealso
#'   \code{\link{ddf.io}},\code{\link{summary.io.fi}},\code{\link{coef.io.fi}},
#'   \code{\link{plot.io.fi}},\code{\link{gof.io.fi}},\code{\link{io.glm}}
#' @references Laake, J.L. and D.L. Borchers. 2004. Methods for incomplete
#'   detection at distance zero. In: Advanced Distance Sampling, eds. S.T.
#'   Buckland, D.R.Anderson, K.P. Burnham, J.L. Laake, D.L. Borchers, and L.
#'   Thomas. Oxford University Press.
#' @keywords Statistical Models
#' @importFrom stats optimHess
ddf.io.fi <- function(model,data,meta.data=list(),control=list(),
                      call="",method){
  # Functions used: assign.default.values, process.data, create.model.frame
  #                 ioglm, predict(predict.io.fi), NCovered (NCovered.io.fi)

  # NOTE: gams are only partially implemented

  # The following are dummy glm and gam functions that are defined here to 
  # provide the list of arguments for use in the real glm/gam functions. 
  # These dummy functions are removed after they are used so the real ones can 
  # be used in the model fitting.
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
    return(list(fct="glm",formula=formula,link=substitute(link)))
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
    return(list(fct="gam",formula=formula,link=substitute(link)))
  }

  # Save current user options and then set design contrasts to treatment style
  save.options<-options()
  options(contrasts=c("contr.treatment","contr.poly"))

  # Set up meta data values
  meta.data=assign.default.values(meta.data, left=0, width=NA, binned=FALSE,
                                  int.range=NA, point=FALSE)

  # Set up control values
  control <- assign.default.values(control,showit = 0,
                                   estimate=TRUE,refit=TRUE,nrefits=25,
                                   initial = NA, lowerbounds = NA,
                                   upperbounds = NA, mono.points=20)

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
                              min(as.numeric(levels(as.factor(xmat$distbegin))))
                             ),
                       as.numeric(levels(as.factor(xmat$distend))))
  }

  # Create result list with some arguments
  result <- list(call=call,data=data,model=model,meta.data=meta.data,
                 control=control,method="io.fi")
  class(result) <- c("io.fi","ddf")

  # Create formula and model frame (if not GAM)
  xmat$offsetvalue <- rep(0,dim(xmat)[1])
  model.formula <- paste("detected",modelvalues$formula)
  p.formula <- as.formula(model.formula)
  xmat2 <- xmat[xmat$observer==2,]
  xmat1 <- xmat[xmat$observer==1,]
  npar <- ncol(model.matrix(p.formula,xmat))

#   fit=optim(par=rep(0,npar),lnl.io,x1=xmat1,x2=xmat2,models=list(p.formula=p.formula),
#     hessian=TRUE,control=list(maxit=5000))
#   fit$hessian=hessian(lnl.removal,x=fit$par,method="Richardson",x1=xmat1,x2=xmat2,models=list(p.formula=p.formula))
  GAM <- FALSE
  #if(modelvalues$fct=="gam") 
  #{
  #   GAM <- TRUE
  #} else
  xmat <- create.model.frame(xmat,as.formula(model.formula),meta.data)
  model.formula <- as.formula(paste(model.formula,"+offset(offsetvalue)"))

  # Fit the conditional detection functions using io.glm 
  suppressWarnings(result$mr <- io.glm (xmat,model.formula,GAM=GAM))

  # if(GAM)result$mr$data=xmat

  if(result$mr$converged){
    # if the glm did converge, then do one quick round of BFGS to 
    #  compute the hessian
    result$hessian<- optimHess(result$mr$coefficients, lnl.io, 
                               x1=xmat1, x2=xmat2, 
                               models=list(p.formula=p.formula))
  }else{
    # only resort to optim if there was not convergence in the glm!

    # Now use optimx with starting values perturbed by 5%
    fit <-try(optimx(result$mr$coefficients*1.05,
                      lnl.io, method="L-BFGS-B", hessian=TRUE,x1=xmat1,x2=xmat2,
                      models=list(p.formula=p.formula)))
    # did this model converge?
    topfit.par <- coef(fit, order="value")[1, ]
    details <- attr(fit,"details")[1,]
    fit <- as.list(summary(fit, order="value")[1, ])
    fit$par <- topfit.par
    fit$message <- ""
    names(fit)[names(fit)=="convcode"] <- "conv"
    fit$hessian<-details$nhatend

    if(fit$conv!=0 | class(fit)=="try-error"){
      # first try the old way of just setting the starting values to zero,
      # this seems (with the crabbie data) to converge back to the values in
      # result$mr$coefficients but without having the convergence issues
      # NEED TO THINK ABOUT THIS MORE...
      fit <- try(optimx(result$mr$coefficients*0,lnl.io, method="L-BFGS-B", 
              hessian=TRUE,x1=xmat1,x2=xmat2,models=list(p.formula=p.formula)))
      topfit.par <- coef(fit, order="value")[1, ]
      details <- attr(fit,"details")[1,]
      fit <- as.list(summary(fit, order="value")[1, ])
      fit$par <- topfit.par
      fit$message <- ""
      names(fit)[names(fit)=="convcode"] <- "conv"
      fit$hessian<-details$nhatend
    }
    # if nothing worked...
    if(fit$conv!=0 | class(fit)=="try-error"){
      stop("No convergence in ddf.io.fi()")
    }else{
      result$mr$mr$coefficients <- fit$par
      result$hessian <- fit$hessian
    }
  }
  # Compute the L_omega portion of the likelihood value, AIC and hessian
  cond.det <- predict(result)
  result$par <- coef(result$mr)  
  npar <- length(result$par)
  p1 <- cond.det$p1
  p2 <- cond.det$p2
  p.c.omega <- p1^result$mr$data$detected[result$mr$data$observer==1]*
       (1-p1)^(1-result$mr$data$detected[result$mr$data$observer==1])*
       p2^result$mr$data$detected[result$mr$data$observer==2]*
       (1-p2)^(1-result$mr$data$detected[result$mr$data$observer==2])

  result$lnl <- sum(log(p.c.omega)) - sum(log(cond.det$fitted))

  #if(GAM)
  #   result$hessian <- result$mr$Vp
  #else
  #   result$hessian <- solve(summary(result$mr)$cov.unscaled)

  # If this is method=io.fi then compute lnl for L_y and add to L_omega before
  # computing AIC. Note the code in predict is currently specific to logit link
  if(method=="io.fi"){
    # Compute fitted values - integral of p(x)/width
    result$fitted <- predict(result,integrate=TRUE)$fitted
    result$Nhat <- NCovered(result$par,result)
    distances <- result$data$distance[result$data$observer==1]
    if(!meta.data$binned){
      if(meta.data$point){
        result$lnl <- result$lnl +
            sum(log(cond.det$fitted*2*distances/meta.data$width^2)) -
            sum(log(result$fitted))
      }else{
        result$lnl<- result$lnl +
            sum(log(cond.det$fitted/meta.data$width)) -
            sum(log(result$fitted))
      }
    }else{
      for(i in 1:(nrow(xmat)/2)){
        int.val <- predict(result,newdata=xmat[(2*(i-1)+1):(2*i),],
                           int.range=as.vector(as.matrix(xmat[(2*(i-1)+1),
                                               c("distbegin","distend")])),
                           integrate=TRUE)$fitted
        result$lnl <- result$lnl + log(int.val)
      }
      result$lnl<- result$lnl- sum(log(result$fitted))
    }
  }

  result$criterion <- -2*result$lnl + 2*npar

  # Restore user options
  options(save.options)

  # Return result
  return(result)
}

io.p01 <- function(xmat1,xmat2,beta){
  p1 <- plogis(xmat1%*%beta)
  p2 <- plogis(xmat2%*%beta)
  return(p2*(1-p1))
}
io.p10 <- function(xmat1,xmat2,beta){
  p1 <- plogis(xmat1%*%beta)
  p2 <- plogis(xmat2%*%beta)
  return(p1*(1-p2))
}

io.p11 <- function(xmat1,xmat2,beta){
  p1 <- plogis(xmat1%*%beta)
  p2 <- plogis(xmat2%*%beta)
  return(p1*p2)
}

io.pdot<-function(xmat1,xmat2,beta){
  p1 <- plogis(xmat1%*%beta)
  p2 <- plogis(xmat2%*%beta)
  return(p1+p2-p1*p2)
}

p.io <- function(par,x1,x2,models){
  # create design matrix for p1,p2 and delta
  x <- rbind(x1,x2)
  dmrows <- nrow(x1)
  xmat <- model.matrix(models$p.formula,x)
  xmat1 <- xmat[1:dmrows,,drop=FALSE]
  xmat2 <- xmat[(dmrows+1):(2*dmrows),,drop=FALSE]
  p01 <- io.p01(xmat1,xmat2,beta=par)
  p10 <- io.p10(xmat1,xmat2,beta=par)
  p11 <- io.p11(xmat1,xmat2,beta=par)
  pdot <- io.pdot(xmat1,xmat2,beta=par)
  return(list(p11=as.vector(p11),
              p01=as.vector(p01),
              p10=as.vector(p10),
              pdot=as.vector(pdot)))
}

lnl.io <- function(par,x1,x2,models){
  # Compute probabilities
  p.list <- p.io(par,x1,x2,models)
  p11 <- p.list$p11
  p01 <- p.list$p01
  p10 <- p.list$p10
  p01[p01==0] <- 1e-6
  p10[p10==0] <- 1e-6
  # Compute log-likelihood value
  lnl <- sum((1-x1$detected)*x2$detected*log(p01)) +
          sum((1-x2$detected)*x1$detected*log(p10))+
          sum(x1$detected*x2$detected*log(p11)) -
          sum(log(p.list$pdot))
  return(-lnl)
}
