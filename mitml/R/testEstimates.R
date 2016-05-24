testEstimates <- function(model, qhat, uhat, var.comp=FALSE, df.com=NULL){
# combine scalar estimates from the analysis of multiply imputed data

  if(missing(model)==(missing(qhat)|missing(uhat))) stop("Either 'model' or both 'qhat' and 'uhat' must be supplied.")

  # ***
  # combine from matrix or list
  #
  if(!missing(qhat)){

    coef.method <- "default"
    if(missing(uhat)) stop("Either 'model' or both 'qhat' and 'uhat' must be supplied.")

    if(is.list(qhat)){
      Qhat <- sapply(qhat, identity)
      Uhat <- sapply(uhat, identity)
    }else{
      Qhat <- qhat
      Uhat <- uhat
    }
    if(is.null(dim(Qhat))){ 
      dim(Qhat) <- c(1,length(qhat))
      nms <- if(is.list(qhat)) names(qhat[[1]]) else if(is.matrix(qhat)) rownames(qhat) else NULL
      dimnames(Qhat) <- list(nms, NULL)
    }
    if(is.null(dim(Uhat))) dim(Uhat) <- dim(Qhat)
    m <- ncol(Qhat)
    if(is.null(rownames(Qhat))) rownames(Qhat) <- paste0("Parameter.",1:nrow(Qhat))

    Qbar <- apply(Qhat,1,mean)
    Ubar <- apply(Uhat,1,mean)
    B <- apply(Qhat,1,var)
    T <- Ubar + (1+m^(-1)) * B

    r <- (1+m^(-1))*B/Ubar 
    v <- vm <- (m-1)*(1+r^(-1))^2
    fmi <- (r+2/(v+3))/(r+1)

    se <- sqrt(T)
    t <- Qbar/se

    if(!is.null(df.com)){
      lam <- r/(r+1)
      vobs <- (1-lam)*((df.com+1)/(df.com+3))*df.com
      v <- (vm^(-1)+vobs^(-1))^(-1)
    }

    p <- 1-pt(abs(t),df=v)
    out <- matrix(c(Qbar,se,t,v,p,r,fmi),ncol=7)
    colnames(out) <- c("Estimate","Std.Error","t.value","df","p.value","RIV","FMI")
    rownames(out) <- names(Qbar)

    # print vout for missing U
    uind <- is.na(Ubar)
    vout <- if(any(uind)) out[uind,1,drop=FALSE] else NULL
    out <- if(all(uind)) NULL else out[!uind,,drop=FALSE]

  }

  # ***
  # combine through model recognition
  #
  if(!missing(model)){

    # * identify procedures for model class
    cls <- class(model[[1]])
    coef.method <- vc.method <- "default"

    if(cls[1]=="lm") vc.method <- "lm"
    if(any(grepl("merMod",cls)) & coef.method=="default"){
      if(!requireNamespace("lme4", quietly=TRUE)) stop("The 'lme4' package must be installed in order to handle 'merMod' class objects.")
      coef.method <- vc.method <- "lmer"
    }
    if(any(grepl("^.?lme$",cls)) & coef.method=="default"){
      if(!requireNamespace("nlme", quietly=TRUE)) stop("The 'nlme' package must be installed in order to handle 'lme' class objects.")
      coef.method <- vc.method <- "nlme"
    }

    # * combine fixed coefficients
    fe <- switch(coef.method,
      lmer=.getCOEF.lmer(model,diagonal=TRUE),
      nlme=.getCOEF.nlme(model,diagonal=TRUE),
      default=.getCOEF.default(model,diagonal=TRUE)
    )

    m <- length(model)
    Qhat <- fe$Qhat
    Uhat <- fe$Uhat
    if(is.null(dim(Qhat))){ 
      dim(Qhat) <- c(1,m)
      dim(Uhat) <- c(1,m)
      dimnames(Qhat) <- dimnames(Uhat) <- list(fe$nms, NULL)
    }

    Qbar <- apply(Qhat,1,mean)
    Ubar <- apply(Uhat,1,mean)
    B <- apply(Qhat,1,var)
    T <- Ubar + (1+m^(-1)) * B

    r <- (1+m^(-1))*B/Ubar 
    v <- vm <- (m-1)*(1+r^(-1))^2
    fmi <- (r+2/(v+3))/(r+1)

    se <- sqrt(T)
    t <- Qbar/se

    if(!is.null(df.com)){
      lam <- r/(r+1)
      vobs <- (1-lam)*((df.com+1)/(df.com+3))*df.com
      v <- (vm^(-1)+vobs^(-1))^(-1)
    }

    p <- 1-pt(abs(t),df=v)
    out <- matrix(c(Qbar,se,t,v,p,r,fmi),ncol=7)
    colnames(out) <- c("Estimate","Std.Error","t.value","df","p.value","RIV","FMI")
    rownames(out) <- names(Qbar)

    # * combine variance components
    vout <- NULL
    if(var.comp){

      vc <- switch(vc.method,
        lmer=.getVC.lmer(model),
        nlme=.getVC.nlme(model),
        lm=.getVC.lm(model),
        default=list(vlist=NULL,addp=NULL)
      )
      if(vc.method=="default") warning("Computation of variance components not supported for objects of class '", cls[1], "' (see ?with.mitml.list for manual calculation).")
   
      vlist <- vc$vlist
      addp <- vc$addp

      if(!is.null(vlist)){
        vlist <- lapply(vlist, function(z) apply(z,1:2,mean) )
        ln <- names(vlist)
        nms <- vout <- c()
        for(vv in 1:length(vlist)){
          vc <- vlist[[vv]]
          rn <- rownames(vc)
          cn <- colnames(vc)
          for(rr in 1:nrow(vc)){
          for(cc in 1:ncol(vc)){
            if(cc>=rr){
              vout <- c(vout, vc[rr,cc])
              nms <- c(nms, paste(rn[rr],"~~",cn[cc],ln[vv],sep=""))
            }
          }}
        }
      }
      if(!is.null(vout)){
        vout <- matrix(vout,ncol=1)
        colnames(vout) <- "Estimate"
        rownames(vout) <- nms
        if(!is.null(addp)) vout <- rbind(vout, as.matrix(addp))
      }
    }

  }
  
  out <- list(
    call=match.call(),
    estimates=out,
    var.comp=vout,
    m=m,
    adj.df=!is.null(df.com),
    df.com=df.com,
    cls.method=coef.method
  )
  class(out) <- "mitml.testEstimates"
  out

}

