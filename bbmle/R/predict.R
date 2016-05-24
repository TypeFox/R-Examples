setGeneric("simulate", function(object, nsim=1, seed=NULL, ...) standardGeneric("simulate"))
setMethod("simulate", "mle2",
          function(object, nsim=1, seed, newdata=NULL,
                   newparams=NULL, ...) {
            if (!is.null(seed)) set.seed(seed)
            if (!is.null(newparams)) {
              object@fullcoef <- newparams
            }
            g <- gfun(object,newdata=newdata, nsim=nsim,op="simulate")
            if (nsim>1) {
              g <- matrix(g,ncol=nsim)
            }
            g
          })

setGeneric("predict", function(object, ...) standardGeneric("predict"))
setMethod("predict", "mle2",
          function(object,newdata=NULL,
                   location="mean",newparams=NULL, ...) {
            if (!is.null(newparams)) {
              object@fullcoef <- newparams
            }
            gfun(object,newdata=newdata,location=location,op="predict")
          })

setGeneric("residuals", function(object, ...) standardGeneric("residuals"))
setMethod("residuals", "mle2",
          function(object,
                   type=c("pearson","response"),
                   location="mean",
                   ...) {
            type <- match.arg(type)
            location <- match.arg(location)
            pred <- predict(object,location)
            ## not sure this will work ...
            obs <- with(object@data,
                        get(gsub("~.+","",object@formula)))
            res <- obs-pred
            if (type=="response") return(res)
            vars <- predict(object,location="variance")
            return(res/sqrt(vars))
          })

## general-purpose function for simulation and
##  prediction (the hard part is evaluating the parameters etc.)
##
gfun <- function(object,newdata=NULL,location=c("mean","median","variance"),
                 nsim,
                 op=c("predict","simulate")) {
  ## notes: should operate on formula
  ## pull out call$formula (not character)
  location <- match.arg(location)
  if (class(try(form <- as.formula(object@call$minuslogl)))!="formula")
    stop("can only use predict() if formula specified")
  LHS <- form[[3]]
  ddist = as.character(LHS[[1]])
  spref <- switch(op,predict="s",simulate="r")
  sdist = gsub("^d",spref,ddist)
  arglist = as.list(LHS)[-1]
  if (!exists(sdist) || !is.function(get(sdist)))
    stop("function ",sdist," does not exist")
  ## evaluate parameters
  ## evaluate sdist [ newdata > coef > data ]
##   if (is.null(object@data)) {
##     comb <- newdata
##   } else {
##     nmatch <- match(names(newdata),names(object@data))
##     comb <- object@data
##     comb[na.omit(nmatch)] <- newdata[!is.na(nmatch)]
##     comb <- c(comb,newdata[is.na(nmatch)])
##   }
##   comb <- c(newdata,object@data)
##   comb <- comb[!duplicated(names(comb))]
##   comb <- comb[sapply(comb,length)>0]
##   rvar <- strsplit(object@formula,"~")[[1]][1]
##   comb <- comb[!names(comb)==rvar] ## remove response variable
  parameters <- eval(object@call$parameters)
  if (!is.null(parameters)) {
    vars <- as.character(sapply(parameters,"[[",2))
    models <-  sapply(parameters,function(z) call.to.char(z[[3]]))
    parameters <- parameters[models!="1"]
    npars <- length(parameters)
    if (npars==0) { ## no non-constant parameters
      parameters <- mmats <- vpos <- NULL
    } else {
      mmats <- list()
      vpos <- list()
      for (i in seq(along=parameters)) {
        vname <- vars[i]
        p <- parameters[[i]]
        p[[2]] <- NULL
        mmat <- with(c(newdata,object@data),
                     model.matrix(p,data=environment()))
        ## c(as.list(newdata),as.list(object@data)))
        pnames <- paste(vname,colnames(mmat),sep=".")
        assign(vname,mmat %*% coef(object)[pnames])
      }
    }
  }
  arglist1 <- lapply(arglist,eval,envir=c(newdata,object@data,
                                    as.list(coef(object))),
                     enclos=sys.frame(sys.nframe()))
  ## HACK: need a way to figure out how many data points there
  ##  are, in the *absence* of an explicit data argument
  ## then replicate constant values to the full length
  if (op=="simulate") {
    if (length(object@data)==0)
      stop("need explicit data argument for simulation")
    ndata <- max(sapply(c(newdata,object@data),length)) ## ???
    arglist1 <- c(arglist1,list(n=ndata*nsim))
  }
  vals <- with(as.list(coef(object)),do.call(sdist,arglist1))
  if (op=="predict") return(vals[[location]]) else return(vals)
}
