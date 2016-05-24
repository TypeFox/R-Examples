################################################################################
# BEST.R.SQ: Do a forward selection of a multivariate linear model with        #
# formula considering only the response variables given by var.subset models   #
# (or all if var.subset not specified)                                         #
# at the present stage this function treats interactions just like any         #
# other variable                                                               
#
# Modified by DW 31/7/14 to return R^2 as well, bug found in steps 2 onwards
# where some potential predictors were accidently cut out of the loop
#
################################################################################

best.r.sq <- function (formula, data = parent.frame(), subset, var.subset,
  n.xvars= min(3, length(xn)), R2="h", ...){

	if(missing(formula)) stop("'formula' is missing")

  if(missing(subset)) subset <- NULL
	# Get the model.frame using the data argument.

  m <- match.call(expand.dots = FALSE)

  if (is.matrix(eval(m$data, parent.frame())))
      m$data <- as.data.frame(data)
  dots <- m$...
  dots <- lapply(dots, eval, data, parent.frame())

  m$var.subset <- m$n.xvars <- m$split.x <- NULL
  m$R2 <- NULL
  m$subset <- NULL
  # m$formula <- foo.new
  m[[1]] <- as.name("model.frame")
  mlm <- as.call(c(as.list(m), list(na.action = NULL)))
  mfm <- eval(mlm, parent.frame())
  mfm.list <- as.list(mfm)
	response <- attr(terms(formula),"response")

  intercept   <- attr(terms(formula),"intercept")
  
  new.foo     <- extend.x.formula(formula, return.interaction=TRUE,
       extend.term=TRUE)
  foo.new     <- new.foo$formula
  is.interact <- new.foo$is.interaction

  terms.foo   <- terms(foo.new)

  xn          <- attr(terms.foo,"term.labels")
  xn.no       <- 1:length(xn)
  
  xn <- xn[!is.interact]  # at the moment: leave out interactions!
  xn.no <- xn.no[!is.interact]
  
  l.vars      <- length(xn)
  
  if(n.xvars>l.vars)
    warning("Select 'n.xvars' <= ", l.vars, ". There are only ", l.vars ,
    " influence variables (excluding Interactions) in the model.")
  n.xvars     <- min(n.xvars, l.vars)
  if (n.xvars == 0) return(NULL)

  if(response!=1) stop("no regular 'formula': response is missing") else
    respname  <- foo.new[[2]]

 # varnames <- names(mfm)
  y  <- as.matrix(mfm[[response]])
    # would be better if the original name could be kept,
    # otherwise something might be overwritten

	if (l.vars == 0) return(NULL)

   # if(missing(var.subset))   new.respname <- deparse(foo.new[[2]])
   #  else new.respname <- paste(deparse(foo.new[[2]]),"[,var.subset]", sep="")
   
   if(!missing(var.subset)) {
      if(is.logical(var.subset) & any(!is.na(var.subset)))
        var.subset <- which(var.subset[!is.na(var.subset)])
      if(any(var.subset<1) | any(var.subset>ncol(y) ))
        stop("'var.subset' must be a vector with values between 1 and",  ncol(y))
      y  <- y[,var.subset, drop=FALSE]
   }
   mfm.list[[1]] <- y
   new.respname <- names(mfm.list)[1] <- "y" # names(mfm)[1] # deparse(foo.new[[2]])


  basic.foo <-  reformulate( xn , response = new.respname )
  attr(mfm.list, "terms") <- NULL # attr(model.frame(basic.foo),"terms")
  basic.foo <- formula(terms(basic.foo)[-(1:length(xn))])
  if(intercept==0)  basic.foo <- update( basic.foo,  ~ . + 0)

   use.kk <-  use.k <- 1:length(xn)
   use.xn <- xn
   order.r.sq <- numeric(n.xvars)
   r.sq.step  <- numeric(n.xvars)
   r.sq.mat <- matrix(NA,l.vars,n.xvars)
   dimnames(r.sq.mat)[[2]]=paste("Step",1:n.xvars)
   dimnames(r.sq.mat)[[1]]=xn 
 
  tm.labs <- attr(terms(basic.foo), "term.labels" )

  for(xchoose in 1:n.xvars){

	r.squared <- numeric( length(xn) )
	
	for (i in 1:length(use.xn) ) {

      use.foo <- reformulate( termlabels = c(tm.labs,
         use.xn[i] ) , response = new.respname )
      # attr(mfm.list, "terms") <- attr(model.frame(use.foo),"terms")

      if(intercept==0)  use.foo <- update( use.foo,  ~ . + 0)
      
      lmy <- do.call( "manylm", c(list( formula = use.foo,
         subset = subset, data=mfm.list), dots) )

			k = use.k[i]
      r.squared[k] <- suppressWarnings( summary(lmy, R2=R2 )$r.squared )
	}
  r.sq.mat[,xchoose] = r.squared
	r.sq.mat[-use.k,xchoose] = NA
  
  order.r.sq[xchoose] <-
      order(r.squared, decreasing=TRUE, na.last = TRUE)[1]
   r.sq.step[xchoose] = max(r.squared)
  
   use.xn <- xn[ -(order.r.sq[1:(xchoose)]) ]

   basic.foo <- reformulate( termlabels = c(tm.labs,
      xn[order.r.sq[xchoose]] ) , response = new.respname )

   if(intercept==0)  basic.foo <- update( basic.foo,  ~ . + 0)
   
   # Now that basic.foo has a term / terms the update tm.labs.
   tm.labs <- attr(terms(basic.foo), "term.labels" )
   
   use.k     <-  use.kk[-(order.r.sq[1:(xchoose)])]

  }
  step.names =  xn[xn.no[order.r.sq]] #naming r.sq.steps by the term that is entered
  names(r.sq.step) = paste(c("",rep("+",length(step.names)-1)),step.names,sep="") 
   
 
  res = list(xs=xn.no[order.r.sq], r2Step =r.sq.step, r2Matrix = r.sq.mat)  
  return(res)
}

