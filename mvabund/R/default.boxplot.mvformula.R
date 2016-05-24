################################################################################
# BOXPLOT.mvformula: Display the abundances of several species and their       #
# relation to explanatory variables in boxplots                                #
################################################################################
default.boxplot.mvformula <- function (x, data = NULL, subset,
  main="", xlab=NULL, ylab="Abundances", col="white", fg= "grey", 
  pch=20, las=NULL, n.vars=12, overall.main="", var.subset=NA,
  write.plot="show", filename="plot.mvabund", scale.lab=NULL, t.lab="o", 
  mfrow=c(min(4,n.vars),3), mfcol=NULL, border="black", all.labels=FALSE,
  ask=if(mixed & write.plot=="show") TRUE else FALSE, ...)
{

  miss.las <- missing(las) | is.null(las)
  m <- match.call(expand.dots = FALSE)
  dots <- m$...
  if(!is.null(dots$log)){
      dots$log <- NULL
      warning("argument 'log' not implemented in 'boxplot.mvformula'")
  }
  m$... <- NULL
  m[[1]]<-NULL
  if (length(m)>0){
	mf <- eval(m, parent.frame())
	mf$x <- NULL
	mf$las <- NULL
  }
 if(!is.null(dots$horizontal)) mf$all.labels <- TRUE
 
 if (length(dots)>0){
	dots <- lapply(dots, eval, parent.frame())
 }
	
 if ("formula" %in% names(mf)){
	x <- mf$formula	
	mf$formula<- NULL
 } else if (missing(x))
        stop("'formula' x is missing")
	  
    if ((length(x) != 3)) 
        stop("'formula' x is incorrect")

  term.labels <- eval(attr(terms(x), "variables"))

  if (attr(terms(x), "response")==1)  term.labels[[1]] <- NULL
  factor.vars <- rep(FALSE, length(term.labels))
  for(i in 1:length(factor.vars)) {
    factor.vars[i] <- is.factor(term.labels[[i]])
  }

  # Restrict to xvar.subset, if this is specified. Note that if dots$xvar.select,
  # but not xvar.subset specified, this will be ignored! xvar.subset is thus not
  # fully functional here!
  if(!is.null(dots$xvar.subset)){
    if(!is.null(dots$xvar.select)) select <- dots$xvar.select else select <- TRUE
    if(select) xvar.subset <- dots$xvar.subset else xvar.subset <- 1:length(factor.vars)
  } else xvar.subset <- 1:length(factor.vars)

  mixed <- any(factor.vars[xvar.subset]) & any(!factor.vars[xvar.subset])

  # If the form of the plot changes, i.e. a new plot will be drawn, make sure
  # to wait for user input, if wanted or by default.
  if(!is.null(ask)){
       aske <- par("ask")
       aske.new <- ask
  } else if(mixed & write.plot=="show") {
       aske <- par("ask")
       aske.new <- TRUE
  } else aske.new <- par("ask")

  mf$ask <- NULL

  if(any(factor.vars[xvar.subset])){
       if(miss.las) las <- 0
       which.factor <- which(factor.vars[xvar.subset])
  	   do.call( "plot.mvformula", c(list(x=x, type="bx",
        xvar.subset=which.factor, ask=aske.new, las=las),mf,dots) )
  }

  if(any(!factor.vars[xvar.subset])){
     if(miss.las) las <- 1
     which.notfactor <- which(!factor.vars[xvar.subset])
	   do.call( "plot.mvformula", c(list(x=x, type="bx",
        xvar.subset=which.notfactor, ask=aske.new, las=las),mf,dots) )
	}
	 if(mixed & write.plot=="show" | !is.null(ask)) par(ask=aske)

}

