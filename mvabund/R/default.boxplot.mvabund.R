################################################################################
# BOXPLOT.MVABUND: Display the abundances of several species in one plot or    #
# compare them for two objects                                                 #
################################################################################
default.boxplot.mvabund <-
function (x, y, range = 1.5, width = NULL,
    # varwidth = FALSE,
    notch = FALSE, outline = TRUE, names,
    # plot = TRUE,
    border = par("fg"), col = NULL,
    # log = "",
    pars = list(boxwex = 0.8, staplewex = 0.5, outwex = 0.5),
    # horizontal = FALSE, add = FALSE,
    at = NULL, xlab, ylab, main, pch=1, fg= "grey", las=1,
    write.plot="show", filename="plot.mvabund", n.vars=min(12,NCOL(x)),
    overall.main, var.subset=NA, transformation="log", scale.lab="s", t.lab="o",
    mfrow=1, mfcol=NULL,
    # checks = TRUE, line.col=NULL,
    ask=if(write.plot=="show"&(is.list(x)|!is.null(dotsnotnamed))) TRUE else FALSE, 
    ...)
{
    if (missing(x)) stop("The mvabund object 'x' is missing.")
    allargs 	<- match.call(expand.dots = FALSE)
    dots  <- allargs$...
    if(!is.null(dots$log)){
       dots$log <- NULL
       warning("argument 'log' not implemented in 'boxplot.mvabund'")
    } 
    allargs$... <- NULL
    allargs$x <- NULL
    allargs$y <- NULL
    allargs[[1]]<- NULL
    allargs <- lapply(allargs, eval, parent.frame())

    if(!is.null(dots)){
	# If args have names: namedargs are TRUE, else FALSE.
	 namedargs <- if (!is.null(attributes(dots)$names)) attributes(dots)$names != ""
                      else rep(FALSE, length.out = length(dots))
 	 dotsnamed <- dots[namedargs]
   	 dotsnamed <- lapply(dotsnamed, eval, parent.frame())
   	 if (any(!namedargs)) {
	    dotsnotnamed <- dots[!namedargs] 
            dotsnotnamed <- lapply(dotsnotnamed, eval, parent.frame())
         } else dotsnotnamed <- NULL
    } else {
	 dotsnotnamed 	<- NULL
	 dotsnamed    	<- NULL
    }

    aske <- par("ask")
    on.exit(par(ask=aske), add = TRUE)
    par(ask=ask)

    if(is.list(x)){
        # Plot each element of the list using several calls to plot.mvabund,
        # it is assumed that no y and no further objects to plot are passed.
        for (i in 1:length(x))
           do.call("plot.mvabund", c(list(x=x[[i]],type="bx"),dotsnamed,allargs),quote=FALSE)
        if( write.plot == "show" ) par(ask = aske)		
    } else {
       # Plot x and each of the unnamed dots values (which are not allowed to be a list),
       # by using several calls to 'plot.mvabund'.
        if (!is.null(dotsnotnamed)){
           for(i in 1:length(dotsnotnamed)){
              do.call( "plot.mvabund", c(list(x=dotsnotnamed[[i]], type="bx"), dotsnamed,allargs), quote=FALSE)
           }           
           if(write.plot == "show" ) par(ask = aske)
        } else {
           # Plot x and an optional y.
           if(!missing(y)) 
	      do.call( "plot.mvabund", c(list(x=x, y=y, type="bx"),dotsnamed,allargs),quote=FALSE)
	   else
	      do.call( "plot.mvabund", c(list(x=x, type="bx"),dotsnamed,allargs), quote=FALSE)	
       }
   }
} 
	
