#############################################################
#                                                           #
#	wle.t.test function                                 #
#	Author: Claudio Agostinelli                         #
#	E-mail: claudio@unive.it                            #
#	Date: Febraury 9, 2001                              #
#	Version: 0.1                                        #
#                                                           #
#	Copyright (C) 2001 Claudio Agostinelli              #
#                                                           #
#	Based on t.test function in                         #
#       ctest package version 1.2.0                         #
#                                                           #
#############################################################

wle.t.test <- function(x, y=NULL, alternative = c("two.sided", "less", "greater"), mu=0, paired = FALSE, var.equal = FALSE, conf.level = 0.95, boot=30, group, num.sol=1, raf="HD", smooth=0.003, tol=10^(-6), equal=10^(-3), max.iter=500) {

    alternative <- match.arg(alternative)
    wy <- NULL
    x.out <- x
    y.out <- y

    if(!missing(mu) && (length(mu) != 1 || is.na(mu)))
        stop("mu must be a single number")
    if(!missing(conf.level) &&
       (length(conf.level) != 1 || !is.finite(conf.level) ||
        conf.level < 0 || conf.level > 1))
        stop("conf.level must be a single number between 0 and 1")
    if( !is.null(y) ) {
	dname <- paste(deparse(substitute(x)),"and",
		       deparse(substitute(y)))
	if(paired)
	    xok <- yok <- complete.cases(x,y)
	else {
	    yok <- !is.na(y)
	    xok <- !is.na(x)
	}
	y <- y[yok]
    }
    else {
	dname <- deparse(substitute(x))
	if( paired ) stop("y is missing for paired test")
	xok <- !is.na(x)
	yok <- NULL
    }
    x <- x[xok]
    if( paired ) {
	x <- x-y
	y <- NULL
    }

    nx <- length(x)
    if(nx < 2) stop("not enough x observations")
    x.est <- wle.normal(x=x,boot=boot, group=group, num.sol=num.sol, raf=raf, smooth=smooth, tol=tol, equal=equal, max.iter=max.iter)
    x.tot.sol <- x.est$tot.sol
    y.tot.sol <- 1
    y.root <- 1
     
    if( !is.null(y) ) {
    ny <- length(y)
    if(ny < 2) stop("not enough y observations")
    y.est <- wle.normal(x=y,boot=boot, group=group, num.sol=num.sol, raf=raf, smooth=smooth, tol=tol, equal=equal, max.iter=max.iter)
    y.tot.sol <- y.est$tot.sol
    }

    wtstat <- matrix(0,ncol=y.tot.sol,nrow=x.tot.sol)
    wdf <- matrix(0,ncol=y.tot.sol,nrow=x.tot.sol)
    num.col <- 1
    if (!is.null(y) & paired==FALSE) num.col <- 2
    westimate <- matrix(0,nrow=x.tot.sol*y.tot.sol,ncol=num.col) 

    for (x.root in 1:x.tot.sol) {

    if (x.tot.sol==1) {
    mx <- x.est$location
    vx <- x.est$scale^2
    wx <- x.est$weights
    if (is.nan(mx) | is.nan(vx)) stop("no solutions are found for location and scale of 'x'")
    } else {
    mx <- x.est$location[x.root]
    vx <- (x.est$scale^2)[x.root]
    wx <- x.est$weights[x.root,]    
    }

    sumwx <- sum(wx)

    estimate <- mx
    if(is.null(y)) {
	df <- sumwx-1
	stderr <- sqrt(vx/sumwx)
	tstat <- (mx-mu)/stderr
	method <- ifelse(paired,"Paired wt-test for normal distributed data","One Sample wt-test for normal distributed data")
	names(estimate) <- ifelse(paired,"mean of the differences","mean of x")
    wtstat[x.root,y.root] <- tstat
    wdf[x.root,y.root] <- df
    westimate[x.root*y.root,] <- estimate

    } else {

        for (y.root in 1:y.tot.sol) {

        if (y.tot.sol==1) {
            my <- y.est$location
            vy <- y.est$scale^2
            wy <- y.est$weights
            if (is.nan(my) | is.nan(vy)) stop("no solutions are found for location and scale of 'y'")
        } else {
            my <- y.est$location[y.root]
            vy <- (y.est$scale^2)[y.root]
            wy <- y.est$weights[y.root,]    
        }

        sumwy <- sum(wy)

	method <- paste(if(!var.equal) "Welch", "Two Sample wt-test for normal distributed data")
	estimate <- c(mx,my)
	names(estimate) <- c("mean of x","mean of y")
	if(var.equal) {
	    df <- sumwx+sumwy-2
	    v <- ((sumwx-1)*vx + (sumwy-1)*vy)/df
	    stderr <- sqrt(v*(1/sumwx+1/sumwy))
	} else {
	    stderrx <- sqrt(vx/sumwx)
	    stderry <- sqrt(vy/sumwy)
	    stderr <- sqrt(stderrx^2 + stderry^2)
	    df <- stderr^4/(stderrx^4/(sumwx-1) + stderry^4/(sumwy-1))
	}
        tstat <- (mx - my - mu)/stderr
        wtstat[x.root,y.root] <- tstat
        wdf[x.root,y.root] <- df
        westimate[x.root*y.root,] <- estimate
    }
    }
    }

    result <- list()

    for (x.root in 1:x.tot.sol) {

    result.x <- list()
    for (y.root in 1:y.tot.sol) {

    tstat <- wtstat[x.root,y.root]
    df <- wdf[x.root,y.root]
    estimate <- westimate[x.root*y.root,]

    if (alternative == "less") {
	pval <- pt(tstat, df)
	cint <- c(NA, tstat + qt(conf.level, df) )
    }
    else if (alternative == "greater") {
	pval <- pt(tstat, df, lower.tail = FALSE)
	cint <- c(tstat - qt(conf.level, df), NA)
    }
    else {
	pval <- 2 * pt(-abs(tstat), df)
	alpha <- 1 - conf.level
        cint <- qt(1 - alpha/2, df)
	cint <- tstat + c(-cint, cint)
    }
    cint <- mu + cint * stderr
    names(tstat) <- "wt"
    names(df) <- "df"
    names(mu) <- if(paired || !is.null(y)) "difference in means" else "mean"
    attr(cint,"conf.level") <- conf.level
    if (x.tot.sol>1) {
       wx <- x.est$weights[x.root,]
    } else {
       wx <- x.est$weights
    }


    if (!is.null(y)) { 
       if (y.tot.sol>1) {
           wy <- y.est$weights[y.root,]
       } else {
           wy <- y.est$weights
       }
    }

    rval <- list(statistic = tstat, parameter = df, p.value = pval,
	       conf.int=cint, estimate=estimate, null.value = mu,
	       alternative=alternative,
	       method=method, data.name=dname, 
               x.weights=wx, y.weights=wy, 
               x.root=x.root, y.root=y.root)

    class(rval) <- "htest"

    result.x <-  c(result.x,list(rval))
    }
    result$test <- c(result$test,list(result.x))
    result.x <- list()
    }

    result$x.tot.sol <- x.tot.sol
    result$y.tot.sol <- y.tot.sol
    result$call <- match.call()   
    result$paired <- paired
    result$x <- x.out
    result$y <- y.out

    class(result) <- "wle.t.test"

    return(result)
}

#############################################################
#                                                           #
#	print.wle.t.test function                           #
#	Author: Claudio Agostinelli                         #
#	E-mail: claudio@unive.it                            #
#	Date: December 23, 2003                              #
#	Version: 0.1-1                                        #
#                                                           #
#	Copyright (C) 2003 Claudio Agostinelli              #
#
#   The printhtest is a copy of the print.htest function to avoid error until a better fix is used. 
#                                        #
#############################################################

print.wle.t.test <- function(x, x.root="ALL", y.root="ALL", digits = 4, quote = TRUE, prefix = "", ...) {


printhtest <-
function(x, digits = 4, quote = TRUE, prefix = "", ...)
{
    cat("\n")
    writeLines(strwrap(x$method, prefix = "\t"))
    cat("\n")
    cat("data: ", x$data.name, "\n")
    out <- character()
    if(!is.null(x$statistic))
        out <- c(out, paste(names(x$statistic), "=",
                            format(round(x$statistic, 4))))
    if(!is.null(x$parameter))
        out <- c(out, paste(names(x$parameter), "=",
                            format(round(x$parameter, 3))))
    if(!is.null(x$p.value))
        out <- c(out, paste("p-value =",
                            format.pval(x$p.value, digits = digits)))
    writeLines(strwrap(paste(out, collapse = ", ")))
    if(!is.null(x$alternative)) {
        cat("alternative hypothesis: ")
	if(!is.null(x$null.value)) {
	    if(length(x$null.value) == 1) {
                alt.char <-
                    switch(x$alternative,
                           two.sided = "not equal to",
                           less = "less than",
                           greater = "greater than")
		cat("true", names(x$null.value), "is", alt.char,
                    x$null.value, "\n")
	    }
	    else {
		cat(x$alternative, "\nnull values:\n")
		print(x$null.value, ...)
	    }
	}
	else cat(x$alternative, "\n")
    }
    if(!is.null(x$conf.int)) {
	cat(format(100 * attr(x$conf.int, "conf.level")),
	    "percent confidence interval:\n",
            format(c(x$conf.int[1], x$conf.int[2])), "\n")
    }
    if(!is.null(x$estimate)) {
	cat("sample estimates:\n")
	print(x$estimate, ...)
    }
    cat("\n")
    invisible(x)
}
  
x.tot.sol <- x$x.tot.sol
y.tot.sol <- x$y.tot.sol

if (x.root!="ALL" & !is.numeric(x.root)) {
    stop("Please, choose one 'x' root, for print all root ALL")
} else if (x.root=="ALL") {
    x.root <- 1:x.tot.sol
} else if (x.tot.sol<x.root) {
    stop(paste("'x' Root ",x.root," not found"))
}

if (y.root!="ALL" & !is.numeric(y.root)) {
    stop("Please, choose one 'y' root, for print all root ALL")
} else if (y.root=="ALL") {
    y.root <- 1:y.tot.sol
} else if (y.tot.sol<y.root) {
    stop(paste("'y' Root ",y.root," not found"))
}

cat("\nCall:\n")
cat(paste(deparse(x$call), sep="\n", collapse = "\n"), "\n\n", sep="") 
cat("\n\nWeighted t test:\n", sep="")

for (xx.root in x.root) {
for (yy.root in y.root) {
    cat("\n'x' Root ",xx.root)
    if (!is.null(x$y) & x$paired==FALSE) cat (" 'y' Root ",yy.root)    
    printhtest(x$test[[xx.root]][[yy.root]], digits=digits, quote=quote, prefix=prefix, ...)
}
}
 
    cat("\nNumber of 'x' solutions ",x.tot.sol,"\n")
    if (!is.null(x$y) & x$paired==FALSE) cat("\nNumber of 'y' solutions ",y.tot.sol,"\n")
    cat("\n")
    invisible(x)

}


