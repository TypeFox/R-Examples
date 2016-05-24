#uncertMC object
#
#Author: S Ellison
#
#Changes:
# 2013-10-24: Amended references to qqplot.default and hist.default to generic; also 
#             amended refs in formals(x.default) to use formals(getS3method("x", "default"))
#
# 2014-03-04: Amended reference to plot.density to plot
#
print.uncertMC<-function(x, digits=NULL, right=FALSE, ..., simplify=TRUE, minimise=FALSE){
        maxwidth<-12L
        cat("\nUncertainty evaluation\n\n")

        cat("Call:\n  ",deparse(x$call), sep="")
        cat("\n\n")
        cat("Expression: ")
        if(class(x$expr)=="formula" ) {
                cat(paste(x$expr, collapse=""))
        } else if(is.function(x$expr)) {
                cat( deparse(x$expr)[1] )
        } else if(is.expression(x$expr)) {
                cat( deparse(x$expr[[1]]) )
        } else if(is.na(x$expr)) {
                cat("NA")
        }

        cat("\n\n")
        cat(paste("Evaluation method: ",x$method, "\n\n"))

        cat("Budget:\n")
        dp<-x$budget[sapply(x$budget, function(x) !all(is.na(x)))]
        if(!is.null(x[["distrib", exact=TRUE]]) ) {
                distrib.labels<- as.vector(
                                sapply(x$distrib, function(x) if(is.function(x)) deparse(x)[1] else paste(x)) 
                        )
                dp$distrib<-sub(paste("(.{",maxwidth,",",maxwidth,"})(.+)", sep=""), "\\1...",distrib.labels)
        }
        if(!is.null(x$distrib.pars)) {
                dp$distrib.pars <- vector("character", length=nrow(x$budget) )
                for(nn in row.names(dp) ) {
                        dp[nn,"distrib.pars"]<-
                                paste(names(x$distrib.pars[[nn]]), format(x$distrib.pars[[nn]], digits=digits), sep="=", collapse=", ")
                }
                
        }
        print.data.frame(dp,digits=digits, right=right, ...)
        if(!is.null(x$additional) ) print(as.data.frame(x$additional), ...)

        cat("\n   y: ", format(x$y))
        cat("\nu(y): ", format(x$u.y), "\n")

        if(!simplify) {
                cat("\nCovariance matrix:\n")
                print(x$cov)
                cat("\nCorrelation matrix:\n")
                print(x$cor)
                cat("\n")
                if(!is.null(x$cov.xy)) {
                        cat("\nX-Y Covariances:\n")
                        print(x$cov.xy)
                }
                if(!is.null(x$cor.xy)) {
                        cat("\nX-Y Correlations:\n")
                        print(x$cor.xy)
                }
        }
        
        cat(sprintf("\nMonte Carlo evaluation using %d replicates:\n", x$B))
        cat("\n   y:\n")
        if(simplify) {
                print(summary(x$MC$y))
        } else {
                print(x$MC$y)
                if(!is.null(x$MC$x) ) {
                        cat("\n   x:\n")
                        print(summary(x$MC$x))
                }
        }
        invisible(x)
}

summary.uncertMC<-function(object, digits=NULL, right=FALSE, ..., simplify=TRUE, minimise=FALSE){
        print.uncertMC(object, digits=digits, right=right, ..., simplify=simplify, minimise=minimise)
}

plot.uncertMC<-function(x, which=1:2, main=paste("Monte Carlo evaluation -",deparse(substitute(x))),
                ask = prod(par("mfcol")) < length(which) && dev.interactive(), 
                caption=list("Histogram", "Q-Q plot", "Density", "Correlation x-y", "Covariance x-y"), 
                xlab=paste(deparse(substitute(x)), "$y", sep=""),
                ..., cex.caption=1.0, cex.main=1.25, lwd.y=2, col.y=2, lty.y=1,
                col.qqline=NULL, lty.qqline=NULL, lwd.qqline=NULL ) {
        
        arglist <- list(...)
        gpars<-arglist[names(arglist) %in% names( par() )]
        one.fig <- prod(par("mfcol")) == 1
        if(prod(par("mfcol"))>2) cex.caption <- cex.caption*0.8
        
        if (ask) {
                oask <- devAskNewPage(TRUE)
                on.exit(devAskNewPage(oask))
        }
        
        if( class(x)[1]=="uncertMC" ) {
                at<-NULL
                x.names<-row.names(x$budget)
                if(1 %in% which) {
                        histpars<-arglist[names(arglist) %in% names(c(formals(getS3method("hist", "default")),par()))]
                        do.call(hist, c(list(x=x$MC$y, main="", xlab=xlab), histpars))
                        if(lwd.y >= 1) abline(v=x$y, col=col.y, lwd=lwd.y, lty=lty.y)
                        mtext(caption[[1]], side = 3, line=0.25, cex=cex.caption)
                        if(one.fig) title(main=main)
                }

                if(2 %in% which) {
                        qqpars<-arglist[names(arglist) %in% names(c(formals(getS3method("qqnorm", "default")),par()))]
                        if(is.null(qqpars$datax)) qqpars$datax = TRUE
                        do.call(qqnorm, c(list(y=x$MC$y,  main=""), qqpars))
                        
                        qqlpars<-arglist[names(arglist) %in% names(c(formals(qqline),par()))]
                        if(is.null(qqlpars$datax)) qqlpars$datax = TRUE
                        if(!is.null(col.qqline)) qqlpars$col<-col.qqline
                        if(!is.null(lty.qqline)) qqlpars$lty<-lty.qqline
                        qqlpars$y<-x$MC$y
                        do.call(qqline, qqlpars)
                        mtext(caption[[2]], side = 3, line=0.25, cex=cex.caption)
                        if(one.fig) title(main=main)
                }
                
                if(3 %in% which) {
                        dpars<-arglist[names(arglist) %in% names(formals(getS3method("density", "default")))]
                        dpars$x<-x$MC$y
                        d<-do.call(density, dpars)
                        
                        dppars<-arglist[names(arglist) %in% 
                                unique(names(c(formals(getS3method("plot", "default")), 
                                        formals(getS3method("plot", "density")), par())))]
                        #do.call(plot.density, c(list(x=d, main=""), dppars))
                        	#Explicit call removed - SLRE
                        do.call(plot, c(list(x=d, main=""), dppars))
                        if(lwd.y >= 1) abline(v=x$y, col=col.y, lwd=lwd.y, lty=lty.y)
                        mtext(caption[[3]], side = 3, line=0.25, cex=cex.caption)
                        if(one.fig) title(main=main)
                }

                if(4 %in% which) {
                        corpars<-arglist[names(arglist) %in% names(formals(cor))]
                        cor.xy<-function(x, y, p) {
                                    do.call(cor, c(list(x=x, y=y), p))
                        }
                        
                        if(!is.null(corpars$method)) {
                                c.methods<- eval(formals(cor)$method)
                                c.m<-c.methods[pmatch(corpars$method,c.methods)][1]
                                if(is.na(c.m)) {
                                        warning(sprintf("Correlation method %s not found: using pearson",corpars$method) )
                                        corpars$method<-"pearson"
                                } else corpars$method <-c.m
                        } else {
                                corpars$method<-"pearson"
                        }
                        
                        if(corpars$method %in% row.names(x$cor.xy)) {
                                xycor<- unlist(x$cor.xy[corpars$method, ])
                        }  else {
                                if(!is.null(x$MC$x)) 
                                    xycor<-sapply(x$MC$x, cor.xy, y=x$MC$y, p=corpars)
                                else xycor<-NULL
                        } 
                        
                        if( !is.null(xycor) ) {
                                names.arg<-if(is.null(names(xycor))) x.names else names(xycor)

                                barpars<-arglist[names(arglist) %in% names(c(formals(getS3method("barplot", "default")), par()))]
                                at<-do.call(barplot, c(list(height=as.vector(xycor), names.arg=names.arg), barpars))
                                corMethod<-paste(toupper(substring(corpars$method, 1,1)), 
                                        substring(corpars$method, 2), sep="", collapse=" ")

                                c4<-if( caption[[4]] == eval(formals(plot.uncertMC)$caption)[[4]] )
                                        paste(corMethod,caption[[4]]) 
                                      else caption[[4]]
                                mtext(c4, side = 3, line=0.25, cex=cex.caption)
                                if(one.fig) title(main=main)
                        } else {
                                warning(sprintf("Missing %s$MC$x and $cor.xy; correlation plot not made", deparse(substitute(x))), call.=TRUE)
                        }
                }

                if(5 %in% which) {
                        covpars<-arglist[names(arglist) %in% names(formals(cov))]
                        cov.xy<-function(x, y, p) {
                                    do.call(cov, c(list(x=x, y=y), p))
                        }
                        
                        if(!is.null(covpars$method)) {
                                c.methods<- eval(formals(cov)$method)
                                c.m<-c.methods[pmatch(covpars$method,c.methods)][1]
                                if(is.na(c.m)) {
                                        warning(sprintf("Covariance method %s not found: using pearson",covpars$method) )
                                        covpars$method<-"pearson"
                                } else covpars$method <-c.m
                        } else {
                                covpars$method<-"pearson"
                        }
                        
                        if(covpars$method %in% row.names(x$cov.xy)) {
                                xycov<- unlist(x$cov.xy[covpars$method, ])
                        }  else {
                                if(!is.null(x$MC$x)) 
                                    xycov<-sapply(x$MC$x, cov.xy, y=x$MC$y, p=covpars)
                                else xycov<-NULL
                        } 
                        
                        if( !is.null(xycov) ) {
                                names.arg<-if(is.null(names(xycov))) x.names else names(xycov)

                                barpars<-arglist[names(arglist) %in% 
                                        names(c(formals(getS3method("barplot", "default")), par()))]
                                at<-do.call(barplot, c(list(height=as.vector(xycov), names.arg=names.arg), barpars))
                                covMethod<-paste(toupper(substring(covpars$method, 1,1)), 
                                        substring(covpars$method, 2), sep="", collapse=" ")

                                c5<-if( caption[[5]] == eval(formals(plot.uncertMC)$caption)[[5]] )
                                        paste(covMethod,caption[[5]]) 
                                      else caption[[5]]
                                mtext(c5, side = 3, line=0.25, cex=cex.caption)
                                if(one.fig) title(main=main)
                        } else {
                                warning(sprintf("Missing %s$MC$x and $cov.xy; covariance plot not made", deparse(substitute(x))), call.=TRUE)
                        }
                }

                if (!one.fig )  {         
                        mtext(main, outer = TRUE, cex = cex.main, 
                                line=if(par("oma")[3L] >= 1) 0 else -1.5)
                }
                
                invisible(NULL)
                
        } else {        
                stop(paste(deparse(substitute(x))), "is not an 'uncertMC' object")
        }       
        
}

