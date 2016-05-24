#Functions for Monte Carlo 
#
#Author: S Ellison
#
#Changes:
#
# 2014-03-04  Removed "require(MASS)" (deprecated) - SLRE
#
uncertMC<-function(expr, x, u, method="MC", df, cor, cov, distrib, distrib.pars, B=200, keep.x=TRUE,  vectorized=TRUE, ...) {

        if(method != "MC") {
                #This mostly covers calls from update()
                rv<-uncert(expr, x=x, u=u, method=method, cor=cor, cov=cov, distrib=distrib, 
                     distrib.pars=distrib.pars, B=B, keep.x=keep.x, ...)
                rv$call <- match.call()
                return(rv)
        }
        
        L <- length(x) 
        
        names.match <- function( names1, names2) {
                all( (names1 %in% names2) ) && all( (names2 %in% names1) )      
        }

        #Make sure everything is a list, 'cos we'll need that later
        x <- as.list(x)
        if(is.null(names(x))) names(x)<-paste("X", 1:L )

        
        #Check for mismatches in variable and parameter names
        var.names<-if(class(expr) == "function" ) 
                                names(formals(expr)) 
                        else all.vars(expr)

        par.names <- names( c(x, ...) )
        
        if( !names.match(var.names, par.names) ) {
                stop("Variables in expr do not match arguments in x and '...'", call.=TRUE )
        }
        
        if(missing(u) && missing(cov)) stop("Either u or cov must be present", call.=TRUE)

        if(!missing(u) && !missing(cov)) warning("Only one of u and cov should be specified: using cov")

        if(!missing(u)) {
                if(length(u) != L) stop("Lengths of x and u do not match", call.=TRUE)
                u<-as.list(u)
                if(is.null(names(u))) names(u)<-names(x)
                   else if( !names.match( names(x), names(u) ) ) {
                        stop("Names in x and u do not match", call.=TRUE )
                }
        } else {
                u<-sqrt(diag(cov))
                names(u)<-names(x)
                u<-as.list(u)
        }

        
        if( missing(df) ) {
                df<-as.list(rep(NA, L))
                names(df)<-names(x)
        } else  {
                df<-as.list(df)
                if(is.null(names(df))) {
                        if( length(df)!=L ) stop("Lengths of x and df do not match and df does not have names", call.=TRUE)
                        names(df) <- names(x)
                } else {
                        if( any( !( names(df) %in% names(x) ) ) ) 
                            stop("df contains names not in names(x)", call.=TRUE)
                }
        }
        
        if(missing(distrib) ) distrib<-NULL
        if(missing(distrib.pars) ) distrib.pars<-NULL
        
        if( is.null(distrib)) {
                distrib<-as.list(rep("norm", L))
                names(distrib)<-names(x)
        } else {
                distrib<-as.list(distrib)
                if(is.null(names(distrib))){
                        if(length(distrib)==L ) {
                                names(distrib)<-names(x)
                        } else {
                                stop("Names missing from distrib and cannot be allocated", .call=TRUE)
                        }
                } else {
                 if( any( !(names(distrib) %in% names(x)) ) )
                            stop("distrib contains names not in names(x)", call.=TRUE)
                }
        }

        if(is.null(distrib.pars)) {
                distrib.pars <- list()
        } else {
                if(is.null(names(distrib.pars)) && length(distrib.pars)==L) {
                        names(distrib.pars)<-names(x)
                } else {
                        stop("Names missing from distrib.pars and cannot be allocated", call.=TRUE)
                }
                 if( any( !(names(distrib.pars) %in% names(x)) ) )
                            stop("distrib.pars contains names not in names(x)", call.=TRUE)
        }
        

                
        
        
        #Finally, it's time to actually do the work...
        if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) 
                runif(1)
        seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)

        if( missing(cor) && missing(cov) ) {
                use.cov <- FALSE
                cov <- NULL
                cor <- diag(1, L)
        } else {
                if(missing(cor)) {
                        if(missing(cov)) cor <- diag(1, L)
                        else cor <- cov / outer( sqrt(diag(cov)), sqrt(diag(cov)),"*" )
                } 

                if(any(abs(cor)>1)) stop("cor contains values outside (-1,1)", call.=TRUE)


                if(missing(cov)) {
                        u <- unlist(u)
                        cov<-outer(u,u,"*")*cor
                } 
                
                if(sum(abs(cor)) - sum(diag(abs(cor))) > .Machine$double.eps*nrow(cor)^2) 
                        use.cov=TRUE
                else
                        use.cov=FALSE
        }
        
        #Back-calculate parameters for distrib[i] to give mean x and var u
        
        if(length(distrib.pars) < L) {
                for( nn in names(x) ) {
                        if(is.null(distrib.pars[[nn]]) || is.na(distrib.pars[[nn]])) 
                               distrib.pars[[nn]]<-.get.pars(distrib[[nn]], x[[nn]], u[[nn]], df[[nn]])
                }
        }
        
        if(use.cov) {
                #Pull out correlated variables
                cor.vars<-which(rowSums(cor)> (1+.Machine$double.eps*nrow(cor)))
                cor.vars.names<-names(x)[cor.vars]
                #Check distribs of these are Normal
                if( any(sapply(distrib, function(d) d!="norm"))) 
                        stop("Correlation is only supported for normal distributions.", call.=TRUE)
                
                #use mvrnorm from MASS to handle correlations
                dfx<-as.data.frame( mvrnorm(n=B, mu=unlist(x[cor.vars]), Sigma=cov[cor.vars,cor.vars]) )
                names(dfx)<-cor.vars.names
                
                #Add any other variables
                uncor.vars<- (1:L)[!(1:L %in% cor.vars)] 
                if(length(uncor.vars)>0) {
                        for(nn in names(x)[uncor.vars] ) {
                                dfx[[nn]]<-do.call(what=paste("r",distrib[[nn]], sep=""), c(list(n=B),distrib.pars[[nn]] ) )
                        }
                }
                dfx<-dfx[, order( c(cor.vars, uncor.vars) )] 
                        #Works even with one empty vector because c(vec, integer(0)) works

        } else {
                dfx<-list()
                for(nn in names(x) ) {
                        dfx[[nn]]<-do.call(what=paste("r",distrib[[nn]], sep=""), c(list(n=B),distrib.pars[[nn]] ) )
                }
                dfx<-as.data.frame(dfx)
        }
        
        constants<-list(...)
        if( is.expression(expr) ) {     
                y0<-eval(expr, c(x, constants))
                if(!vectorized) {
                        y<-apply(dfx, 1, function(Row,const) eval(expr, c(as.list(Row), const)), const=constants) 
                } else {
                        y<-eval(expr, c(dfx, constants))
                }
        } else if( class(expr)=="formula" ) {
            if ((le <- length(expr)) > 1) {
                y0<-eval(expr[[2]], c(x, constants))
                if(!vectorized) {
                        y<-apply(dfx, 1, function(Row,const) eval(expr[[2]], c(as.list(Row), const)), const=constants) 
                } else {
                        y<-eval(expr[[2]], c(dfx, constants))
                }
            } else stop("Invalid formula in uncertMC")
                
        } else if( is.function(expr) ) {
                y0<-do.call(expr, c(x,...))
                if(!vectorized) {
                        y<-apply(dfx, 1, function(Row, additional) do.call(expr, c(as.list(Row), additional)), additional=constants) 
                } else {
                        y <- do.call(expr, c(dfx,constants))
                }
        }
        
        c.est <- coef( lm(y~., data=dfx) )[-1]
        
        mc<-list(seed=seed, y=y)
        if(keep.x) mc$x <- dfx
        
        #cor.xy<-list()
        #for(cor.method in eval(formals(stats::cor)$method)) {
        #        cor.xy[[cor.method]]<-sapply(dfx, stats::cor, y=mc$y, method=cor.method)
        #}
        #cor.xy<-as.data.frame(t(as.data.frame(cor.xy)))
        cor.xy<-as.data.frame(t(stats::cor(dfx, mc$y, method="pearson")))
        row.names(cor.xy) <- "pearson"
        
        #cov.xy<-list()
        #for(cov.method in eval(formals(stats::cor)$method)) {
        #        cov.xy[[cov.method]]<-sapply(dfx, stats::cov, y=mc$y, method=cov.method)
        #}
        #cov.xy<-as.data.frame(t(as.data.frame(cov.xy)))
        cov.xy<-as.data.frame(t(stats::cov(dfx, mc$y, method="pearson")))
        row.names(cov.xy) <- "pearson"

        rv<-.construct.uncert(call.uncert=match.call(), expr=expr, y=y0, u.y=sd(y), method="MC", 
                x=x, df=df[names(x)], u=u[names(x)], ci=c.est, u.c=NULL, x.names=names(x), 
                cov=cov, cor=cor, deriv=NULL, distrib=distrib, distrib.pars=distrib.pars, 
                cov.xy=cov.xy, cor.xy=cor.xy, delta=NULL, MC=mc, B=B, keep.x=keep.x, ...)

        class(rv) <- c("uncertMC", class(rv))
        
        return(rv)
        
}

.apply.expr<-function(expr, x, ...) {

        constants<-list(...)
        if( is.expression(expr) ) {     
                y<-apply(x, 1, function(Row,const) eval(expr, c(as.list(Row), const)), const=constants) 
        } else if( class(expr)=="formula" ) {
                if ((le <- length(expr)) > 1) {
                y<-apply(x, 1, function(Row,const) eval(expr[[2]], c(as.list(Row), const)), const=constants) 
                } else stop("Invalid formula in uncertMC")
                
        } else if( is.function(expr) ) {
                y<-apply(x, 1, function(Row, additional) do.call(expr, as.list(Row, additional))) 
        }
        return(y)

}
