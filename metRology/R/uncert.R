
uncert<-function(obj, ...){
        UseMethod("uncert")
}

uncert.default<-function(obj, c, method=c("GUM", "MC"), cor, cov, distrib=NULL, 
                distrib.pars=NULL, B=200, x=NULL, keep.x=TRUE, u=obj, ...) {

        if(missing(u) && missing(cov)) stop("Either u or cov must be present", call.=TRUE)

        if(missing(cor)) {
                if(missing(cov)) cor <- diag(1, length(u))
                else cor <- cov / outer( sqrt(diag(cov)), sqrt(diag(cov)),"*" )
        } 
        
        if(any(abs(cor)>1)) stop("cor contains values outside (-1,1)", call.=TRUE)


        if(missing(cov)) {
                u <- unlist(u)
                cov<-outer(u,u,"*")*cor
        } else {
                if(!missing(u)) warning("Only one of u and cov should be specified: using cov")
        }

        method <- match.arg(method, several.ok=TRUE)[1]

        if(method=="GUM") {
                cc <- unlist(c)
                v <- ((t(cc) %*% cov ) %*% cc)[1,1]
                cov.xy<-as.vector( cov %*% cc)
                cor.xy<-cov.xy/(u*sqrt(v))
                names(cov.xy)<-names(cor.xy)<-names(x)
                rv<-.construct.uncert(x=x, u=u, ci=cc, y=NA, u.y=sqrt(v), method=method, 
                                call.uncert=match.call(), cor=cor, cov=cov, cov.xy=cov.xy, 
                                cor.xy=cor.xy, ...)
        } else if(method=="MC") {
                xx<-if(is.null(x)) rep(0, length(u)) else x
                if(is.null(names(xx))) names(xx)<-paste("X",1:length(xx), sep="")
                cx<-c
                names(cx)<-paste("c.",names(xx), sep="")
                expr <- parse(text=paste(names(xx), names(cx), sep="*", collapse="+"))[1]
                rv<-do.call(uncertMC, c(list(expr=expr, x=xx, cor=cor, cov=cov, 
                                distrib=distrib, distrib.pars=distrib.pars, B=B), as.list(cx)))
                        #Does not need u as cov is supplied
                if(is.null(x)) {
                        rv$budget$x<-rep(NA, length(u))
                }
                rv$additional <- NULL
                rv$y<-NA
        }
        

        rv$call<-match.call()
        rv$expr<-NA
        
        return(rv)

}

uncert.function<-function(obj, x, u, method=c("NUM", "kragten", "k2", "MC"), cor, cov, 
                        distrib=NULL, distrib.pars=NULL, B=200, delta=0.01, keep.x=TRUE, ...) {

        if(missing(u) && missing(cov)) stop("Either u or cov must be present", call.=TRUE)

        if(missing(cor)) {
                if(missing(cov)) cor <- diag(1, length(u))
                else cor <- cov / outer( sqrt(diag(cov)), sqrt(diag(cov)),"*" )
        } 

        if(any(abs(cor)>1)) stop("cor contains values outside [-1,1]", call.=TRUE)
        
        if(!is.null(names(x)) && !missing(u) ) {
                if(!is.null(names(u)))  u <- u[names(x)]
                        else names(u) <- names(x)
        }
        
        uv<-unlist(u)
        
        #Construct covariance matrix cov if not supplied
        if(missing(cov)) {
                cov<-outer(uv,uv,"*")*cor
        } else {
                if(!missing(u)) warning("Only one of u and cov should be specified: using cov")
                #but we need u for the differentials, so...
                u<-as.list( sqrt( diag(cov) ) )
                names(u) <- names(x)
        }
        
        
        #Obtain differentials ci
        
        method <- match.arg(method, several.ok=TRUE)[1]
        y0<-do.call(obj, c(x,...))
        
        if(method %in% c("NUM", "kragten", "k2") ) {
                
                fxud<-function(x,u,d) x+u*d
                
                if(method == "NUM") {
                        x.plus<-mapply(fxud, x, u, d=delta, SIMPLIFY=FALSE)
                        x.minus<-mapply(fxud, x, u, d=-delta, SIMPLIFY=FALSE)
                        ci <- list()
                        for(i in 1:length(x)) {
                                ci[[names(x)[i]]] <- 
                                        (do.call(obj, c(x[-i],x.plus[i],...)) - do.call(obj, c(x[-i],x.minus[i],...)))/(2*u[[i]]*delta)
                        }
                }

                if(method == "k2") {
                        x.plus<-mapply("+", x, u,SIMPLIFY=FALSE)
                        x.minus<-mapply("-", x, u,SIMPLIFY=FALSE)
                        ci <- list()
                        for(i in 1:length(x)) {
                                ci[[names(x)[i]]] <- 
                                        (do.call(obj, c(x[-i],x.plus[i],...)) - do.call(obj, c(x[-i],x.minus[i],...)))/(2*u[[i]])
                        }
                }

                if(method == "kragten") {
                        x.plus<-mapply(fxud, x, u, d=sign(delta), SIMPLIFY=FALSE)
                        ci <- list()
                        for(i in 1:length(x)) {
                                ci[[names(x)[i]]] <- 
                                        (do.call(obj, c(x[-i],x.plus[i],...)) - y0)/(u[[i]]*sign(delta))
                        }
                }
                
                cc <- unlist(ci)
                v <- ((t(cc) %*% cov ) %*% cc)[1,1]
                cov.xy<-as.vector( cov %*% cc)
                cor.xy<-cov.xy/(uv*sqrt(v))
                names(cov.xy)<-names(cor.xy)<-names(x)
                cov.xy<-as.data.frame(t(cov.xy))
                cor.xy<-as.data.frame(t(cor.xy))
                row.names(cor.xy)<-row.names(cov.xy)<-"theoretical"
                rv <- .construct.uncert(expr=obj, y=y0[1], u.y=sqrt(v), x=x, u=u, ci=ci, 
                        method=method, call.uncert=match.call(),cov=cov, cor=cor, delta=delta,
                        cov.xy=cov.xy, cor.xy=cor.xy, ...)
        
                
        } else if(method=="MC") {
                rv<-uncertMC(obj, x=x, method="MC", 
                        cor=cor, cov=cov, distrib=distrib, distrib.pars=distrib.pars, B=B, keep.x=keep.x,  ...)
                        #Does not need u as cov is supplied
                rv$call <- match.call()
        } else {
                stop(gettextf("method = '%s' is not supported for functions.", 
                        method), domain = NA)
        }
        
        return(rv)

}

uncert.expression<-function(obj, x, u, method=c("GUM", "NUM", "kragten", "k2", "MC"), cor, cov, 
                                distrib=NULL, distrib.pars=NULL, B=200, delta=0.01, keep.x=TRUE, ...) {
        
        if(missing(u) && missing(cov)) stop("Either u or cov must be present", call.=TRUE)

        if(missing(cor)) {
                if(missing(cov)) cor <- diag(1, length(u))
                else cor <- cov / outer( sqrt(diag(cov)), sqrt(diag(cov)),"*" )
        } 

        if(any(abs(cor)>1)) stop("cor contains values outside [-1,1]", call.=TRUE)

        if(!is.null(names(x)) && !missing(u) ) {
                if(!is.null(names(u)))  u <- u[names(x)]
                        else names(u) <- names(x)
        }
        
        obj.names<-all.vars(obj)
        
        if( any( !(obj.names %in% names(c(x,...))) ) | 
                any( !(names(c(x,...)) %in% obj.names) )) {
                stop("Variables in obj do not match arguments in x and '...'", call.=TRUE )
        }

        method <- match.arg(method, several.ok=TRUE)[1]
        if(method %in% c("NUM", "kragten", "k2") ) {
                f.obj<-function(){}
                body(f.obj)<-obj
                formals(f.obj) <- c(x,...)
                rv<-uncert.function(obj=f.obj, x=x, u=u, method=method, 
                                cor=cor, cov=cov, B=B, delta=delta, ...)
                rv$expr<-obj
        } else if(method=="GUM") {
                d.obj<-deriv(obj, obj.names)
                y<-eval(d.obj, c(x,...)) 
                ci<-attr(y, "gradient")[1,names(x)] #Ensures correct order for gradient
                uv<-unlist(u)
                #Construct covariance matrix cov if not supplied
                if(missing(cov)) {
                        cov<-outer(uv,uv,"*")*cor
                }
                v<- ((t(ci) %*% cov) %*% ci)[1,1]
                cov.xy<-as.vector( cov %*% ci)
                cor.xy<-cov.xy/(uv*sqrt(v))
                names(cov.xy)<-names(cor.xy)<-names(x)
                rv <- .construct.uncert( expr=obj, y=y[1], u.y=sqrt(v), x=x, u=u, ci=ci, method="GUM", 
                                cor=cor, cov=cov, cov.xy=cov.xy, cor.xy=cor.xy, deriv=d.obj, ...)
        
        } else if(method=="MC") {
                uv<-unlist(u)
                #Construct covariance matrix cov if not supplied
                if(missing(cov)) {
                        cov<-outer(uv,uv,"*")*cor
                }
                rv<-uncertMC(obj, x=x, method="MC", 
                        cor=cor, cov=cov, distrib=distrib, distrib.pars=distrib.pars, B=B, keep.x=keep.x,  ...)
                        #Does not need u as cov is supplied
                rv$call <- match.call()
        } else {
                stop(gettextf("method = '%s' is not supported for expressions.", 
                        method), domain = NA)
        }
        rv$method=method
        rv$call<-match.call()
        
        return(rv)
}

uncert.formula<-function(obj, x, u, method=c("GUM", "NUM", "kragten", "k2", "MC"), cor, cov, 
                        distrib=NULL, distrib.pars=NULL, B=200, delta=0.01, keep.x=TRUE, ...) {

        if(missing(u) && missing(cov)) stop("Either u or cov must be present", call.=TRUE)

        if(missing(cor)) {
                if(missing(cov)) cor <- diag(1, length(u))
                else cor <- cov / outer( sqrt(diag(cov)), sqrt(diag(cov)),"*" )
        } 

        if(any(abs(cor)>1)) stop("cor contains values outside [-1,1]", call.=TRUE)

        if(!is.null(names(x)) && !missing(u) ) {
                if(!is.null(names(u)))  u <- u[names(x)]
                        else names(u) <- names(x)
        }
        
        obj.names<-all.vars(obj)

        if( any( !(obj.names %in% names(c(x,...))) ) | 
                any( !(names(c(x,...)) %in% obj.names) )) {
                stop("Arguments in obj do not match arguments in x and '...'", call.=TRUE )
        }

        method <- match.arg(method, several.ok=TRUE)[1]

        if(method %in% c("NUM", "kragten", "k2") ) {
                f.obj<-function() {}
                if ((le <- length(obj)) > 1) {
                        body(f.obj)<-obj[[2]]
                } else stop("Invalid formula in uncertainty call")
                formals(f.obj) <- c(x,...)
                rv<-uncert.function(obj=f.obj, x=x, u=u, method=method, 
                                cor=cor, cov=cov, B=B, delta=delta, ...)
        } else if(method=="GUM") {
                d.obj<-deriv(obj, obj.names)
                y<-eval(d.obj, c(x,...)) 
                ci<-attr(y, "gradient")[1,names(x)] #Ensures correct order for gradient
                uv<-unlist(u)
                #Construct covariance matrix cov if not supplied
                if(missing(cov)) {
                        cov<-outer(uv,uv,"*")*cor
                }
                v<- ((t(ci) %*% cov) %*% ci)[1,1]
                cov.xy<-as.vector( cov %*% ci)
                cor.xy<-cov.xy/(uv*sqrt(v))
                names(cov.xy)<-names(cor.xy)<-names(x)
                rv <- .construct.uncert( y=y[1], u.y=sqrt(v), x=x, u=u, ci=ci, method="GUM", 
                                cor=cor, cov=cov, cov.xy=cov.xy, cor.xy=cor.xy, deriv=d.obj, ...)
        
        } else if(method=="MC") {
                uv<-unlist(u)
                #Construct covariance matrix cov if not supplied
                if(missing(cov)) {
                        cov<-outer(uv,uv,"*")*cor
                }
                rv<-uncertMC(obj, x=x, method="MC", 
                        cor=cor, cov=cov, distrib=distrib, distrib.pars=distrib.pars, B=B, keep.x=keep.x, ...)
                        #Does not need u as cov is supplied
                rv$call <- match.call()
        } else {
                stop(gettextf("method = '%s' is not supported for the formula method.", 
                        method), domain = NA)
        }
        rv$expr<-obj
        rv$method<-method
        rv$call<-match.call()
        
        return(rv)
}
