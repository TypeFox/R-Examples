setMethod(f = "DP2EM",signature = "CGHdata",
          definition = function(.Object,mu,theta=NULL){
            if (is.null(theta)){
              xk  = unlist(lapply(names(.Object@Y),FUN=function(m){meanRuptR_c(.Object@Y[[m]], mu[[m]][,2], length(mu[[m]][,2]))}))
              x2k = unlist(lapply(names(.Object@Y),FUN=function(m){meansqRuptR_c(.Object@Y[[m]], mu[[m]][,2], length(mu[[m]][,2]))}))
            } else {
              xk  = unlist(lapply(names(.Object@Y),FUN=function(m){meanRuptR_c(.Object@Y[[m]]-theta, mu[[m]][,2], length(mu[[m]][,2]))}))
              x2k = unlist(lapply(names(.Object@Y),FUN=function(m){meansqRuptR_c(.Object@Y[[m]]-theta, mu[[m]][,2], length(mu[[m]][,2]))}))
            }
            nk  = unlist(lapply(mu,FUN=function(x){x$end-x$begin+1}))
            invisible(list(xk=xk,x2k=x2k,nk=nk))
          }
          )


