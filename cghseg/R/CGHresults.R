setClassUnion("listOrNULL", c("list","NULL"))
setClassUnion("numOrNULL", c("numeric","NULL"))
setClassUnion("factorOrNULL", c("factor","NULL"))

setClass("CGHresults",
         representation=list(
           loglik           = 'list',
           mu               = 'listOrNULL',
           theta            = 'listOrNULL',
           nbiter           = 'numeric',
           from             = 'character',
           genomic.position = "numOrNULL",
           probeID          = "factorOrNULL",
           options          = 'CGHoptions'
           )
         )

setMethod("initialize",
          "CGHresults",
          function(.Object,CGHd,CGHo){
            
            .Object@mu      = NULL
            .Object@theta   = NULL
            .Object@loglik  = list()
            .Object@nbiter  = 0
            .Object@from    = ""
            .Object@options = CGHo
            
            if (CGHo@select !="none"){
              .Object@mu = list()
            }
            
            if (length(names(CGHd@Y)>1)){
              .Object@theta = list(waveffect=0,GCeffect=0)
            }         
            
            .Object@options          = CGHo
            .Object@genomic.position = CGHd@genomic.position
            .Object@probeID          = CGHd@probeID
            .Object@loglik           = list()            
            .Object
          })


setMethod(
          f = "[",
          signature  = "CGHresults",
          definition = function(x,i,j,drop){   
            if (i=="mu")               {return(x@mu)}               else {}
            if (i=="theta")            {return(x@theta)}            else {}
            if (i=="nbiter")           {return(x@nbiter)}           else {}
            if (i=="loglik")           {return(x@loglik)}           else {}
            if (i=="from")             {return(x@from)}             else {}
            if (i=="genomic.position") {return(x@genomic.position)} else {}
            if (i=="probeID")          {return(x@probeID)}          else {}
            if (i=="options")          {return(x@options)}          else {}
          }          
          )



setGeneric("mu<-",function(object,value){standardGeneric("mu<-")})
setReplaceMethod(
                 f="mu",
                 signature="CGHresults",
                 definition=function(object,value){
                   object@mu =value
                   return (object)
                 }
                 )

setGeneric("loglik<-",function(object,value){standardGeneric("loglik<-")})
setReplaceMethod(
                 f="loglik",
                 signature="CGHresults",
                 definition=function(object,value){
                   object@loglik <- value
                   return (object)
                 }
                 )
setGeneric("theta<-",function(object,value){standardGeneric("theta<-")})
setReplaceMethod(
                 f="theta",
                 signature="CGHresults",
                 definition=function(object,value){
                   object@theta =value
                   return (object)
                 }
                 )
setGeneric("nbiter<-",function(object,value){standardGeneric("nbiter<-")})
setReplaceMethod(
                 f="nbiter",
                 signature="CGHresults",
                 definition=function(object,value){
                   object@nbiter =value
                   return (object)
                 }
                 )
setGeneric("from<-",function(object,value){standardGeneric("from<-")})
setReplaceMethod(
                 f="from",
                 signature="CGHresults",
                 definition=function(object,value){
                   object@from =value
                   return (object)
                 }
                 )

setGeneric( name = "getbp"          ,def = function(.Object){standardGeneric("getbp")})
setGeneric( name = "getsegprofiles" ,def = function(.Object){standardGeneric("getsegprofiles")})
setGeneric( name = "getlevels"      ,def = function(.Object){standardGeneric("getlevels")})
setGeneric( name = "getbgoutliers"  ,def = function(.Object,CGHr.smooth,fdr){standardGeneric("getbgoutliers")})

setMethod(
          f = "summary",
          signature = "CGHresults",
          definition = function(object){
            cat("****** Summary of CGHr object ******\n")
            cat("\n")
            cat("[CGHr summary] Patients IDs\n");
            cat(names(object@mu),fill=TRUE,"\n")
            Ktot = sum(unlist(lapply(object@mu,FUN = function(x){dim(x)[1]})))
            cat("[CGHr summary] Total Number of segments: ",Ktot,"\n")
            cat("[CGHr summary] Individual Number of segments:\n")
            print(sapply(object@mu,FUN = function(x){dim(x)[1]}))
            
          }          
          )


setMethod(f = "getbp",signature = "CGHresults",
          definition = function(.Object){

            if (is.null(.Object@genomic.position)){
              Km = dim(.Object@mu[[1]])[1]              
              x  = c(1:.Object@mu[[1]]$end[Km])
            } else {
              x = .Object@genomic.position
            }
            
            bp = lapply(names(.Object@mu),FUN = function(ell){
              cat(ell,"\n")
              bpi = tabulate(.Object@mu[[ell]]$end)
              if (!is.null(.Object@probeID)){
                out = data.frame(probeID = .Object@probeID,position = x,bp=bpi)
              } else {
                out = data.frame(position = x,bp=bpi)
              }
              invisible(out)
            })
            names(bp) = names(.Object@mu)                        
            return(bp)
          } 
          )

setMethod(f = "getsegprofiles",signature = "CGHresults",
          definition = function(.Object){
            
            if (is.null(.Object@genomic.position)){            
              Km      = dim(.Object@mu[[1]])[1]              
              x       = c(1:.Object@mu[[1]]$end[Km])
              predict = sapply(names(.Object@mu), FUN = function(ell){
                nk    = .Object@mu[[ell]]$end -  .Object@mu[[ell]]$begin + 1
                predict = rep(.Object@mu[[ell]]$mean,nk)
                invisible(predict)                    
              })
            } else {
              x       = .Object@genomic.position
              predict = sapply(names(.Object@mu), FUN = function(ell){
                nk = match(.Object@mu[[ell]]$end,x) -  match(.Object@mu[[ell]]$begin,x) + 1
                predict = rep(.Object@mu[[ell]]$mean,nk)
                invisible(predict)                    
              })
            }
            names(predict) = names(.Object@mu)                         
            return(predict)
          })





setMethod(f = "getlevels",signature = "CGHresults",
          definition = function(.Object){
            
            if (.Object@from == "uniseg"){
              cat("[warnings] calling has been performed with uniseg \n")
              cat("[warnings] calling IDs do not correspond and have a meaning for each profile only (not globally)\n",fill=TRUE)
            } else {
              
              if (.Object@options@calling){                
                if (is.null(.Object@genomic.position)){
                  Km    = dim(.Object@mu[[1]])[1]              
                  x     = c(1:.Object@mu[[1]]$end[Km])
                  Z = lapply(.Object@mu, FUN = function(y){
                    nk    = y$end -  y$begin + 1
                    levels = rep(y$levels,nk)
                    if (!is.null(.Object@probeID)){
                      return(data.frame(probeID = .Object@probeID, position=x,data.frame(levels)))
                    } else {
                       return(data.frame(position=x,data.frame(levels)))
                    }
                  })                            
                } else {
                  x = .Object@genomic.position
                  Z = lapply(.Object@mu, FUN = function(y){
                    nk = match(y$end,x) -  match(y$begin,x) + 1 
                    levels = rep(y$levels,nk)
                    if (!is.null(.Object@probeID)){
                      return(data.frame(probeID = .Object@probeID, position=x,data.frame(levels)))
                    } else {
                      return(data.frame(position=x,data.frame(levels)))
                    }
                  })        
                }
                return(Z)
              } else {
                cat("[output warning] getlevels only works when CGHo[\"calling\"]==TRUE \n")
              }              
            }              
          })




setMethod(f = "getbgoutliers",signature = "CGHresults",
          definition = function(.Object,CGHr.smooth,fdr){

            if ( (.Object@options@wavenorm!="position") | (CGHr.smooth@options@wavenorm =="position")){
              cat("[warning] the entry must be CGHr.position,CGHr.spline \n")
              stop()
            } else {
              n       = length(.Object@theta$waveffect)
              bgnoise = .Object@theta$waveffect - CGHr.smooth["theta"]$waveffect
              tau2    = sum(bgnoise^2)/(length(bgnoise)-1)
              Pv      = 2*(1-pnorm(abs(bgnoise),sd=sqrt(tau2)))
              Padj    = p.adjust(Pv,"BH")
              pos     = which(Padj<fdr)
              return(tabulate(pos,n))
            }
          })

