setMethod(f = "getbias",signature = "CGHdata",
          definition = function(.Object,CGHo,mu,bias,phi=c(),tau=c()){

            P        = length(phi)/3
            nplust   = apply(Reduce("rbind",.Object@Y),2,FUN  = function(x){sum(!is.na(x))})
            M        = length(names(.Object@Y))
            GCeffect = 0
            
            if (CGHo["calling"] == FALSE){

              if (CGHo["GCnorm"] =="linear"){
                z = mapply(names(.Object@Y), FUN = function(m){
                  nk  = mu[[m]]$end-mu[[m]]$begin+1
                  muk = mu[[m]]$mean
                  return(.Object@Y[[m]]-rep(muk,nk)-bias$waveffect)
                })  
                y = apply(z,1,FUN=function(x){sum(x,na.rm=TRUE)})/nplust                
                GCeffect = lm(y ~ poly(.Object@GCcontent,2))$fitted.values
              }
              
              z = mapply(names(.Object@Y), FUN = function(m){
                nk  = mu[[m]]$end-mu[[m]]$begin+1
                muk = mu[[m]]$mean
                return(.Object@Y[[m]]-rep(muk,nk)-GCeffect)
              })    
              
            }  else if (CGHo["calling"] ==TRUE){
              
              n.com  = length(.Object@Y[[1]])
              nk     = lapply(mu,FUN=function(x){length(x$mean)})
              end    = cumsum(nk)
              begin  = c(1,end[1:(length(end)-1)]+1)
              rupt   = cbind(begin,end)
              row.names(rupt) = names(nk)
              taulist         = apply(rupt,1,FUN = function(y){
                tmp           = data.frame(matrix(tau[y[1]:y[2],],ncol=P))
                colnames(tmp) = paste(rep("tau",P),c(1:P),sep="")
                return(tmp)
              })
              A        = lapply(1:M,FUN = function(y){cbind(mu[[y]],taulist[[y]])})
              names(A) = names(mu)    
              
              if (CGHo["GCnorm"] =="linear"){               
                z= mapply(1:P,FUN = function(p){
                  D= mapply(names(.Object@Y),FUN=function(m){
                    nk   = A[[m]]$end-A[[m]]$begin+1
                    taup =  A[[m]][,p+3]
                    u = (.Object@Y[[m]]-phi[p]-bias$waveffect)*rep(taup,nk)
                    return(u)
                  })
                  return(apply(D,1,FUN=function(x){sum(x,na.rm=TRUE)}))
                })
                y = apply(z,1,FUN=function(x){sum(x,na.rm=TRUE)})/nplust                            
                GCeffect = lm(y ~ poly(.Object@GCcontent,2))$fitted.values
              }
              
              z= mapply(1:P,FUN = function(p){
                D= mapply(names(.Object@Y),FUN=function(m){
                  nk   = A[[m]]$end-A[[m]]$begin+1
                  taup =  A[[m]][,p+3]
                  u = (.Object@Y[[m]]-phi[p]-GCeffect)*rep(taup,nk)
                  return(u)
                })
                return(apply(D,1,FUN = function(x){sum(x,na.rm=TRUE)}))
              })
            }

            waveffect = apply(z,1,FUN=function(x){sum(x,na.rm=TRUE)})/nplust

            if (CGHo["wavenorm"]=="spline"){
              waveffect = smooth.spline(c(1:length(waveffect)),waveffect)$y
            }
                       
            waveffect = waveffect - mean(waveffect,na.rm=TRUE)
            invisible(list(waveffect=waveffect,GCeffect=GCeffect))
          }
          )

