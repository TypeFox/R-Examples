setMethod(f = "multisegout",signature = "CGHdata",
          definition = function(.Object,seg.rep,Res,Kselect){
            
            M  = length(names(.Object@Y))
            
            out = lapply(names(.Object@Y),FUN = function(m,y){
              i        = which(row.names(seg.rep) == m)
              k        = seg.rep[i,Kselect-M+1]
              rupt     = matrix(Inf,ncol = 2 , nrow= k)
              rupt[,2] = Res[[i]]$t.est[k,1:k]    
              if (k==1){
                rupt[1,1] = 1
              } else {
                rupt[,1] = c(1,rupt[1:(k-1),2]+1)
              }			  
	      resmean = meanRuptR_c(y[[m]], rupt[,2], k)
              mu      = data.frame(begin = rupt[,1],
              end     = rupt[,2],
              mean    = resmean)
              invisible(mu)
            },.Object@Y)
            names(out) = names(.Object@Y)
            invisible(out)
          })
