setMethod(f = "getuniKmax",signature = "CGHdata",
          definition = function(.Object,CGHo,uniKmax=NULL){
            
            if (is.null(uniKmax)){              
              if (CGHo["select"] == "none"){
                cat("[getuniKmax] if no selection is performed while uniKmax is not specified \n")                
              }
              cat("[getuniKmax] uniKmax initialized by pre-screening\n")
              ## uniKmax = lapply(names(.Object@Y),FUN = function(ell){floor(sum(!is.na(.Object@Y[[ell]]))*CGHo["alpha"])})
              select.tmp    = CGHo["select"]
              calling.tmp   = CGHo["calling"]
              calling(CGHo) = FALSE
			  select(CGHo)  = "mBIC"
			  if (CGHo@nbprocs>1){
				  uniKmax = mclapply(.Object@Y, FUN = function(y){
							  Kmax = floor(sum(!is.na(y))*CGHo["alpha"])
							  Kmax = min(200,Kmax)
							  dim(unisegmean(y,CGHo,Kmax)$mu)[1]
#						  	}) # fun argument instead of FUN
						  }, mc.cores = CGHo@nbprocs)
			  }			  
			  else{
				  uniKmax = lapply(.Object@Y,FUN = function(y){
							  Kmax = floor(sum(!is.na(y))*CGHo["alpha"])
							  Kmax = min(200,Kmax)
							  dim(unisegmean(y,CGHo,Kmax)$mu)[1]
						  })
			  }
			  Kmax    = max(unlist(uniKmax))
              uniKmax = lapply(uniKmax,FUN=function(x){x=2*Kmax})
              names(uniKmax) = names(.Object@Y)
              select(CGHo)   = select.tmp
            } else {              
              if ( sum(names(.Object@Y) != names(uniKmax))>0){                
                cat("[getuniKmax] uniKmax should be a list with the same names as CGHd \n")
                cat("[getuniKmax] the names of uniKmax should be patients names \n")      
                stop()
              }              
              lapply(names(.Object@Y), FUN = function(m){
                if (CGHo["calling"]){
                  if ((CGHo["nblevels"]>uniKmax[[m]]) ){
                    cat("[getuniKmax] Error in uniKmax \n")
                    cat("[getuniKmax] Error for profile: ", m, "\n")
                    cat("[getuniKmax] The number of clusters must be lower than the number of segments in this profile","\n")
                    cat("[getuniKmax] Check CGHo[\"nblevels\"] and Kmax for this profile","\n")
                    stop()
                  }
                }
                n.com  = length(.Object@Y[[1]])                
              })
            }
            return(uniKmax)  
          })

















