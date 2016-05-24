setMethod(f = "multisegmean",signature = "CGHdata",
          definition = function(.Object,CGHo,uniKmax,multiKmax){
            
            select.tmp   = CGHo["select"]
            select(CGHo) = "none"
            multiKselect = multiKmax
            M            = length(names(.Object@Y))
            Kseq         = c(M:multiKmax)
            
######   Individual segmentations for patients 
            if (CGHo@nbprocs>1){           
				#cat("multisegmean //               \r")	
					Res = mclapply(names(.Object@Y), FUN = function(m){
								n                           = length(which(!is.na(.Object@Y[[m]])))
								Kmax                        = uniKmax[[m]]
								out                         = unisegmean(.Object@Y[[m]],CGHo,Kmax)
								J.est                       = n*exp(-((2/n)*out$loglik+log(2*pi)+1))
								invisible(list(t.est = out$t.est, loglik = out$loglik,J.est=J.est))
							}, mc.cores = CGHo@nbprocs)
					names(Res) = names(.Object@Y)							
            } 
			else {
              Res = lapply(names(.Object@Y), FUN = function(m){
                n                           = length(which(!is.na(.Object@Y[[m]])))
                Kmax                        = uniKmax[[m]]
                out                         = unisegmean(.Object@Y[[m]],CGHo,Kmax)
                J.est                       = n*exp(-((2/n)*out$loglik+log(2*pi)+1))
                invisible(list(t.est = out$t.est, loglik = out$loglik,J.est=J.est))
              }) 
              names(Res) = names(.Object@Y)
            } 
######   Segment Repartition segments across patients 
            
            J.est              = lapply(Res,FUN = function(x){x$J.est})
            nbdata             = sum(Reduce("c",lapply(.Object@Y,FUN = function(y){length(y[!is.na(y)])})))
            out.ibp            = segibp(data.frame(J.est = unlist(J.est)),unlist(uniKmax),multiKmax)
            multiloglik        = -(nbdata/2)*(log(2*pi*out.ibp[[1]]/nbdata)+1)
            seg.rep            = out.ibp[[2]]
            row.names(seg.rep) = names(.Object@Y)
            
            
######   Model Selection  
            
            if (select.tmp=="none"){
              multiKselect = multiKmax
              dimll        = length(multiloglik)
            } 
			else if (select.tmp=="mBIC"){
				if (CGHo@nbprocs>1){	
                #cat("multisegmean // part 2                  \r")
					mBIC = mclapply(Kseq,FUN=function(K){
								mu      = multisegout(.Object,seg.rep,Res,K)
								getmBIC(K,multiloglik[K-M+1],mu,CGHo)     
							}, mc.cores = CGHo@nbprocs)
				}
				else {
                	mBIC = sapply(Kseq,FUN=function(K){
                  		mu      = multisegout(.Object,seg.rep,Res,K)
                  		getmBIC(K,multiloglik[K-M+1],mu,CGHo)     
                		})
              	}  			
              	multiKselect = Kseq[which.max(mBIC)]
              	dimll        = multiKselect
            }
            
######   Outputs   
            
            mu           = multisegout(.Object,seg.rep,Res,multiKselect)
            select(CGHo) = select.tmp  
#			if (CGHo@nbprocs>1){	
#				cat("multisegmean finished                  \r") 
#			}
            invisible(list(mu=mu,loglik=multiloglik[dimll],nbiter=0))
          } 
          )
