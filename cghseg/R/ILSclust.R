setMethod(f = "ILSclust",signature = "CGHdata",
          definition = function(.Object,CGHo,uniKmax,multiKmax){
            
            P            = CGHo["nblevels"]
            tol          = 1e-2
            select.tmp   = CGHo["select"]
            select(CGHo) = "none"
            options(warn=-1)
            
            command = parse(text = " mu = ILSclust.output(.Object,mu,out.EM$phi,out.EM$tau) \n invisible(list(mu = mu, theta = B,loglik = loglik,nbiter = iter))")    
            
            nbdata   = Reduce("sum",lapply(.Object@Y,FUN = function(y){length(y[!is.na(y)])}) )
            M        = length(names(.Object@Y))
            n.com    = length(.Object@Y[[1]])
            eps      = Inf
            delta    = Inf
            iter     = 0
            
            mu        = multisegmean(.Object,CGHo,uniKmax,multiKmax)$mu
            out.DP2EM = DP2EM(.Object,mu)
            phi       = compactEMinit(out.DP2EM$xk,out.DP2EM$x2k,out.DP2EM$nk,P,CGHo@nbprocs,vh=TRUE)
            out.EM    = compactEMalgo(out.DP2EM$xk,out.DP2EM$x2k,phi,out.DP2EM$nk,P,vh=TRUE)
            mu.test   = ILSclust.output(.Object,mu,out.EM$phi,out.EM$tau)           
            mu.tmp    = mu.test            
            
            while (  (eps > tol) & (iter < CGHo@itermax) & (out.EM$empty==0)){
              iter                = iter+1
              B                   = getbias(.Object,CGHo,mu,B,out.EM$phi,out.EM$tau)
              removebias(.Object) = B$waveffect+B$GCeffect
              mu                  = multisegmixt(.Object,CGHo,uniKmax,multiKmax,out.EM$phi)$mu
              out.DP2EM           = DP2EM(.Object,mu)
	      out.EM              = compactEMalgo(out.DP2EM$xk,out.DP2EM$x2k,phi,out.DP2EM$nk,P,vh=TRUE)
              revertbias(.Object) = B$waveffect+B$GCeffect
              mu.test             = ILSclust.output(.Object,mu,out.EM$phi,out.EM$tau) 
	      eps                 = max(sapply(names(.Object@Y),FUN=function(m,x,y){
                xk = rep(x[[m]]$mean,x[[m]]$end-x[[m]]$begin+1);
                yk = rep(y[[m]]$mean,y[[m]]$end-y[[m]]$begin+1);
                return(max(abs((xk-yk)/xk)))},mu.tmp,mu.test)
                )              
	      mu.tmp              = mu.test              
            }		
            if (out.EM$empty!=0){
              cat("[ILSclust ERROR]: convergence to a solution with empty levels.","\n");
              stop("[ILSclust ERROR]: try a lower nblevels(CGHo)","\n");
            }            
            out.DP2EM    = DP2EM(.Object,mu,theta=Reduce("+",B))
            loglik       = quicklvinc(out.DP2EM$xk,out.DP2EM$x2k,out.EM$phi,out.DP2EM$nk,P,vh=TRUE)$lvinc
            select(CGHo) = select.tmp
            options(warn=0)  
            eval(command)
            
          } #end function
          )

######   auxiliary functions for ILSclust      ########################################

setMethod(f = "ILSclust.output",signature = "CGHdata",
          definition = function(.Object,mu,phi,tau){  
            
            M               = length(names(.Object@Y))
            levels          = apply(tau,1,which.max)
            mutmp           = cbind(mean = phi[levels],levels=levels)
            nk              = lapply(mu,FUN=function(x){length(x$mean)})  
            end             = cumsum(nk)
            begin           = c(1,end[1:(length(end)-1)]+1)
            rupt            = cbind(begin,end)
            rownames(rupt)  = names(.Object@Y)
            mutmp           = apply(rupt,1,FUN = function(y){
              tmp           = data.frame(matrix(mutmp[y[1]:y[2],],ncol=2))
              colnames(tmp) = c("mean","levels")
              invisible(tmp)
            })
            mutmp = lapply(1:M,FUN = function(m){cbind(mu[[m]][,-3],mutmp[[m]])})
            names(mutmp) = names(.Object@Y)  
            invisible(mutmp)
          }
          )
