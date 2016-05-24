setMethod(f = "ILS",signature = "CGHdata",
          definition = function(.Object,CGHo,uniKmax,multiKmax){

            tol          = 1e-2
            select.tmp   = CGHo["select"]
            select(CGHo) = "none"
            command      = parse(text = "invisible(list(mu = mu, theta = B,loglik = loglik,nbiter = iter))")
            
            nbdata   = Reduce("sum",lapply(.Object@Y,FUN = function(y){length(y[!is.na(y)])}) )
            M        = length(names(.Object@Y))
            n.com    = length(.Object@Y[[1]])
            eps      = Inf
            iter     = 0
            
            mu       = multisegmean(.Object,CGHo,uniKmax,multiKmax)$mu
            B        = list(waveffect = rep(0,n.com), GCeffect = rep(0,n.com))
			mu.tmp   = mu 
            
            while ( (eps > tol) & (iter < CGHo@itermax)){
              iter                = iter+1
              B                   = getbias(.Object,CGHo,mu,B)		
              removebias(.Object) = B$waveffect+B$GCeffect
              mu                  = multisegmean(.Object,CGHo,uniKmax,multiKmax)$mu
              revertbias(.Object) = B$waveffect+B$GCeffect
              eps    = max(sapply(names(.Object@Y),FUN=function(m,x,y){
                xk = rep(x[[m]]$mean,x[[m]]$end-x[[m]]$begin+1);
                yk =rep(y[[m]]$mean,y[[m]]$end-y[[m]]$begin+1) ;
                return(max(abs((xk-yk)/xk)))},mu.tmp,mu))
				mu.tmp = mu
            } # end while
            
            out.DP2EM    = DP2EM(.Object,mu,theta=Reduce("+",B))
            RSS          = 0.5*sum(out.DP2EM$nk*(out.DP2EM$x2k-(out.DP2EM$xk)^2))
            n            = Reduce("sum",lapply(.Object@Y,FUN = function(y){length(y[!is.na(y)])}))
            loglik       = -(n/2)*(log(2*pi*RSS/n)+1)
            select(CGHo) = select.tmp
            eval(command)
            
          })

