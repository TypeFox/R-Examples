`bisearch` <-
function(f, interval, lower = min(interval), upper = max(interval),
                     signfu=sign(f(upper)), signfl=sign(f(lower)),
                     tol = .Machine$double.eps^0.25, maxiter = 1000){
            if(signfl*signfu>=0){ stop("f() values at end points not of opposite sign") } else {
               for(i in 1:maxiter){
                    if(f((lower+upper)/2)*signfu>0) upper<-(lower+upper)/2 else lower<-(lower+upper)/2;
                    if(abs(lower-upper)<=tol) break;
                  }
                  if(i==maxiter) warning("maximum number of iterations reached, precision may be insufficient");
                  root <- (lower+upper)/2;
                  list(root=root,f.root=f(root),iter=i,estim.prec=abs(lower-upper))
                 }
            }

