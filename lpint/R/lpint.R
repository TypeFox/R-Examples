lpint <-
  function(jmptimes,jmpsizes=rep(1,length(jmptimes)),
           Y=rep(1,length(jmptimes)), bw=NULL, adjust = 1,
           Tau=max(1.0,jmptimes), p=nu+1, nu=0,
           K = function(x)3/4*(1-x^2)*(x <= 1 & x >= -1),
           n = 101,bw.only=FALSE){
  jmptimes <- jmptimes/Tau;
  ts <- seq(from=0,to=1,length=n);

  Ast.mid <- {
    v <- sapply(0 : (2*p),
                function(i){
                  integrate(function(u)K(u)* u^i, lower=-1.0,
                            upper=1.0)$value
                })
    sapply(0:p,function(i)v[(i+1):(i+p+1)])
  }
  ## browser()
  ek.mid <- function(x){
    ek.s <- function(u){
      solve(Ast.mid,(-u)^(0:p))[nu+1] * K(-u)* gamma(nu+1)
    }
    sapply(x,ek.s)
  };
  if(!is.numeric(bw)){
    const.Knup <- 
      (integrate(function(x)ek.mid(x)^2, -1, 1)$value/
       (integrate(function(x)ek.mid(x)*x^(p+1), -1, 1)$value^2))^(1/(2*p + 3))*
         (gamma(p + 2)^2 * (2*nu+1) /(2*(p+1-nu)) )^(1/(2*p + 3))
    beta.tilde <-
      outer(1:(p+3+1),1:(p+3+1), ## inverse of a Hilbert matrix
            function(i,j)(-1)^(i+j)* (i+j-1)*
            choose(p+3+1 + i - 1,p+3+1 - j) *
            choose(p+3+1 + j - 1,p+3+1 - i) *
            choose(i+j-2,i-1)^2
            )%*%
                t(sapply(0:(p+3),function(i)jmptimes[Y>0]^i))%*%
                    (jmpsizes[Y>0]/Y[Y>0])
    
    ## estimated integrated squared (p+1)st derivative:
    ISD <- sum(outer(p+ 1:3,p+ 1:3,
                     function(i,j){
                       gamma(i+1)*beta.tilde[i+1]/gamma(i-p) *
                         gamma(j+1)*beta.tilde[j+1]/gamma(j-p)/
                           (i+j-2*(p+1)+1)
                     })
               )
    bw <- const.Knup* (sum(jmpsizes[Y>0]/Y[Y>0]^2)/ISD)^(1/(2*p+3))
    if(bw.only)return(bw*Tau)
  }else{bw <- bw/Tau}
  if(length(bw)>1){
    bw <- bw[1]; # relative bandwidth
    warning(paste("length bw=",bw," is larger than 1!\n"))
  }
  bw <- bw*adjust;
  lpint <- function(i){
    ## browser()
    if(ts[i] >= bw && ts[i] <= 1-bw){
      ek <- ek.mid
    }else{##if(i==12)browser()
      Ast <- {
        v <- sapply(0 : (2*p),
                    function(j){
                      integrate(function(u)K(u)* u^j,
                                lower=-min(1.0,ts[i]/bw),
                                upper=min(1.0,(1-ts[i])/bw)
                                )$value
                    })
        sapply(0:p,function(j)v[(j+1):(j+p+1)])
      }
      ##  if(i==12)print(Ast)
      ek <- function(x){
        ek.s <- function(u){
          solve(Ast, (-u)^(0:p))[nu+1] * K(-u) * gamma(nu+1)
        }
        sapply(x,ek.s)
      };
    }
    tmp <- ek((ts[i]-jmptimes[Y>0])/bw)
    c(tmp%*% (jmpsizes[Y>0]/Y[Y>0]) /
      bw^(nu+1) /Tau^(nu+1),
      sqrt((tmp^2)%*% (jmpsizes[Y>0]/Y[Y>0]^2)) /
      bw^(nu+1) /Tau^(nu+1)
      )
  }
  tmp <- sapply(1:n,lpint)
  list(x=ts*Tau,y=tmp[1,],se=tmp[2,],bw=bw*Tau)
}

