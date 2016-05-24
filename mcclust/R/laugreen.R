`laugreen` <-
function(psm, start.cl=NULL, tol= 0.001){
  
            one.constr.val <- function(j,k,Sim){
                if(Sim[j,k]==1) return(c(0,0))
                1    
            }
            one.constr.mat <- function(j,k,Sim){
                if(Sim[j,k]==1){
                    constr.mat <- matrix(0,ncol=ncol(Sim),nrow=2)
                    constr.mat[1:2,c(j,k)] <- c(1,-1,-1,1)
                    return(t(constr.mat))
                } else{constr.mat <- rep(0,ncol(Sim))
                        constr.mat[c(j,k)] <- c(1,1)
                        return(constr.mat)    
                } 
            }

    n <- nrow(psm)
    if(is.null(start.cl)) start.cl <- 1:n
    Sim <- cltoSim(start.cl)
    obj <- sum((Sim*(1-2*psm))[lower.tri(Sim)])
    obj.new <- obj - 2*tol
    iter <- 0
    while(obj.new < obj-tol){
        iter <- iter +1
        obj <- obj.new
        for(i in 1:n){
            f.i <- 1-2*psm[i,]
            all.jk <- combn((1:n)[-i],2)
            constr.val <- unlist(apply(all.jk,2,function(x) one.constr.val(x[1],x[2],Sim=Sim)))
            constr.mat <- matrix(unlist(apply(all.jk,2,function(x) one.constr.mat(x[1],x[2],Sim=Sim))),ncol=n,byrow=TRUE)
            constr.dir <- rep("<=",length(constr.val))
            res.lp <- lp("min", f.i[-i], constr.mat[,-i], constr.dir, constr.val, all.bin=TRUE)
            Sim[i,-i] <- Sim[-i,i] <- res.lp$solution
        }
    obj.new <- sum((Sim*(1-2*psm))[lower.tri(Sim)])
    }
list(cl=Simtocl(Sim), value=sum(abs((Sim-psm)[lower.tri(psm)])), method="laugreen", iter.lg=iter)   
}

