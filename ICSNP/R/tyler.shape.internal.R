# Iteration step for Tyler's shape matrix
    
.tyler.step<-function(V.old,datas,p,n)
        {
        sqrt.V.old<-mat.sqrt(V.old)
        r<-sqrt(rowSums((datas %*% sqrt.V.old)^2))
        M.V.old<-p/n*(t(((1/r)*datas%*%sqrt.V.old))%*%((1/r)*datas%*%sqrt.V.old))
        M.V.old.inv <- solve(M.V.old)
        V.new<-sum(diag(V.old %*% M.V.old.inv))^(-1)*(sqrt.V.old %*% M.V.old.inv %*% sqrt.V.old)
        return(V.new)
        }
