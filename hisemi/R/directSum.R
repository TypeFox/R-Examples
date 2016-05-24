directSum=
function(...){
        mats=list(...)
        for(i in 1:length(mats)){
                if(inherits(mats[[i]],'matrix')) next
                stopifnot(mode(mats[[i]])=='numeric') ## TODO: Sparse matrices
                tmp=as.matrix(mats[[i]])
                if(any(dim(tmp)==0)) dim(tmp)=c(0,0)
                mats[[i]]=tmp
        }
        dims=as.matrix(sapply(mats,dim))
        ans=matrix(0,sum(dims[1,]),sum(dims[2,]))
        ro=0; co=0;
        for(i in 1:length(mats)){
                if(dims[1,i]>0) ans[ro+1:dims[1,i],co+1:dims[2,i]]=mats[[i]]
                ro=ro+dims[1,i]
                co=co+dims[2,i]
        }
        ans
    }

