blockmatrix <-
function(A){
if(any(length(dim(A))!=2,dim(A)[1]!=dim(A)[2])) stop("A must be a square matrix")
order<-dim(A)[1]
if(is.matrix_irreducible(A)) stop("Cannot compute block-permutation for irreducible matrix")
else{
    stage.order<-1:order
    powermatrix<-(diag(1,order)+A)%^%(order-1)
    zeroes1<-numeric(order)
    for(i in 1:order){
        zeroes1[i]<-length(which(powermatrix[i,]==0))
    }
    stage.order1<-order(zeroes1)
    blockpowermatrix<-matrix(0,order,order)
    for(i in 1:order){
        for(j in 1:order){
            blockpowermatrix[i,j]<-powermatrix[stage.order1[i],stage.order1[j]]
        }
    }
    zeroes2<-numeric(order)
    for(i in 1:order){
        if(blockpowermatrix[i,1]==0){
            while(blockpowermatrix[i,1+zeroes2[i]]==0){
                zeroes2[i]<-zeroes2[i]+1
            }
        }
    }
    stage.order2<-order(zeroes2)
    stage.order<-stage.order[stage.order1][stage.order2]
    blockmatrix<-matrix(0,order,order)
    for(i in 1:order){
        for(j in 1:order){
            blockmatrix[i,j]<-A[stage.order[i],stage.order[j]]
        }
    }
    return(list(blockmatrix=blockmatrix,order=stage.order))
}}

