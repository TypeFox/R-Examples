mergeCols <-
function(n,minimum=2){ 
    for (j in length(n[1,]):1){
        if(is.null(dim(n))) return(0)
        if(sum(n[,j])<=minimum) {
            if(j>1) {
                n[,j-1] <- n[,j-1]+n[,j]
                n <- n[,-j]
            }else{
                if(dim(n)[2]>1){
                    n[,2] <- n[,2]+n[,1]
                    n <- n[,-1]
                }else return(0)
            }
        }
    }
    return(n)                
}
