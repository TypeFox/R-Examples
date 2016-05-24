`tri.ineq` <-
function(dist) 
{
    mat<-as.matrix(dist)
    n<-dim(mat)[1]
    ineq <- 0
    for (i in 1:(n-2)) {
        for (j in (i+1):(n-1)) {
            for (k in (j+1):n) {
                sds<-c(mat[j,i], mat[k,i], mat[k,j])
                lng<-max(sds)
                if (lng>(sum(sds)-lng)) ineq<-ineq+1
            }
        }
    }
    return(ineq==0)
}
