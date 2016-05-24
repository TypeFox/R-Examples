get.cnrind <- function(x,y,rw) {
x.crnrs <- rw[c(1,2,2,1)]
y.crnrs <- rw[c(3,3,4,4)]
M1 <- outer(x,x.crnrs,function(a,b){(a-b)^2})
M2 <- outer(y,y.crnrs,function(a,b){(a-b)^2})
MM <- M1 + M2
apply(MM,2,which.min)
}
