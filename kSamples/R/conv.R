conv <-
function(x1,p1,x2,p2){
#----------------------------------------------------
# R routine for calling conv in conv.c
# created 02.05.2012  Fritz Scholz
#----------------------------------------------------
n1 <- length(x1)
n2 <-length(x2)
n <- n1*n2
x <- numeric(n)
p <- numeric(n)

 out <- .C("conv",  x1=as.double(x1), p1=as.double(p1),
			n1=as.integer(n1), x2=as.double(x2), p2=as.double(p2),
			n2=as.integer(n2), x=as.double(x), p=as.double(p),
			n=as.integer(n))
n <- out$n
cbind(out$x[1:n],out$p[1:n])
}
