`design.lsd` <-
function (trt,serie=2,seed=0,kinds="Super-Duper",first=TRUE,randomization=TRUE)
{
number<-10
if(serie>0) number<-10^serie
r <- length(trt)
if (seed == 0) {
genera<-runif(1)
seed <-.Random.seed[3]
}
set.seed(seed,kinds)
parameters<-list(design="lsd",trt=trt,r=r,serie=serie,seed=seed,kinds=kinds,randomization)
a <- 1:(r * r)
dim(a) <- c(r, r)
for (i in 1:r) {
for (j in 1:r) {
k <- i + j - 1
if (k > r)
k <- i + j - r - 1
a[i, j] <- k
}
}
m<-2:r
if(randomization)m<-sample(2:r,r-1)
a<-a[,c(1,m)]
if(randomization){
if (first) {
	m<-sample(1:r,r)
	a<-a[m,]
}}
trat<-trt[a]
columna <- rep(gl(r, 1), r)
fila <- gl(r, r)
fila <- as.character(fila)
fila <- as.numeric(fila)
plots <- fila*number+(1:r)
book <- data.frame(plots, row = as.factor(fila), col = as.factor(columna),
		trat = as.factor(trat))
names(book)[4] <- c(paste(deparse(substitute(trt))))
outdesign<-list(parameters=parameters,sketch=matrix(book[,4], byrow = TRUE, ncol = r),book=book)
return(outdesign)
}

