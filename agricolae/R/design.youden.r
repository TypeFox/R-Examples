`design.youden` <-
function (trt,r,serie=2,seed=0,kinds="Super-Duper",first=TRUE,randomization=TRUE)
{
number<-10
if(serie>0) number<-10^serie
t <- length(trt) # number of treatment
if (seed == 0) {
genera<-runif(1)
seed <-.Random.seed[3]
}
set.seed(seed,kinds)
parameters<-list(design="youden",trt=trt,r=r,serie=serie,seed=seed,kinds=kinds)
a <- matrix(1,nrow=t,ncol=r)
# Generate matrix random
for (j in 1:r){
i<-c(j:t)
if(j>1) i<-c(i,1:(j-1))
a[,j]<-i
}
m<-1:t
if(randomization)m<-sample(1:t,t)
M<-m[a]
a<-trt[M]
dim(a)<-c(t,r)
# Randomize row and column
m<-1:r
if(randomization)m<-sample(1:r,r)
a[,m]-> a
m<-1:t
if(randomization)m<-c(1,sample(2:t,t-1))
a<-a[m,]
if (!first) {
	m<-order(a[1,])
	a<-a[,m]
}
a<-t(a)
trat<-as.vector(a)
columna <- rep(gl(r, 1), t)
fila <- gl(t, r)
fila <- as.character(fila)
fila <- as.numeric(fila)
plots <- fila*number+(1:r)
book <- data.frame(plots, row = as.factor(fila), col = as.factor(columna),
		trat = as.factor(trat))
names(book)[4] <- c(paste(deparse(substitute(trt))))
outdesign<-list(parameters=parameters,sketch=matrix(book[,4], byrow = TRUE, ncol = r),book=book)
return(outdesign)
}
