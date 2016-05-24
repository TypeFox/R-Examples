design.strip<-function (trt1, trt2,r, serie = 2, seed = 0, kinds = "Super-Duper",randomization=TRUE)
{
number<-10
if(serie>0) number<-10^serie
    n1<-length(trt1)
    n2<-length(trt2)
    if (seed == 0) {
    genera<-runif(1)
    seed <-.Random.seed[3]
    }
        set.seed(seed, kinds)
a<-trt1[1:n1]
b<-trt2[1:n2]
if(randomization){
        a<-sample(trt1,n1)
        b<-sample(trt2,n2)
}
        fila<-rep(b,n1)
        columna <- a[gl(n1,n2)]
        block <- rep(1,n1*n2)
    if (r > 1) {
    for (i in 2:r) {
a<-trt1[1:n1]
b<-trt2[1:n2]
if(randomization){
        a<-sample(trt1,n1)
        b<-sample(trt2,n2)
}
        fila<-c(fila,rep(b,n1))
        columna <- c(columna,a[gl(n1,n2)])
        block <- c(block,rep(i,n1*n2))
    }}
    parameters<-list(design="strip",trt1=trt1,trt2=trt2,r=r,serie=serie,seed=seed,kinds=kinds)
    plots <- block*number+1:(n1*n2)
    book <- data.frame(plots, block = as.factor(block), column=as.factor(columna),row = as.factor(fila))
    names(book)[3] <- c(paste(deparse(substitute(trt1))))
    names(book)[4] <- c(paste(deparse(substitute(trt2))))
    outdesign<-list(parameters=parameters,book=book)
    return(outdesign)
}
