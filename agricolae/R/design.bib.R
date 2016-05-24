`design.bib` <-
function (trt, k, r=NULL, serie = 2, seed = 0, kinds = "Super-Duper", maxRep=20,randomization=TRUE)
{
v<-length(trt)
number<-10
if(serie>0) number<-10^serie
    ntr <- length(trt)
    if (seed == 0) {
    genera<-runif(1)
    seed <-.Random.seed[3]
    }
    set.seed(seed, kinds)
#----------
# question r
if(!is.null(r)){
x<- r*(k-1)/(v-1)
b<- v*r/k
y<- ceiling(x)
z<- ceiling(b)
if(!(x==y & b==z)) {
# change r
r<-NULL
for (i in 2:maxRep){
x<- i*(k-1)/(v-1)
b<- v*i/k
y<- ceiling(x)
z<- ceiling(b)
if(x==y & b==z)r<-c(r,i)
}
if(!is.null(r) ) return(cat("\n Change r by ",paste(r, collapse = ", "),"...\n"))
else return(cat("Other k <> ",k,"; 1<k<",v,"\n"))
}
}

if(is.null(r)){
for (i in 2:maxRep){
x<- i*(k-1)/(v-1)
b<- v*i/k
y<- ceiling(x)
z<- ceiling(b)
if(x==y & b==z)r<-c(r,i)
}
r<-r[1]
}
b<- v*r/k
if (requireNamespace("AlgDesign", quietly = TRUE)) {
initial <- AlgDesign::optBlock(~., withinData = factor(1:v), blocksizes = rep(k,b))$row
md <- matrix(initial, byrow = TRUE, ncol = k)
}
#----------
    b<-nrow(md)
    bp<-1:b
    if(randomization)bp<-sample(1:b,b)
    md<- md[bp,]
    for (i in 1:b) {
    bi<-1:k
    if(randomization)bi<-sample(1:k,k)
    md[i,]<- md[i,bi]
    }
mtr<-trt[t(md)]
block <- gl(b,k)
Rep<-as.numeric(block)
plots <- Rep*number+(1:k)
parameters<-list(design="bib",trt=trt,k=k,serie=serie,seed=seed,kinds=kinds)
#plots <- number + 1:(b*k) - 1
book <- data.frame(plots, block = as.factor(block), trt = as.factor(mtr))
names(book)[3] <- c(paste(deparse(substitute(trt))))
r<-as.numeric(table(book[,3])[1])
lambda<-r*(k-1)/(ntr-1)
E<-lambda*ntr/(r*k)
cat("\nParameters BIB\n==============")
cat("\nLambda     :",lambda)
cat("\ntreatmeans :",ntr)
cat("\nBlock size :",k)
cat("\nBlocks     :",b)
cat("\nReplication:",r,"\n")
cat("\nEfficiency factor",E,"\n\n<<< Book >>>\n")
statistics<-data.frame(lambda= lambda,treatmeans=ntr,blockSize=k,blocks=b,r=r,Efficiency=E)
rownames(statistics)<-"values"
outdesign<-list(parameters=parameters,statistics=statistics,
sketch=matrix(book[,3], byrow = TRUE, ncol = k),book=book)
return(outdesign)
}
