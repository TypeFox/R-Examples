design.split<-function (trt1, trt2,r=NULL, design=c("rcbd","crd","lsd"),serie = 2, seed = 0, kinds = "Super-Duper",
first=TRUE,randomization=TRUE )
{
    n1<-length(trt1)
    n2<-length(trt2)
    if (seed == 0) {
    genera<-runif(1)
    seed <-.Random.seed[3]
    }
    set.seed(seed,kinds)
    design <- match.arg(design)
    number<-10^serie +1
    if (design == "crd") {
        plan<-design.crd(trt1,r,serie, seed, kinds,randomization)
        k<-3
        }
    if (design == "rcbd"){
        plan<-design.rcbd(trt1,r,serie, seed, kinds, first,randomization)
        k<-3
        }
    if (design == "lsd") {
        plan<-design.lsd(trt1,serie, seed, kinds, first,randomization)
        r<-n1
        k<-4
        }
book<-plan$book
parameters<-plan$parameters
names(parameters)[2]<-"trt1"
parameters$applied<-parameters$design
parameters$design<-"split"
parameters$trt2<-trt2
j<-0
B<-list()
for(i in c(1,7,2,8,3:6)){
j<-j+1
B[[j]]<-parameters[[i]]
names(B)[j]<-names(parameters)[i]
}
nplot<-nrow(book)
d<-NULL
if(randomization){
for(i in 1:nplot)d<-rbind(d,sample(trt2,n2))
}
else{
d<-rbind(d,trt2[1:n2])
}
aa<-data.frame(book,trt2=d[,1])
for(j in 2:n2) aa<-rbind(aa,data.frame(book,trt2=d[,j]))
aa<-aa[order(aa[,1]),]
splots<-rep(gl(n2,1),nplot)
book <- data.frame(plots=aa[,1],splots,aa[,-1])
rownames(book)<-1:(nrow(book))
    names(book)[k+1] <- c(paste(deparse(substitute(trt1))))
    names(book)[k+2] <- c(paste(deparse(substitute(trt2))))
outdesign<-list(parameters=B,book=book)
return(outdesign)
}
