`design.ab` <-
function(trt, r=NULL,serie=2,design=c("rcbd","crd","lsd"),seed=0,kinds="Super-Duper",
first=TRUE,randomization=TRUE ){
design <- match.arg(design)
if( design=="rcbd" | design=="crd") posicion <- 3
else posicion <- 4
serie<-serie; seed<-seed; kinds<-kinds; first<-first;

# Process to trt to factorial
ntr<-length(trt)
fact<-NULL
tr0<-1:trt[1]
k<-0
a<-trt[1];b<-trt[2]
for(i in 1:a){
for(j in 1:b){
k<-k+1
fact[k]<-paste(tr0[i],j)
}
}

if(ntr >2) {
for(m in 3:ntr){
k<-0
tr0<-fact
fact<-NULL
a<-a*b
b<-trt[m]
for(i in 1:a){
for(j in 1:b){
k<-k+1
fact[k]<-paste(tr0[i],j)
}
}
}
}
#------------------------------
if(design=="rcbd")plan<-design.rcbd(trt=fact, r, serie, seed, kinds, first,randomization )
if(design=="crd")plan<-design.crd(trt=fact, r, serie, seed, kinds,randomization)
if(design=="lsd")plan<-design.lsd(trt=fact, serie, seed, kinds, first,randomization )
parameters<-plan$parameters
parameters$applied<-parameters$design
parameters$design<-"factorial"
plan<-plan$book
trt<-as.character(plan[,posicion])
nplan<-nrow(plan)
A<-rep(" ",nplan*ntr)
dim(A)<-c(nplan,ntr)
colnames(A)<-LETTERS[1:ntr]

for(i in 1:nplan) {
A[i,]<-unlist(strsplit(trt[i], " "))
}
A<-as.data.frame(A)
book<-data.frame(plan[,1:(posicion-1)],A)
outdesign<-list(parameters=parameters,book=book)
return(outdesign)
}
