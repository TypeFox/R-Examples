zigzag<-function(outdesign){
parameters <- outdesign$parameters
design <- parameters$design
book<-outdesign$book
fieldbook<-book
n<-nrow(book)
if(design=="lattice") {
nr<-parameters$r; nc<-1:2
t1<-sqrt(length(parameters$trt))
t2<-t1
}
if(design=="cyclic") {
nr<-nlevels(as.factor(book[,2])); nc<-1:2
t1<-length(parameters$trt)
t2<-n/(nr*t1)
}
if(design=="alpha") {
nr<-parameters$r; nc<-1:2
t2<-parameters$k
t1<-n/(nr*t2)
}
if(design=="strip") {
nr<-parameters$r; nt1<-3; nt2<-4; nc<-1
t2<-length(parameters$trt2)
t1<-length(parameters$trt1)
}
if(design=="split"){
nro<-paste(book[,1],book[,2],sep="-")
nr<-1;  nc<-1:2
t1<-length(parameters$trt1)
t2<-n/t1;
book<-data.frame(nro,book[,3])
}
if(design=="bib") {
t1<-outdesign$statistics$treatmeans
nr<-1;  t2<-n/t1; nc<-1
}
if(design=="youden") {
t2<-parameters$r
nr<-1;  t1<-n/t2; nc<-1
}
if(design=="rcbd" | design=="lsd" |design=="graeco" |design=="factorial"){
t1<-parameters$r
nr<-1;  t2<-n/t1; nc<-1
}
if(design=="dau") {
plots<-book[,1]
ntb<-tapply.stat(book[,3],book[,2],length)[,2]
ntb<-cumsum(ntb)
t1<-nlevels(book[,2])
nc<-1
for(j in seq(2,t1,2)){
x2<-ntb[j]
x1<-ntb[j-1]+1
x3<-plots[x1:x2]
x3<-x3[order(x3,decreasing=TRUE)]
plots[x1:x2]<-x3
}
}
#------------------------------
if(design != "dau"){
r<-nr
X<-array(1:n,c(t2,t1,r))
for(i in 1:r){
for(j in seq(2,t1,2)){
X[,j,i]<-X[order(X[,j,i],decreasing=TRUE),j,i]
}
}
x<-as.numeric(X)
plots<-fieldbook[x,nc]
}
fieldbook[,nc]<-plots
return(fieldbook)
}

