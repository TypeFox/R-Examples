dscor <-
function(data, method=1, option=1){
x=as.data.frame(data)
option=option
f=list(
f1=function(x, option){
l <- combn(ncol(x), 2, function(y) list(data.frame(c(x[, y]))))
cor<- lapply(l, function(li) {
cr=cor.test(li[[1]], li[[2]])$ "estimate"
cr
    })
corp<- lapply(l, function(li) {
cr=cor.test(li[[1]], li[[2]])$ "p.value"
cr
    })
j <- names(x)
    aux <- combn(j, 2)
    w <- apply(aux, 2, paste, collapse = " and ")
f=diag(1,length(x))
f1=matrix(round(as.numeric(cor), 4),nrow=1)
f2=matrix(w, nrow=1)
f3=matrix(round(as.numeric(corp), 4),nrow=1)
f3=ifelse(f3<0.01, "<0.01",f3)
r=rbind(f2,f1,f3)
r=t(r)
r=as.data.frame(r)
colnames(r)=c("pairs","correlation","p-value")
r
d=diag(1,ncol(x))
co=as.numeric(cor)
cp=as.numeric(corp)
ab1=row(d)>col(d)
ab2=row(d)>col(d)
d[ab1]=co
d=t(d)
d[ab2]=cp
d=round(d,4)
colnames(d)=names(x)
rownames(d)=names(x)
res=list(
r,
d)
resp=res[[option]]
return(resp)
}
,
f2=function(x, option){
kw <- function(expr) {
localWarnings <- list()
value <- withCallingHandlers(expr,
warning = function(w) {
localWarnings[[length(localWarnings)+1]] <<- w
invokeRestart("muffleWarning")
})
value=value
}
l <- combn(ncol(x), 2, function(y) list(data.frame(c(x[, y]))))
kw(cor<- lapply(l, function(li) {
cr=cor.test(li[[1]], li[[2]], method =c("spearman"))$ "estimate"
cr
    }))
kw(corp<- lapply(l, function(li) {
cr=cor.test(li[[1]], li[[2]], method =c("spearman"))$ "p.value"
cr
    }))
j <- names(x)
    aux <- combn(j, 2)
    w <- apply(aux, 2, paste, collapse = " and ")
f=diag(1,length(x))
f1=matrix(round(as.numeric(cor), 4),nrow=1)
f2=matrix(w, nrow=1)
f3=matrix(round(as.numeric(corp), 4),nrow=1)
f3=ifelse(f3<0.01, "<0.01",f3)
r=rbind(f2,f1,f3)
r=t(r)
r=as.data.frame(r)
colnames(r)=c("pairs","correlation","p-value")
r
d=diag(1,ncol(x))
co=as.numeric(cor)
cp=as.numeric(corp)
ab1=row(d)>col(d)
ab2=row(d)>col(d)
d[ab1]=co
d=t(d)
d[ab2]=cp
d=round(d,4)
colnames(d)=names(x)
rownames(d)=names(x)
res=list(
r,
d)
resp=res[[option]]
return(resp)
}
)
g=f[[method]]
g(x,option)
}
