na.omit.fdata=function(object,...){
n=nrow(object)
r=apply(object[["data"]],1,count.na)
omit=which(r)
xx=object
if (length(which(r))>0L) {
xx=object[-omit]
temp=setNames(c(1:n)[which(r)],attr(object[["data"]],"dimnames")[[1]][omit])
attr(temp,"class")<-"omit"                                                                                                      
attr(xx,"na.action")<-temp
}
return(invisible(xx))
}
na.fail.fdata=function(object,...){
ok<-complete.cases(object[["data"]])
if (all(ok)) 
invisible(	object)
else stop("missing values in fdata object")
}
