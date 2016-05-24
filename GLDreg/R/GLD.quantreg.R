GLD.quantreg <-
function(q,fit.obj,intercept="",slope="",emp=FALSE){

if(emp==FALSE){
if(tolower(intercept)=="fixed" & slope==""){

r<-sapply(q,function(i,fit.obj) {print(i);
fun.gld.slope.vary.int.fixed(i,fit.obj[[1]],fit.obj[[3]],fun=fit.obj[[1]]$fun,
param=fit.obj[[1]]$param)[[2]]}
,fit.obj)


}

if(tolower(slope)=="fixed" & intercept==""){

r<-sapply(q,function(i,fit.obj) {print(i);
fun.gld.slope.fixed.int.vary(i,fit.obj[[1]],fit.obj[[3]],fun=fit.obj[[1]]$fun,
param=fit.obj[[1]]$param)[[2]]},fit.obj)


}

if(slope=="" & intercept==""){

r<-sapply(q,function(i,fit.obj) {print(i);
fun.gld.all.vary(i,fit.obj[[1]],fit.obj[[3]],fun=fit.obj[[1]]$fun,
param=fit.obj[[1]]$param)[[2]]},fit.obj)

}  }


if(emp==TRUE){
if(tolower(intercept)=="fixed" & slope==""){

r<-sapply(q,function(i,fit.obj) {print(i);
fun.gld.slope.vary.int.fixed.emp(i,fit.obj[[1]],fit.obj[[3]])[[2]]}
,fit.obj)


}

if(tolower(slope)=="fixed" & intercept==""){

r<-sapply(q,function(i,fit.obj) {print(i);
fun.gld.slope.fixed.int.vary.emp(i,fit.obj[[1]],fit.obj[[3]])[[2]]},fit.obj)


}

if(slope=="" & intercept==""){

r<-sapply(q,function(i,fit.obj) {print(i);
fun.gld.all.vary.emp(i,fit.obj[[1]],fit.obj[[3]])[[2]]},fit.obj)

}  }



dimnames(r)[[2]]<-q

print(r)
return(r)


}
