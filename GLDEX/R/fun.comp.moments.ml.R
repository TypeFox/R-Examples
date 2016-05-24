"fun.comp.moments.ml" <-
function(theo.obj,data,name="ML"){
s1<-fun.theo.mv.gld(theo.obj[1,1],theo.obj[2,1],theo.obj[3,1],theo.obj[4,1],
"rs")
s2<-fun.theo.mv.gld(theo.obj[1,2],theo.obj[2,2],theo.obj[3,2],theo.obj[4,2],
"fmkl")
s3<-fun.theo.mv.gld(theo.obj[1,3],theo.obj[2,3],theo.obj[3,3],theo.obj[4,3],
"fmkl")
s4<-unlist(fun.moments(data))
r.mat<-cbind(s4,s1,s2,s3)
dimnames(r.mat)<-list(c("mean","variance","skewness","kurtosis"),
paste(c("DATA","RPRS","RMFMKL","STAR"),name))
eval.mat<-colSums(abs(r.mat[,2:4]-r.mat[,1]))
return(list("r.mat"=r.mat,"eval.mat"=eval.mat))
}

