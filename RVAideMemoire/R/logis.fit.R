logis.fit <-
function(model,int=5,...) {
   x <- model$model[,all.vars(model$call)[2]]
   y <- model$model[,all.vars(model$call)[1]]
   cutr <- cut(x,int)
   probs <- as.vector(tapply(y,cutr,sum)/table(cutr))
   int.moy <- as.vector(tapply(x,cutr,mean))
   std.err <- sqrt(probs*(1-probs)/table(cutr))
   larg <- 0.015*abs(max(x)-min(x))
   points(int.moy,probs,cex=1.2,pch=16,...)
   segments(int.moy,probs-std.err,int.moy,probs+std.err,...)
   segments(int.moy-larg,probs-std.err,int.moy+larg,probs-std.err,...)
   segments(int.moy-larg,probs+std.err,int.moy+larg,probs+std.err,...)
}

