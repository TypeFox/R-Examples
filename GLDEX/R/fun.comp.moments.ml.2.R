"fun.comp.moments.ml.2" <-
function(theo.obj, data, name = "ML")
{
s1 <- fun.theo.mv.gld(theo.obj[1, 1], theo.obj[2, 1], theo.obj[3, 1], 
theo.obj[4, 1], "fmkl")
s2 <- fun.theo.mv.gld(theo.obj[1, 2], theo.obj[2, 2], theo.obj[3, 2], 
theo.obj[4, 2], "fmkl")
s3 <- unlist(fun.moments(data))
r.mat <- cbind(s3, s1, s2)
dimnames(r.mat) <- list(c("mean", "variance", "skewness", "kurtosis"), 
paste(c("DATA", "RMFMKL", "STAR"), name))
eval.mat <- colSums(abs(r.mat[, 2:3] - r.mat[, 1]))
return(list("r.mat"=r.mat, "eval.mat"=eval.mat))
}

