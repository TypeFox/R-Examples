summary.whatif <- function (object, ...)  {

 #Calculate number of counterfactuals
 m <- length(object$in.hull)
 #Calculate number of counterfactuals in convex hull
 m.inhull <- sum(object$in.hull)
 #Calculate average (over all counterfactuals) percent data `nearby'
 mean.near <- mean(object$sum.stat)
 #Create data frame combining results of convex hull test with percent data `nearby'
 sum.df <- data.frame(cfact = seq(1, m, by = 1), in.hull = object$in.hull, 
   per.near = object$sum.stat)

out <- list(call = object$call, m = m, m.inhull = m.inhull, mean.near = mean.near, sum.df = sum.df)
 class(out) <- "summary.whatif"
 return(out)

}
