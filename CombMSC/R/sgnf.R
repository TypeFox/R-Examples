`sgnf` <-
function(obj){
  if(!("msc" %in% class(obj))) stop("Function only defined for objects of class msc")
  temp <- summary(obj)
  temp2 <- as.data.frame(t(sapply(temp$Corners, rbind))[,-(1:length(obj$msc.List))])
  temp3 <- as.data.frame(t(sapply(temp$Minima, function(x) rbind(x[1,])))[,-(1:length(obj$msc.List))])
  res <- apply(temp2, 2, function(x)min(unlist(x))) - apply(temp3, 2, function(x)min(unlist(x)))
  names(res) <- names(obj$summary.Functions)
  res <- rbind(res, res/ apply(temp2, 2, function(x)min(unlist(x))))
  rownames(res) <- c("Absolute", "Proportion")
  zapsmall(res)
}

