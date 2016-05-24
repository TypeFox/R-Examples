goodnessOfFit <- function(obj, iC10=1:10, newdata=NULL,...) {
   UseMethod("goodnessOfFit")
}
goodnessOfFit.iC10 <- function(obj, iC10=1:10, newdata=NULL,...) {
    if (is.null(newdata)) {
       all.data <- obj$fitted
    } else {
       all.data <- rbind(newdata$CN, newdata$Exp)
      }
    td <- matrix(NA, nrow(all.data), 10)
    for (i in iC10) {
        td[,i] <- unclass(apply(all.data[,which(obj$class==i),drop=F], 1, mean))
    }
    cat("\nCorrelation between centroids and predicted profiles:\n")
    x <- rep(NA, length(iC10))
    for (i in iC10) {
        cat("iC", i, ": ", sep="")
	x[i] <- cor(obj$centroids[,i], td[,i],use="pairwise.complete.obs",...)
        cat(round(x[i], 3), "\n")
    }
   total <- weighted.mean(x, w=table(obj$class)[iC10], na.rm=TRUE)
   cat("\nOverall correlation:", round(total, 3), "\n")
   indiv <- rep(NA, ncol(all.data))
   for (i in 1:length(indiv)) {
      indiv[i] <- cor(all.data[,i], obj$centroids[,obj$class[i]],...)
}
   res <- list(cor=x, total=total, indiv=indiv)
   res
}
