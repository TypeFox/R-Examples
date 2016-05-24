print.q.cor <-
function(x, var.content = NULL, initial = NULL, rnd = 2, EXPORT = FALSE, short = FALSE, ...) {
  if(class(x) != "q.cor") {stop("x must be of class 'q.cor'")}
  N.items <- nrow(x$corrs)
  if(is.null(var.content)) {var.content <- c(paste("item", 1:N.items, sep=""))}
  obj2 <- cbind(var.content, x$corrs)
  CombP <- ifelse(x$combined$AbsR[5] < .001, "***", ifelse(x$combined$AbsR[5] < .01, "**", ifelse(x$combined$AbsR[5] < .05, "*", ifelse(x$combined$AbsR[5] < .10, "+", ""))))
  FemP <- ifelse(x$fem.test$AbsR[5] < .001, "***", ifelse(x$fem.test$AbsR[5] < .01, "**", ifelse(x$fem.test$AbsR[5] < .05, "*", ifelse(x$fem.test$AbsR[5] < .10, "+", ""))))
  MaleP <- ifelse(x$male.test$AbsR[5] < .001, "***", ifelse(x$male.test$AbsR[5] < .01, "**", ifelse(x$male.test$AbsR[5] < .05, "*", ifelse(x$male.test$AbsR[5] < .10, "+", ""))))
  if(is.null(initial)) {initial <- "i"}
  if(nchar(N.items)==1) {
    rownames(obj2) <- c(paste(initial, 1:N.items, sep=""))
  }
  if(nchar(N.items)==2) {
    rownames(obj2) <- c(paste(initial, "00", 1:9, sep=""), paste(initial, "0", 10:N.items, sep=""))
  }
  if(nchar(N.items)==3) {
    rownames(obj2) <- c(paste(initial, "00", 1:9, sep=""), paste(initial, "0", 10:99, sep=""), paste(initial, 100:N.items, sep=""))
  }
  if(short==F) {
    res <- obj2[order(obj2[,2], decreasing=T),]
    res.rnd <- data.frame(res[,1], round(res[,2],rnd), res[,3], round(res[,4],rnd), res[,5], round(res[,6],rnd), res[,7], row.names=rownames(res))
    colnames(res.rnd) <- c("Content", "Combined", "p1", "Female", "p2", "Male", "p3")
    print(res.rnd)
  }
  if(short==T) {
    pos <- obj2[obj2[,2] >=0,]
    neg <- obj2[obj2[,2] < 0,]
    pos2 <- pos[order(pos[,2], decreasing=T),]
    neg2 <- neg[order(neg[,2], decreasing=F),]
    comb <- rbind(pos2, neg2)
    short.code <- ifelse(comb[,3]!="   ", 1, ifelse(comb[,5]!="   " & comb[,5]!="+  ", 1, ifelse(comb[,7]!="   " & comb[,7]!="+  ", 1, 0)))
    res <- comb[short.code==1,]
    res.rnd <- data.frame(res[,1], round(res[,2],rnd), res[,3], round(res[,4],rnd), res[,5], round(res[,6],rnd), res[,7], row.names=rownames(res))
    colnames(res.rnd) <- c("Content", "Combined", "p1", "Female", "p2", "Male", "p3")
    print(res.rnd)
  }
  cat("Note. Item content abbreviated. *** = p < .001, ** = p < .01, * = p < .05, + = p < .10. \n")
  cat("Male-Female vector correlation, r = ", round(x$vector.cor,2), ". ", sep="")
  cat("Ns are ", x$N[1], ", ", x$N[2], ", and ", x$N[3], " for Combined, Female, and Male respectively. \n", sep="")
  cat("Average absolute correlations are ", round(x$combined$AbsR[2], 3), CombP, ", ", round(x$fem.test$AbsR[2], 3), FemP, ", and ", round(x$male.test$AbsR[2], 3), MaleP, " for Combined, Female, and Male respectively. \n", sep="")
  if(EXPORT!=F) {
    write.table(res, file=paste(EXPORT, ".csv", sep=""), sep=",", col.names=NA)
  }
}
