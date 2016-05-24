wtd.chi.sq <- function(var1, var2, var3=NULL, weight=NULL, na.rm=TRUE, drop.missing.levels=TRUE, mean1=TRUE){
  if(is.null(weight)){
    weight <- rep(1, length(var1))
  }
  if(mean1==TRUE)
      weight <- weight/mean(weight, na.rm=TRUE)
  if(na.rm==TRUE){
    filt <- (!is.na(var1) & !is.na(var2))
    if(!is.null(var3))
      filt <- (!is.na(var1) & !is.na(var2) & !is.na(var3))
    var1 <- var1[filt]
    var2 <- var2[filt]
    if(!is.null(var3))
      var3 <- var3[filt]
    weight <- weight[filt]
  }
  if(drop.missing.levels==TRUE){
    var1 <- drop.levels(var1)
    var2 <- drop.levels(var2)
    if(!is.null(var3))
      var3 <- drop.levels(var3)
  }
  var12set <- unlist(summary(xtabs(weight~var1+var2))[c("statistic", "parameter", "p.value")])
  names(var12set) <- c("Chisq", "df", "p.value")
  out <- var12set
  if(!is.null(var3)){
    as.numeric(var3)
    var123set <- unlist(summary(xtabs(weight~var1+var2+var3))[c("statistic", "parameter", "p.value")])
    var13set <- unlist(summary(xtabs(weight~var1+var3))[c("statistic", "parameter", "p.value")])
    var23set <- unlist(summary(xtabs(weight~var2+var3))[c("statistic", "parameter", "p.value")])
    v12sep <- unlist(sapply(unique(var3), function(x) unlist(summary(xtabs(weight[var3==x]~var1[var3==x]+var2[var3==x]))[c("statistic", "parameter")])))
    acvar3 <- c(sum(v12sep[1,])-var12set[1], sum(v12sep[2,])-var12set[2], 1-pchisq(sum(v12sep[1,])-var12set[1], sum(v12sep[2,])-var12set[2]))
    v13sep <- unlist(sapply(unique(var2), function(x) unlist(summary(xtabs(weight[var2==x]~var1[var2==x]+var3[var2==x]))[c("statistic", "parameter")])))
    acvar2 <- c(sum(v13sep[1,])-var13set[1], sum(v13sep[2,])-var13set[2], 1-pchisq(sum(v13sep[1,])-var13set[1], sum(v13sep[2,])-var13set[2]))
    v23sep <- unlist(sapply(unique(var1), function(x) unlist(summary(xtabs(weight[var1==x]~var2[var1==x]+var3[var1==x]))[c("statistic", "parameter")])))
    acvar1 <- c(sum(v23sep[1,])-var23set[1], sum(v23sep[2,])-var23set[2], 1-pchisq(sum(v23sep[1,])-var23set[1], sum(v23sep[2,])-var23set[2]))
    out <- rbind(var123set, var12set, var13set, var23set, acvar1, acvar2, acvar3)
    colnames(out) <- c("Chisq", "df", "p.value")
    rownames(out) <- c("Using All Three Variables", "Using Var1 and Var2", "Using Var1 and Var3", "Using Var2 and Var3", "Difference Across Levels of Var1", "Difference Across Levels of Var2", "Difference Across Levels of Var3")
  }
  out
}
