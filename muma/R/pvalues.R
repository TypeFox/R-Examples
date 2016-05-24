pvalues <-
function(file, mtc) {
 pwdfile=paste(getwd(), "/Univariate/DataTable.csv", sep="")
 file=pwdfile
 x <- read.csv(file, sep=",", header=TRUE)
 x.x = x[,3:ncol(x)]
 rownames(x.x) = x[,2]
 k = matrix(x[,1], ncol=1)
 x.n = cbind(k, x.x)
 sorted = x.n[order(x.n[,1]),]
 sorted.x = sorted[,-1]
 g = c()
 for (i in 1:nrow(sorted)) { 
  if (any(g == sorted[i,1])) {g=g} 
  else {g=matrix(c(g,sorted[i,1]), ncol=1)}
 }
NoF=nrow(g)
 dirout.pv = paste(getwd(), "/Univariate/Pvalues/", sep="")
 dir.create(dirout.pv)
 dirout.sign = paste(getwd(), "/Univariate/Significant_Variables/", sep="")
 dir.create(dirout.sign)
for (i in 1:NoF) {
  for (j in 1:NoF) { 
   if (i < j) {
    shapi=paste("ShapiroTest.",i,".csv",sep="")
    shapj=paste("ShapiroTest.",j,".csv",sep="")
    shap.pwdi = paste(getwd(), "/Univariate/Shapiro_Tests/", shapi, sep="")
    shap.pwdj = paste(getwd(), "/Univariate/Shapiro_Tests/", shapj, sep="")
    Si = read.csv(shap.pwdi, header=TRUE)
    Sj = read.csv(shap.pwdj, header=TRUE)
    Si = matrix(Si[,-1], ncol=1)
    Sj = matrix(Sj[,-1], ncol=1)
    welchij = paste("WelchTest_",i,"vs",j, ".csv", sep="")
    welch.pwdij = paste(getwd(), "/Univariate/WelchTest/", welchij, sep="")
    Wij = read.csv(welch.pwdij, header=TRUE)
    Wij = matrix(Wij[,-1], ncol=1)
    wmwp = paste("WMWTest_pvalues_",i,"vs",j,".csv", sep="")
    wmw.pwd = paste(getwd(), "/Univariate/Mann-Whitney_Tests/", wmwp, sep="")
    WMWp = read.csv(wmw.pwd, header=TRUE)
    WMWp = matrix(WMWp[,-1], ncol=1)
    fin=ncol(sorted)-1
    pvalues = matrix(rep(1, fin))
    for (q in 1:fin) {
      if (Si[q,] > 0.05 & Sj[q,] > 0.05) {
	  pvalues[q,]<-Wij[q,]} 
      else {
	pvalues[q,] <- WMWp[q,]}
    }
    if (mtc) {
    pval.corr = matrix(p.adjust(pvalues, method=c("BH"), n = nrow(pvalues)), ncol=1)
    rownames(pval.corr)=colnames(x.x)
    pvalfin = paste("Pvalues_", i,"vs", j, ".csv", sep="")
    assign(pvalfin, pval.corr)
    write.csv(pval.corr, paste(dirout.pv, pvalfin, sep=""))
    sign = c()
    for (q in 1:fin) {
      if (pval.corr[q,]<0.05) {
      sign = matrix(c(sign, colnames(sorted.x)[q]))
      }
    }
    signnam = paste("Significant_Variables_", i, "vs", j, ".csv", sep="")
    assign(signnam, sign)
    write.csv(sign, paste(dirout.sign, signnam, sep=""), row.names=FALSE)
    } else {
    rownames(pvalues)=colnames(x.x)
    pvalfin = paste("Pvalues_", i,"vs", j, ".csv", sep="")
    assign(pvalfin, pvalues)
    write.csv(pvalues, paste(dirout.pv, pvalfin, sep=""))
    sign = c()
    for (q in 1:fin) {
      if (pvalues[q,]<0.05) {
      sign = matrix(c(sign, colnames(sorted.x)[q]))
      }
    }
    signnam = paste("Significant_Variables_", i, "vs", j, ".csv", sep="")
    assign(signnam, sign)
    write.csv(sign, paste(dirout.sign, signnam, sep=""), row.names=FALSE)
    }
   }
  }
 }
}
