col.pvalues <-
function(file) {
 pwdfile=paste(getwd(), "/Univariate/DataTable.csv", sep="")
 file=pwdfile
 x <- read.csv(file, sep=",", header=TRUE)
 x.x = x[,3:ncol(x)]
 rownames(x.x) = x[,2]
 k = matrix(x[,1], ncol=1)
 x.n = cbind(k, x.x)
 sorted = x.n[order(x.n[,1]),]
 g = c()
 for (i in 1:nrow(sorted)) { 
  if (any(g == sorted[i,1])) {g=g} 
  else {g=matrix(c(g,sorted[i,1]), ncol=1)}
 }
NoF=nrow(g)
 all.pvalues = matrix(rep(1, ncol(sorted)-1), ncol=1)
 dirout.col = paste(getwd(), "/Univariate/Pvalues/", sep="")
 fin=ncol(sorted)-1
  for (i in 1:NoF) { 
   for (j in 1:NoF) {
   if (i < j) {
    ni=paste("Pvalues_", i,"vs", j, ".csv", sep="")
    pwdi = paste(getwd(), "/Univariate/Pvalues/", ni, sep="")
    I=read.csv(pwdi, header=TRUE)
    I = matrix(I[,-1])
    for (q in 1:fin) {
     if (I[q,] < 0.05 & all.pvalues[q,] == 1) {all.pvalues[q,] = I[q,]
     } else {all.pvalues[q,] = all.pvalues[q,]}
    }
    }
    }
    }
    colp = matrix(rep(NA, ncol(sorted)-1), ncol=1)
    for (i in 1:fin) {
      if (all.pvalues[i,] <1) {
	colp[i,] = "red"
      } else {
	colp[i,] = "black"
	}
    colnam = "Colors_Pvalues"
    assign(colnam, colp)
    write.csv(colp, paste(dirout.col, colnam, sep=""))
}
}
