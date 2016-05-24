welch <-
function(file){
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
 dirout.w = paste(getwd(), "/Univariate/WelchTest", sep="")
 dir.create(dirout.w)
 # Perform the Welch test for the effective combinations
 for (i in 1:NoF) {
  for (j in 1:NoF) { 
   if (i < j) {
    ni=paste("r.",i,".csv",sep="")
    nj=paste("r.",j,".csv",sep="")
    pwdi = paste(getwd(), "/Univariate/Groups/", ni, sep="")
    pwdj = paste(getwd(), "/Univariate/Groups/", nj, sep="")
    I=read.csv(pwdi, header=TRUE)
    J=read.csv(pwdj, header=TRUE)
    I = I[,-1]
    J = J[,-1]
    fin=ncol(sorted)-1
    we <- matrix(rep(NA, fin))
    for (q in 1:fin) {
    we[q,] <- t.test(I[,q],J[,q],var.equal=F,conf.level=.95,alternative="two.sided")$p.value
    }
    welch.ij=paste("WelchTest_",i,"vs",j, ".csv", sep="")
    assign(welch.ij,we)
    write.csv(we, paste(dirout.w, welch.ij, sep="/"))
   }
  }
 }
}
