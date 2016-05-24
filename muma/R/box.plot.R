box.plot <-
function(file) {
 pwdfile=paste(getwd(), "/Univariate/DataTable.csv", sep="")
 file=pwdfile
 x <- read.csv(file, sep=",", header=TRUE)
 x.x = x[,3:ncol(x)]
 rownames(x.x) = x[,2]
 k = matrix(x[,1], ncol=1)
 colnames(k)[1]="Class"
  x.x = cbind(k, x.x)
 	sorted = x.x[order(x.x[,1]),]
 	sorted.x = as.matrix(sorted[,-1], ncol=ncol(sorted)-1)
 	g = c()
 	for (i in 1:nrow(sorted)) {
  		if (any(g == sorted[i,1])) {
  			g=g} else {g=matrix(c(g,sorted[i,1]), ncol=1)
  				}
 	}
	NoF=nrow(g)
 	dirbox = paste(getwd(), "/Univariate/BoxPlot/", sep="")
 	dir.create(dirbox)
 	for (i in 2:ncol(x.x)) {
 		name = paste(getwd(), "/Univariate/BoxPlot/", colnames(x.x)[i], ".pdf", sep="")
 		pdf(name)
 		boxplot(x.x[,i] ~ x.x[,1], boxfill=c(seq(1,NoF)), ylab=colnames(x.x)[i], xlab ="Groups", border="grey30", main = paste("Boxplot ", colnames(x.x)[i], sep=""))
 		dev.off()
 	}
 	
 }
