corplot <-
function(mat){

	mat.p<-mat[,-ncol(mat)]
	mat.p<-log10(mat.p)

	colnames(mat.p)<-gsub("_.*Pos","",colnames(mat.p))

	par(pch=20)
	par(cex=0.1)
	pairs(as.data.frame(mat.p), lower.panel=panel.smooth.r, upper.panel=panel.cor,gap=0.1,cex.labels=0.8/max(strwidth(colnames(mat.p))))
}
