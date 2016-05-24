runQC <- function(QC, nPath=NULL, B=1e4, pvalCut=.05, pvalAdjust=FALSE, fileForCQCp="c2.all.v3.0.symbols.gmt") {
	QC$RunQC(nPath, B, pvalCut, pvalAdjust, fileForCQCp)
}
