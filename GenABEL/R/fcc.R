"fcc" <-
function(data,pheno) {
	if (is(data,"gwaa.data")) data <- data@gtdata else if (!is(data,"snp.data")) stop ("Wrong argument (should be \"gwaa.data\" or \"snp.data\")")
	out <- .C("fastcc",as.raw(data@gtps),as.integer(pheno),as.integer(data@nids),as.integer(data@nsnps), chi2 = double(7*data@nsnps), PACKAGE="GenABEL")$chi2
	out
}

