getMartDatasets <-
function ()
{
	if (requireNamespace("biomaRt",quietly=TRUE)) {ensemblMart <- biomaRt::useMart("ensembl")} else {stop("biomaRt package must be installed to use this function")}
	if (requireNamespace("biomaRt",quietly=TRUE)) {Datasets <- biomaRt::listDatasets(ensemblMart)} else {stop("biomaRt package must be installed to use this function")}
	return(Datasets)
}
