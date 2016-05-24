## ANISOSTATIS
mpANISOSTATIS <- function(data, anisostatis.option = 'ANISOSTATIS_Type1', column.design, make.columndesign.nominal = TRUE, DESIGN = NULL, make.design.nominal = TRUE, graphs = TRUE)
{
	main <- deparse(substitute(data))
	DESIGN <- designCheck(data, DESIGN, make.design.nominal)
	
	res <- mpSTATIS(data, column.design, make.columndesign.nominal = make.columndesign.nominal, statis.prepro.option = anisostatis.option, graphs = FALSE)
	
	class(res) = c("mpANISOSTATIS","list")

	mpPlotInfo <- mpGraphs(res = res$mexPosition.Data, table = res$mexPosition.Data$Overview$groupmatrix, DESIGN = DESIGN, main = main, graphs = graphs)
	
	return(mpOutputHandler(res=res$mexPosition.Data, mpPlotInfo = mpPlotInfo))
}