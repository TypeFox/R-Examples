.onLoad <- function(libname,pkgname){
#	require(vcd)
#	require(gWidgets)
	options("guiToolkit"="tcltk")
#	xx <- addStockIcons("To_col_vars",getStockIcons()$forward)
#	xx <- addStockIcons("To_row_vars",getStockIcons()$backward)
	enopt <- getOption("ENmisc")
	if (is.null(enopt$mosaicpalette)) enopt$mosaicpalette <- "RdYlGn"
	options(ENmisc=enopt)

#	ENmiscEnvironment <- new.env()
	
#    putENmisc("putENmisc",putENmisc)

#	putENmisc("getENmisc",getENmisc)

#	putENmisc("ENmiscEnv",ENmiscEnv)

#	rm(getENmisc,putENmisc,ENmiscEnv)
	
}
