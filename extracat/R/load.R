#.onLoad <- function(...){
#	if(!"data.table" %in% .packages()){
#		#.GlobalEnv$data.table <- data.frame #function(...) data.frame(...)
#		data.table <- data.frame
#	}
#}
#.onUnload <- function(...){
#	if(!"data.table" %in% .packages()){
#		rm(data.table)
#	}
#}

