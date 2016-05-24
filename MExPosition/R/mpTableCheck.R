mpTableCheck <- function(data,table=NULL,make_table_nominal=TRUE){
	table<-t(designCheck(t(data),table,make_table_nominal))
	if(nrow(table) == 1){
		stop('You only have 1 table. Try a method in ExPosition')
	}
	if(ncol(table) < ncol(data)){
		stop('You have too few column-table assignments.')
	}
	if(ncol(table) > ncol(data)){
		stop('You have too many column-table assignments.')
	}
	return(table)
}