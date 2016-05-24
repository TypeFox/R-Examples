get_associations_points_tw <-
function(point_data, time_window = 180, which_days = NULL, which_locations = NULL) {
	
	#### CHECK INPUTS
	if (length(dim(point_data)) != 2 | dim(point_data)[2] != 4) { stop("Invalid dimensions for point_data") }
	if (is.null(point_data)) { stop("No point_data data!") }
	
	colnames(point_data) <- c("Date","Time","ID","Location")

	# By location (s)
	if (!is.null(which_locations)) {
		point_data <- point_data[which(point_data$Location %in% which_locations),]
	}
	
	# By day(s)
	if (!is.null(which_days)) {
		point_data <- point_data[which(point_data$Date %in% which_days),]
	}


	####  Build Network
	# Look ahead to see if this is a stream from this individual, if it is, remove row from point_data and from fradj (reduces dependency of detections with centrality measures)
	del_rows <- point_data$ID[-nrow(point_data)] == point_data$ID[-1] & point_data$Time[-1] - point_data$Time[-nrow(point_data)] < (180/2)
	point_data <- point_data[-del_rows,]
	
	#fradj <- matrix(0,nrow=nrow(point_data),ncol=length(unique(point_data$ID)))
	#colnames(fradj) <- unique(point_data$ID)
	#rownames(fradj) <- point_data$ID
	
	#for (Row in c(1:(nrow(fradj)-1))) {
		#  Look ahead only so that below we can take rows from both individuals
	#	fradj[Row, which(colnames(fradj) %in% unique(as.character(point_data$ID[which(point_data$Time<=(point_data$Time[Row]+(time_window/2)) & point_data$Time>=(point_data$Time[Row]-(time_window/2)) & point_data$Location == point_data$Location[Row] )])))] <- 1
	#}
	
	get_associates <- function(Row, points, time_window, template) {
    	out <- template
    	out[which(colnames(out) %in% unique(as.character(point_data$ID[which(point_data$Time<=(point_data$Time[Row]+(time_window/2)) & point_data$Time>=(point_data$Time[Row]-(time_window/2)) & point_data$Location == point_data$Location[Row] )])))] <- 1
    	out
	}
	
	template <- matrix(0,nrow=1,ncol=length(unique(point_data$ID)))
	colnames(template) <- unique(point_data$ID)
	fradj <- do.call("rbind",lapply(c(1:nrow(point_data)),FUN=get_associates,points=point_data,time_window=time_window,template=template))
	
	rownames(fradj) <- point_data$ID
	
	list(fradj,point_data$Time,point_data$Date)
}
