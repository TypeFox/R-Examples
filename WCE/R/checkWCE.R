 checkWCE <- function(data, id='Id', event = 'Event',  start='Start', stop='Stop', expos ='dose'){
	Error <- NULL
      if (is.data.frame(data) == F)  {Error = "The data supplied is not in a data.frame format.\n"}
	names(data)[names(data) == id] = 'Id'
	names(data)[names(data) == event] = 'Event'
	names(data)[names(data) == start] = 'Start'
	names(data)[names(data) == stop] = 'Stop'
	names(data)[names(data) == expos] = 'dose'
	if (sum(c("Id","Event","Start","Stop","dose") %in% names(data))!= 5)  {Error <- c(Error,"At least one of id, event, start, stop, or expos is missing from the dataset supplied.\n")} else {
	if (sum(data$Start == data$Stop) > 0)  {Error <- c(Error,"Start and stop values are equal for at least one individual.\n")}
	if (nrow(na.omit(data[,c("Id","Event","Start","Stop","dose")])) < nrow(data[,c("Id","Event","Start","Stop","dose")]))  {Error <- c(Error,"There are missing values in at least one of the id, event, start, stop, or expos columns.\n")}
	if (.gap(data))  {Error <- c(Error,"There is at least one gap or an overlap for at least one individual in the start and stop values.\n")}	
      if (!is.numeric(data[,c("Event")]) | !is.numeric(data[,c("dose")]) | !is.numeric(data[,c("Start")]) |!is.numeric(data[,c("Stop")]))   {Error <- c(Error,"At least one value of start, stop, event or expos is not numeric.\n")}
	if (length(unique(data[,c("Event")]))!=2 | (sum(unique(data[,c("Event")]) %in% c(0,1))!=2))  {Error <- c(Error,"Values supplied for events are not either 0 or 1.\n")}
}
	if (is.null(Error) == F) {Message <- Error} else {Message <- "Data are in the right format for WCE estimation.\n"}
	cat(Message)}