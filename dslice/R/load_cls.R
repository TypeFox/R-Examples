load_cls <- function(file = NULL)
{
	cls_info <- readLines(file)
	label_list <- unlist(strsplit(cls_info[3], " "))
	uniquelist <- sort(unique(label_list))
	y <- vector(length = length(label_list), mode = "integer")
	for(i in 1:length(label_list)){
		y[i] <- which(uniquelist == label_list[i]) - 1L
	}
	return(list(pheotype = uniquelist, value = y))
}