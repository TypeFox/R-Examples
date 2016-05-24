`convertTime` <-
function(tmp.time) {
	tmp <- tail(c(0, 0, rev(as.numeric(rev(unlist(strsplit(sub("sec", "", tmp.time), ":")))))), 3)
	tmp.label <- c("Hour", "Minute", "Second")
	tmp.label[which(tmp!=1)] <- paste(tmp.label, "s", sep="")[which(tmp!=1)]
	return(paste(paste(tmp[tmp!=0], tmp.label[tmp!=0]), collapse=", "))
} ### END convertTime
