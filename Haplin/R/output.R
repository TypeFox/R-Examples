output <- function(object, dirname = ".", ask = T){
##
## WRITES SUMMARY RESULTS OF haplin object TO A DIRECTORY SPECIFIED BY DIRNAME
##
if(class(object) != "haplin") stop("Intended only for haplin objects", call. = F)
#
##
.table <- haptable(object)
.summary <- capture.output(summary(object))
#
## FILE NAMES
.file.table <- paste(dirname, "/haplin_table.txt",  sep = "")
.file.summary <- paste(dirname, "/haplin_summary.txt",  sep = "")
.file.plot <- paste(dirname, "/haplin_plot.jpg",  sep = "")
#
## CHECK IF FILES ALREADY EXIST AND IF SHOULD OVERWRITE

.reply.table <- .reply.summary <- .reply.plot <- T
if(ask){
	if(file.exists(.file.table)) .reply.table <- ("y" == readline(paste("File ", .file.table, " already exists. Overwrite(y/n)?", sep = "")))
	if(file.exists(.file.summary)) .reply.summary <- ("y" == readline(paste("File ", .file.summary, " already exists. Overwrite(y/n)?", sep = "")))
	if(file.exists(.file.plot)) .reply.plot <- ("y" == readline(paste("File ", .file.plot, " already exists. Overwrite(y/n)?", sep = "")))
}

f.vis(.reply.table, vis = F)
f.vis(.reply.summary, vis = F)
f.vis(.reply.plot, vis = F)

if(.reply.table){
	dumpTab(.table, file = .file.table)
}
if(.reply.summary){
	cat(.summary, file = .file.summary, sep = "\n")
}
if(.reply.plot){
	jpeg(filename = .file.plot, res = 300, height = 2000, width = 2000)
		plot(object)
	dev.off()
}

return(invisible())
}
