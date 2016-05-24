dumpTab <- function(x, file){
##
## WRITE A TABLE TO FILE, IN print-FORMAT
##
#
#
.old <- options()
options(width = 10000)
on.exit(options(.old))
.x <- as.data.frame(x)
sink(file = file)
print(.x, quote = F, row.names = F)
sink()


}
