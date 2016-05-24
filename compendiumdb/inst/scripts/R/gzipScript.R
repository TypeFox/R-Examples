options("warn"=-1)
msg <- file("message.log", open="wt")
sink(msg,type="message",append=TRUE)

library(GEOquery)
filename <- commandArgs(trailingOnly = T)
#gds <- gzfile(filename[1])
#write(readLines(gds),paste(filename[2],".soft",sep=""))
gunzip(filename[1],paste(filename[2],".soft",sep=""),remove=FALSE)