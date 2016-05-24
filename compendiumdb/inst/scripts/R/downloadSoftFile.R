msg <- file("message.log", open="wt")
sink(msg,type="message")
sink(msg,type="output")
options(warn = -1, showWarnCalls=FALSE)
args <- commandArgs(trailingOnly = T)

#print(args)
url1 <- "http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?targ=self&acc="
url2 <- "&form=text&view=full"

urlink <- paste(url1,args[1],url2,sep="")
download.file(urlink,paste(args[2],".soft",sep=""),quiet=TRUE)
