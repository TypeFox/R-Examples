getGeneSeq <- function(chr, start, end, organism){
  # Options:
  # Organism: susScr3
  reply <- scan(paste("http://genome.ucsc.edu/cgi-bin/das/",organism,"/dna?segment=chr",chr,":",start,",",end,sep=""),what="raw")
  getStart <- which(grepl("length=",reply))+1
  getEnd <- which(grepl("/DNA",reply))-1
  paste(reply[getStart:getEnd],collapse="")
}