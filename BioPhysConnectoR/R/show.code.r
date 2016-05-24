#############################################
#   This code is subject to the license as stated in DESCRIPTION.
#   Using this software implies acceptance of the license terms:
#    - GPL 2
#
#   (C) by F. Hoffgaard, P. Weil, and K. Hamacher in 2009.
#
#  keul(AT)bio.tu-darmstadt.de
#
#
#  http://www.kay-hamacher.de
#############################################



show.code<-function(code=0:19,offset=0){
	seqs<-c("CYS","MET","PHE","ILE","LEU","VAL","TRP","TYR","ALA",
	  "GLY","THR","SER","ASN","GLN","ASP","GLU","HIS","ARG","LYS","PRO")
	cat("CODING:\n")
	for(i in 1:5){
	  cat(c(seqs[i],code[i]+offset),sep=": ")
	  cat("   ")}
	  cat("\n")
	for(i in 6:10){
	  cat(c(seqs[i],code[i]+offset),sep=": ")
	  cat("   ")}
	  cat("\n")
	for(i in 11:15){
	  cat(c(seqs[i],code[i]+offset),sep=": ")
	  cat("  ")}
	  cat("\n")
	for(i in 16:20){
	  cat(c(seqs[i],code[i]+offset),sep=": ")
	  cat("  ")}
	  cat("\n")
	}
