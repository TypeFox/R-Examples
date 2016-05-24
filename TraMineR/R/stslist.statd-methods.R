## ===========================
## Methods for stsstatd objects
## ===========================

print.stslist.statd <- function(x, digits=2, ...) {

	## computing max length of state label
	## to align valid states and entropy tables
	rl <- max(nchar(rownames(x$Frequencies)))
	## width <- max(nchar(colnames(x$Frequencies)), 4)	

	ident1 <- rep(" ",rl)
	ident2 <- paste(rep(" ",rl-1), collapse="")

	cat(ident1,"[State frequencies]\n")
	print(x$Frequencies, digits=digits)

	VS <- t(as.matrix(x$ValidStates))
	rownames(VS) <- paste("N",ident2,sep="")
	cat("\n", ident1,"[Valid states]\n")
	print(VS, digits=digits)
	
	H <- t(as.matrix(x$Entropy))	
	rownames(H) <- paste("H",ident2,sep="")
	cat("\n", ident1,"[Entropy index]\n")
	print(H, digits=digits)
}

"[.stslist.statd" <- function(...) {
	stop(" [!] Operation not allowed", call.=FALSE)
} 
