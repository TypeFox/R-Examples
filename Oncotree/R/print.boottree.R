"print.boottree" <-
function(x,...){
	orig.string <- paste(x$original$parent.num, collapse=".")
	boot.idx <- match(orig.string, as.character(x$tree.list$Tree))
	cat("Out of the", sum(x$tree.list$Freq), "replicates",
	  "there are", nrow(x$tree.list), "unique trees with frequencies from", 
	    max(x$tree.list$Freq), "down to", min(x$tree.list$Freq),"\n")
	cat("The bootstrap process found the original tree", x$tree.list$Freq[boot.idx] , "times\n")
	invisible(x)
}

