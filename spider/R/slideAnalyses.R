slideAnalyses <-
function(DNAbin, sppVector, width, interval = 1, distMeasures = TRUE, treeMeasures = FALSE){
	#Produce distance matrix, number of zero cells of full sequence
	boxplot_out <-FALSE
	dat <- as.matrix(DNAbin)
	dimnames(dat)[[1]] <- sppVector
	datd <- dist.dna(dat, pairwise.deletion = TRUE)
	dat_zero_out <- sum(as.numeric(datd == 0))/length(datd)
	#Create the windows
	win <- slidingWindow(dat, width, interval = interval)
	pos_out <- sapply(win, function(x) attr(x, "window")[1])
	win_dist <- lapply(win, function(x) dist.dna(x, pairwise.deletion = TRUE))
	#Distance metrics
	if(distMeasures){
		#Mean of distance matrix
		dist_mean_out <- sapply(win_dist, function(x) mean(x, na.rm=TRUE)) 
		#Number of zero cells
		zero_out <- sapply(win_dist, function(y) sum(as.numeric(y == 0), na.rm=TRUE)/length(y))
		##################
		#Threshold measures REMOVED
		#thres_above_out <- sapply(win_dist, function(x) sum( as.numeric(x >= thresA) ) )
		#thres_below_out <- sapply(win_dist, function(x) sum( as.numeric(x <= thresB) ) )
		##################
		#Diagnostic nucleotides 
		nd_out <- slideNucDiag(DNAbin, sppVector, width, interval)
		nd_out <- colSums(nd_out)
		#Nearest non-conspecific distance
		noncon_out <- sapply(win_dist, function(x) nonConDist(x, propZero = TRUE, rmNA=TRUE))
		
	}
	if(treeMeasures){
		#Produce NJ tree of full sequence
		dat_tr <- nj(datd)
		depth <- which(node.depth(dat_tr)[node.depth(dat_tr) > 1] <= median(node.depth(dat_tr)[node.depth(dat_tr) > 1]))
		#Remove windows with NA distances
		dist_na <- sapply(win_dist, function(x)  sum( as.numeric( is.na(x) ) ) )
		pos_tr_out <- pos_out[ dist_na < 1 ]
		win_tr <- lapply(win_dist[ dist_na < 1 ], nj)
		#Comparing clade composition with full alignment
		comp_out <- sapply(win_tr, function(x) tree.comp(dat_tr, x))
		comp_depth_out <- sapply(win_tr, function(x) tree.comp(dat_tr, x, method="shallow"))
		#Monophyly
		win_mono <- lapply(win_tr, function(x) monophyly(x, sppVector))
		win_mono_out <- sapply(win_mono, function(x) length(which(x))/length(x))
		}
rm(list = ls()[!ls() %in% c(ls(pattern="_out"), ls(pattern="res"))])
output <- as.list(environment())
class(output) <- "slidWin"
output
}

