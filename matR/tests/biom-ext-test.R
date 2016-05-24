library(matR)
N <- 1:4
List <- mget (paste0 ("xx", N), inherits=TRUE)

#-----------------------------------------------------------------------------------------
#  OK FOR CRAN
#-----------------------------------------------------------------------------------------

f <- function (x, y) {										# merge and check integrity
	z <- merge (x,y)
#	applyBiomMethods (z)									# truncated for CRAN
#	message ("rows:\t", nrow(x), "(x)\t", nrow(y), "(y)\t", nrow(z), "(merge)\t", 
#		length (intersect (rownames(x), rownames(y))), "(in common)\n")
#	message ("cols:\t", ncol(x), "(x)\t", ncol(y), "(y)\t", ncol(z), "(merge)\t", 
#		length (intersect (colnames(x), colnames(y))), "(in common)\n")
	}

for (j in matrix2list (t (combn (N, 2)))) {
#-----------------------------------------------------------------------------------------
#  merge()
#-----------------------------------------------------------------------------------------
	x <- List [[j [1]]]
	y <- List [[j [2]]]
	f(x,y) ; f(y,x)											# merge each pair
	}

for (x in List) {
#-----------------------------------------------------------------------------------------
#  rows() and columns()
#-----------------------------------------------------------------------------------------
	str (rows (x))											# all row/column annotations
	str (columns (x))
	str (rows (x, "a"))										# many matches
	str (columns (x, "a"))
	str (rows (x, "syzygy"))								# no match --- not sure this result is correct
	str (columns (x, "syzygy"))
	}	

for (x in List) {
#-----------------------------------------------------------------------------------------
#  rows<-() and columns<-()
#-----------------------------------------------------------------------------------------
	rows (x, "junk") <- 1:nrow(x)							# assign data of adequate length
	columns (x, "junk") <- 1:ncol(x)
	rows (x, "junk") <- 1:2									# assign too-short data
	columns (x, "junk") <- 1:2
	}

for (x in List) {
#-----------------------------------------------------------------------------------------
#  subselection
#-----------------------------------------------------------------------------------------
	i <- rep (c(TRUE,FALSE), len=nrow(x))
	j <- rep (c(TRUE,FALSE), len=ncol(x))
	x [i, ]													# just test a few random things
	x [,j]
	x [i,j]	
	x [1,]													# but the special case of one index, in particular
	x [,1]
	}

for (x in List) {
#-----------------------------------------------------------------------------------------
#  dimnames<-()
#-----------------------------------------------------------------------------------------
	rownames(x) <- 1:nrow(x)								# assign data of adequate length
	colnames(x) <- 1:ncol(x)
	rownames(x) <- 1										# assign too-short data
	colnames(x) <- 1
	}
