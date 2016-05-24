f.sort.alleles.new <- function(data, xchrom = F, sex)
{
## 
##
if(xchrom & missing(sex)) stop("Sex variable must be specified for X-chromosome
analyses!\n", call. = F)
## PREPARATIONS:
.alleles <- attr(data, "alleles")
.nalleles <- sapply(.alleles, length)
.nmarkers <- dim(data)[2]/6
#
#
## SET IN "CANONICAL" ORDERING, WITH SMALLEST FIRST:
#
for(i in seq(1, .nmarkers * 6, 2)){
	.tmpcol1 <- pmin(data[,i], data[,i+1])
	.tmpcol2 <- pmax(data[,i], data[,i+1])
	data[,i] <- .tmpcol1	
	data[,i+1] <- .tmpcol2		
}
#
## AGGREGATE DATA WITH A FREQUENCY COUNT:
if(!xchrom) .tmpd <- data
if(xchrom) .tmpd <- cbind(data, sex = sex)
#
.data.agg <- f.aggregate(.tmpd)
#
if(xchrom){# STORE SEX VARIABLE SEPARATELY
	.sex <- .data.agg$sex
	.data.agg$sex <- NULL
	#
	if(any(is.na(.sex))) stop(paste(sum(.data.agg$freq[is.na(.sex)]), " missing
    values found in sex variable! Must be removed from file before analysis.\n",
    sep = ""), call. = F)
	   .tmp <- sort(unique(.sex))
	   if(!all(is.element(.tmp, c(1,2)))) stop(paste("The sex variable is coded
       ", paste(.tmp, collapse = " "), ". It should be coded 1 (males) and 2
       (females).", sep = ""), call. = F) #
	#if(verbose) cat("\nNote: The following sex variable coding has been assumed: males = 1, females = 2")
}
.lines <- attr(.data.agg, "orig.lines")
.nlines.unique <- dim(.data.agg)[1]
#
#
## RESHAPE AND PERMUTE COLUMNS (ALL POSSIBLE PERMUTATIONS) WITHIN MOTHER, FATHER 
## AND CHILD, RESULTING
## IN 6 COLUMNS (ROW)SORTED BY MARKER AND BY PERMUTATION WITHIN MARKER, 8 ROWS PR.
## TRIAD PR MARKER.
#
.data.long <- as.matrix(.data.agg[,1:(6 * .nmarkers)])
.data.long <- t(matrix(as.numeric(t(.data.long)), nrow = 6)) # STACK DATA LINE BY LINE
.swap <- c(c(1,2,3,4,5,6), c(2,1,3,4,5,6), c(1,2,4,3,5,6), c(2,1,4,3,5,6), c(1,2,3,4,6,5), c(2,1,3,4,6,5), c(1,2,4,3,6,5), c(2,1,4,3,6,5))
.data.long <- .data.long[,.swap]
.data.long <- t(matrix(as.numeric(t(.data.long)), nrow = 6))
dimnames(.data.long) <- list(NULL, c("m1", "m2", "f1", "f2", "c1", "c2"))
#
## ADD A VECTOR ind.unique.line WHICH COUNTS LINE NUMBER AMONG THE UNIQUE LINES
## (DOES NOT REFER TO THE ORIGINAL LINE NUMBERS), AND A VECTOR ind.unique.line.marker
## WHICH IS A UNIQUE TAG FOR EACH UNIQUE LINE/MARKER COMBINATION:
if(!xchrom) .data.long <- cbind(.data.long, ind.unique.line = rep(1:.nlines.unique, each = 8*.nmarkers), ind.marker = rep(1:.nmarkers, rep(8, .nmarkers))) 
if(xchrom) .data.long <- cbind(.data.long, sex = rep(.sex, each = 8 * .nmarkers), ind.unique.line = rep(1:.nlines.unique, each = 8*.nmarkers), ind.marker = rep(1:.nmarkers, rep(8, .nmarkers))) 
#
.ind.unique.line.marker <- f.pos.in.grid(A = c(.nmarkers, .nlines.unique), comb = .data.long[,c("ind.marker", "ind.unique.line")])	
.data.long <- cbind(.data.long, ind.unique.line.marker = .ind.unique.line.marker)


#
## IDENTIFY MENDELIANLY CONSISTENT COMBINATIONS:
#
if(!xchrom){
	.valid2 <- (.data.long[,2] == .data.long[,5]) # SECOND IN MOTHER EQUALS FIRST IN CHILD
	.valid2[is.na(.valid2)] <- T # IF ONE OR BOTH ARE MISSING THE COMBINATION IS VALID
	.valid4 <- (.data.long[,4] == .data.long[,6]) # SECOND IN FATHER EQUALS SECOND IN CHILD
	.valid4[is.na(.valid4)] <- T # IF ONE OR BOTH ARE MISSING THE COMBINATION IS VALID
#
	.valid <- .valid2 & .valid4
}
if(xchrom){
	## CHECK FATHER HAVE TWO EQUAL ALLELES
	.ia3 <- is.na(.data.long[,3])
	.ia4 <- is.na(.data.long[,4])
	.valid.f <- (.data.long[,3] == .data.long[,4]) # FATHER'S ALLELES EQUAL
	.valid.f[.ia3 & .ia4] <- T # IF BOTH ARE MISSING THE COMBINATION IS VALID
	.valid.f[is.na(.valid.f)] <- F # ELSE IT IS NOT
	#
	## CHECK THAT BOYS HAVE TWO EQUAL ALLELES
	.ia5 <- is.na(.data.long[,5])
	.ia6 <- is.na(.data.long[,6])
	.girl <- (.data.long[,"sex"] == 2)
	.valid.c <- (.girl | (.ia5 & .ia6) | (.data.long[,5] == .data.long[,6])) # IF A BOY, ALLELES SHOULD BE EQUAL OR BOTH MISSING
	.valid.c[is.na(.valid.c)] <- F # ELSE IT IS NOT VALID
	#
	## CHECK THAT SECOND ALLELE OF MOTHER IS INHERITED
	.valid2 <- (.data.long[,2] == .data.long[,5]) # SECOND IN MOTHER EQUALS FIRST IN CHILD
	.valid2[is.na(.valid2)] <- T # IF ONE OR BOTH ARE MISSING THE COMBINATION IS VALID
	#
	## CHECK THAT SECOND ALLELE OF FATHER IS INHERITED BY DAUGHTERS
	.valid4 <- (!.girl | .data.long[,4] == .data.long[,6]) # IF A GIRL, SECOND IN FATHER EQUALS SECOND IN CHILD
	.valid4[is.na(.valid4)] <- T # IF ONE OR BOTH ARE MISSING THE COMBINATION IS VALID
	#
	.valid <- .valid2 & .valid4 & .valid.f & .valid.c
}# END if(xchrom)
#	

.valid.markers <- tapply(.valid, as.dframe(.data.long[,c("ind.unique.line", "ind.marker")]), any)
.valid.unique.lines <- apply(.valid.markers, 1, all) # NOTE: ASSUMES COMPLETE LIST OF INTEGERS FOR ind.unique.line, AND SORTED..
.rows.with.Mendelian.inconsistency <- unlist(.lines[!.valid.unique.lines])	
#
## REMOVE ALL INCONSISTENT LINES (WARNING: COMPLETELY INCONSISTENT WILL HERE DISAPP.!):
#
.keep <- .valid.unique.lines[.data.long[,"ind.unique.line"]] & .valid # KEEP ALL COSISTENT FOR THOSE WITH AT LEAST ONE CONSISTENT
.data.long <- .data.long[.keep,]
#
## REDUCE TO A 4-COLUMN MATRIX (STILL LONG FORMAT) FOR ONLY MOTHER AND FATHER, REPLACE THE NAs WHENEVER
# POSSIBLE, AND LEAVE REAL NAs OPEN:
#
.ia2 <- is.na(.data.long[,2])
.ia4 <- is.na(.data.long[,4])
if(!xchrom){
	.data.long[.ia2, 2] <- .data.long[.ia2, 5]
	.data.long[.ia4, 4] <- .data.long[.ia4, 6]
}
if(xchrom){
	.girl <- (.data.long[,"sex"] == 2)
	.data.long[.ia2, 2] <- .data.long[.ia2, 5]
	.data.long[.ia4 & .girl, 3] <- .data.long[.ia4 & .girl, 6]
	.data.long[.ia4 & .girl, 4] <- .data.long[.ia4 & .girl, 6]
}
#
.data.long <- .data.long[,-c(5,6)]

#
## REMOVE DUPLICATE ROWS:
#
.tag.allcol <- f.create.tag(.data.long)
.ind.unique.allcol <- !duplicated(.tag.allcol)
#
.data.long <- .data.long[.ind.unique.allcol, , drop = F]
#
#
## EXPAND ALL NAs WITH SEQUENCE OF ALL POSSIBLE ALLELES FOR THE CORRESPONDING MARKER:
#
# MERK: SKULLE VEL EGENTLIG UNNGAATT AT GUTTER
# BLE EKSPANDERT DOBBELT I xchrom, SLIK SOM DET
# UNNGAAS AT FEDRE BLIR DET!
if(!xchrom){
	for(i in 4:1){# FOR EACH OF THE 4 ALLELES
		for(j in 1:.nmarkers){
			# HOW MANY LINES TO EXPAND EACH MISSING INTO:
			.nexpand <- rep(1, dim(.data.long)[1])			
			.nexpand[is.na(.data.long[,i]) & (.data.long[,"ind.marker"] == j)] <- .nalleles[j]
			# CREATE THE INDEX FOR EXPANSION AND EXPAND:
			.ind.expand <- rep(seq(along = .nexpand), .nexpand)
			.data.long <- .data.long[.ind.expand, , drop = F]
			# REPLACE THE NAs IN THE EXPANDED DATA WITH SEQUENCE OF ALL POSSIBLE ALLELES AT MARKER:
			.data.long[is.na(.data.long[,i]) & (.data.long[,"ind.marker"] == j),i] <- 1:.nalleles[j]	
		}
	}
}
if(xchrom){
	for(i in 3:1){# FOR EACH OF THE 4 ALLELES, BUT SIMULTANEOUSLY FOR THE TWO ALLELES OF THE FATHER
		for(j in 1:.nmarkers){
			# HOW MANY LINES TO EXPAND EACH MISSING INTO:
			.nexpand <- rep(1, dim(.data.long)[1])			
			.nexpand[is.na(.data.long[,i]) & (.data.long[,"ind.marker"] == j)] <- .nalleles[j]
			# CREATE THE INDEX FOR EXPANSION AND EXPAND:
			.ind.expand <- rep(seq(along = .nexpand), .nexpand)
			.data.long <- .data.long[.ind.expand,]
			# REPLACE THE NAs IN THE EXPANDED DATA WITH SEQUENCE OF ALL POSSIBLE ALLELES AT MARKER:
			if(i == 3) .data.long[is.na(.data.long[,3]) & (.data.long[,"ind.marker"] == j),3:4] <- 1:.nalleles[j]	
			else .data.long[is.na(.data.long[,i]) & (.data.long[,"ind.marker"] == j),i] <- 1:.nalleles[j]	
		}
	}
}	# THIS MATRIX DOES NOT HAVE ANY REDUNDANT ROWS (?), BUT MAYBE AFTER PLACING ON SINGLE LINE...(?)
##
## SET MARKERS SIDE BY SIDE, WITH ALL POSSIBLE COMBINATIONS:
.line.long <- 1:(dim(.data.long)[1])
# MATRIX OF LINE NUMBERS CORRESPONDING TO EACH COMB OF UNIQUE LINES AND MARKERS:
.line.bits <- tapply(.line.long, as.dframe(.data.long[,c("ind.unique.line", "ind.marker"), drop = F]), function(x)x, simplify = F)
# CREATE ALL POSSIBLE COMBINATIONS OF MARKER GENOTYPES WITH OTHER MARKER GENOTYPES:
.line.seq <- apply(.line.bits, 1, function(x) as.numeric(t(as.matrix(do.call("expand.grid", x))))) # WARNING: IF ALL EXPANSIONS HAVE SAME LENGTH, THIS IS A MATRIX, NOT A VECTOR OF MODE LIST.
.line.seq <- as.numeric(unlist(.line.seq)) # as.numeric ALSO HANDLES .line.seq WHEN A MATRIX
# EXPAND DATA WITH ALL POSSIBLE COMBINATIONS:
.data.long <- .data.long[.line.seq,,drop = F]
.navn <- dimnames(.data.long)[[2]] # SAVE NAMES BEFORE RESHAPING
# REFORMAT DATA TO SET SIDE BY SIDE:
.data.long <- t(matrix(as.numeric(t(.data.long)), nrow = dim(.data.long)[2]*.nmarkers))
#
####
####
# SET CORRECT NAMES AFTER RESHAPING:
.navn <- as.character(t(outer(paste("l", 1:.nmarkers, sep = ""), .navn, function(x,y) paste(x,y,sep="."))))
dimnames(.data.long) <- list(NULL, .navn)
#
##
## REORDER COLUMNS AND FIX COLNAMES:
#
if(!xchrom) .indtmp <- c(as.numeric(outer(7*(0:(.nmarkers - 1)), 1:4,  "+")), 5)
if(xchrom) .indtmp <- c(as.numeric(outer(8*(0:(.nmarkers - 1)), 1:4,  "+")), c(6, 5)) # SEX IS ADDED
.data.long <- .data.long[, .indtmp, drop = F]
dimnames(.data.long)[[2]][dimnames(.data.long)[[2]] == "l1.ind.unique.line"] <- "ind.unique.line"
if(xchrom) dimnames(.data.long)[[2]][dimnames(.data.long)[[2]] == "l1.sex"] <- "sex"
#
##
##	COMPUTE HAPLOTYPE NUMBER FOR MARKER COMBINATIONS:
#
.d2 <- ncol(.data.long)
if(!xchrom) .pos.haplo <- f.pos.to.haplocomb(A = .nalleles, comb = .data.long[, 1:(.d2 - 1)])
if(xchrom) .pos.haplo <- f.pos.to.haplocomb(A = .nalleles, comb = .data.long[,1:(.d2 - 2)])
.haplo.comb <- f.pos.to.haplocomb(pos = .pos.haplo, A = prod(.nalleles))

colnames(.haplo.comb) <- c("m1", "m2", "f1", "f2")

#
##
## PREPARE OUTPUT WITH HAPLOTYPE NUMBERS, INDEX OF UNIQUE LINES (THOSE REMAINING AFTER REMOVING MENDLIAN INCONSISTENCIES)
## AND THE CORRESPONDING FREQUENCIES:
## NOTE: ind.unique.line REFERS TO THE SEQUENCE 1:.nlines.unique (THOSE REMAINING) AND NOT THE ORIGINAL LINE NUMBERS
#


### BOER SJEKKE DENNE?:
##
##
if(!xchrom) .ut <- cbind(comb = .haplo.comb, ind.unique.line = .data.long[,"ind.unique.line"], freq = .data.agg$freq[.data.long[,"ind.unique.line"]])
if(xchrom) .ut <- cbind(comb = .haplo.comb, sex = .data.long[,"sex"], ind.unique.line = .data.long[,"ind.unique.line"], freq = .data.agg$freq[.data.long[,"ind.unique.line"]])
.ut <- as.dframe(.ut)
#
attr(.ut, "alleles") <- .alleles
attr(.ut, "rows.with.Mendelian.inconsistency") <- .rows.with.Mendelian.inconsistency
attr(.ut, "orig.lines") <- .lines # NOTE: .lines CAN BE INDEXED BY .tag AND .tag.unique (AS CHARACTER) OR BY ind.unique.line 
#	
return(.ut)

}
