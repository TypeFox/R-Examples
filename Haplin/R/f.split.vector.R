f.split.vector <- function(vec, split, na.strings = "NA"){
##
## EFFECTIVE SPLIT OF GENOTYPE VECTOR vec INTO TWO COLUMN MATRIX, SPLIT BY SEPARATOR split
## CHECKS THAT THERE'S ONLY ONE SPLIT. C;NA ETC IS SPLIT TO NA NA
##
#
## FIND ALL UNIQUE GENETIC VALUES, NA FIRST
.unique <- sort(unique(vec), na.last = F)
#
## TO SAVE TIME, SPLIT ONLY THE UNIQUE VALUES
.unique.split <- strsplit(.unique, split = split, fixed = T)
## JUST TO MAKE SURE, CHECK MISSING
.una <- which(is.na(.unique.split)) # SHOULD BE 1 IF MISSING, OTHERWISE EMPTY
if(length(.una) > 1) stop("Something's wrong with the data reading!\n")
#
## REPLACE THE SINGLE MISSING WITH TWO MISSING
if(length(.una) == 1) .unique.split[[.una]] <- c(NA_character_, NA_character_)
#
## CHECK THAT SPLITTING IS OK BY CHECKING THAT ALL GENETIC VALUES HAVE BEEN SPLIT INTO TWO
.ind <- sapply(.unique.split, length)
.ind <- which(.ind != 2)
if(length(.ind) > 0){
	cat("\n-----------------------\nSomething's wrong with the separator\nin the following data values:\n")
	print(.unique[.ind])
	stop('Correct data vector or check "split" argument!', call. = F)
}
#
## CONVERT TO MATRIX WITH TWO COLUMNS, ONE FOR EACH ALLELE
.unique.split <- t(sapply(.unique.split, function(x)x))
if(dim(.unique.split)[1] != length(.unique)) stop()
#
## MAKE SURE THERE ARE NO LEFTOVER MISSING OF THE TYPE NA;T ETC,
## RECODE THE T (OR WHATEVER) TO MISSING.
## THOSE OF TYPE NA;NA SHOULD ALREADY HAVE BEEN TAKEN CARE OF
.unique.split[.unique.split == na.strings] <- NA
.rna <- rowSums(is.na(.unique.split))
.unique.split[.rna > 0, ] <- NA
#
## EXPAND .unique.split TO MATCH THE vec LENGTH
.mat <- .unique.split[match(vec, .unique), ]
#
##
return(.mat)
}
