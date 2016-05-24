f.split.matrix <- function(mat, split, tag.sep = "_"){
##
## EFFECTIVE SPLIT OF (EACH COLUMN OF) GENOTYPE MATRIX mat INTO TWO 
## COLUMNS SIDE BY SIDE IN MATRIX, SPLIT BY SEPARATOR split
## CHECKS THAT THERE'S ONLY ONE SPLIT. C;NA ETC IS SPLIT TO NA NA
#
## DIMENSION AND DIMNAMES OF ORIGINAL MATRIX
.d <- dim(mat)
.dimnames <- dimnames(mat)
#
## RESHAPE TO CHARACTER VECTOR
.mat <- as.vector(mat)
#
## SPLIT ALLELES, RESULT IS A TWO-COLUMN MATRIX
.mat <- f.split.vector(.mat, split = split)
#
## SPLIT INTO TWO FILES, ONE FOR EACH ALLELE
.mat1 <- matrix(.mat[,1], nrow = .d[1], ncol = .d[2])
.mat2 <- matrix(.mat[,2], nrow = .d[1], ncol = .d[2])
#
## DIMENSION FOR JOINED DATA SET
.d <- dim(.mat1) * c(1,2)
## NEW JOINED DATA SET
.ut <- matrix(NA_character_, nrow = .d[1], ncol = .d[2])
.ut[, seq(1, .d[2], 2)] <- .mat1
.ut[, seq(2, .d[2], 2)] <- .mat2
#
## ADD CORRECT NAMES
row.names(.ut) <- rownames(mat)
.colnames <- outer(colnames(mat), 1:2, paste, sep = tag.sep)
.colnames <- as.vector(t(.colnames))
colnames(.ut) <- .colnames
#
##
return(.ut)
}
