#########################################################################
# File Name: ScoringMatrix.R
# Author: Ge Tan
# mail: gtan@me.com
# Created Time: Tue Jul 29 20:48:49 2014
#########################################################################

### ------------------------------------------------------------------
### CharacterSet
### Not Exported!
.AAOrder <- c("A","R","N","D","C","Q","E","G","H","I",
              "L","K","M","F","P","S","T","W","Y","V")
AAGapCharacterSet <- c("A","R","N","D","C","Q","E","G","H","I",
                        "L","K","M","F","P","S","T","W","Y","V", "-")
AACharacterSet <- .AAOrder
AmbiguousAACharacterSet <- c("A","R","N","D","C","Q","E","G","H","I",
                              "L","K","M","F","P","S","T","W","Y","V",
                              "B","Z","J","X")
AmbiguousAAGapCharacterSet <- c("A","R","N","D","C","Q","E","G","H","I",
                                 "L","K","M","F","P","S","T","W","Y","V",
                                 "B","Z","J","X", "-")
DNACharacterSet <- c("A", "C", "G", "T")
DNAGapCharacterSet <- c("A", "C", "G", "T", "-")
RNACharacterSet <- c("A", "C", "G", "U")
RNAGapCharacterSet <- c("A", "C", "G", "U", "-")

### ------------------------------------------------------------------
### AA to Int index
### Exported!
AAToInt <- function(AA){
  mapping <- 1:length(AACharacterSet)
  names(mapping) <- AACharacterSet
  return(mapping[AA])
}

.validatePAMMatrix <- function(PAM){
  ## First test whether it is a matrix
  if(!is.matrix(PAM))
    stop("The input must be an object of matrix.")
  ## Second test whether it has correct row number and column number
  if(nrow(PAM) != 20L || ncol(PAM) != 20L)
    stop("The input matrix must be a 20 * 20 square matrix")
  ## Third test the AA order.
  if(!identical(rownames(PAM), .AAOrder))
    stop("The row names are not identical with default AA order: ", .AAOrder)
  if(!identical(colnames(PAM), .AAOrder))
    stop("The column names are not identical with default AA order: ", .AAOrder)
  ## Fourth test non-negative value
  if(any(PAM < 0))
    stop("The PAM matrix must be non-negative")
  ## Fifth test the PAM row sum is 1
  if(!all(sapply(unname(rowSums(PAM)), all.equal, 1)))
    stop("The row sum of PAM matrix must be 1")
  return("success")
}

.validateBF <- function(BF){
  ## First test whether it is a vector
  if(!is.vector(BF))
    stop("The input must be an object of vector")
  ## Second test whether it has correct number of elements
  if(length(BF) != 20L)
    stop("The input vector must be a vector with 20 elements")
  ## Third test the AA order
  if(!identical(names(BF), .AAOrder))
    stop("The names of input vector are not identical with default AA order: ",
         .AAOrder)
  ## Fourth test non-negative value
  if(any(BF < 0))
    stop("The background frequency must be non-negative")
  ## Fifth test the sum is 1
  if(!all.equal(sum(BF), 1))
    stop("The sum of background frequency must be 1")
  return("success")
}


### -----------------------------------------------------------------
### Calculate the n-PAM matrices from PAM1 mutation matrix and n.
### To compute n-PAM matrices, we multiply the PAM1 matrix through itself 
### N times, which is most efficiently achieved through 
### n additions in log space.
### Different from Darwin, here PAM1[i,j] means the transition probability
### from i to j.
### Exported!
PAMn <- function(PAM1, n){
  ## Validated by Darwin
  .validatePAMMatrix(PAM1)
  ans <- expm.Higham08(n * logm(PAM1))
  dimnames(ans) <- dimnames(PAM1)
  return(ans)
}
### PAM250 <- PAMn(GONNET, 250)

### -----------------------------------------------------------------
### Computing Dayhoff matrices from PAM mutation matrices and AA frequency.
### In the introduction, we motivated the quest for mutation matrices 
### through the need of a method to quantify the evolutionary relationship 
### of two AA. We now have such matrices. 
### Remember that we are interested in the ratio: 
### P("alignment i and j arose through evolution") /  P("alignment i and j arose by chance")
### Exported!
Dayhoffn <- function(PAM1, BF, n){
  ## Verified by Darwin
  .validatePAMMatrix(PAM1)
  .validateBF(BF)
  pamn <- PAMn(PAM1, n)
  ans <- 10 * log10(sweep(pamn, MARGIN=2, BF, FUN="/"))
  return(ans)
}
### Dayhoff250 <- Dayhoffn(GONNET, GONNETBF, 250)



