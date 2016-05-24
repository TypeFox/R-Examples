#'Random Matrix Generator
#'@description Create a random matrix with specified dimensions, percentage fill, marginal distributions, and matrix type (abundance or binary) to test the behavior of null model randomization algorithms.
#'@param aBetaRow parameter shape1 for beta of row marginals
#'@param bBetaRow parameter shape2 for beta of row marginals
#'@param aBetaCol parameter shape1 for beta of column marginals
#'@param bBetaCol parameter shape2 for beta of column marginals
#'@param numRows number of rows in random matrix
#'@param numCols number of columns in random matrix
#'@param mFill proportion of matrix cells to be filled with non-zero elements
#'@param abun expected total abundance (from Poisson distribution) of individuals in the matrix. If set to the default value of 0, abundances are not shown but are replaced by occurrences (1s)
#'@param emptyRow specifies whether empty rows in random matrix will be retained (TRUE) or stripped out (FALSE)
#'@param emptyCol specified wither empty columns will be retained (TRUE) or stripped out (FALSE)
#'@details Understanding the behavior of null models with artificial data is an essential step before applying them to real data. This function creates stochastic community matrices in which each row is a species, each column is a site or island, and each entry is the occurrence (presence-absence) or abundance of a species in a site. For the analysis of niche overlap, the sites can be treated as unordered niche categories, and the abundances are the utilization values for each species.
#'
#'Row and column marginal distributions are described from a beta distribution, with user supplied coefficients. Marginal distributions are rescaled to one, and the conjoint probability of a species occurring in a site is determined with the outer product of the species (=row) and site(= column) marginals. This simple calculation assumes sites and species are independent, and excludes site x species interactions as well as species x species interactions.
#'
#'The user specifies the percent fill of the matrix, and this number of cells are randomly selected without replacement using the cojoint probabilities calculated from each marginal distribution. If the user has requested a presence-absence matrix (`abun = 0'), these cells are assigned a value of 1. If the user has reqested an abundance matrix ('abun > 0'), then the value of abundance specifies the summed abundance of all individuals in the matrix. The value of `abun` is used to set the lambda parameter for each occupied cell, and then a single draw from a Poisson distribution is used for the abundance in that cell. Small conjoint marginal probabilities can lead to empty rows or columns and the user can specify whether or not to retain empty rows and columns. The matrix rows and columns are sorted in descending order according to the marginal frequencies, and these are returned (matrix$rowMarg and matrix$colMarg) along with the matrix (matrix$m) in list form.
#'
#' `aBetaRow`, `bBetaRow`, `aBetaCol`, and `bBetaCol` specify the two shape parameters for the row and column marginals. The marginal values are created by a single random draw from these beta distributions, and then are rescaled so they sum to 1.0. Thus, the mean parameter value specified by the beta distribution does not matter in the calculation. Instead, it is the size of the variance that determines the amount of heterogeneity among row or column margins. Small values for the two shape parameters generate greater heterogeneity among the rows or columns marginals of the matrix.
#' 
#' Thus, a distribution with `aBetaRow=1000` and `bBetaRow=1000` will generate marginal probabilities that are virtually identical for the different species (=rows), whereas `aBetaRow=1` and `bBetaRow=1` will generate uniform probabilities. These default values applied to both rows and columns will generate a typical presence-absence matrix, with some common and some sparse species, and some species-rich and species-poor sites.
#' 
#' Setting `numRows`, `numCols`, and `mFill` allow the test matrix to be tailored to match the observed matrix. However, it may be necessary to increase `numRows` and `nCols` if the parameters often generate empty rows or columns.
#' 
#' If low values of `abun` are specified, some occupied cells may be set to 0 because of a random draw from the Poisson distribution for that matrix cell.
#' 
#' Once the test matrix is created, it can be used to explore any of the combinations of algorithm and matrix that are available in EcoSimR.
#' @examples \dontrun{
#' ## Create a null matrix similar to MacArthur's warblers
#' testMatrix <- ranMatGen(aBetaCol=1000,bBetaCol=1000,
#'                         aBetaRow=1,bBetaRow=1,
#'                         numRows=5,numCols=16,
#'                         mFill=0.75, abun=1000,
#'                         emptyRow=FALSE,emptyCol=TRUE)$m
#' 
#' ## Run the null model
#' testMod <- niche_null_model(testMatrix)
#' 
#' ## Summary and niche utilization plot
#' summary(testAnalysis)
#' plot(testMod,type="niche") 
#' 
#' 
#' ## Create a null matrix similar to West Indies Finches
#' testMatrix <- ranMatGen(aBetaCol=0.5,bBetaCol=0.5,
#'                         aBetaRow=0.5,bBetaRow=0.5,
#'                         numRows=30,numCols=30,
#'                         mFill=0.25,abun=0,emptyRow=FALSE,
#'                         emptyCol=FALSE)$m 
#' 
#' ## Run the null model
#' testMod <- cooc_null_model(testMatrix$m, 
#'                            algo="simFast",
#'                            burnin=10000,n.reps=1000)
#' 
#' ## Summary and matrix, burn-in plots
#' summary(testMod)
#' plot(testMod,type="cooc")
#' plot(testMod,type="burnin") 
#'}
#'
#'@export


# function ranMatGen
#
# NJG
# 2 August 2014
# Creates a random matrix (binary or integer)
# The user specifies the matrix dimensions and fill
# and gives parameters to specify beta distributions
# for row and column sums.
#
# The function returns a random matrix,
# and the marginal weights from which it was constructed

ranMatGen <- function (aBetaRow = 1, bBetaRow = 1,
                   aBetaCol = 1, bBetaCol = 1,
                   numRows = 20, numCols = 5,
                   mFill = 0.5, abun = 0,
                   emptyRow=FALSE, emptyCol=FALSE)
{

  
# set the marginal distributions, sort them, and rescale
rowMarg <- sort(rbeta(n=numRows,shape1=aBetaRow,shape2=bBetaRow),decreasing=TRUE)
rowMarg <- rowMarg/(sum(rowMarg))

colMarg <- sort(rbeta(n=numCols,shape1=aBetaCol,shape2=bBetaCol),decreasing=TRUE)
colMarg <- colMarg/(sum(colMarg))

# create the heat map matrix using independent marginals
heatMatrix <- outer(rowMarg,colMarg)

# draw specified number of cells
drawCells <- sample(x=heatMatrix,size=round(numRows*numCols*mFill),replace=FALSE,prob=heatMatrix)


# create array to hold coordinates of chosen cells
matCoor <- matrix(nrow=length(drawCells),ncol=2)

# determine coordinates of chosen cells
for (i in 1:length(drawCells)) {
matCoor[i,c(1,2)] <- which(heatMatrix==drawCells[i],arr.ind=TRUE)[c(1,2)]
}

# create empty random matrix
ranMat <- matrix(0,nrow=numRows,ncol=numCols)

# fill in for binary presence absence matrix or abundance

for (i in 1:nrow(matCoor)) {

  # if abun==0, fill in with occurrences
  if(abun==0)ranMat[matCoor[i,1],matCoor[i,2]] <- 1 else

    # if abun > 0 use poisson function for abundances  
    ranMat[matCoor[i,1],matCoor[i,2]] <- rpois(1,lambda=abun*(drawCells[i]/sum(drawCells)))
}


# optionally remove empty rows and/or columns
emptyRows <- which(rowSums(ranMat)==0)
emptyCols <- which(colSums(ranMat)==0)

if (!emptyRow & length(emptyRows)>0) ranMat <- ranMat[-(emptyRows),] 

if (!emptyCol & length(emptyCols)>0) ranMat <- ranMat[,-(emptyCols)] 

return(list(m=ranMat,rowMarg=rowMarg,colMarg=colMarg))
}
