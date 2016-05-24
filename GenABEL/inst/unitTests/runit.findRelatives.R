### --- Test setup ---

if(FALSE) {
    ## Not really needed, but can be handy when writing tests
    library(RUnit)
    library(GenABEL)
    library(DatABEL)
}

### do not run
#stop("SKIP THIS TEST")
###

### ---- common functions and data -----

#source(paste("../inst/unitTests/shared_functions.R"))
source(paste(path,"/shared_functions.R",sep=""))

### --- Test functions ---

test.findRelatives <- function()
{
    nloci <- 2000
    set.seed(1)
    q <- runif(nloci,min=0.05,max=0.95)
#
# g1---g2
#    |
#    g3----g4
#      __|__
#     |  |  |
#    g5 g6 g7---g8
#            _|_
#           |   |
#          g9  g10---g11
#                __|_
#               |   /\
#             g12 g13 g14
#
    nids <- 14
    sex <- c(1,0,1,0,0,0,1,0,0,1,0,0,0,0)
    names(sex) <- 1:14
    gt <- matrix(ncol=nloci,nrow=nids)
    gt[1,] <- rbinom(nloci,2,q)
    gt[2,] <- rbinom(nloci,2,q)
    gt[4,] <- rbinom(nloci,2,q)
    gt[8,] <- rbinom(nloci,2,q)
    gt[11,] <- rbinom(nloci,2,q)
    gt[3,] <- generateOffspring(gt[1,],gt[2,],q=q)
    gt[5,] <- generateOffspring(gt[3,],gt[4,],q=q)
    gt[6,] <- generateOffspring(gt[3,],gt[4,],q=q)
    gt[7,] <- generateOffspring(gt[3,],gt[4,],q=q)
    gt[9,] <- generateOffspring(gt[7,],gt[8,],q=q)
    gt[10,] <- generateOffspring(gt[7,],gt[8,],q=q)
    gt[12,] <- generateOffspring(gt[10,],gt[11,],q=q)
    gt[13,] <- generateOffspring(gt[10,],gt[11,],q=q)
    gt[14,] <- gt[13,]
    aa<-findRelatives(gt,q=q,nmei=c(1:2))
    checkIdentical(aa$compressedGuess[1,3],"1")
    checkIdentical(aa$compressedGuess[2,3],"1")
    checkIdentical(aa$compressedGuess[1,5],"2")
    checkIdentical(aa$compressedGuess[1,6],"2")
    checkIdentical(aa$compressedGuess[2,5],"2")
    checkIdentical(aa$compressedGuess[2,6],"2")
    checkIdentical(aa$compressedGuess[5,6],"2+2")
    checkIdentical(aa$compressedGuess[9,10],"2+2")
    checkIdentical(aa$compressedGuess[13,14],"0")
	aaPed <- reconstructNPs(aa$guess,sex)
    realPed <- matrix(c(rep("1",14),1:14,
                    0,0,1,0,3,3,3,0,7,7,0,10,10,10,
                    0,0,2,0,4,4,4,0,8,8,0,11,11,11),ncol=4)
    rownames(realPed) <- c(1:nids)
    colnames(realPed) <- c("fid","id","father","mother")
	checkIdentical(aaPed[,colnames(realPed)],realPed)
}