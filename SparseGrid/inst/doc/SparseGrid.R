### R code from vignette source 'SparseGrid.Rnw'

###################################################
### code chunk number 1: setSweaveOptions
###################################################
# have an (invisible) initialization noweb chunk
# to remove the default continuation prompt '>'
options(continue = " ")
options(width = 60)

# eliminate margin space above plots
options(SweaveHooks=list(fig=function()
    par(mar=c(5.1, 4.1, 1.1, 2.1))))


###################################################
### code chunk number 2: NumberOfGridPoints
###################################################
# load library
library('SparseGrid')

# generate number of grid points
num.dims <- 9
accuracy <- 5
tab1 <- sapply( 1:num.dims, 
    function(d) { 
        x <- createSparseGrid( 'GQU', 
                                dimension=d, 
                                k=accuracy )
        return( c( length( x$weights ), accuracy^d ) )
    } )

# write header for LaTeX table
cat( "\\begin{tabular}[", 1 + ncol(tab1), "]{l|", paste(rep("r", ncol(tab1)), collapse=''), "}
\\hline \\hline \n", sep='')

cat( "$D$    & ", paste( 1:ncol(tab1), sep='', collapse=" & " ), "\\\\\\hline \n")
cat( paste( c("SGI", tab1[1,]), collapse=' & ' ), "\\\\ \n")
cat( paste( c("$k^D$", tab1[2,]), collapse=' & ' ), "\\\\ \n")

# write footer
cat("\\hline\\hline \n\\end{tabular} \n")


###################################################
### code chunk number 3: InstallSparseGrid (eval = FALSE)
###################################################
## install.packages('SparseGrid')


###################################################
### code chunk number 4: LoadSparseGrid (eval = FALSE)
###################################################
## library('SparseGrid')


###################################################
### code chunk number 5: HelpSparseGrid (eval = FALSE)
###################################################
## ?SparseGrid


###################################################
### code chunk number 6: CiteSparseGrid
###################################################
citation('SparseGrid')


###################################################
### code chunk number 7: HelpFunctions (eval = FALSE)
###################################################
## ?createSparseGrid
## ?createProductRuleGrid
## ?createMonteCarloGrid
## ?createIntegrationGrid


###################################################
### code chunk number 8: loadLibrary
###################################################
library('SparseGrid')


###################################################
### code chunk number 9: createSparseGrid
###################################################
dimension <- 10
k         <- 2
sgrid     <- createSparseGrid( type='KPU', dimension=dimension, k=k )


###################################################
### code chunk number 10: numberSparseGridNodes
###################################################
length( sgrid$weights ) 


###################################################
### code chunk number 11: defineFunction
###################################################
g <- function( x, mu=0, sigma=1 ) {
    return( prod( exp(-.5*((x-mu)/sigma)^2)/sqrt(2*pi*sigma^2) ) )
}


###################################################
### code chunk number 12: defineApproximationFunction
###################################################
approximate.integral <- function( func, sgrid, ... ) {
    gx <- apply( sgrid$nodes, 1, function(x) { func(x, ...) } )
    return( sum( gx * sgrid$weights )  )
}


###################################################
### code chunk number 13: approximateSGIntegral
###################################################
sigma <- 2
approximate.integral( g, sgrid, mu=0, sigma=sigma )


###################################################
### code chunk number 14: defineTrueValue
###################################################
trueval <-
    ( pnorm( 1, mean=0, sd=2 ) - pnorm( 0, mean=0, sd=sigma ) )^dimension


###################################################
### code chunk number 15: createMCGrid
###################################################
num.sim <- 1000
set.seed( 3141 )
mcgrid <- createMonteCarloGrid( 
                runif, dimension=dimension, num.sim=num.sim )


###################################################
### code chunk number 16: approximateMCIntegral
###################################################
approximate.integral( g, mcgrid, mu=0, sigma=sigma )


###################################################
### code chunk number 17: defineParameters
###################################################
set.seed( 3141 )
dimension   <- 10   # dimension of integral
maxk        <- 4    # max. accuracy level (pol. exactness wil be 2k-1)


###################################################
### code chunk number 18: initResultMatrix
###################################################
# create matrix to hold results
res <- matrix( NA, nrow=maxk-1, ncol=5 )
colnames( res ) <- c("D", "k", "nodes", "SG error", "MC error")
rownames( res ) <- rep( "", maxk-1 )


###################################################
### code chunk number 19: performComparison
###################################################
# loop over different accuracy levels
for ( k in 2:maxk ) {

    # sparse grid integration
    sgrid   <- createSparseGrid('KPU', dimension, k)
    SGappr  <- approximate.integral( g, sgrid, mu=0, sigma=sigma )
    SGerror <- sqrt((SGappr - trueval)^2) / trueval
    
    # Monte Carlo integration with the same number of nodes
    # 1000 simulation repetitions
    num.nodes <- length( sgrid$weights )
    MCappr    <- rep(0, 1000)
    for (r in 1:1000) {
        mcgrid      <- createMonteCarloGrid( 
                            runif, 
                            dimension=dimension, 
                            num.sim=num.nodes )
        MCappr[ r ] <- approximate.integral( 
                            g, mcgrid, mu=0, sigma=sigma )
    }
    MCerror = sqrt(mean((MCappr-trueval)^2)) / trueval
    
    # save results in row of matrix res
    res[k-1,] <- c(dimension, k, num.nodes, SGerror, MCerror)
}


###################################################
### code chunk number 20: showResults
###################################################
res


