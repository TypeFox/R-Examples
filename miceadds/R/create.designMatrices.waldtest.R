

#################################################
# create design matrices for waldtest
create.designMatrices.waldtest <- function( pars , k ){
        NP <- length(pars)
        Cdes <- matrix( 0 , nrow=k , ncol=NP)
        colnames(Cdes) <- pars
        rdes <- rep(0,k)
        res <- list( Cdes = Cdes , rdes=rdes )
        return(res)
                }
######################################################