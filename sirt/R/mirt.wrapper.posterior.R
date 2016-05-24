
#############################################################
# calculation of posterior
mirt.wrapper.posterior <- function( mirt.obj , weights=NULL ){
	# adapt function for multiple groups
    # ****
	mobj <- mirt.obj
	# extract theta
	Theta <- mobj@Theta
	# extract theta distribution
	# pi.k <- mobj@bfactor$Prior[[1]]
	pi.k <- mobj@Prior[[1]]
	# extract full data
	fulldata <- mobj@Data$fulldata[[1]]
	#load items
	I <- ncol( mobj@Data$data )
	items <- vector('list', I)
	for(ii in 1:I){
		 items[[ii]] <- mirt::extract.item(mobj, ii)
				}
				
	# check whether prodlist exists
    prodlist <- attr(mobj@pars, "prodlist") 
	# Theta <- mobj@Theta 
    Theta1 <- Theta
	if ( ! is.null(prodlist) ){
				Theta1 <- mirt_prodterms(Theta, prodlist)	
							}				
	# item-wise probabilities for each Theta
	traces <- Probtrace_sirt(items, Theta1)	
	# log-Likelihood
	f.yi.qk <- exp( fulldata %*% t(log(traces))  )
	# compute individual posterior
	N <- nrow( fulldata )
	TP <- length(pi.k)
	piM <- matrix( pi.k , nrow=N , ncol=TP , byrow=TRUE )
	f.qk.yi <- f.yi.qk * piM 
	f.qk.yi <- f.qk.yi / matrix( rowSums( f.qk.yi ) , nrow=N , ncol=TP , byrow=FALSE )
	# maximum category
	maxK <- apply( mobj@Data$data , 2 , max , na.rm=TRUE)+1
	resp.ind <- 1- is.na(mobj@Data$data)
	resp <- mobj@Data$data
	resp[ resp.ind == 0 ] <- 0
	# calc counts
	group <- NULL	# only applies to single groups for now
	if (is.null(weights) ){ 
			pweights <- rep(1,N) } else {
			pweights <- weights 
					}
    # Theta is only used for calculating dimension size					
	n.ik <- mirt.wrapper.calc.counts( resp, theta=Theta , resp.ind=resp.ind , 
				group=group , maxK=max(maxK) , pweights=pweights , hwt=f.qk.yi )
    probs <- traces		
    probs <- array( probs , dim = c(TP,max(maxK),I) )	
	probs <- aperm( probs , c(3,2,1) )		
	# result list
	res <- list( "theta.k" = Theta , "pi.k" = pi.k ,
			"f.yi.qk" = f.yi.qk , "f.qk.yi" = f.qk.yi ,
			"n.ik"=n.ik , "probs" = probs ,
			"N"= N , "TP"=TP , "I" = I , "data" = mobj@Data$data ,
			"maxK" = maxK )
	class(res) <- "mirt"
	return(res)
	}
######################################################################
# auxiliary function
# trace function for all items
Probtrace_sirt <- function(items, Theta){
		 traces <- lapply(items, mirt::probtrace, Theta=Theta)
		 ret <- do.call(cbind, traces)
		 ret
	}
##########################################################################
# prodterms function from mirt package
# This function is not exported and hence redefined in sirt
mirt_prodterms <- function (theta0, prodlist) 
{
    products <- matrix(1, ncol = length(prodlist), nrow = nrow(theta0))
    for (i in 1L:length(prodlist)) {
        tmp <- prodlist[[i]]
        for (j in 1L:length(tmp)) products[, i] <- products[, 
            i] * theta0[, tmp[j]]
    }
    ret <- cbind(theta0, products)
    ret
}
