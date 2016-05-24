 


#----------------------------------------------------------------------------------
# person fit with logistic regression
pf.logist.regression <- function( data , itemdiff , perc = seq(5,100,5) ){
    P <- nrow( data )
    fits <- rep(0,P)
    convdis <- .display.conv( maxiter=P , perc = perc , P=P )
    for (pp in 1:P){
        # pp <- 40
        y1 <-  t(data[pp,])[,1]
        mod.pp <- stats::glm( y1 ~ itemdiff , family = "binomial")
        fits[pp] <- stats::coef(mod.pp)[2]
        ind.pp <- which( convdis[,1] %in% pp )
        if (length(ind.pp) > 0){ cat("\n" ,paste(convdis[ind.pp,2], "%" , sep="") ) ; 
		utils::flush.console()} 
            }
        cat("\n")
    return( fits)
    }
#...........................................................................
# display convergence progress
.display.conv <- function( maxiter , perc =  seq( 5 , 100 , 5 ) , P=P ){
    h1 <- 100 * seq( 1 , P ) / P 
    h2 <- sapply( perc , FUN = function(perc.pp){
            # perc.pp <- perc[1]
            which.min( abs( h1 - perc.pp ) )[1]
            } )
    cbind( "iter" = h2 , "perc" = perc )
        }
#-----------------------------------------------------------------------------------------
# Reise fit
##NS export(pf.reisefit)
# pf.reisefit <- function( data , itemdiff ){
#    dat1 <- data
#    a1 <- data.frame( "y" = matrix(  t(dat1) , ncol= 1 , byrow=T) , 
#                    "person" = rep(  seq( 1 , nrow(dat1) ) , each=I ) ,
#                    "itemdiff" = rep( itemdiff , nrow(dat1) )
#                        )
#    mod <- lmer( y ~ itemdiff + ( 1 + itemdiff | person ) , data = a1 , family = "binomial" )
#    l1 <- list( "reise.fit" = fixef(mod)[2] + ranef(mod)$person[,2]   , "summary" = summary(mod) , "model" = mod )
#    return(l1)
#    }
#-----------------------------------------------------------------------------------------
# calculate l0 and lz likelihood based statistic
##NS export(pf.l0)
pf.l0 <- function( data , wle , itemdiff ){ 
    dat1 <- data 
	dat1.resp <- 1 - is.na(dat1)
	dat1[ is.na(dat1) ] <- 0		
    # I <- ncol(dat1)
	Iresp <- rowSums( dat1.resp )
    pmatrix <- stats::plogis( outer( wle , itemdiff , "-" ) )
    # l0
    l0 <- rowSums( ( dat1 * log( pmatrix ) + ( 1 - dat1) * log( 1 - pmatrix ) )*dat1.resp / 
				Iresp )
    # standardized l0
    E.l0 <- rowSums( dat1.resp*( pmatrix * log( pmatrix ) + ( 1 - pmatrix) * log( 1 - pmatrix ) )  )
    Var.l0 <- rowSums( dat1.resp*( pmatrix * ( 1 - pmatrix ) * ( stats::qlogis( pmatrix ) )^2 )  )
    # lz <- ( I*l0 - E.l0 ) / sqrt( Var.l0 )
	lz <- ( Iresp*l0 - E.l0 ) / sqrt( Var.l0 )
    list( "l0" = l0 , "lz" = lz )
    }
#-----------------------------------------------------------------------------------
# Infit and Outfit (Rasch model)
##NS export(pf.outfit.infit)
pf.outfit.infit <- function( data , wle , itemdiff ){ 
    dat1 <- data 
    # I <- ncol(dat1)
    pmatrix <- stats::plogis( outer( wle , itemdiff , "-" ) )
	
	dat1.resp <- 1 - is.na(dat1)
	dat1[ is.na(dat1) ] <- 0
	
	Iresp <- rowSums(dat1.resp)
	
    # Outfit (unweighted fit)
    U <- rowSums(  ( ( dat1 - pmatrix )^2 / Iresp / pmatrix / ( 1 - pmatrix ) ) * dat1.resp )
    # Infit (weighted fit)
    W <- rowSums(  dat1.resp * ( dat1 - pmatrix )^2  ) / 
					rowSums(  dat1.resp * pmatrix * ( 1 - pmatrix ) )
    list( "Outfit" = U , "Infit" = W )
        }
########################################################################		
		
		
#.......................................
# U3 statistic                      
##NS export(pf.U3) 
pf.U3 <- function( data , pval ){
    dat1 <- data
    # data processing
    data.proc <- .guttman.data.frame( data = dat1 , pval = pval )
    data.proc2 <- .inverse.guttman.data.frame( data = dat1 , pval = pval )
    P <- nrow(dat1)
    I <- ncol(dat1)
    # calculate A_i
    pmatrix <- outer( rep(1,P) , data.proc$pval ) 
    pmatrix <- log( pmatrix / ( 1 - pmatrix ) )
    T1 <- rowSums( data.proc$data.guttman * pmatrix)
    T2 <- rowSums( data.proc$data * pmatrix) 
    T3 <- rowSums( data.proc2$data.inverse.guttman  * pmatrix )
    GI <- ( T1 - T2 ) / ( T1 - T3 )
    GI 
        }
#.......................................
# Caution Index            
##NS export(pf.caution)             
pf.caution <- function( data , pval ){
    dat1 <- data
    # data processing
    data.proc <- .guttman.data.frame( data = dat1 , pval = pval )
    data.proc2 <- .inverse.guttman.data.frame( data = dat1 , pval = pval )
    P <- nrow(dat1)
    I <- ncol(dat1)
    # calculate caution index
    pmatrix <- outer( rep(1,P) , data.proc$pval ) 
    T1 <- rowSums( data.proc$data.guttman * pmatrix)
    T2 <- rowSums( data.proc$data * pmatrix) 
    T3 <- rowSums( data.proc2$data.inverse.guttman  * pmatrix )
    GI <- ( T1 - T2 ) / ( T1 - T3 )
    GI 
        }
#----------------------------------------------------------------
# Correlation Index  
##NS export(pf.rbis)                                            
pf.rpbis <- function( data , pval ){ 
    I <- ncol(data)
    P <- nrow(data)
    cor.rowwise( data , outer( rep(1,P) , pval ) )
    }  
# Correlation Index 
##NS export(pf.rpbis.itemdiff)                                             
pf.rpbis.itemdiff <- function( data , itemdiff ){ 
    I <- ncol(data)
    P <- nrow(data)
    cor.rowwise( data , outer( rep(1,P) , itemdiff ) )
    }  
#----------------------------------------------------------------
# Disagreement index
##NS export(pf.dependability)
pf.dependability <- function( data , pval ){ 
    dat1 <- data
    # data processing
    data.proc <- .guttman.data.frame( data = dat1 , pval = pval )
    P <- nrow(dat1)
    I <- ncol(dat1)
    # calculate A_i
    pmatrix <- outer( rep(1,P) , data.proc$pval ) 
    A.max <- rowSums( data.proc$data.guttman * pmatrix) / I
    A.i <- rowSums( data.proc$data * pmatrix) / I
    # disagreement index
    DI <- A.max - A.i
    DI 
    }
#...................................................................
# ECI indexes
##NS export(pf.ECI)
pf.eci <- function( data , pval , wle , itemdiff){ 
    dat1 <- data
    # data processing
    data.proc <- .guttman.data.frame( data = dat1 , pval = pval )
    itemdiff2 <- itemdiff[ data.proc$index ]
    P <- nrow(dat1)
    I <- ncol(dat1)
    # calculate matrices
    P.matr <- stats::plogis( outer( wle , itemdiff2 , "-")  )
    X.matr <- data.proc$data 
    n.matr <- outer( rep( 1 , P) , data.proc$pval )
    G.matr <- outer( rep(1,P) , colMeans( P.matr )  )
    # indices
    ECI1 <- 1 - cov.rowwise( X.matr , n.matr ) / cov.rowwise( P.matr , n.matr )
    ECI2 <- 1 - cov.rowwise( X.matr , G.matr ) / cov.rowwise( P.matr , G.matr )
    ECI3 <- 1 - cor.rowwise( X.matr , G.matr ) / cor.rowwise( P.matr , G.matr )
    ECI4 <- 1 - cov.rowwise( X.matr , P.matr ) / cov.rowwise( G.matr , P.matr )
    ECI5 <- 1 - cor.rowwise( X.matr , P.matr ) / cor.rowwise( G.matr , P.matr )
    ECI6 <- 1 - cov.rowwise( X.matr , P.matr ) / cov.rowwise( P.matr , P.matr )
    ECI <- data.frame( ECI1 , ECI2 , ECI3 , ECI4 , ECI5 , ECI6 )
    colnames(ECI) <- paste("ECI" , 1:6 , sep="")
    list( "ECI" = ECI , "ECI1" = ECI1 ,"ECI2" = ECI2 , "ECI3" = ECI3 , "ECI4" = ECI4 , "ECI5" = ECI5 , "ECI6" = ECI6 )
    }
#..................................................................
# calculate covariance rowwise for matrices
cov.rowwise <- function( m1 , m2 ){ 
    I <- ncol(m1)
    ( rowSums( m1 * m2 ) - I * rowMeans(m1) * rowMeans(m2) ) / ( I - 1 )
    }
#...............................................................
# rowwise SD's
sd.rowwise <- function(m1){
    sqrt( cov.rowwise( m1 , m1 ) )
    }
#...............................................................
# rowwise correlations
cor.rowwise <- function(m1,m2){
    cov.rowwise(m1,m2) / sd.rowwise( m1 ) / sd.rowwise( m2 )
    }
#................................................................
# create Guttman data frame
.guttman.data.frame <- function( data , pval ){
    data.sort <- .sort.items( data = data , pval = pval )
    data1 <- data.sort$data
    pval1 <- data.sort$pval
    I <- ncol(data1)
    # create Guttman pattern data matrix data.guttman
    score <- rowSums( data1 )
    data.guttman <- 0*data1
    for (ii in 1:I){ 
			data.guttman[,ii] <- 1*( score >= ii ) 
					}
	res <- list( data = data1 , pval = pval1 , data.guttman = data.guttman , index = data.sort$index )					
    return( res )
        }
#.........................................................................
# create inverse Guttman data frame
.inverse.guttman.data.frame <- function( data , pval ){
    data.sort <- .sort.items( data = data , pval = pval )
    data1 <- data.sort$data
    pval1 <- data.sort$pval
    I <- ncol(data1)
    # create Guttman pattern data matrix data.guttman
    score <- rowSums( data1 )
    data.guttman <- 0*data1
    for (ii in 1:I){ data.guttman[,I + 1 - ii] <- 1*( score >= ii ) }
    return( list( data = data1 , pval = pval1 , data.inverse.guttman = data.guttman ) )
        }
#.........................................................................
#........................................................................
# auxiliary function for sorting items according to p values
.sort.items <- function( data , pval ){
    l1 <- sort( pval , decreasing=TRUE , index.return=TRUE	)
    dat1 <- data[ , l1$ix ]
    pval1 <- pval[ l1$ix ]    
    return( list( data = dat1 , pval = pval1 , index = l1$ix) )
        }
#........................................................................


