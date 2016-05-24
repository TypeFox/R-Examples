
#--------------------------------------------------------------------------
# calculation of PRMSE for Subscores according to Haberman (2007)
prmse.subscores <- function( data.X , data.Z){
    ind <- which( rowSums( is.na( data.Z ) ) > 0 )
#    require(psy)
    if (length(ind) > 0){
        data.X <- data.X[ - ind , ]
        data.Z <- data.Z[ - ind , ]
                    }    
    aX <- .cronbach.alpha( data.X )
    aZ <- .cronbach.alpha( data.Z )
    res <- list( "N" = aX$sample.size ,  "nX" = aX$number.of.items )
    score.X <- rowSums( data.X )
    score.Z <- rowSums( data.Z )
    res$M.X <- mean( score.X )
    res$Var.X <- stats::var( score.X )
    res$SD.X <- sqrt( res$Var.X )
    res$alpha.X <- aX$alpha
    res$Var.TX <- res$Var.X * res$alpha.X
    res$Var.EX <- res$Var.X * ( 1 - res$alpha.X )
    res$nZ <- aZ$number.of.items
    res$M.Z <- mean( score.Z )
    res$Var.Z <- stats::var( score.Z )
    res$SD.Z <- sqrt( res$Var.Z )
    res$alpha.Z <- aZ$alpha
    res$Var.TZ <- res$Var.Z * res$alpha.Z
    res$Var.EZ <- res$Var.Z * ( 1 - res$alpha.Z )
    res$cor.X_Z <- stats::cor( score.X , score.Z )
    res$cor.X_Y <- stats::cor( score.X , score.Z - score.X )
    res$cor.TX_TY <- res$cor.X_Y / sqrt( res$alpha.X ) / 
                sqrt( .cronbach.alpha( data.Z[ , setdiff( colnames(data.Z) , colnames(data.X) ) ] )$alpha )
    res$Var.TX <- res$Var.X - res$Var.EX
    res$Var.TZ <- res$Var.Z - res$Var.EZ
    
    res$cor.TX_TZ <- res$cor.X_Z / sqrt( res$alpha.X * res$alpha.Z ) - 
                        res$Var.EX / sqrt( res$Var.TX * res$Var.Z )
    res$cor.TX_Z <- res$cor.TX_TZ * sqrt( res$alpha.Z )
    # RMSE basierend auf Subscores (Kelleyformel)
    res$rmse.X <- sqrt( res$Var.TX * ( 1 - res$alpha.X ) )
    # RMSE basierend auf Total Scores
    res$rmse.Z <- sqrt( res$Var.TX * ( 1 - res$cor.TX_Z^2 ) )
    # calculation of regression coefficients
    regr <- matrix( 0 , nrow=3 , ncol=3 )
    colnames(regr) <- c("Int" , "beta.X" , "beta.Z" )
    rownames(regr) <- c("Mod.X" , "Mod.Z" , "Mod.XZ" )
    regr[ "Mod.X" , "beta.X" ] <- res$alpha.X
    regr[ "Mod.X" , "Int" ] <- res$M.X - res$alpha.X * res$M.X
    regr[ "Mod.Z" , "beta.Z" ] <- res$cor.TX_Z * sqrt(res$Var.TX) / sqrt( res$Var.Z )
    regr[ "Mod.Z" , "Int" ] <- res$M.X - regr[ "Mod.Z" , "beta.Z" ] * res$M.Z
    regr[ "Mod.XZ" , "beta.X" ] <- ( sqrt( res$Var.TX ) * ( sqrt( res$alpha.X ) - res$cor.TX_Z * res$cor.X_Z ) ) /
                                        ( res$SD.X * ( 1 - res$cor.X_Z^2 ) )
    regr[ "Mod.XZ" , "beta.Z" ] <- ( sqrt( res$Var.TX ) * ( res$cor.TX_Z - sqrt(res$alpha.X) * res$cor.X_Z ) ) /
                                        ( res$SD.Z * ( 1 - res$cor.X_Z^2 ) )
    regr[ "Mod.XZ" , "Int" ] <- res$M.X - regr[ "Mod.X" , "beta.X" ] * res$M.X - regr[ "Mod.Z" , "beta.Z" ] * res$M.Z
    # calculation of RMSE of the regression on both subscores and total score
    pcor <- (( res$cor.TX_Z - sqrt( res$alpha.X) * res$cor.X_Z ) /
                    sqrt( 1 - res$alpha.X ) / sqrt( 1 - res$cor.X_Z^2 ) ) 
    res$rmse.XZ <- res$rmse.X * sqrt( 1 - pcor^2)
    res$prmse.X <- res$alpha.X
    res$prmse.Z <- res$cor.TX_Z^2
    res$prmse.XZ <- 1 - ( 1 - res$alpha.X )*( 1 -pcor^2 )
    list( "res" = res , "regr" = regr )
        }
#..................................................................................



#--------------------------------------------------------------------------------
# calculation of PRMSE for a number of subscales
prmse.subscores.scales <- function( data , subscale ){ 
        # data ... original data frame
        # scales ... classification of items in data into scales
        scales <- sort( unique( subscale ) )
        dfr <- NULL
        for ( ss in scales ){
            # ss <- scales[1]
            res1 <- prmse.subscores(  data.X = data[ , subscale == ss  ] , data.Z = data )
            dfr <- cbind( dfr , unlist( res1$res ) )
                }
        colnames(dfr) <- scales
        return(dfr)
            }
#--------------------------------------------------------------------------------




# aux function for Cronbach's Alpha
    .cronbach.alpha <- function( data ){ 
        # covariance
        c1 <- stats::cov( data , use = "pairwise.complete.obs" )
        # mean covariance
        c1a <- c1 ; diag(c1a) <- 0
        I <- ncol(data)
        mc <- sum(c1a) / ( I^2 - I )
        # mean variance
        mv <- mean( diag(c1) )
        alpha <- I * mc / ( mv + (I-1) *mc )
        mean.tot <- mean( rowSums(data) )
        var.tot <- var( rowSums( data ) )
         res <- list( "I" = I , "alpha" = alpha , 
				"mean.tot" = mean.tot ,  "var.tot" = var.tot ,
				"sample.size" = nrow(data) , 
				"number.of.items" = I )
        return(res)
            }

