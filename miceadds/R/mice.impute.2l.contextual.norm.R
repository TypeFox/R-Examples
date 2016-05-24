mice.impute.2l.contextual.norm <- function (y, ry, x, type , ridge = 10^(-5) , imputationWeights = NULL , 
                interactions=NULL , quadratics = NULL ,  ...){
   vname <- get("vname", pos = parent.frame()) # get variable name            
   newstate <- get( "newstate" , pos = parent.frame() )  
   # data preparation
   xcov <- .a2l.contextual.auxiliary( y = y  , ry=ry , x=x , type=type , ...)     
   #------
   # norm imputation at level 2
   ximp <- mice.impute.weighted.norm( y= y , ry=ry, x = xcov , ridge = ridge , imputationWeights = imputationWeights , 
                  interactions= interactions , quadratics = quadratics ,  ... ) 
   return(ximp)
    }

	
	

#......................................
.a2l.contextual.auxiliary <- function( y , ry , x , type , ...){
   # extract cluster index
   clusterx <- x[,type == -2 ]
   #x1 <- cbind(1, as.matrix(x[,type %in% c(1,2) ]))   
   x1 <-  as.matrix(x[,type %in% c(1,2) ])
   if ( sum( type==2)  > 0 ){   
        z <-  as.matrix(x[,type == 2 ]) 
        # calculate aggregated values
        a1 <- stats::aggregate( z , list( clusterx ) , mean , na.rm=F)
        colnames(a1)[-1] <- paste( "M." , colnames(z) , sep="")
                            } 
   # calculate aggregated value for y
   a21 <- stats::aggregate( y , list( clusterx ) , sum , na.rm=F)
   a22 <- stats::aggregate( 1+0*y , list( clusterx ) , sum , na.rm=F)
   ic <- match( clusterx , a21[,1] )
   y2 <- ( a21[ ic , 2] - y ) / ( a22[ ic , 2 ] - 1 )
   y2[ is.na(y2) ] <- mean(y2,na.rm=T)                   
   if ( sum( type==2)  > 0 ){   
      xcov <- as.matrix( cbind(  x1 , a1[ ic , -1 ] , y2 ) )
                            } else {
      xcov <- as.matrix( cbind(  x1 ,  y2 ) )
                               }
   vname <- get("vname", pos = parent.frame()) 
   colnames(xcov)[ ncol(xcov) ] <- paste("M1." , vname , sep="")
   return(xcov)     
        }
###########################################################################################




#-----------------------------------------------------------------------------
# function for inclusion of group level predictors
.include.2l.predictors <- function( y, x , ry , type , ... ){
         X <- as.matrix( x[ , type %in% c(1,2)] )
         X <- cbind( 1 , X )        
        # group level predictors
        if ( sum( type == -2 ) > 0 ){
            cluster <- x[ , type == -2 ]
          if ( sum( type == 2 ) > 0 ){          
            x1a <- as.matrix( cbind( x[ , type== - 2 ] , x[ , type== 2 ]  ) )
            colnames(x1a) <- c( colnames(x)[ type == -2 ] , colnames(x)[ type == 2 ] )
            gm0 <- mice.impute.2l.groupmean(y = y , ry = FALSE * ry , x = x1a, 
                        type = c(-2,rep(1,ncol(x1a)-1 ) ) , grmeanwarning=FALSE , ... )
            gm0 <- as.matrix(gm0)
            colnames(gm0) <- paste("M." , colnames(x1a)[-1] , sep="" )
            X <- as.matrix(cbind( X , gm0 ))
                        }
                     }   else { cluster <- NULL }
          res <- list( "cluster" = cluster , "X" = X )
          return(res)
          }
#-----------------------------------------------------------------------------
          