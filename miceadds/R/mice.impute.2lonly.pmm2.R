
mice.impute.2lonly.pmm2 <- function (y, ry, x, type , ...){
    imp <- .imputation.level2a( y = y , ry = ry , x = x , type = type , 
                               imputationMethod = "pmm" , ... )
}

mice.impute.2lonly.norm2 <- function (y, ry, x, type , ...){
    imp <- .imputation.level2a( y = y , ry = ry , x = x , type = type , 
                               imputationMethod = "norm" , ... )
    return(imp)
}

#******************************************
# imputation function at level 2
# can be done with norm and pmm
.imputation.level2a <- function( y , ry , x , type , imputationMethod , ... ){
    if ( sum(type==-2 ) != 1 ){
        stop( "No class variable")
    }
    # extract cluster index
    clusterx <- x[,type == -2 ]
    x <- cbind(1, as.matrix(x[,type %in% c(1,2)]))      # calculate aggregated values	
    # change ARb 2013-02-12
    # a1 <- aggregate( (cbind(x,y)) , list( clusterx ) , mean , na.rm=F)
    a2 <- rowsum( cbind(x,y) , clusterx , na.rm=FALSE)
	#~~~~~
	# change ARb 2014-02-18
	clusterx0 <- as.numeric(paste0(rownames(a2)))
    a2 <- a2 / rowsum( 1+0*y , clusterx , na.rm=FALSE )[,1] 
    a1 <- cbind( clusterx0  , a2 )
	#~~~~~
    #*****			
    N1 <- ncol(a1)
    cly2 <- unique( clusterx[  ry ] )  # clusters without missings on y
    ry2 <- a1[,1] %in% cly2  
    x1 <- as.matrix(a1[, -c(1,N1)])
    # norm imputation at level 2
    if ( imputationMethod == "norm" ){ 
        ximp2 <- mice::mice.impute.norm( y= as.matrix(a1[,N1]), ry=ry2, x = x1[,-1] , ...) 
    }
    # pmm imputation at level 2
    if ( imputationMethod == "pmm" ){ 
        ximp2 <- mice::mice.impute.pmm( y= as.matrix(a1[,N1]), ry=ry2, x = x1[,-1] , ...) 
    }
    # data postprocessing
    cly2 <- a1[ ! ry2 , 1] 
    i1 <- match( clusterx, cly2 )
    ximp <- ( ximp2[i1] )[ ! ry ]
    return(ximp)	
}
