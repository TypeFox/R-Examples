
#############################################################
print.designMatrices <-
  function( X , ... ){
    x <- X
    BB <- x$flatB
    colnames(BB) <- paste("B_", colnames(BB), sep ="")
    out <- cbind( x$flatA, BB )
    
    NAs <- apply( x$flatA , 1 , function(fA) all(is.na(fA)) )
    out <- out[!NAs, ]
    
    print(out)
    invisible( out )
  }


rownames.design <- function(X){
  Y <- apply(X, 2, as.numeric )
  Y <- sapply(1:ncol(Y), function(vv) 
    paste( colnames(Y)[vv], add.lead(Y[,vv], ceiling(log( max(as.numeric(Y[,vv])), 10)) ), sep ="" )
  )
  
  rownames(X) <- apply(Y, 1, paste, collapse = "-")
  return(X)
} 

rownames.design2 <- function(X){
  Y <- apply(X, 2, as.numeric )
  Y <- sapply(1:ncol(Y), function(vv) 
    # paste( colnames(Y)[vv], add.lead(Y[,vv], ceiling(log( max(as.numeric(Y[,vv])), 10)) ), sep ="" )
    paste( colnames(Y)[vv], add.lead(Y[,vv], 1) , sep ="" )
  )
  
  rownames(X) <- apply(Y, 1, paste, collapse = "-")
  return(X)
} 


###########################################################
## function .A.matrix
.A.matrix <-
  function( resp, formulaA = ~ item + item*step, facets = NULL,  
            constraint = c("cases", "items") , progress=FALSE ,
            maxKi = NULL ){
    z0 <- Sys.time()			
    ### redefine facets matrix
    facets0 <- facets
    NF <- length(facets)
    facet.list <- as.list( 1:NF )
    names(facet.list) <- colnames(facets)	
    if (NF==0){ facet.list <- NULL }
    if (NF>0){
      for (ff in 1:NF){
        #		ff <- 2
        uff <- sort( unique( facets[,ff] ) )
        facets[,ff] <- match( facets[,ff] , uff )
        facet.list[[ff]] <- data.frame(
          "facet.label" = paste0( colnames(facets)[ff] , uff ) , 
          "facet.index" = paste0( colnames(facets)[ff] , seq(1,length(uff) ) ) )
      }
    }
# cat(" +++  v62" ) ; z1 <- Sys.time() ; print(z1-z0) ; z0 <- z1     						
    ### Basic Information and Initializations
    constraint <- match.arg(constraint)
    if ( is.null(maxKi) ){
      maxKi <- apply( resp , 2 , max , na.rm=TRUE )
    }
    maxK <- max( maxKi )
    nI <- ncol( resp )
# cat(" +++  v70" ) ; z1 <- Sys.time() ; print(z1-z0) ; z0 <- z1     						
    
    # stop processing if there are items with a maximum score of 0
    i11 <- names(maxKi)[ maxKi == 0 ]
    if ( length(i11) > 0 ){
      stop( cat( "Items with maximum score of 0:" , paste(i11 , collapse=" " ) ) )
    }
    
    tf <- terms( formulaA )	
    fvars <- as.vector( attr(tf,"variables"), mode = "character" )[-1]
    #cat("fvars 212") ; print(fvars)	
    otherFacets <- setdiff( fvars, c("item", "step") )
    contr.list <- as.list( rep( "contr.sum", length(fvars) ) )
    names( contr.list ) <- fvars
# cat(" +++  v80" ) ; z1 <- Sys.time() ; print(z1-z0) ; z0 <-z1
    #****
    # ARb: 2013-03-27
    # no contrasts for items
    nitems <- ncol(resp)
    #	contr.list[["item"]] <- diag(1,nitems)
    #******	
    ### prepare data-Object for model.matrix()
    # expand.list <- 
    #   as.vector( c( list( if( "item" %in% fvars ) factor(1:nI),
    #                       if( "step" %in% fvars ) factor(1:maxK) ),
    #                 if( length( otherFacets ) == 1){
    #                   list( factor( 1:max(facets[, otherFacets]) ) )
    #                   #    				list( factor( unique( facets[ , otherFacets] ) ) )
    #                 } else if( length( otherFacets ) > 1 ){
    #                   #                      apply( as.matrix( facets[, otherFacets] ), 2, 
    #                   #							function(ff){ as.factor(1:max(ff)) }
    #                   #									)
    #                   # Bug for equal numbers of levels within facets
    #                   # Correction 2013-09-03
    #                   sapply( otherFacets , FUN = function(ff){
    #                     fff <- facets[, ff]
    #                     as.factor(1:max(fff)) 
    #                   } , simplify=FALSE )
    #                 }                     
    #   ) )
    # TK: 2014-03-12
    # Consider long vector response matrix and "item" in facets input
    expand.list <- vector(mode = "list", length = 0)
    if( "item" %in% fvars )         expand.list <- c(expand.list, if("item" %in% names(facet.list)) list(as.factor(sort(unique(facets[,"item"])))) else list(factor(1:nI)) )
    if( "step" %in% fvars )         expand.list <- c(expand.list, if("step" %in% names(facet.list)) list(as.factor(sort(unique(facets[,"step"])))) else list(factor(1:maxK)) )
    if( length( otherFacets ) == 1) expand.list <- c(expand.list, list(factor(1:max(facets[, otherFacets]))) )     
    if( length( otherFacets ) > 1 ) expand.list <- c(expand.list, sapply( otherFacets , FUN = function(ff) as.factor(1:max(facets[, ff])) , simplify=FALSE ))
    
    names( expand.list ) <- fvars
    
# cat(" +++  v100" ) ; z1 <- Sys.time() ; print(z1-z0) ; z0 <- z1     							  	  
    #     expand.list <- expand.list[ !unlist( lapply(expand.list, is.null) ) ]
    
    for (vv in seq(1 , length(expand.list) ) ){
      expand.list[[vv]] <- paste( expand.list[[vv]] ) 
    }
# cat(" +++  v110" ) ; z1 <- Sys.time() ; print(z1-z0) ; z0 <- z1     					
    X <- rownames.design2( expand.grid(expand.list) )
    ### constraints and formulaA
    if( constraint == "cases" ) formulaA <- update.formula(formulaA, ~0+.)
    NX <- ncol(X)
    for (ff in 1:NX){
      uff <- length( unique(X[,ff] ) )
      if (uff==1){ cat(paste0("          - facet " ,
                              colnames(X)[ff] , " does only have one level!" ) , "\n") }
    }
# cat(" +++  v120" ) ; z1 <- Sys.time() ; print(z1-z0) ; z0 <- z1     						
    mm <- - stats::model.matrix(formulaA, X, contrasts = contr.list)
    #    mm <- - model.matrix(formulaA, X )
    if( constraint == "items" ) mm <- mm[,-1]
    
    ############################################################
    ###*** ARb 2013-03-28
    ### generate all interactions	
    xsi.constr <- .generate.interactions(X , facets , formulaA , mm )									
    ###############################################################
# cat(" +++  v130" ) ; z1 <- Sys.time() ; print(z1-z0) ; z0 <- z1     						  	  		
    ### Postprocessing
    # model.matrix _ case: step in fvars
    if( "step" %in% fvars ){
      if( ncol( attr(tf, "factors") ) == 1 ){
        return( warning("Can't proceed the estimation: 
                        Factor of order 1 other than step must be specified.") )
      } 
      if( all( attr(tf, "factors")["step",] != 1 ) ){
        return( warning("Can't proceed the estimation: 
                        Lower-order term is missing.") )
      } 
      
      A <- NULL
      
      stepgroups <- unique( gsub( "(^|-)+step([[:digit:]])*", "\\1step([[:digit:]])*", rownames(X) ) )
      X.out <- data.frame(as.matrix(X), stringsAsFactors = FALSE)
# cat(" +++  v150" ) ; z1 <- Sys.time() ; print(z1-z0) ; z0 <- z1     						  	  		
      if (progress){
        cat("        o Create design matrix A\n")
        ip <- length(stepgroups)
        VP <- min( ip , 10 )
        cat(paste0("          |",paste0( rep("*" , VP) , collapse="") , "|\n"))
        cat("          |") ; flush.console()
        if (VP<10){ disp_progress <- 1:ip } else {
          disp_progress <- 100* ( 1:ip ) / (ip+1)
          disp_progress <- sapply( seq(5,95,10) , FUN = function(pp){ # pp <- 5
            which.min( abs( disp_progress - pp ) )[1] }
          )
        }						
      }
#cat(" +++  v155" ) ; z1 <- Sys.time() ; print(z1-z0) ; z0 <- z1     						  
      ii <- 0 ; vv <- 1
##!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	  
##!!!! This loop is time consuming!	
## changed ARb 2015-03-28 

	NRX <- length( rownames(X) )
	rownames_X_matr <- strsplit( rownames(X) , split="-")
	rownames_X_matr <- matrix( unlist( rownames_X_matr ) , nrow=NRX , byrow=TRUE )
	step_col <- 0
	for (ff in 1:( ncol( rownames_X_matr ) ) ){
		if ( length( grep( "step1" , rownames_X_matr[,ff] ) ) > 0 ){
			step_col <- ff
					} 
				}
	rownames_X_matr2 <- rownames_X_matr[ , - step_col , drop=FALSE ]
	N2 <- ncol( rownames_X_matr2 )
	rownames_X_matr2_collapse <- rownames_X_matr2[,1]
	if (N2>1){
		for (nn in 2:N2){ 
		rownames_X_matr2_collapse <- paste0( rownames_X_matr2_collapse , "-" ,
						rownames_X_matr2[,nn] )
							}
						}
	stepgroups2 <- unique(rownames_X_matr2_collapse)
	match_stepgroups <- match( rownames_X_matr2_collapse , stepgroups2 )


	index_matr <- cbind( match_stepgroups , 1:NRX)
	index_matr <- index_matr[ order( index_matr[ , 1] ) , ]

	SG <- length(stepgroups2)	
	res <- .Call( "a_matrix_cumsum" , as.matrix(index_matr)-1 , as.matrix(mm) , SG ,
				PACKAGE="TAM")
	mm.sg.temp <- res$cumsum_mm
	rownames(mm.sg.temp) <- paste0("I", seq(1,nrow(mm.sg.temp) ) )
	ind2 <- seq( 1 , NRX+SG , maxK+1 )
	rownames(mm.sg.temp)[ind2] <- gsub("step([[:digit:]])*", "step0", stepgroups, fixed=T) 
	rownames(mm.sg.temp)[setdiff( seq(1,NRX+SG) , ind2) ] <- rownames(mm)[ index_matr[,2] ]
	colnames(mm.sg.temp) <- colnames(mm)
	A1 <- rbind(A , mm.sg.temp)
	
# cat(" +++  v157" ) ; z1 <- Sys.time() ; print(z1-z0) ; z0 <- z1     						  

	index_matr2 <- index_matr
	index_matr2 <- index_matr2[ index_matr2[,1] != c(0 , index_matr2[ -NRX , 1] ) , ]
	x.sg.temp <- X.out[ index_matr2[,2] , ]	
	x.sg.temp[,"step"] <- 0
	rownames(x.sg.temp) <- gsub("step([[:digit:]])*", "step0", stepgroups, fixed=T) 
	X.out1 <- rbind( X.out , x.sg.temp )
if (TRUE){
	X.out <- X.out1
	A <- A1
		}
if (FALSE){	
      for( sg in stepgroups ){
#        # sg <- stepgroups[1]
		ind.mm <- grep(sg, rownames(mm))
        mm.sg.temp <- rbind( 0, apply( mm[ ind.mm ,], 2, cumsum ) )        
#		mm.sg.temp <- rbind( 0 , colCumsums.sirt( mm[ ind.mm ,] ) )				
        # substitute the following line later if the sirt function
        # colCumsums.sirt is available at CRAN
        #		mm.sg.temp <- rbind( 0 , colCumsums.sirt( mm[ grep(sg, rownames(mm)) ,] ) )
        rownames(mm.sg.temp)[1] <- gsub("step([[:digit:]])*", "step0", sg, fixed=T)
        A <- rbind(A, mm.sg.temp)
		isg <- grep(sg, rownames(X.out))[1]	
        x.sg.temp <- X.out[ isg , ]
        x.sg.temp[,"step"] <- 0
        rownames(x.sg.temp) <- gsub("step([[:digit:]])*", "step0", sg, fixed = TRUE)
        X.out <- rbind(X.out, x.sg.temp)		 	

	if ( progress ){
          ii <- ii+1
          if (( ii == disp_progress[vv] ) & (vv<=10) ){
            cat("-") ; flush.console()
            vv <- vv+1			
          }
        }
      }   # end for (sg in stepgroups) ...
  }
	  

	  
      if ( progress ){
        cat("|\n") ; flush.console()
      }
#cat(" +++  v160" ) ; z1 <- Sys.time() ; print(z1-z0) ; z0 <- z1     				  	  		      
    } else {
      # model.matrix _ case: step not in fvars  
      rownames(mm) <- paste( rownames(X) , "-step1", sep = "")
      A <- mm
      
      for( kk in setdiff(0:maxK, 1) ){
        mm.k.temp <- mm*kk
        rownames(mm.k.temp) <- paste( rownames(X) , "-step", kk , sep ="")
        A <- rbind(A, mm.k.temp)
      }
# cat(" +++  v170" ) ; z1 <- Sys.time() ; print(z1-z0) ; z0 <- z1     				  	  		      
      X.out <- expand.grid( c( expand.list, list("step"=factor(0:maxK)) ) )
      X.out <- rownames.design2( data.frame(as.matrix(X.out), stringsAsFactors = FALSE) )
      
    }# end step in fvars
    
#cat(" +++  v180" ) ; z1 <- Sys.time() ; print(z1-z0) ; z0 <- z1     

    # facet design
    facet.design <- list( "facets" = facets , "facets.orig" = facets0 , 
                          "facet.list" = facet.list[otherFacets])
    A <- A[ ! duplicated( rownames(A) ) , ]
    A <- A[order(rownames(A)), ,drop = FALSE]      
    X.out <- X.out[order(rownames(X.out)), ,drop = FALSE]

	#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ ARb 2015-10-16	
	#@@@@ clean xsi.constr	
	xsi1 <- xsi.constr$xsi.constraints	
	xsi.constr$intercept_included <- FALSE
	ind <- grep("(Intercept" , rownames(xsi1) , fixed=TRUE)
	if ( length(ind) > 0 ){
		xsi1 <- xsi1[ - ind , ]
		xsi.constr$xsi.constraints <- xsi1
		xsi.constr$intercept_included <- TRUE	
							}
	xsi1 <- xsi.constr$xsi.table	
	ind <- grep("(Intercept" , paste(xsi1$parameter) , fixed=TRUE)
	if ( length(ind) > 0 ){
		xsi1 <- xsi1[ - ind , ]
		xsi.constr$xsi.table <- xsi1	
							}													
	#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
	
#cat(" +++  out .A.matrix" ) ; z1 <- Sys.time() ; print(z1-z0) ; z0 <- z1     						
    return(list( "A"=A, "X"=X.out, "otherFacets"=otherFacets , "xsi.constr"=xsi.constr ,
                 "facet.design" = facet.design ) )
  }
## end .A.matrix
#####################################################


####################################################
# create ConQuest parametrization for 
# partial credit model
.A.PCM2 <- function( resp , Kitem=NULL , constraint = "cases" , Q=NULL ){
  if ( is.null(Kitem) ){ 
    Kitem <- apply( resp , 2 , max , na.rm=T ) + 1
  }
  maxK <- max(Kitem)
  I <- ncol(resp)
  Nxsi <- sum(Kitem) - I
  A <- array( 0 , dim=c( I , maxK , Nxsi ) )
  vv <- 1
  for (ii in 1:I){
    A[ ii , 2:Kitem[ii] , vv ] <- - ( 2:Kitem[ii] - 1 )
    if ( Kitem[ii] < maxK ){
      A[ ii , ( Kitem[ii] + 1 ):maxK , ] <- NA                
    }
    vv <- vv+1
  }
  for (ii in 1:I){
    if ( Kitem[ii] > 2 ){
      for (kk in 2:(Kitem[ii] - 1) ){
        A[ ii , kk:(Kitem[ii]-1) , vv ] <- - 1
        vv <- vv + 1
      }
    }
  }
  dimnames(A)[[1]] <- colnames(resp)
  vars <- colnames(resp)
  # constraints
  unidim <- TRUE
  if ( ! is.null(Q) ){
	unidim <- ncol(Q) == 1 
				}
  if ( constraint == "items" ){  
     if ( unidim ){	
		I <- ncol(resp)
		x1 <- matrix( - A[I,,I] , nrow=dim(A)[2] , ncol=I-1 , byrow=FALSE )
		A[ I ,, seq(1,I-1) ] <- x1
		A <- A[,,-I]
		vars <- vars[ - I ]
					}
	 if (!unidim){	
          rem.pars <- NULL	 
			D <- ncol(Q)
			for (dd in 1:D){
				ind.dd <- which( Q[,dd] != 0 )
				I <- ind.dd[ length(ind.dd) ]				
				x1 <- matrix( - A[I,,I] , nrow=dim(A)[2] , ncol= length(ind.dd)-1 , byrow=FALSE )
				A[ I ,, ind.dd[ - length(ind.dd) ] ] <- x1
	
                rem.pars <- c(rem.pars , I )
						}
        vars <- vars[ - rem.pars ]
		A <- A[,, - rem.pars ]		
					}
					}					
  vars <- c(vars , unlist( sapply( (1:I)[Kitem>2] , FUN = function(ii){
    paste0( colnames(resp)[ii] , "_step" , 1:(Kitem[ii] - 2) ) } )	) )
  dimnames(A)[[3]] <- vars			
  return(A) 
}
#############################################################



####################################################
# dispersion model
.A.PCM3 <- function( resp , Kitem =NULL ){
  if ( is.null(Kitem) ){
    Kitem <- apply( resp , 2 , max , na.rm=T ) + 1
  }
  maxK <- max(Kitem)
  I <- ncol(resp)	
  Nxsi <- I + sum( Kitem > 2 )
  A <- array( 0 , dim=c( I , maxK , Nxsi ) )
  vv <- 1
  for (ii in 1:I){
    A[ ii , 2:Kitem[ii] , vv ] <- - ( 2:Kitem[ii] - 1 )
    if ( Kitem[ii] < maxK ){
      A[ ii , ( Kitem[ii] + 1 ):maxK , ] <- NA                
    }
    vv <- vv+1
  }
  for (ii in 1:I){
    if ( Kitem[ii] > 2 ){
      Kii <- Kitem[ii]-1
      A[ ii , 1:(Kii+1) , vv ] <- ( 0:Kii ) * ( Kii - ( 0:Kii) )						
      vv <- vv + 1
    }
  }
  dimnames(A)[[1]] <- colnames(resp)
  vars <- colnames(resp)
  #	vars <- c(vars , unlist( sapply( 1:I , FUN = function(ii){
  #		paste0( colnames(resp)[ii] , "_step" , 1:(Kitem[ii] - 2) ) } )	) )
  vars1 <- paste0( vars[ Kitem > 2 ] , "_disp" )
  vars <- c( vars , vars1 )
  dimnames(A)[[3]] <- vars			
  return(A) 
}
#############################################################
