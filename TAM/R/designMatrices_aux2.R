

###########################################################
## function .A.matrix
.A.matrix2 <-
  function( resp, formulaA = ~ item + item*step, facets = NULL,  
            constraint = c("cases", "items") , progress=FALSE ,
            maxKi = NULL , Q=Q ){
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
    #cat(" +++  v62" ) ; z1 <- Sys.time() ; print(z1-z0) ; z0 <- z1     						
    ### Basic Information and Initializations
    constraint <- match.arg(constraint)
    if ( is.null(maxKi) ){
      maxKi <- apply( resp , 2 , max , na.rm=TRUE )
    }
    maxK <- max( maxKi )
    nI <- ncol( resp )
    #cat(" +++  v70" ) ; z1 <- Sys.time() ; print(z1-z0) ; z0 <- z1     						
    
    # stop processing if there are items with a maximum score of 0
    i11 <- names(maxKi)[ maxKi == 0 ]
    if ( length(i11) > 0 ){
      stop( cat( "Items with maximum score of 0:" , paste(i11 , collapse=" " ) ) )
    }
    
    tf <- stats::terms( formulaA )	
    fvars <- as.vector( attr(tf,"variables"), mode = "character" )[-1]
    #cat("fvars 212") ; print(fvars)	
    otherFacets <- setdiff( fvars, c("item", "step") )
    contr.list <- as.list( rep( "contr.sum", length(fvars) ) )
    names( contr.list ) <- fvars
    #cat(" +++  v80" ) ; z1 <- Sys.time() ; print(z1-z0) ; z0 <- z1     							
    #****
    # ARb: 2013-03-27
    # no contrasts for items
    nitems <- ncol(resp)
    #	contr.list[["item"]] <- diag(1,nitems)    
    # TK: 2014-03-12
    # Consider long vector response matrix and "item" in facets input
    expand.list <- vector(mode = "list", length = 0)
    if( "item" %in% fvars )         expand.list <- c(expand.list, if("item" %in% names(facet.list)) list(as.factor(sort(unique(facets[,"item"])))) else list(factor(1:nI)) )
    if( "step" %in% fvars )         expand.list <- c(expand.list, if("step" %in% names(facet.list)) list(as.factor(sort(unique(facets[,"step"])))) else list(factor(1:maxK)) )
    if( length( otherFacets ) == 1) expand.list <- c(expand.list, list(factor(1:max(facets[, otherFacets]))) )     
    if( length( otherFacets ) > 1 ) expand.list <- c(expand.list, sapply( otherFacets , FUN = function(ff) as.factor(1:max(facets[, ff])) , simplify=FALSE ))
    
    names( expand.list ) <- fvars
    
	
    #cat(" +++  v100" ) ; z1 <- Sys.time() ; print(z1-z0) ; z0 <- z1     							  	  
    #     expand.list <- expand.list[ !unlist( lapply(expand.list, is.null) ) ]
    
    for (vv in seq(1 , length(expand.list) ) ){
      expand.list[[vv]] <- paste( expand.list[[vv]] ) 
    }
	
	
    # cat(" +++  v110" ) ; z1 <- Sys.time() ; print(z1-z0) ; z0 <- z1     								  	
	g2 <- g1 <- expand.grid(expand.list)
	diffK <- ( stats::sd( maxKi) > 0 )
	# diffK <- FALSE
	diffK <- TRUE
	# reduced combinations of items
	if (diffK){	
		I <- length(maxKi)
		# g1 <- as.data.frame(g1)
		for (ii in 1:I){
			ind <- which( ( ( as.numeric(paste0(g1$item)) == ii ) & 
					( as.numeric(paste0(g1$step)) > maxKi[ii]  ) )	 )
			if ( length(ind) > 0 ){
				g1 <- g1[ - ind  , ]		
								}
					}
			}
			
			
			
    X <- rownames.design2( g1 )
	if (diffK){ 
		X2 <- rownames.design2( g2 )
				}
	
	#**** include maxKi (2014-05-30)
#	X$maxKi <- maxKi[ X$item ]
#	X$par.exist <- 1 * ( as.numeric(paste(X$step)) <= X$maxKi )	
	#****
    ### constraints and formulaA
    if( constraint == "cases" ) formulaA <- stats::update.formula(formulaA, ~0+.)
    NX <- ncol(X)
    for (ff in 1:NX){
      uff <- length( unique(X[,ff] ) )
      if (uff==1){ cat(paste0("          - facet " ,
                              colnames(X)[ff] , " does only have one level!" ) , "\n") }
    }
    # cat(" +++  v120" ) ; z1 <- Sys.time() ; print(z1-z0) ; z0 <- z1     								  	 
    mm <- - stats::model.matrix(formulaA, X, contrasts = contr.list)
	
	if (diffK){
		mm2 <- - stats::model.matrix(formulaA, X2, contrasts = contr.list)
				}


				
    #    mm <- - model.matrix(formulaA, X )
    if( constraint == "items" ){ mm <- mm[,-1] }


    
    ############################################################
    ###*** ARb 2013-03-28
    ### generate all interactions	
    xsi.constr <- .generate.interactions2(X , facets , formulaA , mm )										
    ###############################################################
    #cat(" +++  v130" ) ; z1 <- Sys.time() ; print(z1-z0) ; z0 <- z1     						  	  		
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
      #cat(" +++  v150" ) ; z1 <- Sys.time() ; print(z1-z0) ; z0 <- z1     						  	  		
      if (progress){
        cat("        o Create design matrix A\n")
        ip <- length(stepgroups)
        VP <- min( ip , 10 )
        cat(paste0("          |",paste0( rep("*" , VP) , collapse="") , "|\n"))
        cat("          |") ; utils::flush.console()
        if (VP<10){ disp_progress <- 1:ip } else {
          disp_progress <- 100* ( 1:ip ) / (ip+1)
          disp_progress <- sapply( seq(5,95,10) , FUN = function(pp){ # pp <- 5
            which.min( abs( disp_progress - pp ) )[1] }
          )
        }						
      }
	  #******
	  # collect xsi parameters to be excluded
	  xsi.elim.index <- xsi.elim <- NULL 
      ii <- 0 ; vv <- 1
	  
	  
      for( sg in stepgroups ){
#         sg <- stepgroups[2]
	#	mm1 <- mm[ grep(sg, rownames(mm)) ,]
	   mm1 <- grep(paste0("(", sg, ")+$"), rownames(mm))		
	   # ind2 <- grep(sg, rownames(mm))	   
	   ind2 <- grep( paste0("(", sg, ")+$") , rownames(mm))
#	   if (length(ind2)>0){
          mm.sg.temp <- rbind( 0, apply( mm[ ind2 ,,drop=FALSE], 2, cumsum ) )
#						}	
	if ( is.null(rownames(mm.sg.temp)) ){
		rownames(mm.sg.temp) <- paste0("rn" , seq(0,nrow(mm.sg.temp)-1) )
								}
        # substitute the following line later if ...
        rownames(mm.sg.temp)[1] <- gsub("step([[:digit:]])*", "step0", sg, fixed=T)
		rownames(mm.sg.temp)[-1] <- rownames(mm[ind2,,drop=FALSE])
		#****
		# set entries to zero if there are no categories in data
		sg1 <- strsplit( sg , split= "-")[[1]]
		ii <- as.numeric( gsub("item" , "" , sg1[1] ) )
		if ( maxKi[ii] < maxK ){
		for (kk in (maxKi[ii]+1):maxK){
#			kk <- 2
			# set rows in A matrix to zero
			mm.sg.temp[ grep( paste0( "-step" , kk ) , rownames(mm.sg.temp) ) ,  ] <- NA				
#			mm.sg.temp[ grep( paste0( "-step" , kk ) , rownames(mm.sg.temp) ) ,  ] <- 0					
			i1 <- grep( paste0(sg1[1] ,"\\:" ) , colnames(mm.sg.temp) , value=TRUE)
			i2 <- grep( paste0(":step" , kk-1) ,  colnames(mm.sg.temp) , value=TRUE) 			
			i3 <- intersect( i1 , i2 )			
			if ( length(i3) > 0 ){
				xsi.elim <- c( xsi.elim , i3 )
				xsi.elim.index <- c( xsi.elim.index , 
							which( colnames(mm.sg.temp ) %in% i3 ) )
#				mm.sg.temp[  , i3 ] <- NA
#				mm.sg.temp[ ! ( is.na( mm.sg.temp[  , i3 ] ) ) , i3 ] <- 0

				#@@@ suggestion Michal Modzelewski 				
				 mm.sg.temp[ , i3 ][ ! ( is.na( mm.sg.temp[ , i3 ] ) )] <- 0				
		 
							}  #**** end i3 
			#*******************				
							
							
						}
					}				
			A <- rbind(A, mm.sg.temp)
		if ( maxKi[ii] < maxK ){
		for (kk in (maxKi[ii]+1):maxK){
			vv <- paste0( sg1[1] , "-step" , kk ) 
			X.out[ grep( vv , rownames(X.out) ) , 2  ] <- 0							
						}
					}						
			x.sg.temp <- X.out[grep(sg, rownames(X.out))[1], ]					
			x.sg.temp[,"step"] <- 0
			rownames(x.sg.temp) <- gsub("step([[:digit:]])*", "step0", sg, fixed = TRUE)
			X.out <- rbind(X.out, x.sg.temp)
        if ( progress ){
          ii <- ii+1
          if (( ii == disp_progress[vv] ) & (vv<=10) ){
            cat("-") ; utils::flush.console()
            vv <- vv+1			
              }
            }   # end progress
	
      }  # end stepgroups
      if ( progress ){
        cat("|\n") ; utils::flush.console()
      }
      # cat(" +++  v160" ) ; z1 <- Sys.time() ; print(z1-z0) ; z0 <- z1     				  	  		      
    } else {
      # model.matrix _ case: step not in fvars  
      rownames(mm) <- paste( rownames(X) , "-step1", sep = "")
      A <- mm
      
      for( kk in setdiff(0:maxK, 1) ){
        mm.k.temp <- mm*kk
        rownames(mm.k.temp) <- paste( rownames(X) , "-step", kk , sep ="")
        A <- rbind(A, mm.k.temp)
      }
      #cat(" +++  v170" ) ; z1 <- Sys.time() ; print(z1-z0) ; z0 <- z1     				  	  		      
      X.out <- expand.grid( c( expand.list, list("step"=factor(0:maxK)) ) )
      X.out <- rownames.design2( data.frame(as.matrix(X.out), stringsAsFactors = FALSE) )
      
    }# end step in fvars

	
	#***
	# set entries in A to zero for constraints
	A[ rowMeans(is.na(A)) < 1 , xsi.elim ] <- 0

    #cat(" +++  v180" ) ; z1 <- Sys.time() ; print(z1-z0) ; z0 <- z1     								  	 
    # facet design
    facet.design <- list( "facets" = facets , "facets.orig" = facets0 , 
                          "facet.list" = facet.list[otherFacets])
    A <- A[ ! duplicated( rownames(A) ) , ]
    A <- A[order(rownames(A)), ,drop = FALSE]      
    X.out <- X.out[order(rownames(X.out)), ,drop = FALSE]

	
	
	
	#*** elimination
	if ( ! is.null(xsi.elim) ){
		xsi.elim <- data.frame( xsi.elim , xsi.elim.index )
		xsi.elim <- xsi.elim[ ! duplicated( xsi.elim[,2] ) , ]
		xsi.elim <- xsi.elim[ order( xsi.elim[,2] ) , ]
#		A <- A[,-xsi.elim[,2] ]		 						
				}
				
				
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
                 "facet.design" = facet.design , "xsi.elim"=xsi.elim ) )
  }
## end .A.matrix
#####################################################

