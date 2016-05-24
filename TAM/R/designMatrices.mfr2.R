

#########################################################################
designMatrices.mfr2 <-
  function( resp, formulaA = ~ item + item:step, facets = NULL,  
            constraint = c("cases", "items"), ndim = 1,
            Q=NULL, A=NULL, B=NULL , progress=FALSE ){
    z0 <- Sys.time()
	# sort items according alpha numeric sorting
	cn <- colnames(resp)
	if ( ! is.null(cn) ){
		cns <- sort( cn , index.return=TRUE )$ix
		resp <- resp[ , cns ]
		if ( !is.null(A) ){ A <- A[ cns ,, ] }
		if ( !is.null(B) ){ B <- B[ cns ,, ] }		
		if ( !is.null(Q) ){ Q <- Q[ cns , ] }				
	
		}


    ### Basic Information and Initializations
    constraint <- match.arg(constraint)
    ## restructure formulaA
    t1 <- attr( stats::terms( formulaA ) , "term.labels" )
    t2 <- intersect( c("item" , "step" , "item:step") , t1 )
        
# cat(" ---  z20" ) ; z1 <- Sys.time() ; print(z1-z0) ; z0 <- z1     	
    formulaA <- paste(  paste( c(t2 , setdiff(t1 , t2 ) ) , collapse= " + " ) )
    formulaA <- stats::as.formula( paste( " ~ " , formulaA ) )	
    
    #********************************
    # change formate in facets
    FF <- ncol(facets)
    NFF <- nrow(facets)
    if (progress){ 
      cat( "        o Check facets (" , paste(Sys.time()) , ")\n") ; 
	  utils::flush.console();
    }
    if ( is.null(FF) ){ FF <- 0 }
    if (FF>0){	
      for (ff in 1:FF){			
        # ff <- 1
        #**** inclusion ARb 2013-09-07
        #		cff <- nchar(facets[,ff] )
        cff <- nchar(paste( facets[,ff] ) )
        Mff <- max(cff)
        sff <- paste( rep("_" , Mff ) , collapse="" )
        if( min(cff) < Mff ){
          facets.ff0 <- facets[,ff]
          #			facets[,ff] <- paste0( facets[,ff] , substring( sff , 1 , Mff - cff ) )
          facets[,ff] <- paste0( "_" , facets[,ff] , substring( sff , 1 , Mff - cff ) )
          if (progress){
            u1 <- unique( setdiff( paste(facets[,ff]) , paste( facets.ff0 ) ) )			
            p1 <- paste0( "          * Changed levels of facet ", colnames(facets)[ff], ":" )
            p1 <- paste( p1 , paste( paste0("'",u1,"'") , collapse= " " ) )
            cat(p1, "\n")
          }
        }
      }
    }
# cat(" ---  z50" ) ; z1 <- Sys.time() ; print(z1-z0) ; z0 <- z1     				
    
    #********************************	
    #    resp[ is.na(resp) ] <- 0
    maxKi <- apply( resp , 2 , max , na.rm=TRUE )
    maxK <- max( maxKi )
    I <- nI <- ncol( resp )
    item <- rep( 1:nI , maxKi+1 )
	
    if ( is.null( colnames(resp) ) ){
      colnames(resp) <- paste0( "item" , 1:nI )
    }
    
# cat(" ---  before .A.matrix" ) ; z1 <- Sys.time() ; print(z1-z0) ; z0 <- z1     						


    
    # A Matrix
    if( is.null(A) ){
      AX <- .A.matrix2( resp, formulaA = formulaA, facets = facets, constraint = constraint ,
                       progress=progress , Q = Q)			
					   
      A <- AX$A; X <- AX$X; otherFacets <- AX$otherFacets
	  xsi.elim <- AX$xsi.elim  
      xsi.constr <- AX$xsi.constr
      facet.design <- AX$facet.design
      facet.list <- facet.design$facet.list
      facets <- facet.design$facets	  
      X.noStep <- unique(X[,- grep("step", colnames(X)), drop = FALSE ])   
	## bug
		# Fehler in `row.names<-.data.frame`(`*tmp*`, value = value) : 
		#   duplicate 'row.names' are not allowed  
#      rownames(X.noStep) <- gsub("-step([[:digit:]])*", "", rownames(X.noStep))
		Xnames.noStep <- gsub("-step([[:digit:]])*", "", rownames(X.noStep))

    }

	
#Revalprstr("X.noStep") 	
# cat(" ---  created A matrix (I) " ) ; z1 <- Sys.time() ; print(z1-z0) ; z0 <- z1     					
    # TK: 2014-03-12
    # Consider long vector response matrix
    if( "item" %in% colnames(facets) & "item" %in% t2 ) otherFacets <-  c("item", otherFacets)
    X.ind <-  if( "item" %in% colnames(facets) | nI==1) rep(1,nrow(X)) else as.numeric(X[,"item"])
    X.noStep.ind <-  if( "item" %in% colnames(facets) | nI==1) 
			rep(1,nrow(X.noStep)) else as.numeric(X.noStep[,"item"])
    
    if (progress){ 
      cat( "        o Created A Matrix (" , paste(Sys.time()) , ")\n") ; 
	  utils::flush.console();
    }
# Revalpr("X.noStep.ind")

    
# cat(" ---  created A matrix (II) " ) ; z1 <- Sys.time() ; print(z1-z0) ; z0 <- z1     					
    # Q Matrix
    
    if( is.null(Q) ){
      if( ndim > 1 ) warning("random q matrix")
      Q <- matrix( 0 , nrow = nI , ncol = ndim )
      Q[ cbind( 1:nI , sample(1:ndim, nI, replace=T) ) ] <- 1 
    } 
    Q <- Q[X.ind,,drop=FALSE]
    dimnames(Q) <- list( rownames(X), 
                         paste( "Dim" , add.lead( 1:ndim, ceiling(log(ndim, 10)) ), sep ="") 
    )
    # ndim
    ndim <- dim(Q)[2]

    # cat(" ---  after Q" ) ; z1 <- Sys.time() ; print(z1-z0) ; z0 <- z1        
    # B Matrix
    if( is.null(B) ){ 
      B.store.in <- NULL
      B <- Q * as.numeric(X$step)
    }else{
      B.store.in <- B
      B <- Q * as.numeric(X$step)
    }
    if (progress){ 
      cat( "        o Created B Matrix (" , paste(Sys.time()) , ")\n") ; 
	  utils::flush.console();
    }	

# cat(" ---  created B matrix " ) ; z1 <- Sys.time() ; print(z1-z0) ; z0 <- z1     						

	
    # gresp
    ind.resp.cols <- as.numeric(X.ind)
# cat(" ---  before gresp   " ) ; z1 <- Sys.time() ; print(z1-z0) ; z0 <- z1    
    gresp <- resp[,ind.resp.cols]	
# cat(" ---  after gresp selection   " ) ; z1 <- Sys.time() ; print(z1-z0) ; z0 <- z1    	
    #    gresp <- 1* ( gresp == t(array( X[, "step"], dim = dim(t(gresp)) )) )	
#    gresp <- 1*(gresp==matrix( as.numeric(X[,"step"]) , nrow(gresp) , ncol(gresp) , byrow=TRUE ))
#     cat(" ---  after gresp   " ) ; z1 <- Sys.time() ; print(z1-z0) ; z0 <- z1    
    # This step is time-consuming!!	
	#**** ARb 2014-05-30
	gresp <- .Call("gresp_extend" ,  as.matrix(gresp) , as.numeric( X[,"step"] ) ,
				PACKAGE="TAM")
    # cat(" ---  after gresp   " ) ; z1 <- Sys.time() ; print(z1-z0) ; z0 <- z1    
    # This step is time-consuming!!
# cat(" ---  after gresp (2) " ) ; z1 <- Sys.time() ; print(z1-z0) ; z0 <- z1     					    
  
    ind.resp.cols <- as.numeric(X.noStep.ind)
    gresp.noStep <- resp[,ind.resp.cols]
# cat(" ---  gresp no step    " ) ; z1 <- Sys.time() ; print(z1-z0) ; z0 <- z1
    if( length(otherFacets) > 0 ){
      rnFacets <- rownames( rownames.design2( as.matrix(facets[,otherFacets]) ))
      rnX <-      rownames( rownames.design2( as.matrix(X[,otherFacets]) ))
      rnX.noStep <-      rownames( rownames.design2( as.matrix(X.noStep[,otherFacets]) ))   
# cat("rownames.design2  X" ) ; z1 <- Sys.time() ; print(z1-z0) ; z0 <- z1		  
      #*** ARb 2013-03-26:
      #*** Set all entries in gresp and gresp.noStep to missing
      #*** if they are not observed.
      #      gresp <- gresp * (1* outer(rnFacets, rnX, "=="))
#      gresp[ outer(rnFacets, rnX, "!=") ] <- NA
      #      gresp.noStep <- gresp.noStep * (1* outer(rnFacets, rnX.noStep, "=="))
#      gresp.noStep[ outer(rnFacets, rnX.noStep, "!=") ] <- NA
# cat("gresp NA " ) ; z1 <- Sys.time() ; print(z1-z0) ; z0 <- z1	
		#**** ARb 2014-05-30
	  gresp <- .Call("gresp_na_facets" , as.matrix(gresp) , rnFacets , rnX ,
					PACKAGE="TAM")
	  gresp.noStep <- .Call("gresp_na_facets" , as.matrix(gresp.noStep ) , rnFacets , rnX.noStep ,
					PACKAGE="TAM")	  
#cat("gresp NA2 " ) ; z1 <- Sys.time() ; print(z1-z0) ; z0 <- z1	
    }

#cat(" ---  after other facets" ) ; z1 <- Sys.time() ; print(z1-z0) ; z0 <- z1    	
    colnames(gresp) <- rownames(X)
#    X$empty <- 1* (colSums( gresp, na.rm=TRUE ) == 0)
	X$empty <- .Call( "colsums_gresp" , gresp , PACKAGE="TAM")
	
    colnames(gresp.noStep) <- rownames(X.noStep)
#    X.noStep$empty <- 1* (colSums( gresp.noStep, na.rm=TRUE ) == 0)
	X.noStep$empty <- .Call( "colsums_gresp" , gresp.noStep , PACKAGE="TAM")

    ### output
    ind <- X[,"empty"] == 1
    nStep <- maxK+1
    nGenit <- nrow(X)
    
    .generate.3d <- function(x){      
      return( aperm( array( as.matrix(x), dim= c( nStep, nGenit/nStep, ncol(x) ), 
                            dimnames= list( paste("_step",0:maxK, sep= ""), 
                                            unique(gsub("-step([[:digit:]])*", "", rownames(x))),
                                            colnames(x) ) )
                     , c(2,1,3) 
      ) )
    }
    # generate B
    .generateB.3d <- function(x , diffK=FALSE){  
	  #@@ ARb 2014-10-21
		dim2 <- nGenit/nStep
		dimnames2 <- unique(gsub("(^|-)+step([[:digit:]])*", "", rownames(x)))
		dim2 <- length(dimnames2)
  
      x2 <- array( 0 , c(nStep , dim2 , ncol(x) ) )
      dimnames(x2) <- list( paste("_step",0:maxK, sep= ""), 
                            dimnames2 ,
                            colnames(x) )
		#@@@@@						
	for (ss in 0:(nStep-1)){
			  str.ss <- paste0("step",ss )
			  iss <- grep(  paste0(str.ss,"+(-|$)") , rownames(x) )# , fixed=TRUE )
			  str.ss2 <- gsub( paste0("(^|-)+",str.ss) , "" , rownames(x)[iss] )
			  x2[ss+1,str.ss2,] <- x[ iss , ]
			}	  
	  
      x2 <- aperm( x2 , c(2,1,3) )
	  	  	  
	  #****
      # handle differing number of categories
	  maxK <- max( maxKi )
      items.diffK <- which( maxKi < maxK )
	  if ( length( items.diffK) == 0 ){
	     diffK <- FALSE }
	  if ( diffK ){
		  for (ii in items.diffK ){
		#	  ii <- 2	
			  item.ii <- colnames(resp)[ii]
			  Ki.ii <- maxKi[ii]
			  ind1 <- grep( paste0( item.ii , "-") , dimnames(x2)[[1]] )
			  ind2 <- grep( paste0( "step" ,  Ki.ii ) , dimnames(x2)[[2]] )
			  ind3 <- grep( paste0( item.ii , ":") , dimnames(x2)[[3]] )
			  ind3b <- grep( paste0( "step")  , dimnames(x2)[[3]] )
			  ind3 <- intersect( ind3 , ind3b )
			  if ( length(ind1)*length(ind2)*length(ind3) > 0 ){
				   x2[ ind1 , ind2 , ind3 ] <- 0
									}
			 for (vv in seq( Ki.ii+1 , maxK)){ 
				 # vv <- Ki.ii + 1 
				 ind1 <- grep( paste0( item.ii , "-") , dimnames(x2)[[1]] )
				 if ( length(grep("-",dimnames(x2)[[1]] )) == 0 ){
#					ind1 <- grep( paste0( item.ii ) , dimnames(x2)[[1]] )					
					#***********************
					#### fix Michal Modzelewski 2014-12-09
					#### quick fix for the problem with general/specific names ####
					ind1 <- grep( paste0( item.ii, "$" ) , dimnames(x2)[[1]] )										
								}
				 ind2 <- grep( paste0( "step" ,  vv ) , dimnames(x2)[[2]] )
				 if ( length(ind1)*length(ind2) > 0 ){
					   x2[ ind1 , ind2 , ] <- NA
										}
							}
				}
			}
      return(x2)
    }    
# cat(" ---  before item rename" ) ; z1 <- Sys.time() ; print(z1-z0) ; z0 <- z1    
    #***
    # debugging ind manually
    ind <- FALSE * ind
    #*************
    # rename items
    # TK 2014-03-12: consider "item" in facets or long vector response MX
    if( "item" %in% colnames(facets) ){
      itemren <- data.frame( "item" = unique(facet.design$facets.orig[,"item"]) , "itemren" = paste0( "item" , unique(facet.design$facets[,"item"]) ) )
    } else {
#      itemren <- data.frame( "item" =  colnames(resp) , "itemren" = paste0( "item" , 1:nI ) )
	#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ ARb 2015-10-16		
	  n1 <- nI
		if ( xsi.constr$intercept_included ){ n1 <- nI - 1  }
      itemren <- data.frame( "item" =  colnames(resp)[1:n1] , 
			"itemren" = paste0( "item" , 1:n1 ) )
	#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@		

    }
# cat(".....\nbefore rename A" ) ; z1 <- Sys.time() ; print(z1-z0) ; z0 <- z1	
    # print("g100")	
    A <- .rename.items( matr=A , itemren )

    dimnames(A)[[1]] <- .rename.items2aa( vec=dimnames(A)[[1]] ,
                                          facet.list=facet.list , I=I )
    
    # cat(".rename.items (A)" ) ; z1 <- Sys.time() ; print(z1-z0) ; z0 <- z1		
    xsi.table <- xsi.constr$xsi.table
    #	A <- .rename.items3( matr=A , facet.list , I )	
    #cat(".rename.items3 (A)" ) ; z1 <- Sys.time() ; print(z1-z0) ; z0 <- z1			
    A <- .rename.items3a( matr=A , facet.list , I , cols=TRUE , xsi.table )	
    #cat(".rename.items3a (A)" ) ; z1 <- Sys.time() ; print(z1-z0) ; z0 <- z1			
    B <- .rename.items( matr=B , itemren )	
    # cat(".rename.items (B)" ) ; z1 <- Sys.time() ; print(z1-z0) ; z0 <- z1			
    dimnames(B)[[1]] <- dimnames(A)[[1]]	
    #	B <- .rename.items3( matr=B , facet.list )	
    gresp <- t( .rename.items( matr=t(gresp) , itemren , cols=FALSE)	)
    # cat(".rename.items (gresp)" ) ; z1 <- Sys.time() ; print(z1-z0) ; z0 <- z1				
    #	gresp <- t( .rename.items3( matr=t(gresp) , facet.list , cols=FALSE)	)	
    dimnames(gresp)[[2]] <- dimnames(A)[[1]]	
    gresp.noStep <- t( .rename.items( matr=t(gresp.noStep) , itemren , cols=FALSE)	)	
#cat("h2" ) ; z1 <- Sys.time() ; print(z1-z0) ; z0 <- z1					
    #	gresp.noStep <- t( .rename.items3( matr=t(gresp.noStep) , facet.list , I , cols=FALSE)	)	
    gresp.noStep <- t( .rename.items3a( matr=t(gresp.noStep) , facet.list , I , cols=FALSE ,
                                        xsi.table )	)		
#cat(".rename.items (gresp.noStep)" ) ; z1 <- Sys.time() ; print(z1-z0) ; z0 <- z1					
    Q <- .rename.items( matr=Q , itemren , cols=FALSE)

    dimnames(Q)[[1]] <- dimnames(A)[[1]]
    # cat(".rename.items (Q)" ) ; z1 <- Sys.time() ; print(z1-z0) ; z0 <- z1						
    # 	Q <- .rename.items3( matr=Q , facet.list , cols=FALSE)
    X <- .rename.items( matr=X , itemren , cols=FALSE)
    dimnames(X)[[1]] <- dimnames(A)[[1]]		
    #	X <- .rename.items3( matr=X , facet.list , cols=FALSE)	
    # cat(".rename.items (Q,X)" ) ; z1 <- Sys.time() ; print(z1-z0) ; z0 <- z1					
    X.noStep <- .rename.items( matr=X.noStep , itemren , cols=FALSE)
    #	X.noStep <- .rename.items3( matr=X.noStep , facet.list , cols=FALSE)	
    #***
    G1 <- xsi.constr$xsi.table 	
    G1$parameter <- .rename.items2( paste( G1$parameter) , itemren) 	
    # cat(".rename.items2 (G1$parameter)  " ) ; z1 <- Sys.time() ; print(z1-z0) ; z0 <- z1					
    #	G1$parameter <- .rename.items2a( paste( G1$parameter) , facet.list , I) 	
    #cat(".rename.items2a (G1$parameter)  " ) ; z1 <- Sys.time() ; print(z1-z0) ; z0 <- z1					
    G1$parameter <- .rename.items2b( paste( G1$parameter) , facet.list , I , xsi.table ) 	
    xsi.constr$xsi.table <- G1	
 # cat(".rename.items2b (G1$parameter)  " ) ; z1 <- Sys.time() ; print(z1-z0) ; z0 <- z1	
    #***
    G1 <- xsi.constr$xsi.constraints
    rownames(G1) <- .rename.items2( rownames(G1) , itemren) 	
    colnames(G1) <- .rename.items2( colnames(G1) , itemren) 
    # cat(".rename.items2 (colnames(G1))  " ) ; z1 <- Sys.time() ; print(z1-z0) ; z0 <- z1					
    colnames(G1) <- .rename.items2b( colnames(G1) , facet.list , I , xsi.table , sel1=1) 
    rownames(G1) <- .rename.items2b( rownames(G1) , facet.list , I , xsi.table , sel1=2) 		
    #cat(".rename.items2a (colnames(G1))  " ) ; z1 <- Sys.time() ; print(z1-z0) ; z0 <- z1					
    G1 -> xsi.constr$xsi.constraints
#cat("rename items" ) ; z1 <- Sys.time() ; print(z1-z0) ; z0 <- z1
	    
    if (progress){ 
      cat( "        o Relabeled Variable Names (" , paste(Sys.time()) , ")\n") ; 
	  utils::flush.console();
    }
    
# cat(" ---  after all item renames" ) ; z1 <- Sys.time() ; print(z1-z0) ; z0 <- z1    			
    
    # A
    A.flat.0 <- A.flat <- A; A.flat.0[ind,] <- 0
# cat("d300\n")
# Revalpr("A.flat")	
    A.3d <- .generateB.3d( A.flat , diffK=TRUE )		# here occurs the problem if one works with a reduced expand.grid	
# cat("d400\n")	
    A.flat <- A.flat[!ind,]
    A.3d.0 <- .generateB.3d( A.flat.0 , diffK=TRUE)
    # B                  
    B.flat.0 <- B.flat <- B; B.flat.0[ind,] <- 0
    B.3d <- .generateB.3d( B.flat )
    B.flat <- B.flat[!ind,]
    B.3d.0 <- .generateB.3d( B.flat.0 )
    if(!is.null(B.store.in)) B.3d.0[] <- B.store.in
    
    # Q                  
    Q.flat.0 <- Q.flat <- Q; Q.flat.0[ind,] <- 0
    Q.3d <- .generateB.3d( Q.flat )
    Q.flat <- Q.flat[!ind,]
    Q.3d.0 <- .generateB.3d( Q.flat.0 )     
# cat(" ---  output mfr" ) ; z1 <- Sys.time() ; print(z1-z0) ; z0 <- z1    			
    # out
    out <- list( "gresp" = list("gresp"=gresp, "gresp.noStep"=gresp.noStep), 
                 "A" = list("A.flat"=A.flat, "A.flat.0"=A.flat.0, 
                            "A.3d"=A.3d, "A.3d.0"=A.3d.0), 
                 "B" = list("B.flat"=B.flat, "B.flat.0"=B.flat.0, 
                            "B.3d"=B.3d, "B.3d.0"=B.3d.0), 
                 "Q" = list("Q.flat"=Q.flat, "Q.flat.0"=Q.flat.0,
                            "Q.3d"=Q.3d, "Q.3d.0"=Q.3d.0), 
                 "X" = list("X"=X, "X.noStep"=X.noStep) ,
                 "xsi.constr" = xsi.constr  , "xsi.elim"=xsi.elim
    )
    class(out) <- "designMatrices.mfr"
    return(out)
  }


