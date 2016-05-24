

#########################################################################
designMatrices.mfr <-
  function( resp, formulaA = ~ item + item:step, facets = NULL,  
            constraint = c("cases", "items"), ndim = 1,
            Q=NULL, A=NULL, B=NULL , progress=FALSE ){
			
tamcat_active <- TRUE
tamcat_active <- FALSE			
z0 <- Sys.time()	

    ### Basic Information and Initializations
    constraint <- match.arg(constraint)
    ## restructure formulaA
    t1 <- attr( stats::terms( formulaA ) , "term.labels" )
    t2 <- intersect( c("item" , "step" , "item:step") , t1 )
   
z0 <- tamcat( " ---  z20" , z0 , tamcat_active )    

    formulaA <- paste(  paste( c(t2 , setdiff(t1 , t2 ) ) , collapse= " + " ) )
    formulaA <- stats::as.formula( paste( " ~ " , formulaA ) )	
    
    #********************************
    # change formate in facets
    FF <- ncol(facets)
    NFF <- nrow(facets)
    if (progress){ 
      cat( "        o Check facets (" , paste(Sys.time()) , ")\n") ; flush.console();
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

z0 <- tamcat( " ---  z50" , z0 , tamcat_active )    
   
    #********************************	
    #    resp[ is.na(resp) ] <- 0
    maxKi <- apply( resp , 2 , max , na.rm=TRUE )
    maxK <- max( maxKi )
    I <- nI <- ncol( resp )
    item <- rep( 1:nI , maxKi+1 )
    if ( is.null( colnames(resp) ) ){
      colnames(resp) <- paste0( "item" , 1:nI )
    }
    
z0 <- tamcat( " ---  before .A.matrix" , z0 , tamcat_active )    
    
    # A Matrix
    if( is.null(A) ){
#!!!! work on speeding this step!!	
      AX <- .A.matrix( resp, formulaA = formulaA, facets = facets, constraint = constraint ,
                       progress=progress)
z0 <- tamcat( " ---  after .A.matrix" , z0 , tamcat_active )    

      A <- AX$A; X <- AX$X; otherFacets <- AX$otherFacets
      xsi.constr <- AX$xsi.constr
      facet.design <- AX$facet.design
      facet.list <- facet.design$facet.list
      facets <- facet.design$facets
      X.noStep <- unique(X[,- grep("step", colnames(X)), drop = FALSE ])  
      rownames(X.noStep) <- gsub("-step([[:digit:]])*", "", rownames(X.noStep))      
    } 	

z0 <- tamcat( " ---  A matrix (is.null(A)" , z0 , tamcat_active )    

    # TK: 2014-03-12
    # Consider long vector response matrix
    if( "item" %in% colnames(facets) & "item" %in% t2 ) otherFacets <-  c("item", otherFacets)
    X.ind <-  if( "item" %in% colnames(facets) | nI==1) rep(1,nrow(X)) else as.numeric(X[,"item"])
    X.noStep.ind <-  if( "item" %in% colnames(facets) | nI==1)
         	rep(1,nrow(X.noStep)) else as.numeric(X.noStep[,"item"])
    
    if (progress){ 
      cat( "        o Created A Matrix (" , paste(Sys.time()) , ")\n") ; flush.console();
    }
   
z0 <- tamcat( " ---  created A matrix" , z0 , tamcat_active )    

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

z0 <- tamcat( " ---  after Q" , z0 , tamcat_active )    

    # B Matrix
    if( is.null(B) ){ 
      B.store.in <- NULL
      B <- Q * as.numeric(X$step)
    }else{
      B.store.in <- B
      B <- Q * as.numeric(X$step)
    }
    if (progress){ 
      cat( "        o Created B Matrix (" , paste(Sys.time()) , ")\n") ; flush.console();
    }	

    # gresp
    ind.resp.cols <- as.numeric(X.ind)
z0 <- tamcat( " ---  before gresp" , z0 , tamcat_active )    	
    gresp <- resp[,ind.resp.cols]		
z0 <- tamcat( " ---  after gresp selection   " , z0 , tamcat_active )    		

#	res <- gresp_selection( as.matrix(resp) , ind.resp.cols-1 )
# no much time gain with Rcpp function
# cat(" ---  after gresp selection (Rcpp)  " ) ; z1 <- Sys.time() ; print(z1-z0) ; z0 <- z1    	
 #    gresp <- 1* ( gresp == t(array( X[, "step"], dim = dim(t(gresp)) )) )	        
#    gresp <- 1*(gresp==matrix( as.numeric(X[,"step"]) , nrow(gresp) , ncol(gresp) , byrow=TRUE ))

    # This step is time-consuming!!
	#**** ARb 2014-05-30
	gresp <- .Call("gresp_extend" , as.matrix(gresp) , as.numeric( X[,"step"] ) ,
					PACKAGE="TAM")    
z0 <- tamcat( " ---  after gresp   " , z0 , tamcat_active )      
        
    ind.resp.cols <- as.numeric(X.noStep.ind)
    gresp.noStep <- resp[,ind.resp.cols]
z0 <- tamcat( " ---  gresp ind.resp.cols   " , z0 , tamcat_active )      	

    if( length(otherFacets) > 0 ){
      rnFacets <- rownames( rownames.design2( as.matrix(facets[,otherFacets]) ))
      rnX <-      rownames( rownames.design2( as.matrix(X[,otherFacets]) ))
      rnX.noStep <-      rownames( rownames.design2( as.matrix(X.noStep[,otherFacets]) ))   
z0 <- tamcat( " ---  rownames.design2   " , z0 , tamcat_active )  

      #*** ARb 2013-03-26:
      #*** Set all entries in gresp and gresp.noStep to missing
      #*** if they are not observed.
      #      gresp <- gresp * (1* outer(rnFacets, rnX, "=="))
#      gresp[ outer(rnFacets, rnX, "!=") ] <- NA
      #      gresp.noStep <- gresp.noStep * (1* outer(rnFacets, rnX.noStep, "=="))
#      gresp.noStep[ outer(rnFacets, rnX.noStep, "!=") ] <- NA
	  gresp <- .Call("gresp_na_facets" ,  as.matrix(gresp) , rnFacets , rnX ,
				PACKAGE="TAM")
z0 <- tamcat( " ---  gresp na facets gresp   " , z0 , tamcat_active )  				
	  gresp.noStep <- .Call("gresp_na_facets", as.matrix(gresp.noStep ) , rnFacets , rnX.noStep ,
				PACKAGE="TAM"	  )	
z0 <- tamcat( " ---  gresp na facets gresp.noStep   " , z0 , tamcat_active )  				
      
    }
#  cat(" ---  after other facets" ) ; z1 <- Sys.time() ; print(z1-z0) ; z0 <- z1    	
    colnames(gresp) <- rownames(X)

#    X$empty <- 1* (colSums( gresp, na.rm=TRUE ) == 0)
# cat(" ---  col sums (gresp) in X " ) ; z1 <- Sys.time() ; print(z1-z0) ; z0 <- z1 	
# X$empty <- colsums_gresp( gresp )
# cat(" ---  col sums (gresp) in X (Rcpp)" ) ; z1 <- Sys.time() ; print(z1-z0) ; z0 <- z1 		
# X$empty <- colsums_gresp2( gresp )
# data communication with Rcpp does need no computation time

	X$empty <- .Call( "colsums_gresp" , gresp , PACKAGE="TAM")
z0 <- tamcat( " ---  col sums (gresp) in X (Rcpp)" , z0 , tamcat_active )  				 


    colnames(gresp.noStep) <- rownames(X.noStep)	
#    X.noStep$empty <- 1* (colSums( gresp.noStep, na.rm=TRUE ) == 0)
	X.noStep$empty <- .Call( "colsums_gresp" , gresp.noStep , PACKAGE="TAM")
z0 <- tamcat( " ---  col sums (gresp noStep) in X (Rcpp)" , z0 , tamcat_active )  				 

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
	
#**************************
    # generate B
    .generateB.3d <- function(x){  
      x2 <- array( 0 , c(nStep , nGenit/nStep , ncol(x) ) )
      dimnames(x2) <- list( paste("_step",0:maxK, sep= ""), 
                            unique(gsub("(^|-)+step([[:digit:]])*", "", rownames(x))),
                            colnames(x) )							
	for (ss in 0:(nStep-1)){
			  str.ss <- paste0("step",ss )
			  iss <- grep(  paste0(str.ss,"+(-|$)") , rownames(x) )# , fixed=TRUE )
			  str.ss2 <- gsub( paste0("(^|-)+",str.ss) , "" , rownames(x)[iss] )
			  x2[ss+1,str.ss2,] <- x[ iss , ]
			}
      x2 <- aperm( x2 , c(2,1,3) )
      return(x2)
    }    
#**************************	
	
# cat(" ---  before item rename" ) ; z1 <- Sys.time() ; print(z1-z0) ; z0 <- z1    
    #***
    # debugging ind manually
    ind <- FALSE * ind
    #*************
    # rename items
    # TK 2014-03-12: consider "item" in facets or long vector response MX
	
    if( "item" %in% colnames(facets) ){
      itemren <- data.frame( 
			"item" = unique(facet.design$facets.orig[,"item"]) , 
			"itemren" = paste0( "item" , unique(facet.design$facets[,"item"]) ) )		
    } else {
	#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ ARb 2015-10-16		
	  n1 <- nI
		if ( xsi.constr$intercept_included ){ n1 <- nI - 1  }
      itemren <- data.frame( "item" =  colnames(resp)[1:n1] , 
			"itemren" = paste0( "item" , 1:n1 ) )
	#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@			
    }

# cat(" --- .....before rename A" ) ; z1 <- Sys.time() ; print(z1-z0) ; z0 <- z1	
    A <- .rename.items( matr=A , itemren )
z0 <- tamcat( " --- .rename.items (A)" , z0 , tamcat_active )  				 


    # print( dimnames(A) )
    # print(facet.list)
    dimnames(A)[[1]] <- .rename.items2aa( vec=dimnames(A)[[1]] ,
                                          facet.list=facet.list , I=I )
z0 <- tamcat( " --- .rename.items2aa (A)" , z0 , tamcat_active )  				     
    xsi.table <- xsi.constr$xsi.table
    A <- .rename.items3a( matr=A , facet.list , I , cols=TRUE , xsi.table )	
z0 <- tamcat( " --- .rename.items3a (A)" , z0 , tamcat_active )

#    B <- .rename.items( matr=B , itemren )		
# z0 <- tamcat( " --- .rename.items (B)" , z0 , tamcat_active )  				     	
  
	dimnames(B)[[1]] <- dimnames(A)[[1]]		
z0 <- tamcat( " ---  dimnames(B)" , z0 , tamcat_active )  				     			
    #	B <- .rename.items3( matr=B , facet.list )		
#    gresp <- t( .rename.items( matr=t(gresp) , itemren , cols=FALSE) )
# z0 <- tamcat( " --- .rename.items (gresp)" , z0 , tamcat_active )  				     	 	
	dimnames(gresp)[[2]] <- dimnames(A)[[1]]
z0 <- tamcat( " ---  colnames(gresp) <- dimnames(A)" , z0 , tamcat_active )  				     	 

#Revalpr("sum(dimnames(A)[[1]] != colnames(gresp))")

 
    # gresp <- t( .rename.items3( matr=t(gresp) , facet.list , cols=FALSE)	)	
    dimnames(gresp)[[2]] <- dimnames(A)[[1]]	
    gresp.noStep <- t( .rename.items( matr=t(gresp.noStep) , itemren , cols=FALSE)	)	
z0 <- tamcat( " --- .rename.items (gresp.noStep)" , z0 , tamcat_active )

  				     	 
    #	gresp.noStep <- t( .rename.items3( matr=t(gresp.noStep) , facet.list , I , cols=FALSE)	)	
    gresp.noStep <- t( .rename.items3a( matr=t(gresp.noStep) , facet.list , I , cols=FALSE ,
                                        xsi.table )	)		
z0 <- tamcat( " --- .rename.items (gresp.noStep) facet list" , z0 , tamcat_active ) 


    Q <- .rename.items( matr=Q , itemren , cols=FALSE)
    dimnames(Q)[[1]] <- dimnames(A)[[1]]
z0 <- tamcat( " --- .rename.items (Q)" , z0 , tamcat_active ) 					
    # 	Q <- .rename.items3( matr=Q , facet.list , cols=FALSE)
    X <- .rename.items( matr=X , itemren , cols=FALSE)
    dimnames(X)[[1]] <- dimnames(A)[[1]]
    #	X <- .rename.items3( matr=X , facet.list , cols=FALSE)	
z0 <- tamcat( " --- .rename.items (Q,X)" , z0 , tamcat_active )	
    X.noStep <- .rename.items( matr=X.noStep , itemren , cols=FALSE)
z0 <- tamcat( " --- .rename.items (X.noStep)" , z0 , tamcat_active )	
    #	X.noStep <- .rename.items3( matr=X.noStep , facet.list , cols=FALSE)	
    #***
    G1 <- xsi.constr$xsi.table 	
    G1$parameter <- .rename.items2( paste( G1$parameter) , itemren) 	
z0 <- tamcat( " --- .rename.items2 (G1$parameter)" , z0 , tamcat_active )	
    #	G1$parameter <- .rename.items2a( paste( G1$parameter) , facet.list , I) 	
# cat(" --- .rename.items2a (G1$parameter)  " ) ; z1 <- Sys.time() ; print(z1-z0) ; z0 <- z1	
     G1$parameter <- .rename.items2b( paste( G1$parameter) , facet.list , I , xsi.table ) 	
    xsi.constr$xsi.table <- G1	
z0 <- tamcat( " --- .rename.items2b (G1$parameter)" , z0 , tamcat_active )	
    #***
    G1 <- xsi.constr$xsi.constraints
    rownames(G1) <- .rename.items2( rownames(G1) , itemren) 	
    colnames(G1) <- .rename.items2( colnames(G1) , itemren) 
z0 <- tamcat( " --- .rename.items2 (colnames(G1))" , z0 , tamcat_active )	
    colnames(G1) <- .rename.items2b( colnames(G1) , facet.list , I , xsi.table , sel1=1) 
    rownames(G1) <- .rename.items2b( rownames(G1) , facet.list , I , xsi.table , sel1=2) 		
z0 <- tamcat( " --- .rename.items2b (dimnames(G1))" , z0 , tamcat_active )	
    G1 -> xsi.constr$xsi.constraints
 #cat(" --- rename items" ) ; z1 <- Sys.time() ; print(z1-z0) ; z0 <- z1
    
    if (progress){ 
      cat( "        o Relabeled Variable Names (" , paste(Sys.time()) , ")\n") ; flush.console();
    }   
z0 <- tamcat( " --- after all renames" , z0 , tamcat_active )	
    # A
    A.flat.0 <- A.flat <- A; A.flat.0[ind,] <- 0
    A.3d <- .generateB.3d( A.flat )
    A.flat <- A.flat[!ind,]
    A.3d.0 <- .generateB.3d( A.flat.0 )
z0 <- tamcat( " --- A generate 3d" , z0 , tamcat_active )	
 
    # B                  
    B.flat.0 <- B.flat <- B; B.flat.0[ind,] <- 0
    B.3d <- .generateB.3d( B.flat )
    B.flat <- B.flat[!ind,]
    B.3d.0 <- .generateB.3d( B.flat.0 )
    if(!is.null(B.store.in)) B.3d.0[] <- B.store.in
z0 <- tamcat( " --- B generate 3d" , z0 , tamcat_active )	
    
    # Q                  
    Q.flat.0 <- Q.flat <- Q; Q.flat.0[ind,] <- 0
    Q.3d <- .generateB.3d( Q.flat )
    Q.flat <- Q.flat[!ind,]
    Q.3d.0 <- .generateB.3d( Q.flat.0 )     
z0 <- tamcat( " --- Q generate 3d" , z0 , tamcat_active )	

    # out
    out <- list( "gresp" = list("gresp"=gresp, "gresp.noStep"=gresp.noStep), 
                 "A" = list("A.flat"=A.flat, "A.flat.0"=A.flat.0, 
                            "A.3d"=A.3d, "A.3d.0"=A.3d.0), 
                 "B" = list("B.flat"=B.flat, "B.flat.0"=B.flat.0, 
                            "B.3d"=B.3d, "B.3d.0"=B.3d.0), 
                 "Q" = list("Q.flat"=Q.flat, "Q.flat.0"=Q.flat.0,
                            "Q.3d"=Q.3d, "Q.3d.0"=Q.3d.0), 
                 "X" = list("X"=X, "X.noStep"=X.noStep) ,
                 "xsi.constr" = xsi.constr 
    )
    class(out) <- "designMatrices.mfr"
    return(out)
  }


