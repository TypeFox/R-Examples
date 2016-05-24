
#**********************************************************
# calculates a latent response under the din function
din.latent.response <- function( q.matrix , S , rule="DINA"){
  Q <- as.matrix(q.matrix)
  S <- as.matrix(S)
  L <- matrix(nrow = nrow(Q), ncol = nrow(S))
  SQ <- S %*% t(Q)
  nums <- rowSums(Q)
  nums <- ifelse( rule=="DINO" , 1 , nums )
  nums <- matrix( nums , nrow=nrow(SQ) , ncol=ncol(SQ) , byrow=TRUE )
  SQ <- 1 * ( SQ >= nums  )
  L <- t(SQ)
  colnames(L) <- rownames(S)
  return(L)
}

###################################
# calculation of equivalent skill classes
din.equivalent.class <-function( q.matrix , rule="DINA"){
  Q <- q.matrix
  # Matrix with all skill classes
  S<-expand.grid( as.data.frame( t( matrix( rep( c(0,1), each = ncol(Q) ), ncol=2 ) ) ) )
  J <- nrow(Q)
  if ( length(rule) == 1){ rule <- rep( rule , J ) }
  rownames(S) <- paste0("Skills_" , apply( S , 1 , 
			FUN = function(ll){ paste(ll , collapse="" ) } )  )
    
  # Calculation of latent response of every skill class  
  A<-din.latent.response(Q,S,rule=rule)
  A<-t(A)
  I<-nrow(A)
  # calculate latent responses
  latent.response <- paste0("LatResp_" , 
		sapply( 1:I, FUN = function(ii){ paste( A[ ii, ], collapse="" )  } ) )
  
  skillclasses <- data.frame( "skillclass" = rownames(S) )
  skillclasses$latent.response <- latent.response
  
  # define distinguishable skill classes
  skillclasses$distinguish.class <- match( latent.response , unique( latent.response ) )
  # calculates how many skill classes correspond to the same latent response
   latent.response <- table( latent.response )
  six <- sort( latent.response, index.return=FALSE, decreasing=TRUE)  
	gini.CDM <- .gini( as.numeric(six) )
    res <- list( "latent.responseM" = A , "latent.response" = latent.response ,
				"S" = S , "gini" = gini.CDM , "skillclasses" = skillclasses )
  cat( nrow(S) , "Skill classes |" , max( skillclasses$distinguish.class ) , 
		" distinguishable skill classes |" ,
		"Gini coefficient =" , round( gini.CDM ,3 ) , "\n")
  return(res)
		}

#********************
# Gini coefficient , simply copied from the R ineq package
.gini <- function (x) 
{
    n <- length(x)
    x <- sort(x)
    G <- sum(x * 1:n)
    G <- 2 * G/(n * sum(x))
    G - 1 - (1/n)
}





