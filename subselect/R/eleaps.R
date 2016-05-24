eleaps<-function(mat,kmin=length(include)+1,kmax=ncol(mat)-length(exclude)-1,nsol=1,exclude=NULL,include=NULL,
		  criterion="default",pcindices="first_k",timelimit=15,H=NULL,r=0,
		  tolval=1000*.Machine$double.eps,tolsym=1000*.Machine$double.eps,maxaperr=1E-4) 
{

###############################
# auxiliary  variables        #
###############################

	p <- ncol(mat)    				# Number of original variables
	nexclude <- length(exclude)     		# Number of excluded variables
	ninclude <- length(include)     		# Number of included variables
	if (pcindices[1]!="first_k") esp <- TRUE	# The user has specified the set of Principal Components to be used with the GCD criterion
	else esp <- FALSE				# The user has not specified the set of Principal Components to be used with the GCD criterion

        
###############################
# general validation of input #
###############################

         initialization(mat,criterion,r)
         if (validation(mat,kmin,kmax,exclude,include,criterion,pcindices,tolval,tolsym,"eleaps") == "singularmat")  singularmat <- TRUE  
         else singularmat <- FALSE         


#############################################################
# validation specific for the input to the eleaps function  #
#############################################################

	if (timelimit <= 0) stop("\n The time limit argument must be a positive real number\n")
        if ((criterion =="CCR1_2") && (r > 3)) stop("\n The 'eleaps' function does not accept, for the CCR1_2 criterion, \n an effects matrix (H) with  rank greater than 3\n")
        
######################################################################
# Parameter validation if the criterion is one of "TAU_2", "XI_2",   #
# "ZETA_2", "CCR1_2" or "WALD"                                       #
######################################################################

	if (criterion =="TAU_2" || criterion =="XI_2" || criterion =="ZETA_2" || criterion =="CCR1_2" || criterion == "WALD" )
	{
	   validnovcrit(mat,criterion,H,r,p,tolval,tolsym)
	   if (r==1 && criterion!="WALD")  criterion <- "XI_2"	# In multivariate linear problems with r=1 all four criteria are equivalent
	}	    						# but searches based on XI_2 or ZETA_2 are computationally more efficient 
        
##########################################
# Initializations for the C++ subroutine #
##########################################

       # Normalize matrices in order to improve numerical accuracy
	
	if (r==0) mat <- mat/max(mat)
	else {
		scl <- sqrt(diag(mat)) 
        	sclotprd <- outer(scl,scl)
		mat <- mat/sclotprd
        	H <- H/sclotprd
	}  

	# Matrix initializations
	
	if (singularmat==FALSE) Si <- solve(mat)
	else Si <- NULL					

	if (criterion == "RV")  S2 <- mat %*% mat
	else S2 <- NULL

	if (criterion == "GCD")  {
	  SSpectd <- eigen(mat,symmetric=TRUE)
	  Segval <- SSpectd$values[pcindices]
	  Segvct <- SSpectd$vectors[,pcindices]
        }
        else  {
           Segval <- NULL
           Segvct <- NULL
        }
 
	if (criterion == "TAU_2" || criterion == "ZETA_2" || (criterion == "CCR1_2" && r > 1) ) {
	   E <- mat - H
	   if (singularmat==FALSE) Ei <- solve(E)	
	   else Ei <- NULL				
        }
	else  E <- Ei <- NULL

	if (criterion == "XI_2" || criterion == "WALD" || criterion == "ZETA_2" || criterion == "CCR1_2" ) 	{
	  HSpectd <- eigen(H,symmetric=TRUE)
	  if (r > 1) Hegvct <- HSpectd$vectors[,1:r] %*% sqrt(diag(HSpectd$values[1:r])) 
	  else Hegvct <- HSpectd$vectors[,1] * sqrt(HSpectd$values[1]) 
        }
	else  Hegvct <- NULL
	if ( (criterion == "XI_2" || criterion == "WALD" || criterion == "CCR1_2" ) && singularmat==FALSE) 
		HegvctTinv <- solve(mat,Hegvct)
	else HegvctTinv <- NULL
	if ( (criterion == "ZETA_2" || (criterion == "CCR1_2" && r == 3)) && singularmat==FALSE)
		HegvctEinv <- solve(E,Hegvct) 
	else HegvctEinv <- NULL

	if ( ( criterion == "TAU_2" || (criterion == "CCR1_2" && r > 1) ) && singularmat==FALSE)
		Wilksval <- det(E) / det(mat)
	else Wilksval <- NULL
	if ( ( criterion == "XI_2" || criterion == "WALD" || criterion == "CCR1_2" ) && singularmat==FALSE)
		HSi <- t(solve(mat,H))
	if ( ( criterion == "XI_2" || (criterion == "CCR1_2" && r > 1) ) && singularmat==FALSE)
		BartPival <- sum(diag(HSi))
	else BartPival <- NULL
        if ( criterion == "WALD") Waldval <- sum(diag(HSi))
        else Waldval <- NULL
	if ( ( criterion == "ZETA_2" || (criterion == "CCR1_2" && r == 3) ) && singularmat==FALSE) 
		LawHotval <- sum(diag(solve(E,H)))	
	else LawHotval <- NULL
	if ( ( criterion == "CCR1_2") && singularmat==FALSE) 
		CCR12val <- as.numeric(eigen(HSi,symmetric=FALSE,only.values=TRUE)$values[1])	
	else CCR12val <- NULL
        
	if (criterion == "WALD")   {

###       Convert a min Wald problem into an (artificial) equivalent max XI2 problem
        
		Hegvct <- Hegvct / sqrt(Waldval)
		HegvctTinv <- HegvctTinv / sqrt(Waldval)
		criterion <- "XI_2"
		criterio <- 5
		BartPival <- 1.
	}

############################
# Call to the C subroutine #
############################

	 Cout <- .Call("eleaps",
           as.double(mat),
           as.double(S2),
           as.double(Si),
           as.double(Segval),
           as.double(Segvct),
           as.double(E),
           as.double(Ei),
           as.double(Hegvct),
           as.double(HegvctTinv),
           as.double(HegvctEinv),
           as.double(Wilksval),	
           as.double(BartPival),
           as.double(LawHotval),
           as.double(CCR12val),
	   as.integer(r),
           as.integer(kmin),
           as.integer(kmax),
           as.integer(nsol),
           as.integer(exclude),
           as.integer(include),
           as.integer(nexclude),
           as.integer(ninclude),
           as.character(criterion),
           as.logical(esp),
           as.integer(pcindices),
           as.integer(length(pcindices)),
           as.integer(p),
	   as.double(timelimit),
	   as.double(maxaperr), 	    
	   as.logical(singularmat), 	    
           PACKAGE="subselect"   
        ) 
	if (Cout$nomemory == TRUE) return(NULL)

#######################################
# Preparing and returning the output  #
#######################################

	 output <- prepRetOtp(c(Cout[1:4],call=match.call()),kmin,kmax,nsol,Waldval,"eleaps",Optimal=Cout$found,tl=timelimit) 
	 if (is.null(output)) return(invisible(NA))
	 output   # return(output)
}

