improve<-function(mat, kmin, kmax=kmin, nsol=1, exclude=NULL,
include=NULL, setseed = FALSE, criterion="default", pcindices="first_k",
initialsol=NULL, force=FALSE, H=NULL,r=0, tolval=1000*.Machine$double.eps,
tolsym=1000*.Machine$double.eps)
{

#########################################################################################
# auxiliary  variables (includes declarations to avoid "no visible binding" NOTE)       #
#########################################################################################
  
	p <- ncol(mat)    				# Number of original variables
	nexclude <- length(exclude)     		# Number of excluded variables
	ninclude <- length(include)     		# Number of included variables				 
	if (pcindices[1]!="first_k") esp <- TRUE		# The user has specified the set of Principal Components to be used with the GCD criterion
	else esp <- FALSE				# The user has not specified the set of Principal Components to be used with the GCD criterion
	if (!is.null(initialsol)) silog <- TRUE		# The user has specified a set of initial solutions
	else silog <- FALSE				# The user has not specified a set of initial solutions
	exc <- exclude					# Name for vector of excluded variables after 'validation' routine 
	inc <- include					# Name for vector of included variables after 'validation' routine 

        
###############################
# general validation of input #
###############################

  
        initialization(mat, criterion, r)
         if (validation(mat, kmin, kmax, exclude, include, criterion, pcindices, tolval, tolsym, "improve") == "singularmat")  singularmat = TRUE  
         else singularmat = FALSE         
        maxnovar = 400
        if ((p > maxnovar) & (force==FALSE)) stop("\n For very large data sets, memory problems may crash the R session. \n To proceed anyways, repeat the function call with \n the argument 'force' set to 'TRUE' (after saving anything important \n from the current session)\n")


######################################################################
# Parameter validation if the criterion is one of "TAU_2", "XI_2",   #
# "ZETA_2" or "CCR1_2"  or "WALD"                                    #
######################################################################

if (criterion == "TAU_2" || criterion == "XI_2" || criterion ==
"ZETA_2" || criterion == "CCR1_2" || criterion == "WALD") validnovcrit(mat,criterion,H,r,p,tolval,tolsym)


##########################################################################
# more specific validation of input for the anneal and improve functions #
##########################################################################

        implog<-TRUE
        validannimp(kmin, kmax, nsol, exclude, nexclude, include, ninclude, initialsol, implog)


###########################################
# initializations for Fortran subroutine  #
###########################################

       # Normalize matrices in order to improve numerical accuracy
	
	if (r==0) mat <- mat/max(mat)
	else {
		scl <- sqrt(diag(mat)) 
        	sclotprd <- outer(scl,scl)
		mat <- mat/sclotprd
        	H <- H/sclotprd
	}  

        if (setseed == TRUE) set.seed(2,kind="default")
        valores<-rep(0.0,(kmax-kmin+1)*nsol)
        vars<-rep(0,nsol*length(kmin:kmax)*kmax)
        bestval<-rep(0.0,length(kmin:kmax))
        bestvar<-rep(0,kmax*length(kmin:kmax))
        if (criterio == 3) {
          decespectral<-eigen(mat,symmetric=TRUE)
          valp<-decespectral$values
          vecp<-decespectral$vectors
               }
           else {
              valp<-rep(0,p)
              vecp<-matrix(nrow=p,ncol=p,rep(0,p*p))
                 }

###############################
# call to Fortran subroutine  #
###############################

 	if (criterion == "WALD")   {

###       Convert a min Wald problem into an (artificial) equivalent max XI2 problem
        
        	Waldval <- sum(diag(solve(mat,H)))
            	H <- H / Waldval     
		criterion <- "XI_2"
		criterio <- 5
	}
        else Waldval <- NULL

       Fortout<-.Fortran("improve",as.integer(criterio),as.integer(p),
          as.double(as.vector(mat)),
          as.integer(kmin),as.integer(kmax),as.double(valores),
          as.integer(vars),as.double(bestval),as.integer(bestvar),
          as.integer(nexclude),as.integer(exc),as.integer(ninclude),
          as.integer(inc),as.integer(nsol),
          as.integer(length(pcindices)),as.integer(pcindices),
          as.logical(esp),as.logical(silog),
          as.integer(as.vector(initialsol)),as.double(valp),
          as.double(as.vector(vecp)),as.double(as.vector(H)),as.integer(r),
	  as.logical(singularmat),as.double(tolval),logical(1),
	  PACKAGE="subselect")

#######################################
# Preparing and returning the output  #
#######################################

        output <- prepRetOtp(c(Fortout[6:9],call=match.call()),kmin,kmax,nsol,Waldval,"improve",Numprb=Fortout[[26]]) 
	if (is.null(output)) return(invisible(NA))
	output   # return(output)
}


