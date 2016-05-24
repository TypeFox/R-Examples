validation<-function(mat, kmin, kmax, exclude, include, criterion, pcindices, tolval, tolsym, algorithm)
{

##########################################################
#  general validation of input for all search functions  #
##########################################################

         if (tolval < 0) stop("\n The 'tolval' argument must be non-negative.")

#################################################
# checking acceptability of criterion requested #
#################################################
         
 	labelsrm<-c("RM","Rm","rm","1",1)
       	labelsrv<-c("RV","Rv","rv","2",2)
       	labelsgcd<-c("GCD","Gcd","gcd","3",3)
	labelstau<-c("TAU_2","Tau_2","tau_2","TAU2","Tau2","tau2","TAU","Tau","tau","WILKS","Wilks","wilks")
	labelsxi<-c("XI_2","Xi_2","xi_2","XI2","Xi2","xi2","XI","Xi","xi","BARTLLET","Bartllet","bartllet")
	labelszeta<-c("ZETA_2","Zeta_2","zeta_2","ZETA2","Zeta2","zeta2","ZETA","Zeta","zeta","HOTELLING","Hotelling","hotelling")
       	labelsccr1<-c("CCR1_2","Ccr1_2","ccr1_2","CCR12","Ccr12","ccr12","CCR2","Ccr2","ccr2","CCR1","Ccr1","ccr1","ROY","Roy","roy")
	labelswald<-c("WALD","Wald","wald")
         	
	if (sum(criterion == c(labelsrm,labelsrv,labelsgcd,labelstau,labelsxi,labelszeta,labelsccr1,labelswald)) == 0) 
		stop("criterion requested is not catered for, or has been misspecified\n")

       	if (sum(criterion == labelsrm) > 0) criterio<-1
       	if (sum(criterion == labelsrv) > 0) criterio<-2
       	if (sum(criterion == labelsgcd) > 0) criterio<-3
       	if (sum(criterion == labelstau) > 0) criterio<-4
       	if (sum(criterion == labelsxi) > 0) criterio<-5
       	if (sum(criterion == labelszeta) > 0) criterio<-6
       	if (sum(criterion == labelsccr1) > 0) criterio<-7
       	if (sum(criterion == labelswald) > 0) criterio<-8
      		
	if (criterio == 1)  criterion <- "RM"
	if (criterio == 2)  criterion <- "RV"
	if (criterio == 3)  criterion <- "GCD"
	if (criterio == 4)  criterion <- "TAU_2"
	if (criterio == 5)  criterion <- "XI_2"
	if (criterio == 6)  criterion <- "ZETA_2"
	if (criterio == 7)  criterion <- "CCR1_2"
	if (criterio == 8)  criterion <- "WALD"

####################################################################
# checking for an input matrix that must be square, of full rank,  #
#    symmetric, positive definite  - calls the 'validmat' routine  #
####################################################################

         p <- dim(mat)[2]

        if (criterion=="TAU_2" || criterion=="XI_2" || criterion=="ZETA_2" || criterion=="CCR1_2")
        {   		
         	if (validmat(mat,p,tolval,tolsym,allowsingular=TRUE,algorithm) == "singularmat")  singularmat <- TRUE  	
         	else singularmat <- FALSE						 			
	}                                                				 			
        else {						 							
		validmat(mat,p,tolval,tolsym,allowsingular=FALSE,algorithm)   						
         	singularmat <- FALSE						 				
	}					 								
        
######################################################################
#  checking consistency of the requested values of kmin, kmax        #
#  (extreme sizes of variable subsets). The validation that kmax<=p  #
#   is made later, considering the no. of excluded variables         #
######################################################################         

         if (!is.numeric(kmin) || !is.numeric(kmax)) stop("\n Arguments kmin and kmax must be numeric.\n Perhaps unnamed arguments in the wrong order?")
         if (kmax < kmin) {
          aux<-kmin
          kmin<-kmax
          kmax<-aux
          warning("the argument kmin should precede the argument kmax.\n Since the value of kmin exceeded that of kmax, they have been swapped \n")}
          if (kmin != floor(kmin)) {
                  kmin<-floor(kmin)
                  warning("\n The value of 'kmin' requested was not an integer. It has been truncated. \n")
            }                  
          if (kmin >= p) {
             kmin<-p-1
             warning("\n The value of kmin requested is equal to or exceeds the number \n of variables. It has been set at p-1. \n")
            }
          if (kmin < 1) {
             kmin <- 1
             warning("\n The value of kmin requested is less than 1. It has been set at 1. \n")
            }
          if (kmax != floor(kmax)) {
                  kmax<-floor(kmax)
                  warning("\n The value of 'kmax' requested was not an integer. It has been truncated. \n")
            }   
         if (kmax >= p) {
         	kmax<-p-1
         	warning("\n The value of kmax requested is equal to or exceeds the number \n of variables. It has been set at p-1. \n")
            }

#############################################################
# checking for consistency of requests to exclude/include   #
#  certain variables.                                       #
#############################################################

         if (sum(duplicated(c(exclude,include))) > 0)
          stop("\n You have requested that the same variable be both included in, and excluded \n from, the subset.\n")
         nexclude<-length(exclude)
         if (nexclude !=0) {
         if (kmax >= p-nexclude) {
                                  kmax<-p-nexclude-1
                                  warning("\n Cardinalities requested are too large for the requested number of excluded \n variables, and kmax has been set at p-nexclude-1 \n")}
         exclude<-sort(exclude)}  
         exc<-c(0,exclude)
         ninclude<-length(include)
         if (ninclude !=0){
         if (kmin <= ninclude) {kmin<-ninclude+1
                                warning("\n Cardinalities requested are too small for the requested number of included \n variables, and kmin has been set to ninclude + 1 \n")}
         include<-sort(include)}
         inc<-c(0,include)
         if (kmax<kmin) stop("\n After trying to adapt to the requests for exclusion and inclusion of \n variables, kmax is now smaller than kmin. There must be a mistake.\n")


######################
# Checking pcindices #
######################


         
# Value always assigned to esp and to non-numeric pcindices to enable passing to parent.frame and to Fortran or C++ routines, for criteria other than GCD

         esp<-FALSE
         if (criterio == 3) { 
                            if (is.null(pcindices)) stop("\n For criterion GCD, argument pcindices must be non-NULL. \n")
                            if (is.numeric(pcindices))  {esp<-TRUE
                                    if  (sum(!(as.integer(pcindices) == pcindices)) > 0) stop("\n The PC indices must be integers.\n")
                                    if (max(pcindices)  >  p) stop("\n PCs of rank larger than the data set were requested. \n")
                                                        }
                            else {if (pcindices[1] != "first_k")
                                    {stop("\n unrecognized value for 'pcindices' argument \n")}}}
         if ((pcindices[1] == "first_k") || is.null(pcindices)) {pcindices<-1:kmax}

#################################
# assigning any changed values  #
#################################
         
         assign("mat",mat,pos=parent.frame())         
         assign("kmax",kmax,pos=parent.frame())         
         assign("kmin",kmin,pos=parent.frame())         
         assign("exclude",exclude,pos=parent.frame())         
         assign("include",include,pos=parent.frame())         
         assign("nexclude",nexclude,pos=parent.frame())         
         assign("ninclude",ninclude,pos=parent.frame())         
         assign("criterio",criterio,pos=parent.frame())   
         assign("criterion",criterion,pos=parent.frame())    
         assign("pcindices",pcindices,pos=parent.frame())         
         assign("p",p,pos=parent.frame())                  
         assign("exc",exc,pos=parent.frame())         
         assign("inc",inc,pos=parent.frame())         
         assign("esp",esp,pos=parent.frame())

         if (singularmat==TRUE) return("singularmat")   
         else return("validmat")			
       }


