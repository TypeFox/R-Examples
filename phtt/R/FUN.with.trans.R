FUN.with.trans <- function(z, N, T, is.intercept, effect = c("none", "individual", "time", "twoways")){
	with.trans <- match.arg(effect)
        const <- ifelse(is.intercept, mean(z), 0)
        
	switch(with.trans	
           , none = {
			Z  	<- z - const
             	liste <- list(
			   "Tr"  = "none", # Name of *Tr*ansformation
                           "I"   = ifelse(is.intercept, TRUE, FALSE), 
                           "ODM" = z,        # *O*rig. *D*ata *M*atrix
                           "TDM" = Z,        # *T*ransformed *D*ata *M*atrix
                           "TDV" = c(Z),     # *T*ransformed *D*ata *V*ector
                           "TRm" = list(     # *TR*ansformation *m*eans. 
						"OVc"  = const,	# *OV*erall *c*onstant
						"InC"  = rep.int(const, N),	# *In*dividual *C*onstants
						"TiVC" = rep.int(const, T))	# *Ti*m V*arying *C*onstants
				   )    
                  return(liste)
           }
           , individual = {
			InC	<- colMeans(z)
			Z  	<- z - matrix(InC, T, N, byrow = TRUE) 
             	liste <- list(
			   "Tr"  = "individual", # Name of *Tr*ansformation
                           "I"   = ifelse(is.intercept, TRUE, FALSE), 
                           "ODM" = z,            # *O*rig. *D*ata *M*atrix
                           "TDM" = Z,            # *T*ransformed *D*ata *M*atrix
                           "TDV" = c(Z),         # *T*ransformed *D*ata *V*ector
                           "TRm" = list(	  	 # *TR*ansformation *m*eans. 
						"OVc"  = const,	# *OV*erall *c*onstant
						"InC"  = InC,	# *In*dividual *C*onstants
						"TiVC" = rep.int(const, T))	# *Ti*m V*arying *C*onstants
				   )  
                  return(liste)								
           }
           , time = {
			TiVC <- rowMeans(z)
	           	Z  	<- z - TiVC
             	liste <- list(
				   "Tr"  = "time",       # Name of *Tr*ansformation
                           "I"   = ifelse(is.intercept, TRUE, FALSE), 
                           "ODM" = z,            # *O*rig. *D*ata *M*atrix
                           "TDM" = Z,            # *T*ransformed *D*ata *M*atrix
                           "TDV" = c(Z),         # *T*ransformed *D*ata *V*ector
                           "TRm" = list(	  	 # *TR*ansformation *m*eans. 
						"OVc"  = const,	# *OV*erall *c*onstant
						"InC"  = rep.int(const, N), 	# *In*dividual *C*onstants
						"TiVC" = TiVC)	# *Ti*m V*arying *C*onstants
				   )  
                  return(liste)								
           }
           , twoways = {
			InC   <- colMeans(z)
			TiVC  <- rowMeans(z)
			Cons  <- mean(z)
             	Z  	<- z - matrix(InC, T, N, byrow = TRUE) - TiVC + Cons
             	liste <- list(
                              "Tr"  = "twoway",   # Name of *Tr*ansformation
                              "I"   = ifelse(is.intercept, TRUE, FALSE), 
                              "ODM" = z,          # *O*rig. *D*ata *M*atrix
                              "TDM" = Z,          # *T*ransformed *D*ata *M*atrix
                              "TDV" = c(Z),       # *T*ransformed *D*ata *V*ector
                              "TRm" = list(	  	 # *TR*ansformation *m*eans. 
						"OVc"  = Cons,	# *OV*erall *c*onstant
						"InC"  = InC, 	# *In*dividual *C*onstants
						"TiVC" = TiVC)	# *Ti*m V*arying *C*onstants
				   )  
                  return(liste)								
           }
    )		 		
  }

## Wrapper-Fct =================================================================================================
WithTrans <- function(z, intercept= TRUE, heto.effect = c("none", "individual", "time", "twoways"))
  {
	is.regular.panel(z, stopper = TRUE)
	nr <- nrow(z)
	nc <- ncol(z)
	with.trans <- match.arg(heto.effect)
	FUN.with.trans(z, N = nc, T = nr, is.intercept = intercept, effect = with.trans)
  }
##========================================================================================================
standardize <- function(z, MARGIN = c("rows", "columns")){
	is.regular.panel(z, stopper = TRUE)
	margc <-  match.arg(MARGIN)
	marg <-  switch(margc, rows = 1, columns = 2) 
	R <- apply(z, MARGIN = marg, function(x) (x- mean(x))/sqrt(var(x)))
	R <- t(R)
}

