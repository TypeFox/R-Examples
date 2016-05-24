#################################################################
#
# stmodelGLM.R
#
#######################
# stepp model: GLM     #
#######################
setClass("stmodelGLM",
	   representation(coltrt   = "numeric",	# treatment
				colY     = "numeric",   # outcome 
				trts     = "numeric",	# trt encoding
				MM	   = "ANY",   	# model matrix
				glm	   = "character", # glm type
				link	   = "character", # link function
				debug    = "numeric"    # debug flag
				),
	   prototype(coltrt = 0,
			 colY   = 0,
			 trts   = 0,
			 MM     = NULL,
			 glm	  = "",
			 link   = "",
			 debug  = 0
			)
	   )


setMethod("estimate",	
	    signature="stmodelGLM",
	    definition=function(.Object, sp, ...){
		nsubpop   <- sp@nsubpop
		subpop    <- sp@subpop
		coltrt    <- .Object@coltrt
		trts	    <- .Object@trts
		OuTcOmE   <- .Object@colY
		covariate <- .Object@MM
		glm	    <- .Object@glm
		link	    <- .Object@link
		debug	    <- .Object@debug != 0
 
		txassign  <- rep(NA, length(coltrt))
		txassign[which(coltrt == trts[1])] <- 1
		txassign[which(coltrt == trts[2])] <- 0

		sglmObs1  <- rep(0, nsubpop)
    		sglmObs2  <- rep(0, nsubpop)
    		sglmSE1   <- rep(0, nsubpop)
    		sglmSE2   <- rep(0, nsubpop)

		## Create a formula for the model:
		if (!is.null(covariate)){
		  has.intercept <- colnames(covariate)[1]=="(Intercept)"
		  if (has.intercept){
		    cstart <- 2
		    xnam <- c("txassign", colnames(covariate)[cstart:dim(covariate)[2]])
		    fmla <- as.formula(paste("OuTcOmE ~ ", paste(xnam, collapse= "+")))
		  } else {
		    cstart <- 1
		    xnam <- c("txassign", colnames(covariate)[cstart:dim(covariate)[2]])
		    fmla <- as.formula(paste("OuTcOmE ~ 0+", paste(xnam, collapse= "+")))  
		  }		  
		} else {
		  cstart <- 2	# default model has an intercept
		  fmla <- as.formula("OuTcOmE ~ txassign")
		}

	  	DIFF       <- rep(0, nsubpop)
    	  	DIFFSE     <- rep(0, nsubpop)
		logRATIO   <- rep(0, nsubpop)
		logRATIOSE <- rep(0, nsubpop)
  
    	  	for (i in 1:nsubpop) {
		  seli    <- subpop[,i]==1

		  if (glm == "gaussian"){
		      m1     <- glm(fmla, subset=seli, family=gaussian(link="identity"))
		  }
		  else
		  if (glm == "binomial"){
			m1     <- glm(fmla, subset=seli, family=binomial(link="logit"))
		  }
		  else
		  if (glm == "poisson") {
			m1     <- glm(fmla, subset=seli, family=poisson(link="log"))
		  }
		  else {
		      stop("Unknown Model !")
		  }
		  # use predict 
		  if (!is.null(covariate)){
		    mm1  <- covariate[seli,,drop=FALSE]
		    mm1  <- cbind(txassign[seli], mm1)
		    cm   <- apply(mm1, 2, mean)
		    if (has.intercept) cm   <- cm[-2]
		    names(cm)[1] <- "txassign"
		    cm[1]<- 1	# treatment indicator is 1
		    p1   <- predict(m1, newdata=data.frame(t(cm)), se.fit=TRUE)
		    cm[1]<- 0	# treatment indicator is 0
		    p2   <- predict(m1, newdata=data.frame(t(cm)), se.fit=TRUE)
		  } else {
		    p1   <- predict(m1, newdata=data.frame(txassign=1), se.fit=TRUE)
		    p2   <- predict(m1, newdata=data.frame(txassign=0), se.fit=TRUE)
		  }

		  if (glm == "gaussian"){
                sglmObs1[i]   <- p1$fit
	          sglmSE1[i]    <- p1$se

                sglmObs2[i]   <- p2$fit
                sglmSE2[i]    <- p2$se

	          DIFF[i]		<- coef(m1)["txassign"]
		    DIFFSE[i]	<- sqrt(diag(vcov(m1))[cstart])
		    #print(vcov(m1))
		    #print(c(sqrt(diag(vcov(m1))[1]), p1$se, p2$se)) 
		    logRATIO[i]	<- log(sglmObs1[i])-log(sglmObs2[i])
		    # use delta method
		    logRATIOSE[i]	<- sqrt((sglmSE1[i]^2)/(sglmObs1[i]^2) + (sglmSE2[i]^2)/(sglmObs2[i]^2))
		  } else
		  if (glm == "binomial"){
		    exp1 		<- exp(p1$fit)
                sglmObs1[i]   <- exp1/(1+exp1)
	          sglmSE1[i]    <- p1$se*(exp1/((1+exp1)^2))	# delta method

		    exp2 		<- exp(p2$fit)
                sglmObs2[i]   <- exp2/(1+exp2)
                sglmSE2[i]    <- p2$se*(exp2/((1+exp2)^2))	# delta method
	
	          DIFF[i]		<- sglmObs1[i] - sglmObs2[i]
		    DIFFSE[i]	<- sqrt(sglmSE1[i]^2 + sglmSE2[i]^2)
		    
		    logRATIO[i]	<- coef(m1)["txassign"]	# log OR from the model coefficient
		    logRATIOSE[i]	<- sqrt(diag(vcov(m1))[cstart])
		  } else
		  if (glm == "poisson"){
                sglmObs1[i]   <- exp(p1$fit)
	          sglmSE1[i]    <- p1$se*exp(p1$fit)	# delta method

                sglmObs2[i]   <- exp(p2$fit)
                sglmSE2[i]    <- p2$se*exp(p2$fit)	# delta method

	          DIFF[i]		<- sglmObs1[i] - sglmObs2[i]
	          DIFFSE[i]	<- sqrt(sglmSE1[i]^2 + sglmSE2[i]^2)
		    logRATIO[i]	<- coef(m1)["txassign"] 
		    logRATIOSE[i]	<- sqrt(diag(vcov(m1))[cstart])	# log RR from the model coefficient
		  }
	      }  # for each subpopulation

		if (glm == "gaussian"){ 
	  	  sdiffw    <- sum(DIFF/DIFFSE)
		  logRATIOw <- sum(logRATIO/logRATIOSE)
	  	  m1all     <- glm(fmla, family=gaussian(link="identity"))
		} else 
		if (glm == "binomial"){
	  	  sdiffw   	<- sum((sglmObs1-sglmObs2)/sqrt(sglmSE1^2 + sglmSE2^2))
		  slogRRw 	<- exp(sum(logRATIO/logRATIOSE))		
	  	  m1all     <- glm(fmla, family=binomial(link="logit"))
		} else 
		if (glm == "poisson"){
	  	  sdiffw   	<- sum((sglmObs1-sglmObs2)/sqrt(sglmSE1^2 + sglmSE2^2))
		  slogRRw 	<- exp(sum(logRATIO/logRATIOSE))
	  	  m1all     <- glm(fmla, family=poisson(link="log"))
		}
		if (!is.null(covariate)){
		  om1    <- covariate
		  om1    <- cbind(txassign, om1)
		  onv1   <- apply(om1,2,mean)
		  if (has.intercept) onv1 <- onv1[-2]
		  names(onv1)[1] <- "txassign"
		  onv1[1]<- 1	# treatment indicator is 1
		  op1    <- predict(m1all, newdata=data.frame(t(onv1)), se.fit=TRUE)
		  onv1[1]<- 0	# treatment indicator is 0
		  op2    <- predict(m1all, newdata=data.frame(t(onv1)), se.fit=TRUE)
		} else {
		  op1    <- predict(m1all, newdata=data.frame(txassign=1), se.fit=TRUE)
		  op2    <- predict(m1all, newdata=data.frame(txassign=0), se.fit=TRUE)
		}

		if (glm == "gaussian"){
	  	  overallSglmObs1   <- op1$fit
		  overallSglmSE1    <- op1$se

	  	  overallSglmObs2   <- op2$fit
		  overallSglmSE2    <- op2$se

	  	  overallDIFF       <- coef(m1all)["txassign"]
	  	  overallDIFFSE	  <- sqrt(diag(vcov(m1all))[cstart])
		  overalllogRATIO	  <- log(overallSglmObs1)-log(overallSglmObs2)
		  # delta method
		  overalllogRATIOSE <- sqrt((overallSglmSE1^2)/(overallSglmObs1^2)
							+(overallSglmSE2^2)/(overallSglmObs2^2)
						  )	
		} else 
		if (glm == "binomial"){
		  oexp1 		<- exp(op1$fit)
	  	  overallSglmObs1 <- oexp1/(1+oexp1)
		  overallSglmSE1  <- op1$se*(exp1/((1+exp1)^2))	# delta method

		  oexp2 		<- exp(op2$fit)
	  	  overallSglmObs2 <- oexp2/(1+oexp2)
		  overallSglmSE2  <- op2$se*(exp2/((1+exp2)^2))	# delta method

	  	  oRDIFF      	<- overallSglmObs1 - overallSglmObs2
	  	  oRDIFFSE	  	<- sqrt(overallSglmSE1^2 + overallSglmSE2^2)
		    
		  ologRR		<- coef(m1all)["txassign"]
		  ologRRSE		<- sqrt(diag(vcov(m1all))[cstart])	
		} else
		if (glm == "poisson") {
	  	  overallSglmObs1 <- exp(op1$fit)
		  overallSglmSE1  <- op1$se*exp(op1$fit)	# delta method

	  	  overallSglmObs2 <- exp(op2$fit)
		  overallSglmSE2  <- op2$se*exp(op2$fit)	# delta method

	  	  oRDIFF      	<- overallSglmObs1 - overallSglmObs2
	  	  oRDIFFSE	  	<- sqrt(overallSglmSE1^2 + overallSglmSE2^2)
		  ologRR		<- coef(m1all)["txassign"]
		  ologRRSE 	  	<- sqrt(diag(vcov(m1all))[cstart])
    	 	}
    	
    	  	if (sum(is.na(sglmObs1)) != 0 | sum(is.na(sglmObs2)) != 0 | 
                    is.na(overallSglmObs1) != FALSE | is.na(overallSglmObs2) != FALSE) {
          	  cat("\n")
        	  print(paste("Unable to estimate the effect because there are too few events within one or more subpopulation(s)."))
        	  print(paste("The problem may be avoided by constructing larger subpopulations."))
        	  stop()
    	  	}
		  
		if (glm == "gaussian"){ 
	  	  #
	  	  #   create the estimates (GLM) object - to be saved
	   	  #
    	  	  estimate <- list( model	  = "GLMGe",
		  			  sObs1       = sglmObs1,
        	              	  sSE1        = sglmSE1,
                          	  oObs1       = overallSglmObs1,
                          	  oSE1        = overallSglmSE1,
                          	  sObs2       = sglmObs2,
                          	  sSE2        = sglmSE2,
                          	  oObs2       = overallSglmObs2,
                          	  oSE2        = overallSglmSE2,
				  	  sglmw	  = sdiffw,
				  	  RD		  = DIFF,
				  	  RDSE	  = DIFFSE,
				  	  oRD		  = overallDIFF,
				  	  oRDSE	  = overallDIFFSE,
					  logR	  = logRATIO,
					  logRSE	  = logRATIOSE,
					  ologR	  = overalllogRATIO,
					  ologRSE	  = overalllogRATIOSE,
					  sglmlogrw	  = logRATIOw
				      )
		  } else
		  if (glm=="binomial"){
	  	    #
	  	    #   create the estimates (GLM) object - to be saved
	  	    #
    	  	    estimate <- list( model	  = "GLMBe",
					    sObs1     = sglmObs1,
        	              	    sSE1  	  = sglmSE1,
                          	    oObs1     = overallSglmObs1,
                          	    oSE1  	  = overallSglmSE1,
                          	    sObs2     = sglmObs2,
                          	    sSE2  	  = sglmSE2,
                          	    oObs2     = overallSglmObs2,
                          	    oSE2  	  = overallSglmSE2,
				  	    sglmw	  = sdiffw,
				  	    RD	  = DIFF,
				  	    RDSE	  = DIFFSE,
				  	    oRD	  = oRDIFF,
				  	    oRDSE	  = oRDIFFSE,
					    logR	  = logRATIO,
					    logRSE	  = logRATIOSE,
					    ologR	  = ologRR,
					    ologRSE	  = ologRRSE,
					    sglmlogrw = slogRRw
				        )
 		  } else
		  if (glm=="poisson"){
		    #
	  	    #   create the estimates (GLM) object - to be saved
	  	    #
    	  		estimate <- list( model	  = "GLMPe",
					    sObs1     = sglmObs1,
        	              	    sSE1      = sglmSE1,
                          	    oObs1     = overallSglmObs1,
                          	    oSE1      = overallSglmSE1,
                          	    sObs2     = sglmObs2,
                          	    sSE2      = sglmSE2,
                          	    oObs2     = overallSglmObs2,
                          	    oSE2      = overallSglmSE2,
				  	    sglmw	  = sdiffw,
				  	    RD	  = DIFF,
				  	    RDSE	  = DIFFSE,
				  	    oRD	  = oRDIFF,
				  	    oRDSE	  = oRDIFFSE,
					    logR	  = logRATIO,
					    logRSE	  = logRATIOSE,
					    ologR	  = ologRR,
					    ologRSE	  = ologRRSE,
					    sglmlogrw = slogRRw
				        )
 		  }
	    return(estimate)
	    }
)


setMethod("test",
	    signature="stmodelGLM",
	    definition=function(.Object, nperm, sp, effect, showstatus=TRUE){

		nsubpop	<- sp@nsubpop
		subpop	<- sp@subpop
		coltrt	<- .Object@coltrt
		trts		<- .Object@trts
		OuTcOmE   	<- .Object@colY
		covariate 	<- .Object@MM
		glm	      <- .Object@glm
		link		<- .Object@link
		debug		<- .Object@debug != 0

		test 		<- NULL

		if (nperm > 0){

		  txassign  <- rep(NA, length(coltrt))
		  txassign[which(coltrt == trts[1])] <- 1
		  txassign[which(coltrt == trts[2])] <- 0
		  #
		  #   do the permutations
		  #
    		  pvalue      <- NA
    		  differences <- matrix(rep(0, (nperm * nsubpop)), ncol = nsubpop)
		  diffha      <- matrix(rep(0, (nperm * nsubpop)), ncol = nsubpop)
    		  logratios   <- matrix(rep(0, (nperm * nsubpop)), ncol = nsubpop)
      	  logratioha  <- matrix(rep(0, (nperm * nsubpop)), ncol = nsubpop)

    		  tPerm       <- rep(0, nperm)
    		  no          <- 0
    		  p           <- 0
	    	  terminate   <- 0
    		  Ntemp       <- nrow(subpop)
    		  IndexSet1   <- (1:Ntemp)[txassign == 1]
    		  IndexSet2   <- (1:Ntemp)[txassign == 0]

		  ## Create a formula for the model:
		  ## Cannot handle a no intercept model here
		  if (!is.null(covariate)){
		    has.intercept <- colnames(covariate)[1]=="(Intercept)"
		    if (has.intercept){
			cstart <- 2
		      xnam <- c("txassign", colnames(covariate)[cstart:dim(covariate)[2]])
		      fmla <- as.formula(paste("OuTcOmE ~ ", paste(xnam, collapse= "+")))
		    } else {
			cstart <- 1
		      xnam <- c("txassign", colnames(covariate)[cstart:dim(covariate)[2]])
		      fmla <- as.formula(paste("OuTcOmE ~ 0+", paste(xnam, collapse= "+")))  
		    }	
		  } else {
		    cstart <- 2	# default model has an intercept
		    fmla <- as.formula("OuTcOmE ~ txassign")
		  }
	  	  if (showstatus){
		    title <- paste("\nComputing the pvalue with ", nperm)
		    title <- paste(title, "number of permutations\n")
		    cat(title)
	  	    pb <- txtProgressBar(min=0, max=nperm-1, style=3)
		  }
        	  while (no < nperm) {
		    ## PERMUTE THE VECTOR OF SUBPOPULATIONS WITHIN TREATMENTS ##
      	    if (showstatus) setTxtProgressBar(pb, no)

      	    Subpop        <- as.matrix(subpop)
      	    permuteSubpop <- matrix(0, nrow = Ntemp, ncol = nsubpop)
      	    permuteSubpop[txassign == 1] <- Subpop[sample(IndexSet1),]
      	    permuteSubpop[txassign == 0] <- Subpop[sample(IndexSet2),]
      	    subpop        <- permuteSubpop
      	    sglm1         <- rep(0, nsubpop)
      	    sglm2         <- rep(0, nsubpop)
		    sglmSE1	      <- rep(0, nsubpop)
		    sglmSE2	      <- rep(0, nsubpop)
      	    slogratios    <- rep(0, nsubpop)
		    slogratioSEs  <- rep(0, nsubpop)

      	    for (i in 1:nsubpop) {
			seli <- subpop[,i]==1

			if (!is.null(covariate)){
		        mm1  <- covariate[seli,,drop=FALSE]
			  mm1  <- cbind(txassign[seli], mm1)
		        nv1  <- apply(mm1,2,mean)
			  if (has.intercept) nv1 <- nv1[-2]
		        names(nv1)[1] <- "txassign"
			  nv1[1] <- 1	# treatment indicator is 1
		        nv2  <- nv1
			  nv2[1] <- 0	# treatment indicator is 0
			}

			if (glm == "gaussian"){
			  m1 <- glm(fmla, subset=seli, family=gaussian(link="identity"))

			  if (!is.null(covariate)){
		          p1  <- predict(m1, newdata=data.frame(t(nv1)), se.fit=TRUE)
			    p2  <- predict(m1, newdata=data.frame(t(nv2)), se.fit=TRUE)
		  	  } else {
		          p1  <- predict(m1, newdata=data.frame(txassign=1), se.fit=TRUE)
			    p2  <- predict(m1, newdata=data.frame(txassign=0), se.fit=TRUE)
			  }
                    sglm1[i] 	 <- p1$fit
	              sglmSE1[i] <- p1$se	        
                    sglm2[i] 	 <- p2$fit
                    sglmSE2[i] <- p2$se
		    	  slogratios[i]	<- log(sglm1[i])-log(sglm2[i])
		        slogratioSEs[i]	<- sqrt((sglmSE1[i]^2)/(sglm1[i]^2) + (sglmSE2[i]^2)/(sglm2[i]^2))
			} else
			if (glm == "binomial"){
	      	  m1 <- glm(fmla, subset=seli, family=binomial(link="logit"))

			  if (!is.null(covariate)){
		          p1   <- predict(m1, newdata=data.frame(t(nv1)), se.fit=TRUE)
			    p2   <- predict(m1, newdata=data.frame(t(nv2)), se.fit=TRUE)
			  } else {
		          p1   <- predict(m1, newdata=data.frame(txassign=1), se.fit=TRUE)
			    p2   <- predict(m1, newdata=data.frame(txassign=0), se.fit=TRUE)
			  }
		        exp1 <- exp(p1$fit)
                    sglm1[i]   <- exp1/(1+exp1)
	              sglmSE1[i] <- p1$se*(exp1/((1+exp1)^2))	# delta method
		        exp2 <- exp(p2$fit)
                    sglm2[i]   <- exp2/(1+exp2)
                    sglmSE2[i] <- p2$se*(exp2/((1+exp2)^2))	# delta method

		    	  slogratios[i]	<- coef(m1)["txassign"]
		        slogratioSEs[i]	<- sqrt(diag(vcov(m1))[cstart])
			} else
			if (glm == "poisson"){
	      	  m1 <- glm(fmla, subset=seli, family=poisson(link="log"))

			  if (!is.null(covariate)){
		          p1   <- predict(m1, newdata=data.frame(t(nv1)), se.fit=TRUE)
		          p2   <- predict(m1, newdata=data.frame(t(nv2)), se.fit=TRUE)
		  	  } else {
		          p1   <- predict(m1, newdata=data.frame(txassign=1), se.fit=TRUE)
			    p2   <- predict(m1, newdata=data.frame(txassign=0), se.fit=TRUE)
			  }
                    sglm1[i]   <- exp(p1$fit)
	              sglmSE1[i] <- p1$se*exp(p1$fit)	# delta method
                    sglm2[i]   <- exp(p2$fit)
                    sglmSE2[i] <- p2$se*exp(p2$fit)	# delta method
		    	  slogratios[i]	<- coef(m1)["txassign"]
		        slogratioSEs[i]	<- sqrt(diag(vcov(m1))[cstart])
			} else
			stop ("Unsupported GLM !")

        	  }
	  	  sglmwha         <- sum((sglm1-sglm2)/sqrt(sglmSE1^2 + sglmSE2^2))
		  if (glm == "gaussian")
		    slogrglmha	<- sum(slogratios/slogratioSEs)
		  else if (glm == "binomial" | glm == "poisson")
		    slogrglmha	<- exp(sum(slogratios/slogratioSEs))	

		  if (!is.null(covariate)){
		    mm1all <- covariate
		    mm1all <- cbind(txassign, mm1all)
		    nv1    <- apply(mm1all,2,mean)
		    if (has.intercept) nv1 <- nv1[-2]
		    names(nv1)[1] <- "txassign"
		    nv1[1] <- 1
		    nv2    <- nv1
		    nv2[1] <- 0
		  } 

		  if (glm == "gaussian"){
	  	    m1all         <- glm(fmla, family=gaussian(link="identity"))

		    if (!is.null(covariate)){
		      p1  	     	<- predict(m1all, newdata=data.frame(t(nv1)), se.fit=TRUE)
			p2  		<- predict(m1all, newdata=data.frame(t(nv2)), se.fit=TRUE)
		    } else {
		      p1   		<- predict(m1all, newdata=data.frame(txassign=1), se.fit=TRUE)
			p2  		<- predict(m1all, newdata=data.frame(txassign=0), se.fit=TRUE)
		    }
                overallSglm1 	<- p1$fit		    
                overallSglm2 	<- p2$fit
		    overalllogRatio <- log(overallSglm1)-log(overallSglm2)
		  } else
		  if (glm == "binomial"){
	  	    m1all         <- glm(fmla, family=binomial(link="logit"))
			
		    if (!is.null(covariate)){
		      p1   		<- predict(m1all, newdata=data.frame(t(nv1)), se.fit=TRUE)
			p2   		<- predict(m1all, newdata=data.frame(t(nv2)), se.fit=TRUE)
		    } else {
		      p1  		<- predict(m1all, newdata=data.frame(txassign=1), se.fit=TRUE)
			p2  		<- predict(m1all, newdata=data.frame(txassign=0), se.fit=TRUE)
		    }
		    exp1 		<- exp(p1$fit)
                overallSglm1  <- exp1/(1+exp1)		    
		    exp2 		<- exp(p2$fit)
                overallSglm2  <- exp2/(1+exp2)
		    overalllogRatio <- coef(m1all)["txassign"]
		  } else
		  if (glm == "poisson"){
	  	    m1all         <- glm(fmla, family=poisson(link="log"))

		    if (!is.null(covariate)){
		      p1   		<- predict(m1all, newdata=data.frame(t(nv1)), se.fit=TRUE)
		      p2   		<- predict(m1all, newdata=data.frame(t(nv2)), se.fit=TRUE)
		    } else {
		      p1  		<- predict(m1all, newdata=data.frame(txassign=1), se.fit=TRUE)
			p2  		<- predict(m1all, newdata=data.frame(txassign=0), se.fit=TRUE)
		    }
	  	    overallSglm1  <- exp(p1$fit)
	  	    overallSglm2  <- exp(p2$fit)
		    overalllogRatio <- coef(m1all)["txassign"]
		  }

      	  if (sum(is.na(sglm1)) == 0 & sum(is.na(sglm2)) == 0 & is.na(overallSglm1) == FALSE & is.na(overallSglm2) == FALSE) {
        		no <- no + 1
        		p <- p + 1
        		for (s in 1:nsubpop) {
          	  		differences[p, s] <- (sglm1[s] - sglm2[s]) - (overallSglm1 - overallSglm2)
				diffha[p,s]		<- (sglm1[s] - sglm2[s]) - sglmwha
          	  	  	logratios[p,s]    <- slogratios[s] - overalllogRatio
			  	logratioha[p,s]   <- slogratios[s] - slogrglmha
          	  	}
      	  }

      	  terminate <- terminate + 1
      	  if (terminate >= nperm + 10000) {
        		print(paste("After permuting ", nperm, "plus 10000, or ", 
        	      nperm + 10000, " times, the program is unable to generate the permutation distribution based on ", 
        	      nperm, "permutations of the data"))
        		print(paste("Consider creating larger subpopulations or selecting a different timepoint for estimation"))
        		stop()
      	  }
    	      }
	      if (showstatus) close(pb)

	        # generating the sigmas and pvalues
	        sigmas       <- ssigma(differences)
	        sigma1	   <- sigmas$sigma
	        chi2pvalue   <- ppv (differences, sigmas$sigmainv, effect$sObs1-effect$sObs2, effect$oObs1-effect$oObs2, nperm)
	        pvalue       <- ppv2(differences,                  effect$sObs1-effect$sObs2, effect$oObs1-effect$oObs2, nperm)

	        sigmas       <- ssigma(diffha)
	        hasigma      <- sigmas$sigma
  	        hapvalue     <- ppv (diffha,      sigmas$sigmainv, effect$sObs1-effect$sObs2, effect$sglmw, nperm)

	        sigmas       <- ssigma(logratios)
	        logRsigma	   <- sigmas$sigma
	        logRpvalue   <- ppv2(logratios,                    effect$logR, 			  effect$ologR, nperm, debug=FALSE)

	        sigmas       <- ssigma(logratioha)
	        halogRsigma  <- sigmas$sigma
  	        halogRpvalue <- ppv (logratioha,  sigmas$sigmainv, effect$logR, 			  effect$sglmlogrw, nperm)

	        test = list( model	     = "GLMt",
				   sigma	     = sigma1,
				   hasigma	     = hasigma,
				   logRsigma     = logRsigma,
				   halogRsigma   = halogRsigma,
				   pvalue	     = pvalue,
				   chi2pvalue    = chi2pvalue,
				   hapvalue	     = hapvalue,
				   logRpvalue    = logRpvalue,
				   halogRpvalue  = halogRpvalue
				  )
		}

		return(test)

	  }
)

print.estimate.GLM <- function(stobj, family){
	  model <- stobj@model
	  sp	  <- stobj@subpop
	  est	  <- stobj@effect

	  if (family == "gaussian"){
          cat("\n")
          write(paste("Effect estimates for treatment group", model@trts[1]), file = "")
          temp <- matrix(c(1:sp@nsubpop, round(est$sObs1, digits = 4), round(est$sSE1, 
            digits = 4)), ncol = 3)
          write("     Subpopulation     Effect		 Std. Err.", 
            file = "")
          for (i in 1:sp@nsubpop) {
            write(paste(format(temp[i, 1], width = 12),
				format(temp[i, 2], width = 19, nsmall = 4),
				format(temp[i, 3], width = 15, nsmall = 4)),
			 	file = "")
          }
          write(paste("        Overall",
				format(round(est$oObs1, digits = 4), nsmall = 4, width = 16),
				format(round(est$oSE1,  digits = 4), nsmall = 4, width = 15)),
				file = "")
          cat("\n")
          write(paste("Effect estimates for treatment group", model@trts[2]), file = "")
          temp <- matrix(c(1:sp@nsubpop, round(est$sObs2, digits = 4),
					round(est$sSE2, digits = 4)), ncol = 3)
           write("     Subpopulation     Effect            Std. Err.", 
            file = "")
          for (i in 1:sp@nsubpop) {
            write(paste(format(temp[i, 1], width = 12),
				format(temp[i, 2], width = 19, nsmall = 4),
				format(temp[i, 3], width = 15, nsmall = 4)),
				file = "")
          }
          write(paste("        Overall",
				format(round(est$oObs2, digits = 4), nsmall = 4, width = 16), 
				format(round(est$oSE2,  digits = 4), nsmall = 4, width = 15)),
				file = "")
          cat("\n")

          write(paste("Effect differences"), file = "")
          temp <- matrix(c(1:sp@nsubpop, round(est$sObs1 - est$sObs2, digits = 4), 
            round(sqrt(est$sSE1^2 + est$sSE2^2), digits = 4)), ncol = 3)
          write("                          Effect", file = "")
          write("     Subpopulation      Difference      Std. Err.", 
            file = "")
          for (i in 1:sp@nsubpop) {
            write(paste(format(temp[i, 1], width = 12),
				format(temp[i, 2], width = 19, nsmall = 4),
				format(temp[i, 3], width = 15, nsmall = 4)),
				file = "")
          }
          write(paste("        Overall",
			format(round(est$oObs1 - est$oObs2, digits = 4), nsmall = 4, width = 16), 
            	format(round(sqrt(est$oSE1^2 + est$oSE2^2), digits = 4), nsmall = 4, width = 15)),
			file = "")
          cat("\n")

          write(paste("Effect ratios"), file = "")
          temp <- matrix(c(1:sp@nsubpop, round(est$logR, digits = 4), 
            		round(est$logRSE, digits = 4), round(exp(est$logR), digits=4)), ncol = 4)	
          write("                          Effect", file = "")
          write("     Subpopulation      log Effect Ratio  Std. Err.	    Effect Ratio", file = "")
          for (i in 1:sp@nsubpop) {
            write(paste(format(temp[i, 1], width = 12, ),
				format(round(temp[i, 2], digits = 4), width = 19, nsmall = 4),
				format(round(temp[i, 3], digits = 4), width = 15, nsmall = 4),
				format(round(temp[i, 4], digits = 4), width = 19, nsmall = 4)),
				file = "")
          }
          write(paste("        Overall",
		format(round(est$ologR,      digits = 4),  nsmall = 4, width = 16), 
            format(round(est$ologRSE,    digits = 4),  nsmall = 4, width = 15),
		format(round(exp(est$ologR), digits = 4),  nsmall = 4, width = 19)),
		file = "")
          cat("\n")

        }
	  else
	  if (family == "binomial"){
          cat("\n")
          write(paste("Risk estimates for treatment group", model@trts[1]), file = "")
          temp <- matrix(c(1:sp@nsubpop, round(est$sObs1, digits = 4), round(est$sSE1, digits = 4)), ncol = 3)
          write("     Subpopulation     Risk		 Std. Err.", file = "")
          for (i in 1:sp@nsubpop) {
            write(paste(format(temp[i, 1], width = 12),
				format(temp[i, 2], width = 19, nsmall = 4),
				format(temp[i, 3], width = 15, nsmall = 4)),
			file = "")
          }
          write(paste("        Overall", 
				format(round(est$oObs1, digits = 4), nsmall = 4, width = 16),
				format(round(est$oSE1,  digits = 4), nsmall = 4, width = 15)),
			file = "")
          cat("\n")
          write(paste("Risk estimates for treatment group", model@trts[2]), file = "")
          temp <- matrix(c(1:sp@nsubpop, round(est$sObs2, digits = 4), round(est$sSE2, digits = 4)), ncol = 3)
          write("     Subpopulation     Risk           Std. Err.", file = "")
          for (i in 1:sp@nsubpop) {
            write(paste(format(temp[i, 1], width = 12),
				format(temp[i, 2], width = 19, nsmall = 4), 
				format(temp[i, 3], width = 15, nsmall = 4)),
			file = "")
          }
          write(paste("        Overall",
				format(round(est$oObs2, digits = 4), nsmall = 4, width = 16),
				format(round(est$oSE2,  digits = 4), nsmall = 4, width = 15)),
			file = "")
          cat("\n")

          write(paste("Risk differences"), file = "")
          temp <- matrix(c(1:sp@nsubpop, round(est$sObs1 - est$sObs2, digits = 4), 
            		round(sqrt(est$sSE1^2 + est$sSE2^2), digits = 4)), ncol = 3)
          write("                          Risk", file = "")
          write("     Subpopulation      Difference      Std. Err.", 
            file = "")
          for (i in 1:sp@nsubpop) {
            write(paste(format(temp[i, 1], width = 12), 
				format(temp[i, 2], width = 19, nsmall = 4),
				format(temp[i, 3], width = 15, nsmall = 4)),
			file = "")
          }
          write(paste("        Overall",
				format(round(est$oObs1 - est$oObs2, digits = 4), nsmall = 4, width = 16), 
            		format(round(sqrt(est$oSE1^2 + est$oSE2^2), digits = 4), nsmall = 4, width = 15)),
			file = "")
          cat("\n")

          write(paste("Odds Ratio"), file = "")
          temp <- matrix(c(1:sp@nsubpop, 
					round(est$logR,      digits = 4), 
            			round(est$logRSE,    digits = 4), 
					round(exp(est$logR), digits = 4)), 
				ncol = 4)
          write("                        log Odds", file = "")
          write("     Subpopulation        Ratio          Std. Err.	    Odds Ratio", file = "")
          for (i in 1:sp@nsubpop) {
            write(paste(format(temp[i, 1], width = 12), 
				format(round(temp[i, 2], digits = 4), width = 19, nsmall = 4), 
				format(round(temp[i, 3], digits = 4), width = 15, nsmall = 4),
				format(round(temp[i, 4], digits = 4), width = 19, nsmall = 4)), 
			file = "")
          }
          write(paste("        Overall",
				format(round(est$ologR,     digits = 4), nsmall = 4, width = 16), 
            		format(round(est$ologRSE,   digits = 4), nsmall = 4, width = 15),
				format(round(exp(est$ologR), digits = 4), nsmall = 4, width = 19)), 
			file = "")
          cat("\n")

        } else
	  if (family == "poisson"){
          cat("\n")
          write(paste("    Effect estimates for treatment group", model@trts[1]), file = "")
          temp <- matrix(c(1:sp@nsubpop, round(est$sObs1, digits = 4), round(est$sSE1, digits = 4)), ncol = 3)
          write("     Subpopulation        Effect	Std. Err.", file = "")
          for (i in 1:sp@nsubpop) {
            write(paste(format(temp[i, 1], width = 12), 
				format(temp[i, 2], width = 19, nsmall = 4), 
				format(temp[i, 3], width = 15, nsmall = 4)), 
			file = "")
          }
          write(paste("        Overall", 
				format(round(est$oObs1, digits = 4), nsmall = 4, width = 16),
			 	format(round(est$oSE1,  digits = 4), nsmall = 4, width = 15)), 
			file = "")
          cat("\n")
          write(paste("    Effect estimates for treatment group", model@trts[2]), file = "")
          temp <- matrix(c(1:sp@nsubpop, round(est$sObs2, digits = 4), round(est$sSE2, digits = 4)), ncol = 3)
           write("     Subpopulation        Effect      Std. Err.", file = "")
          for (i in 1:sp@nsubpop) {
            write(paste(format(temp[i, 1], width = 12), 
				format(temp[i, 2], width = 19, nsmall = 4), 
				format(temp[i, 3], width = 15, nsmall = 4)), 
			file = "")
          }
          write(paste("        Overall", 
				format(round(est$oObs2, digits = 4), nsmall = 4, width = 16), 
				format(round(est$oSE2,  digits = 4), nsmall = 4, width = 15)), 
			file = "")
          cat("\n")

          write(paste("Effect differences"), file = "")
          temp <- matrix(c(1:sp@nsubpop, 
				round(est$sObs1 - est$sObs2, digits = 4), 
            		round(sqrt(est$sSE1^2 + est$sSE2^2), digits = 4)), ncol = 3)
          write("                         Effect ", file = "")
          write("     Subpopulation      Difference      Std. Err.", 
            file = "")
          for (i in 1:sp@nsubpop) {
            write(paste(format(temp[i, 1], width = 12), 
				format(temp[i, 2], width = 19, nsmall = 4), 
				format(temp[i, 3], width = 15, nsmall = 4)), 
			file = "")
          }
          write(paste("        Overall", 
				format(round(est$oObs1 - est$oObs2, digits = 4), nsmall = 4, width = 16), 
            		format(round(sqrt(est$oSE1^2 + est$oSE2^2), digits = 4), nsmall = 4, width = 15)), 
			file = "")
          cat("\n")

          write(paste("Relative    Effect "), file = "")
          temp <- matrix(c(1:sp@nsubpop,  
				round(est$logR,      digits = 4), 
				round(est$logRSE,    digits = 4),
				round(exp(est$logR), digits = 4)), 
				ncol = 4)
          write("                           log", file = "")
          write("     Subpopulation        Effect         Std. Err.        Relative Effect", 
            file = "")
          for (i in 1:sp@nsubpop) {
            write(paste(format(temp[i, 1], width = 12), 
				format(temp[i, 2], width = 19, nsmall = 4), 
				format(temp[i, 3], width = 15, nsmall = 4),
				format(temp[i, 4], width = 19, nsmall = 4)),
			 file = "")
          }
          write(paste("        Overall", 
				format(round(est$ologR,      digits = 4), nsmall = 4, width = 16), 
            		format(round(est$ologRSE,    digits = 4), nsmall = 4, width = 15),
				format(round(exp(est$ologR), digits = 4), nsmall = 4, width = 19)), 
			file = "")
          cat("\n")
	  }
}

print.cov.GLM <- function(stobj){
	write(paste("The covariance matrix of the effect differences estimates for the",
			stobj@subpop@nsubpop, "subpopulations is:"), file = "")
      print(stobj@result$sigma)

	cat("\n")
	write(paste("The covariance matrix of the log effect ratios for the", 
          		stobj@subpop@nsubpop, "subpopulations is:"), file = "")
      print(stobj@result$logRsigma)
      cat("\n")

	write(paste("The covariance matrix (based on homogeneous association) of the effect differences for the", 
           		stobj@subpop@nsubpop, "subpopulations is:"), file = "")
	print(stobj@result$hasigma)
      
	cat("\n")
	write(paste("The covariance matrix (based on homogeneous association) of the log effect ratios for the", 
          		stobj@subpop@nsubpop, "subpopulations is:"), file = "")
      print(stobj@result$halogRsigma)
      cat("\n")

}

print.stat.GLM <- function(pvalue, chi2pvalue, hapvalue, Rpvalue, haRpvalue){

	write(paste("Interaction P-value based on effect differences estimates :", pvalue), file = "")

      cat("\n")
	write(paste("Chisquare interaction P-value based on effect differences estimates :", chi2pvalue), file = "")

      cat("\n")
      write(paste("Homogeneous association interaction P-value based on effect differences estimates :", hapvalue), file = "")

      cat("\n")
	write(paste("Interaction P-value based on effect ratio estimates :", Rpvalue), file = "")

      cat("\n")
      write(paste("Homogeneous association interaction P-value based on effect ratio estimates :", haRpvalue), file = "")

	cat("\n")


}


setMethod("print",
	    signature="stmodelGLM",
	    definition=function(x, stobj, estimate=TRUE, cov=TRUE, test=TRUE, ...){

		#
		#  1. estimates
		#
	      if (estimate){
		  print.estimate.GLM(stobj, x@glm)
      	}

		if (!is.null(stobj@result)){
  		  #
		  #   2. covariance matrices
		  #
    		  if (cov){
		    print.cov.GLM(stobj)
 		  }
  
		  #
		  #   3. Supremum test and Chi-square test results
		  #
		  if (test){
 	          t 	  <- stobj@result
		    print.stat.GLM(t$pvalue, t$chi2pvalue, t$hapvalue, t$logRpvalue, t$halogRpvalue)
		  }
		}
 	 }
)

# constructor function for stmodelGLM
stepp.GLM <- function(coltrt, trts, colY, MM=NULL, glm, debug=0){
	model <- new("stmodelGLM", coltrt=coltrt, trts=trts, colY=colY, MM=MM, glm=glm, debug=debug)
	return(model)
}


