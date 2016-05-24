

## FUNCTIONS FOR TURNING DATA INTO GROWTH AND SURVIVAL OBJECTS ############################
## uses linear and logistic regressions, with various polynomials
# dataf must have columns "size", "sizeNext", "surv",
# and facultatively, "covariate" and "covariatel" (if single discrete covariate, like light)
# or covariate1, covariate2, covariate3, covariate4,... (if continuous and discrete)
# and "fec", "age";  age is for picking out seedling sizes for full IPM


# =============================================================================
# =============================================================================
## 1. Growth  models  #############################################


#More general function for growth than older versions listed below
#
# Parameters - dataf - the data-frame (which must contain the headings found in the formulas)
#              formula= - a model formula that requires
#                    "sizeNext" or "incr" as a reponse variable
#                    "size" as a possible covariate possibly with
#                        other combinations including size2 (size^2), size3(size^3), logsize(log(size)), expsize(exp(size))
#                        and potentially a discrete covariate (called covariate)
#              responseType - the response variable desired, crucial for building
#                             the right kind of object. Possible levels are "sizeNext", "incr", "logincr"
#              regType - options are "constantVar" (i.e. use lm), changingVar (use gls with model for variance)
# 
# Returns - a growth object                  
#
#

makeGrowthObj <- function(dataf=NULL,
		Formula=sizeNext~size,
		regType="constantVar",
		Family="gaussian", link=NULL, coeff=NULL, sd=NULL) {
	
	if (!is.null(dataf)) { 
	
	dataf <- subset(dataf, is.na(dataf$size) == FALSE)
	
	if (length(dataf$sizeNext)>0) dataf <- subset(dataf, is.na(dataf$sizeNext) == FALSE)
	
	
	if (length(dataf$offspringNext) > 0) 
		dataf <- subset(dataf, !dataf$offspringNext %in% c("sexual", 
						"clonal"))
	
	if (length(grep("incr", as.character(Formula))) > 0 & length(dataf$incr) == 0) {
		print("building incr as sizeNext - size")
		dataf$incr <- dataf$sizeNext - dataf$size
	}
	if (length(grep("logincr", as.character(Formula))) > 0 & length(dataf$logincr) == 0) {
		print("building logincr as log(sizeNext - size) - pre-build if this is not appropriate")
		dataf$logincr <- log(dataf$sizeNext - dataf$size)
	}
	
	#eliminate growth in dead individual
	if (length(grep("incr", as.character(Formula))) > 0) {
		if (sum(!is.na(dataf$incr) & dataf$surv==0,na.rm=TRUE)>0) {
		print("measures of growth exist where individual has died (surv==0); replacing these with NA")
		dataf$incr[dataf$surv==0] <- NA
	}}
	if (length(grep("sizeNext", as.character(Formula))) > 0) { 
		if(sum(!is.na(dataf$sizeNext) & dataf$surv==0,na.rm=TRUE)>0) {
		print("measures of growth exist where individual has died (surv==0); replacing these with NA")		
		dataf$sizeNext[dataf$surv==0] <- NA
	}}
	if (length(grep("logincr", as.character(Formula))) > 0) {
		if (sum(!is.na(dataf$logincr) & dataf$surv==0,na.rm=TRUE)>0) {
		print("measures of growth exist where individual has died (surv==0); replacing these with NA")		
		dataf$logincr[dataf$surv==0] <- NA
	}}
	
	#create appropriate size based covariates
	dataf$size2 <- dataf$size ^ 2
	dataf$size3 <- dataf$size ^ 3
	if (length(grep("expsize", as.character(Formula))) > 0) dataf$expsize <- exp(dataf$size)
	if (length(grep("logsize", as.character(Formula))) > 0) dataf$logsize <- log(dataf$size)
	
	#setup for discrete covariates if data suggests may be implemented by the
	#presence of "covariate" and "covariateNext"
	if ("covariate" %in% unlist(strsplit(as.character(Formula), "[+-\\* ]")) & length(dataf$covariate) > 0) { 
		dataf$covariate <- as.factor(dataf$covariate)
		levels(dataf$covariate) <- 1:length(unique(dataf$covariate))
	}
	if ("covariateNext" %in% unlist(strsplit(as.character(Formula), "[+-\\* ]")) & length(dataf$covariateNext) > 0) { 
		dataf$covariateNext <- as.factor(dataf$covariateNext)
		levels(dataf$covariateNext) <- 1:length(unique(dataf$covariateNext))
	}
	
	if (length(intersect(all.vars(Formula),colnames(dataf)))<length(all.vars(Formula))) print("warning: not all variables in the formula are present in dataf; model cannot be fit")
	
	#eval fit and make the objects
	if (Family=="gaussian") { 
		if (regType == "constantVar")  {
			fit <- lm(Formula, data=dataf)
		} else { 
			if (regType == "changingVar"){
				fit.here <- gls(Formula, na.action = na.omit, weights = varExp(form =  ~fitted(.)), data = dataf)
				fit <- list(coefficients = fit.here$coefficients,
						sigmax2 = summary(fit.here)$sigma^2,
						var.exp.coef = as.numeric(fit.here$modelStruct$varStruct[1]), 
						fit = fit.here)
				#print(class(fit.here))
			}
		}
	} else {
		if (regType != "constantVar") print("Warning: your regType is ignored because a non-gaussian model is fitted using glm")
		if (Family=="negbin"){
			if (is.null(link)) { print("setting link to identity"); link <- "identity" }
			if (link=="identity") fit <- glm.nb(Formula, data=dataf,link="identity")
			if (link=="log") fit <- glm.nb(Formula, data=dataf,link="log")
			#if (link=="sqrt") fit <- glm.nb(Formula, data=dataf,link="sqrt")
			if (link!="identity" & link!="log"){stop("unknown link specified for negative binomial")}
			fit.here <- list()
			if (link=="log") fit.here[[1]] <- glm.convert(fit)
			if (link=="identity") { 
				fit.dummy <-  glm(Formula, data=dataf)
				fit.dummy$coefficients <- fit$coefficients
				fit.here[[1]] <- fit.dummy
			}
			fit.here[[2]] <- fit$theta
			fit.here[[3]] <- fit  
		} else {
			fit <- glm(Formula, data=dataf, family=Family)
			fit.here <- fit
			#print("here")
		}           
	}
	
	#make the objects
	#with sizeNext as response
	if (length(grep("sizeNext", as.character(Formula))) > 0) { 
		if (class(fit)[1] == "lm") { 
			gr1 <- new("growthObj")
			gr1@fit <- fit
			gr1@sd <- summary(fit)$sigma
		} else {
			if (class(fit.here)[1] == "gls") { 
				gr1 <- new("growthObjDeclineVar")
				gr1@fit <- fit
			} else {
				if (class(fit)[1] == "glm") { 
					if (Family=="poisson") { gr1 <- new("growthObjPois"); gr1@fit <- fit } else {print("unidentified object class")}
				} else {
					if (class(fit)[1] == "negbin") {
						gr1 <- new("growthObjNegBin"); gr1@fit <- fit.here
					} 
				}
			}
		}    
	} else {
		if (length(grep("incr", as.character(Formula))) > 0 & 
				length(grep("logincr", as.character(Formula))) == 0) { 
			
			if (class(fit)[1] == "lm") { 
				gr1 <- new("growthObjIncr")
				gr1@fit <- fit
				gr1@sd <- summary(fit)$sigma
			} else {
				if (class(fit.here)[1] == "gls") { 
					gr1 <- new("growthObjIncrDeclineVar")
					gr1@fit <- fit
				} else {print("undefined object class")}
			}
		} else {
			if (length(grep("logincr", as.character(Formula))) > 0) { 					
				if (class(fit)[1] == "lm") { 
					gr1 <- new("growthObjLogIncr")
					gr1@fit <- fit
					gr1@sd <- summary(fit)$sigma
				} else {
					if (class(fit.here)[1] == "gls") { 
						gr1 <- new("growthObjLogIncrDeclineVar")
						gr1@fit <- fit
					} else {print("undefined object class")}
				}
			}
		}    
	}
	
	} else {
		
		if (is.null(coeff) | is.null(sd)) stop("require coefficients and standard devaition if data is not supplied")
		gr1 <- .createGrowthObj(Formula=Formula, coeff=coeff, sd=sd)
		
	}
	
	
	
	return(gr1)
}


# =============================================================================
# =============================================================================
### IDENTICAL TO makeGrowthObj - JUST COPIED
##### ONLY DIFFERENCE IS DATA - don't throw away all NAs in size - and 
##### focus on ones where offspringNext is sexual or clonal

makeOffspringObj <- function(dataf=NULL,
		Formula=sizeNext~size,
		regType="constantVar",
		Family="gaussian", link=NULL, coeff=NULL, sd=NULL, reproductionType="sexual") {
	
	if (!is.null(dataf)) { 
			
		if (length(dataf$stage) > 0) 
			dataf <- subset(dataf,is.na(dataf$stage)&dataf$stageNext=="continuous")
		
		if (length(dataf$offspringNext) > 0) 
			dataf <- subset(dataf, dataf$offspringNext %in% reproductionType)
		
		#eliminate growth in dead individual
		if (length(grep("sizeNext", as.character(Formula))) > 0) { 
			if(sum(!is.na(dataf$sizeNext) & dataf$surv==0,na.rm=TRUE)>0) {
				print("measures of growth exist where individual has died (surv==0); replacing these with NA")		
				dataf$sizeNext[dataf$surv==0] <- NA
			}}
		
		#create appropriate size based covariates
		dataf$size2 <- dataf$size ^ 2
		dataf$size3 <- dataf$size ^ 3
		if (length(grep("expsize", as.character(Formula))) > 0) dataf$expsize <- exp(dataf$size)
		if (length(grep("logsize", as.character(Formula))) > 0) dataf$logsize <- log(dataf$size)
		
		#setup for discrete covariates if data suggests may be implemented by the
		#presence of "covariate" and "covariateNext"
		if ("covariate" %in% unlist(strsplit(as.character(Formula), "[+-\\* ]")) & length(dataf$covariate) > 0) { 
			dataf$covariate <- as.factor(dataf$covariate)
			levels(dataf$covariate) <- 1:length(unique(dataf$covariate))
		}
		if ("covariateNext" %in% unlist(strsplit(as.character(Formula), "[+-\\* ]")) & length(dataf$covariateNext) > 0) { 
			dataf$covariateNext <- as.factor(dataf$covariateNext)
			levels(dataf$covariateNext) <- 1:length(unique(dataf$covariateNext))
		}
		
		if (length(intersect(all.vars(Formula),colnames(dataf)))<length(all.vars(Formula))) print("warning: not all variables in the formula are present in dataf; model cannot be fit")
		
		#eval fit and make the objects
		if (Family=="gaussian") { 
			if (regType == "constantVar")  {
				fit <- lm(Formula, data=dataf)
			} else { 
				if (regType == "changingVar"){
					fit.here <- gls(Formula, na.action = na.omit, weights = varExp(form =  ~fitted(.)), data = dataf)
					fit <- list(coefficients = fit.here$coefficients,
							sigmax2 = summary(fit.here)$sigma^2,
							var.exp.coef = as.numeric(fit.here$modelStruct$varStruct[1]), 
							fit = fit.here)
					#print(class(fit.here))
				}
			}
		} else {
			if (regType != "constantVar") print("Warning: your regType is ignored because a non-gaussian model is fitted using glm")
			if (Family=="negbin"){
				if (is.null(link)) { print("setting link to identity"); link <- "identity" }
				if (link=="identity") fit <- glm.nb(Formula, data=dataf,link="identity")
				if (link=="log") fit <- glm.nb(Formula, data=dataf,link="log")
				#if (link=="sqrt") fit <- glm.nb(Formula, data=dataf,link="sqrt")
				if (link!="identity" & link!="log"){stop("unknown link specified for negative binomial")}
				fit.here <- list()
				if (link=="log") fit.here[[1]] <- glm.convert(fit)
				if (link=="identity") { 
					fit.dummy <-  glm(Formula, data=dataf)
					fit.dummy$coefficients <- fit$coefficients
					fit.here[[1]] <- fit.dummy
				}
				fit.here[[2]] <- fit$theta
				fit.here[[3]] <- fit  
			} else {
				fit <- glm(Formula, data=dataf, family=Family)
				fit.here <- fit
				#print("here")
			}           
		}
		
		#make the objects
		#with sizeNext as response
		if (class(fit)[1] == "lm") { 
				gr1 <- new("growthObj")
				gr1@fit <- fit
				gr1@sd <- summary(fit)$sigma
			} else {
				if (class(fit.here)[1] == "gls") { 
					gr1 <- new("growthObjDeclineVar")
					gr1@fit <- fit
				} else {
					if (class(fit)[1] == "glm") { 
						if (Family=="poisson") { gr1 <- new("growthObjPois"); gr1@fit <- fit } else {print("unidentified object class")}
					} else {
						if (class(fit)[1] == "negbin") {
							gr1 <- new("growthObjNegBin"); gr1@fit <- fit.here
						} 
					}
				}
			}    
				
		} else {
			
		if (is.null(coeff) | is.null(sd)) stop("require coefficients and standard devaition if data is not supplied")
		gr1 <- .createGrowthObj(Formula=Formula, coeff=coeff, sd=sd)
	}
	
	return(gr1)
}

### Function to make a new object of class offspringObj

# =============================================================================
# =============================================================================
## Function to create a new Hossfeld growth object
#
# Parameters - dataf - a dataframe
#
# Returns - a Hossfeld growth object

#no covariate, and one polynom, linear regression
makegrowthObjHossfeld <- function(dataf) {  

	#subset data to include only growth of individuals
	dataf<-subset(dataf,is.na(dataf$size)==FALSE&is.na(dataf$sizeNext)==FALSE)
	if (length(dataf$offspringNext)>0) dataf<-subset(dataf,!dataf$offspringNext%in%c("sexual","clonal"))
	
	if (length(dataf$incr)==0) dataf$incr <- dataf$sizeNext-dataf$size
	dataf$incr[dataf$incr<0] <- 0
	tmp <- optim(c(1, 1, 1), wrapHossfeld, dataf = dataf, method = "Nelder-Mead")
	print(tmp$convergence)
	gr1 <- new("growthObjHossfeld")
	gr1@paras <- tmp$par
	resids <- Hossfeld(dataf$size, tmp$par) - dataf$incr 
	gr1@sd <- sd(resids, na.rm = T)
	return(gr1)
}


# 2. Survival models  #########################################################

# =============================================================================
# =============================================================================
#General function for survival 
#
# Parameters - dataf - the data-frame (which must contain the headings found in the formulas)
#              formula - a model formula that requires
#                    "surv" as a reponse variable
#                    "size" as a possible covariate possibly with
#                        other combinations including size2 (size^2), size3(size^3), logsize(log(size)), expsize(exp(size))
#                        and potentially a discrete covariate (called covariate)
# 
# Returns - a survival object                   
#
#
makeSurvObj <- function(dataf=NULL,
		Formula=surv~size+size2, coeff=NULL){
	
	if (!is.null(dataf)){	
	
	#subset data to include only survival status of individuals with continuous size at the beginning of the transition
	dataf<-subset(dataf,is.na(dataf$surv)==FALSE)
	if (length(dataf$offspringNext)>0) dataf<-subset(dataf,!dataf$offspringNext%in%c("sexual","clonal"))
	
	#build appropriate size based covariates
	dataf$size2 <- dataf$size^2
	dataf$size3 <- dataf$size^3
	if (length(grep("expsize",as.character(Formula)))>0) dataf$expsize <- exp(dataf$size)
	if (length(grep("logsize",as.character(Formula)))>0) dataf$logsize <- log(dataf$size)
	

	#setup for discrete covariates if data suggests may be implemented by the
	#presence of "covariate" and "covariateNext"
	if ("covariate" %in% unlist(strsplit(as.character(Formula), "[+-\\* ]")) & length(dataf$covariate) > 0) { 
		dataf$covariate <- as.factor(dataf$covariate)
		levels(dataf$covariate) <- 1:length(unique(dataf$covariate))
	}
	if ("covariateNext" %in% unlist(strsplit(as.character(Formula), "[+-\\* ]")) & length(dataf$covariateNext) > 0) { 
		dataf$covariateNext <- as.factor(dataf$covariateNext)
		levels(dataf$covariateNext) <- 1:length(unique(dataf$covariateNext))
	}
	
	if (length(intersect(all.vars(Formula),colnames(dataf)))<length(all.vars(Formula))) print("warning: not all variables in the formula are present in dataf; model cannot be fit")
	
	
	fit <- glm(Formula,family=binomial,data=dataf)
	sv1 <- new("survObj")
	sv1@fit <- fit
	
	} else {
		
		if (is.null(coeff)) stop("require coefficients if data is not supplied")
		sv1 <- .createSurvObj(Formula=Formula, coeff=coeff)
		
	}
	
	return(sv1)
}



# 3. fecundity models  ########################################################

# =============================================================================
# =============================================================================
makeFecObj <- function(dataf=NULL,
		fecConstants=data.frame(NA),
		Formula=list(fec~size),
		Family="gaussian",
		Transform="none",
		meanOffspringSize=NA,
		sdOffspringSize=NA,
		offspringSplitter=data.frame(continuous=1),
		vitalRatesPerOffspringType=data.frame(NA),
		fecByDiscrete=data.frame(NA),
		offspringSizeExplanatoryVariables="1", 
		coeff=NULL, doOffspring=TRUE, 
		reproductionType="sexual"){
		
	if (!is.null(dataf)) { 
	
	#make sure Formula is a formula or a list of formulas
	if (class(Formula)=="list") {
		if (class(Formula[[1]])!="formula") stop("Error - the entries in your Formula list should be of class 'formula': e.g. fec~size without quotation marks")
	} else {
		if (class(Formula)!="formula") stop("Error - the Formula entry should by of class 'formula' or a list of such entries:  e.g. fec~size without quotation marks")
        Formula <- list(Formula)	
	}
			
	#if stage or stageNext do not exist in dataf, create them assuming that 
	#everything is continuous. 
	if (length(dataf$stage)==0) { 
		dataf$stage <- rep("continuous",nrow(dataf))
		dataf$stage[is.na(dataf$size)] <- NA
		dataf$stage <- as.factor(dataf$stage)
	}
	if (length(dataf$stageNext)==0) {
		dataf$stageNext <- rep("continuous",nrow(dataf))
		dataf$stageNext[dataf$surv==0] <- "dead"
		dataf$stageNext <- as.factor(dataf$stageNext)
	}
	
	#order stage names from discrete to continuous
	stages <- names(tapply(c(levels(dataf$stage),levels(dataf$stageNext)),c(levels(dataf$stage),levels(dataf$stageNext)),length))
	stages <- stages[stages!="dead"] 
	if ((sum(names(offspringSplitter)%in%stages)/length(offspringSplitter))<1) {
		stages <- c(stages,names(offspringSplitter))
		print("Warning - the variable names in your offspringSplitter data.frame are not all part of the levels of stage or stageNext in your data file. This could be because of an mismatch in stage names, or because you included discrete stages in offspringSplitter that are not in the data file but wchich you will introduce in makeDiscreteTrans (in which case you can ignore this warning).")
	}
	stages <- unique(stages)
	stages <- c(stages[stages!="continuous"],"continuous") 
	dummy<-rep(0,length(stages));names(dummy)<-stages;dummy<-as.data.frame(t(as.matrix(dummy)))
	for (i in names(offspringSplitter)) dummy[i]<-offspringSplitter[i]
	offspringSplitter <- dummy
	
	##warnings
	if (ncol(offspringSplitter)>1 & (ncol(offspringSplitter)-1)!=ncol(fecByDiscrete)) {
		print("Warning - offspring splitter indicates more than just continuous stages. No fecundity by the discrete stages supplied in fecByDiscrete; assumed that is 0")
		fecByDiscrete <- offspringSplitter[,1:(ncol(offspringSplitter)-1)]
		fecByDiscrete[] <- 0
	}
	
	if (sum(offspringSplitter)!=1) {
		print("Warning - offspring splitter does not sum to 1. It is now rescaled to sum to 1.")
		offspringSplitter <- offspringSplitter / sum(offspringSplitter) 
	}
	
	if ("covariate" %in% unlist(strsplit(as.character(Formula), "[+-\\* ]")) & length(dataf$covariate) > 0) { 
		dataf$covariate <- as.factor(dataf$covariate)
		levels(dataf$covariate) <- 1:length(unique(dataf$covariate))
	}
	if ("covariateNext" %in% unlist(strsplit(as.character(Formula), "[+-\\* ]")) & length(dataf$covariateNext) > 0) { 
		dataf$covariateNext <- as.factor(dataf$covariateNext)
		levels(dataf$covariateNext) <- 1:length(unique(dataf$covariateNext))
	}
	
	f1 <- new("fecObj")
	
	dataf$size2 <- dataf$size^2
	dataf$size3 <- dataf$size^3
	if (length(grep("expsize",unlist(as.character(Formula))))>0) dataf$expsize <- exp(dataf$size)
	if (length(grep("logsize",unlist(as.character(Formula))))>0) dataf$logsize <- log(dataf$size)
	
	if (length(Formula)>length(Family)) {
		misE <- (length(Family)+1):length(Formula)
		print(c("number of families not the same as the number of Formula supplied, using default of `gaussian' for missing ones which are:",Formula[[misE]],". (which might be exactly what you want)"))
		Family <- c(Family,rep("gaussian",length(Formula)-length(Family)))
	}
	if (length(Formula)>length(Transform)) {
		misE <- (length(Transform)+1):length(Formula)
		print(c("number of transforms not the same as the number of fecundity columns in the data file, using default of `none' for missing ones which are:",Formula[[misE]],". (which might be exactly what you want)"))
		Transform <- c(Transform,rep("none",length(Formula)-length(Transform)))
	}
	
	fecNames <- rep(NA,length(Formula))
	for (i in 1:length(Formula)) {
		
		fecNames[i] <- all.vars(Formula[[i]])[1]
		
		if (Transform[i]=="exp") dataf[,fecNames[i]] <- exp(dataf[,fecNames[i]])
		if (Transform[i]=="log") dataf[,fecNames[i]] <- log(dataf[,fecNames[i]])
		if (Transform[i]=="sqrt") dataf[,fecNames[i]] <- sqrt(dataf[,fecNames[i]])
		if (Transform[i]=="-1") dataf[,fecNames[i]] <- dataf[,fecNames[i]]-1
		dataf[!is.finite(dataf[,fecNames[i]]),fecNames[i]] <- NA
		if (length(intersect(all.vars(Formula[[i]]),colnames(dataf)))<length(all.vars(Formula[[i]]))) print("warning: not all variables in the formula are present in dataf; model cannot be fit")
		f1@fitFec[[i]] <- glm(Formula[[i]],family=Family[i],data=dataf)
	}
	
	if (doOffspring) { 
	if (offspringSplitter$continuous>0) {
		if (is.na(meanOffspringSize[1])|is.na(sdOffspringSize[1])) {
			if (length(dataf$offspringNext)==0) {
				offspringData<-subset(dataf,is.na(dataf$stage)&dataf$stageNext=="continuous")
				if (nrow(offspringData) == 0) stop(paste("Error - no offspring size data are given: these can be given through either the meanOffspringSize and sdOffspringSize slots, or through individual data added to your data file (with stage equals NA, or an offspringNext column indicating ",reproductionType, " offspring)", sep=""))
			} else {
				offspringData<-subset(dataf,dataf$offspringNext==reproductionType &dataf$stageNext=="continuous")
				if (nrow(offspringData) == 0) stop(paste("Error - no offspring size data are given: these can be given through either the meanOffspringSize and sdOffspringSize slots, or through individual data added to your data file (with stage equals NA, or a offspringNext column indicating ", reproductionType, " offspring)",sep=""))
			}
			## relationship defining offspring size - note that the mean value is ALWAYS taken from
			## a lm now (but equivalent to just fitting an intercept if that is desired....)
			## [worth keeping sd separate from lm though (extracted from lm or not) because otherwise is a pain to adjust (as shown in growth model)]
			
		
			f1@offspringRel <- lm(paste('sizeNext~',offspringSizeExplanatoryVariables,sep=''),data=offspringData)
			f1@sdOffspringSize <- summary(f1@offspringRel)$sigma
		} else {
			f1@offspringRel <- lm(rep(meanOffspringSize[1],21)~1)
			f1@sdOffspringSize <- sdOffspringSize
		}
	}	
	} #end of doOffspring loop
	
	if (sum(dim(vitalRatesPerOffspringType)==c(1,1))<2) {
		if ((sum(vitalRatesPerOffspringType==0,na.rm=T)+sum(vitalRatesPerOffspringType==1,na.rm=T))<(ncol(vitalRatesPerOffspringType)*nrow(vitalRatesPerOffspringType))) stop("Error - in vitalRatesPerOffspringType data.frame only 0's and 1's are allowed: a 1 indicates that a fecundity rate applies to that offspring type. ")
		#if (sum(names(vitalRatesPerOffspringType)==names(offspringSplitter))<length(offspringSplitter)) stop("Error - the offspring names in vitalRatesPerOffspringType should match those in offspringSplitter - and in the same order, with continuous last")
		if (sum(rownames(vitalRatesPerOffspringType)==c(fecNames,names(fecConstants)))<(length(Formula)+length(fecConstants))) stop ("Error - the row names in vitalRatesPerOffspringType should consist of (in order) the names of the fec columns in the dataset and then the names of the fecConstants.")
	} else {
		vitalRatesPerOffspringType <- as.data.frame(matrix(1,ncol=length(offspringSplitter),nrow=length(Formula)+length(fecConstants)),row.names=c(fecNames,names(fecConstants)))
		vitalRatesPerOffspringType <- subset(vitalRatesPerOffspringType,dimnames(vitalRatesPerOffspringType)[[1]]!="NA.")
		names(vitalRatesPerOffspringType) <- names(offspringSplitter)
	}
	
	if (doOffspring & length(f1@sdOffspringSize)>0) {
		if (is.na(f1@sdOffspringSize)) {
			print("Warning - could not estimate parameters for the distribution of offspring size; defaults must be supplied for meanOffspringSize and sdOffspringSize; you will not be able to construct an IPM without these values.")
		}
	}
	
	f1@fecNames <- fecNames
	f1@fecConstants <- fecConstants
	f1@offspringSplitter <- offspringSplitter 
	f1@vitalRatesPerOffspringType <- vitalRatesPerOffspringType 
	f1@fecByDiscrete <- fecByDiscrete
	f1@Transform <- Transform
	
	} else {
		
	if (is.null(coeff)) stop("require coefficients if data is not supplied")
	if (doOffspring) if (is.na(meanOffspringSize) | is.na(sdOffspringSize)) stop("require meanOffspringSize and sdOffspringSize if data is not supplied")
	
	
	f1 <- .createFecObj(Formula=Formula, 
							coeff=coeff,
							Family = Family,
							Transform = Transform,
							meanOffspringSize = meanOffspringSize, 
							sdOffspringSize = sdOffspringSize, 
							offspringSplitter = offspringSplitter, 
							vitalRatesPerOffspringType = vitalRatesPerOffspringType, 
							fecByDiscrete = fecByDiscrete, 
							offspringSizeExplanatoryVariables = offspringSizeExplanatoryVariables,
							fecConstants = fecConstants, doOffspring=doOffspring, 
							reproductionType=reproductionType) 
	
		
	}
	
	return(f1)
}








# 3a. clonality models  #######################################################################################################

# =============================================================================
# =============================================================================
#now as a wrapper...
makeClonalObj <- function(dataf=NULL,
		fecConstants=data.frame(NA),
		Formula=list(fec~size),
		Family="gaussian",
		Transform="none",
		meanOffspringSize=NA,
		sdOffspringSize=NA,
		offspringSplitter=data.frame(continuous=1),
		vitalRatesPerOffspringType=data.frame(NA),
		fecByDiscrete=data.frame(NA),
		offspringSizeExplanatoryVariables="1", 
		coeff=NULL, doOffspring=TRUE) { 

	f1 <-makeFecObj(dataf=dataf,
				fecConstants=fecConstants,
				Formula=Formula,
				Family=Family,
				Transform=Transform,
				meanOffspringSize=meanOffspringSize,
				sdOffspringSize=sdOffspringSize,
				offspringSplitter=offspringSplitter,
				vitalRatesPerOffspringType=vitalRatesPerOffspringType,
				fecByDiscrete=fecByDiscrete,
				offspringSizeExplanatoryVariables=offspringSizeExplanatoryVariables, 
				coeff=coeff, doOffspring=doOffspring,
				reproductionType="clonal")
	
	return(f1)
}




# 4. Discrete Transition models  #######################################################################################################

# =============================================================================
# =============================================================================
## Function to take a data-frame and make a discrete transition object
## for combining with a continuous P matrix
#
# Parameters - dataf - dataframe with headings of at least
#                      size, sizeNext, surv, fec, stage, stageNext, number
#
# Returns - an object of class discreteTrans
#

makeDiscreteTrans <- function(dataf, 
		stages = NA,
		discreteTrans = NA,
		meanToCont = NA,
		sdToCont = NA,
		continuousToDiscreteExplanatoryVariables = "size") {
	
	#order stage names from discrete to continuous
	if (is.na(stages[1])) {
		stages <- names(tapply(c(levels(dataf$stage),levels(dataf$stageNext)),c(levels(dataf$stage),levels(dataf$stageNext)),length))
		if (!is.na(discreteTrans[1])) stages<-c(stages,dimnames(discreteTrans)[[2]])
	}
	stages <- unique(stages)
	stages <- c(stages[!stages%in%c("continuous","dead")],"continuous","dead") 
	if (length(stages)==2) stop("Error - no discrete stages found. If no discrete stages are included in your data file, please specify them in the discreteTrans argument of the makeDiscreteTrans function.")
	#if no number of instances are not specified, assume each row represents one individual
	if (("number"%in%names(dataf)) == FALSE) dataf$number <- 1
	#define the number of discrete classes
	nDiscreteClasses <- length(stages)-2
	#define the transition between all classes
	if (is.na(discreteTrans[1])&length(discreteTrans)==1) {
		discreteTrans <- matrix(0,nrow=nDiscreteClasses+2,ncol=nDiscreteClasses+1, dimnames=list(stages,stages[1:(length(stages)-1)]))
		for (j in stages[1:(length(stages)-1)]) {
			for (i in stages) discreteTrans[i,j] <- sum(dataf[dataf$stage==j & dataf$stageNext==i,]$number,na.rm=TRUE)
		}
	}
	if (class(discreteTrans)!="matrix") stop("Error - the discreteTrans you entered should be a matrix")
	if (nrow(discreteTrans)!=length(stages)|ncol(discreteTrans)!=(length(stages)-1)) stop("Error - the discreteTrans matrix you entered should be a square matrix with dimensions equal to the number of stages (including continuous)")
	if (sum(dimnames(discreteTrans)[[1]]==stages)<length(stages)) stop("Error - the row names of your discreteTrans matrix should be in alphabetical order, with continuous being the last one")
	if (sum(dimnames(discreteTrans)[[2]]==stages[1:(length(stages)-1)])<(length(stages)-1)) stop("Error - the column names of your discreteTrans matrix should be in alphabetical order, with continuous being the last one")
	for (j in stages[1:(length(stages)-1)]) discreteTrans[,j] <- discreteTrans[,j] / sum(discreteTrans[,j], na.rm = TRUE)
	#define the mean size of individuals coming from discrete stages to the continuous stage
	if (is.na(meanToCont[1])&length(meanToCont)==1) {
		meanToCont <- matrix(NA,nrow=1,ncol=nDiscreteClasses,dimnames=list(1,stages[1:nDiscreteClasses]))
		for (j in stages[which(as.numeric(discreteTrans["continuous",1:nDiscreteClasses])>0)]) {
			meanToCont[,j] <- mean(dataf[dataf$stage == j & dataf$stageNext == "continuous",]$sizeNext, na.rm = TRUE)
		}
	}
	if (class(meanToCont)!="matrix") stop("Error - the meanToCont matrix you entered should be a matrix")
	if (nrow(meanToCont)!=1) stop("Error - the meanToCont matrix you entered should contain just 1 row with means (or NA's for those discrete stages from which no individuals move to the continuous class")
	if (sum(dimnames(meanToCont)[[2]]==stages[1:nDiscreteClasses])<nDiscreteClasses) stop("Error - the column names of the meanToCont matrix you entered should be in alphabetical order and match the column names of the discrete classes in discreteTrans (so without continuous). If some of the discete stages are not mentioned in your data file, this error can be fixed by adding those stages first: levels(dataf$stage)<-c(levels(dataf$stage),<unmentioned_discrete_stages>)")
	#define the sd size of individuals coming from discrete stages to the continuous stage
	if (is.na(sdToCont[1])&length(sdToCont)==1) {
		sdToCont <- matrix(NA,nrow=1,ncol=nDiscreteClasses,dimnames=list(1,stages[1:nDiscreteClasses]))
		for (j in stages[which(as.numeric(discreteTrans["continuous",1:nDiscreteClasses])>0)]) {
			sdToCont[,j] <- sd(dataf[dataf$stage == j & dataf$stageNext == "continuous",]$sizeNext, na.rm = TRUE)
		}
	}
	if (class(sdToCont)!="matrix") stop("Error - the sdToCont matrix you entered should be a matrix")
	if (nrow(sdToCont)!=1) stop("Error - the sdToCont matrix you entered should contain just 1 row with means (or NA's for those discrete stages from which no individuals move to the continuous class")
	if (sum(dimnames(sdToCont)[[2]]==stages[1:nDiscreteClasses])<nDiscreteClasses) stop("Error - the column names of the sdToCont matrix you entered should be in alphabetical order and match the column names of the discrete classes in discreteTrans (so without continuous). If some of the discete stages are not mentioned in your data file, this error can be fixed by adding those stages first: levels(dataf$stage)<-c(levels(dataf$stage),<unmentioned_discrete_stages>)")
	# make the regression to relate the probability of individuals moving to any of the discrete stages as a function of their size 
	if (sum(discreteTrans[stages[1:nDiscreteClasses],"continuous"])==0) {
		moveToDiscrete <- glm(rep(0,21)~1, family = binomial)
	} else {
		subData <- subset(dataf, dataf$stage == "continuous" & dataf$surv == 1)
		subData$contToDiscrete <- 1
		subData$contToDiscrete[subData$stageNext == "continuous"] <- 0
		subData$size2 <- subData$size ^ 2
		subData$size3 <- subData$size ^ 3
		if (length(grep("expsize", as.character(continuousToDiscreteExplanatoryVariables))) > 0) subData$expsize <- exp(subData$size)
		if (length(grep("logsize", as.character(continuousToDiscreteExplanatoryVariables))) > 0) subData$logsize <- log(subData$size)
		moveToDiscrete <- glm(paste('contToDiscrete~',continuousToDiscreteExplanatoryVariables,sep=''), family = binomial, data = subData)
	}
	
	#define new object
	disTrans <- new("discreteTrans")
	disTrans@discreteTrans <- discreteTrans
	disTrans@meanToCont <- meanToCont
	disTrans@sdToCont <- sdToCont
	disTrans@moveToDiscrete <- moveToDiscrete
	return(disTrans)
}


## 5. Models including Data Augmentation models #############################################################


# =============================================================================
# =============================================================================
### integer related functions

makeFecObjInteger <- function(dataf,
		fecConstants=data.frame(NA),
		Formula=list(fec~size),
		Family="gaussian",
		Transform="none",
		meanOffspringSize=NA,
		thetaOffspringSize=NA,
		offspringSplitter=data.frame(continuous=1),
		vitalRatesPerOffspringType=data.frame(NA),
		fecByDiscrete=data.frame(NA),
		offspringSizeExplanatoryVariables="1",
		distOffspring = "poisson", 
		coeff=NULL, doOffspring=TRUE, 
		reproductionType="sexual"){
	
	
	#make sure Formula is a formula or a list of formulas
	if (class(Formula)=="list") {
		if (class(Formula[[1]])!="formula") stop("Error - the entries in your Formula list should be of class 'formula': e.g. fec~size without quotation marks")
	} else {
		if (class(Formula)!="formula") stop("Error - the Formula entry should by of class 'formula' or a list of such entries:  e.g. fec~size without quotation marks")
		Formula <- list(Formula)	
	}
	
	#if stage or stageNext do not exist in dataf, create them assuming that 
	#everything is continuous. 
	if (length(dataf$stage)==0) { 
		dataf$stage <- rep("continuous",nrow(dataf))
		dataf$stage[is.na(dataf$size)] <- NA
		dataf$stage <- as.factor(dataf$stage)
	}
	if (length(dataf$stageNext)==0) {
		dataf$stageNext <- rep("continuous",nrow(dataf))
		dataf$stageNext[dataf$surv==0] <- "dead"
		dataf$stageNext <- as.factor(dataf$stageNext)
	}
	
	#order stage names from discrete to continuous
	stages <- names(tapply(c(levels(dataf$stage),levels(dataf$stageNext)),c(levels(dataf$stage),levels(dataf$stageNext)),length))
	stages <- stages[stages!="dead"] 
	if ((sum(names(offspringSplitter)%in%stages)/length(offspringSplitter))<1) {
		stages <- c(stages,names(offspringSplitter))
		print("Warning - the variable names in your offspringSplitter data.frame are not all part of the levels of stage or stageNext in your data file. This could be because of an mismatch in stage names, or because you included discrete stages in offspringSplitter that are not in the data file but wchich you will introduce in makeDiscreteTrans (in which case you can ignore this warning).")
	}
	stages <- unique(stages)
	stages <- c(stages[stages!="continuous"],"continuous") 
	dummy<-rep(0,length(stages));names(dummy)<-stages;dummy<-as.data.frame(t(as.matrix(dummy)))
	for (i in names(offspringSplitter)) dummy[i]<-offspringSplitter[i]
	offspringSplitter <- dummy
	
	##warnings
	if (ncol(offspringSplitter)>1 & (ncol(offspringSplitter)-1)!=ncol(fecByDiscrete)) {
		print("Warning - offspring splitter indicates more than just continuous stages. No fecundity by the discrete stages supplied in fecByDiscrete; assumed that is 0")
		fecByDiscrete <- offspringSplitter[,1:(ncol(offspringSplitter)-1)]
		fecByDiscrete[] <- 0
	}
	
	if (sum(offspringSplitter)!=1) {
		print("Warning - offspring splitter does not sum to 1. It is now rescaled to sum to 1.")
		offspringSplitter <- offspringSplitter / sum(offspringSplitter) 
	}
	
	if ("covariate" %in% unlist(strsplit(as.character(Formula), "[+-\\* ]")) & length(dataf$covariate) > 0) { 
		dataf$covariate <- as.factor(dataf$covariate)
		levels(dataf$covariate) <- 1:length(unique(dataf$covariate))
	}
	if ("covariateNext" %in% unlist(strsplit(as.character(Formula), "[+-\\* ]")) & length(dataf$covariateNext) > 0) { 
		dataf$covariateNext <- as.factor(dataf$covariateNext)
		levels(dataf$covariateNext) <- 1:length(unique(dataf$covariateNext))
	}
	
	f1 <- new("fecObjInteger")
	
	dataf$size2 <- dataf$size^2
	dataf$size3 <- dataf$size^3
	if (length(grep("expsize",unlist(as.character(Formula))))>0) dataf$expsize <- exp(dataf$size)
	if (length(grep("logsize",unlist(as.character(Formula))))>0) dataf$logsize <- log(dataf$size)
	
	if (length(Formula)>length(Family)) {
		misE <- (length(Family)+1):length(Formula)
		print(c("number of families not the same as the number of Formula supplied, using default of `gaussian' for missing ones which are:",Formula[[misE]],". (which might be exactly what you want)"))
		Family <- c(Family,rep("gaussian",length(Formula)-length(Family)))
	}
	if (length(Formula)>length(Transform)) {
		misE <- (length(Transform)+1):length(Formula)
		print(c("number of transforms not the same as the number of fecundity columns in the data file, using default of `none' for missing ones which are:",Formula[[misE]],". (which might be exactly what you want)"))
		Transform <- c(Transform,rep("none",length(Formula)-length(Transform)))
	}
	
	fecNames <- rep(NA,length(Formula))
	for (i in 1:length(Formula)) {
		
		fecNames[i] <- all.vars(Formula[[i]])[1]
		
		if (Transform[i]=="exp") dataf[,fecNames[i]] <- exp(dataf[,fecNames[i]])
		if (Transform[i]=="log") dataf[,fecNames[i]] <- log(dataf[,fecNames[i]])
		if (Transform[i]=="sqrt") dataf[,fecNames[i]] <- sqrt(dataf[,fecNames[i]])
		if (Transform[i]=="-1") dataf[,fecNames[i]] <- dataf[,fecNames[i]]-1
		dataf[!is.finite(dataf[,fecNames[i]]),fecNames[i]] <- NA
		if (length(intersect(all.vars(Formula[[i]]),colnames(dataf)))<length(all.vars(Formula[[i]]))) print("warning: not all variables in the formula are present in dataf; model cannot be fit")
		f1@fitFec[[i]] <- glm(Formula[[i]],family=Family[i],data=dataf)
	}
	
	if (doOffspring) {
	if (offspringSplitter$continuous>0) {
		if (is.na(meanOffspringSize[1])) {
			if (length(dataf$offspringNext)==0) {
				offspringData<-subset(dataf,is.na(dataf$stage)&dataf$stageNext=="continuous")
				if (nrow(offspringData) == 0) stop ("Error - no offspring size data are given: these can be given through either the meanOffspringSize and thetaOffspringSize slots, or through individual data added to your data file (with stage equals NA, or a offspringNext column indicating 'sexual/clonal' offspring)")
			} else {
				offspringData<-subset(dataf,dataf$offspringNext==reproductionType&dataf$stageNext=="continuous")
				if (nrow(offspringData) == 0) stop ("Error - no offspring size data are given: these can be given through either the meanOffspringSize  slot, or through individual data added to your data file (with stage equals NA, or a offspringNext column indicating 'sexual/clonal' offspring)")
			}
			## relationship defining offspring size  
			if (distOffspring=="poisson") f1@offspringRel <- glm(paste('sizeNext~',offspringSizeExplanatoryVariables,sep=''),data=offspringData, family="poisson")
			if (distOffspring=="negBin") {
				f1@offspringRel <- glm.nb(paste('sizeNext~',offspringSizeExplanatoryVariables,sep=''),data=offspringData)
				f1@thetaOffspringSize <- f1@offspringRel$theta;
				f1@offspringRel <- glm.convert(f1@offspringRel)}
		} else {
			if (distOffspring=="poisson") 
				f1@offspringRel <- glm(rep(meanOffspringSize[1],21)~1, family="poisson")
			if (distOffspring=="negBin") { 
				f1@offspringRel <- glm.nb(rep(meanOffspringSize[1],21)~1)
				f1@thetaOffspringSize <- f1@offspringRel$theta
				f1@offspringRel <- glm.convert(f1@offspringRel)}
		}
	}
	} ## end doOffspring

	if (sum(dim(vitalRatesPerOffspringType)==c(1,1))<2) {
		if ((sum(vitalRatesPerOffspringType==0,na.rm=T)+sum(vitalRatesPerOffspringType==1,na.rm=T))<(ncol(vitalRatesPerOffspringType)*nrow(vitalRatesPerOffspringType))) stop("Error - in vitalRatesPerOffspringType data.frame only 0's and 1's are allowed: a 1 indicates that a fecundity rate applies to that offspring type. ")
		#if (sum(names(vitalRatesPerOffspringType)==names(offspringSplitter))<length(offspringSplitter)) stop("Error - the offspring names in vitalRatesPerOffspringType should match those in offspringSplitter - and in the same order, with continuous last")
		if (sum(rownames(vitalRatesPerOffspringType)==c(fecNames,names(fecConstants)))<(length(Formula)+length(fecConstants))) stop ("Error - the row names in vitalRatesPerOffspringType should consist of (in order) the names of the fec columns in the dataset and then the names of the fecConstants.")
	} else {
		vitalRatesPerOffspringType <- as.data.frame(matrix(1,ncol=length(offspringSplitter),nrow=length(Formula)+length(fecConstants)),row.names=c(fecNames,names(fecConstants)))
		vitalRatesPerOffspringType <- subset(vitalRatesPerOffspringType,dimnames(vitalRatesPerOffspringType)[[1]]!="NA.")
		names(vitalRatesPerOffspringType) <- names(offspringSplitter)
	}
	
	if (doOffspring & length(f1@thetaOffspringSize)>0 & distOffspring=="negBin") {
		if (is.na(f1@thetaOffspringSize)) {
			print("Warning - could not estimate parameters for the distribution of offspring size; defaults must be supplied for meanOffspringSize and thetaOffspringSize; you will not be able to construct an IPM without these values.")
		}
	}
	
	f1@fecNames <- fecNames
	f1@fecConstants <- fecConstants
	f1@offspringSplitter <- offspringSplitter 
	f1@vitalRatesPerOffspringType <- vitalRatesPerOffspringType 
	f1@fecByDiscrete <- fecByDiscrete
	f1@Transform <- Transform
	f1@distOffspring <- distOffspring
	
	return(f1)
}

# =============================================================================
# =============================================================================
makeClonalObjInteger <- function(dataf,
		fecConstants=data.frame(NA),
		Formula=list(fec~size),
		Family="gaussian",
		Transform="none",
		meanOffspringSize=NA,
		thetaOffspringSize=NA,
		offspringSplitter=data.frame(continuous=1),
		vitalRatesPerOffspringType=data.frame(NA),
		fecByDiscrete=data.frame(NA),
		offspringSizeExplanatoryVariables="1",
		distOffspring = "poisson",
		coeff=NULL, doOffspring=TRUE){
	
	
	f1 <- makeFecObjInteger(dataf=dataf,
			fecConstants=fecConstants,
			Formula=Formula,
			Family=Family,
			Transform=Transform,
			meanOffspringSize=meanOffspringSize,
			thetaOffspringSize=thetaOffspringSize,
			offspringSplitter=offspringSplitter,
			vitalRatesPerOffspringType=vitalRatesPerOffspringType,
			fecByDiscrete=fecByDiscrete,
			offspringSizeExplanatoryVariables=offspringSizeExplanatoryVariables,
			distOffspring = distOffspring , 
			coeff=coeff, doOffspring=doOffspring, 
			reproductionType="clonal")
	
	return(f1)
}


# =============================================================================
# =============================================================================
makeDiscreteTransInteger <- function(dataf, 
		stages = NA,
		discreteTrans = NA,
		meanToCont = NA,
		thetaToCont = NA,
		continuousToDiscreteExplanatoryVariables = "size",
		distToCont="poisson") {
	
	#order stage names from discrete to continuous
	if (is.na(stages[1])) {
		stages <- names(tapply(c(levels(dataf$stage),levels(dataf$stageNext)),c(levels(dataf$stage),levels(dataf$stageNext)),length))
		if (!is.na(discreteTrans[1])) stages<-c(stages,dimnames(discreteTrans)[[2]])
	}
	stages <- unique(stages)
	stages <- c(stages[!stages%in%c("continuous","dead")],"continuous","dead") 
	if (length(stages)==2) stop("Error - no discrete stages found. If no discrete stages are included in your data file, please specify them in the discreteTrans argument of the makeDiscreteTrans function.")
	#if no number of instances are not specified, assume each row represents one individual
	if (("number"%in%names(dataf)) == FALSE) dataf$number <- 1
	#define the number of discrete classes
	nDiscreteClasses <- length(stages)-2
	#define the transition between all classes
	if (is.na(discreteTrans[1])&length(discreteTrans)==1) {
		discreteTrans <- matrix(0,nrow=nDiscreteClasses+2,ncol=nDiscreteClasses+1, dimnames=list(stages,stages[1:(length(stages)-1)]))
		for (j in stages[1:(length(stages)-1)]) {
			for (i in stages) discreteTrans[i,j] <- sum(dataf[dataf$stage==j & dataf$stageNext==i,]$number,na.rm=TRUE)
		}
	}
	if (class(discreteTrans)!="matrix") stop("Error - the discreteTrans you entered should be a matrix")
	if (nrow(discreteTrans)!=length(stages)|ncol(discreteTrans)!=(length(stages)-1)) stop("Error - the discreteTrans matrix you entered should be a square matrix with dimensions equal to the number of stages (including continuous)")
	if (sum(dimnames(discreteTrans)[[1]]==stages)<length(stages)) stop("Error - the row names of your discreteTrans matrix should be in alphabetical order, with continuous being the last one")
	if (sum(dimnames(discreteTrans)[[2]]==stages[1:(length(stages)-1)])<(length(stages)-1)) stop("Error - the column names of your discreteTrans matrix should be in alphabetical order, with continuous being the last one")
	for (j in stages[1:(length(stages)-1)]) discreteTrans[,j] <- discreteTrans[,j] / sum(discreteTrans[,j], na.rm = TRUE)
	#define the mean size of individuals coming from discrete stages to the continuous stage
	if (is.na(meanToCont[1])&length(meanToCont)==1) {
		meanToCont <- matrix(NA,nrow=1,ncol=nDiscreteClasses,dimnames=list(1,stages[1:nDiscreteClasses]))
		#define theta also so as can use it in the loop....
		if (is.na(thetaToCont[1])&length(thetaToCont)==1) 
			thetaToCont <- matrix(NA,nrow=1,ncol=nDiscreteClasses,dimnames=list(1,stages[1:nDiscreteClasses]))
			
		
		for (j in stages[which(as.numeric(discreteTrans["continuous",1:nDiscreteClasses])>0)]) {
		if (distToCont=="poisson")
			tmp <- glm(dataf[dataf$stage == j & dataf$stageNext == "continuous",]$sizeNext~1,
					family="poisson")
		if (distToCont=="negbin") {
			tmp <- glm.nb(dataf[dataf$stage == j & dataf$stageNext == "continuous",]$sizeNext~1)
			thetaToCont[,j] <- tmp$theta
		}
		
		meanToCont[,j] <- exp(tmp$coefficient[1])
			
			
	}
	}
	if (class(meanToCont)!="matrix") stop("Error - the meanToCont matrix you entered should be a matrix")
	if (nrow(meanToCont)!=1) stop("Error - the meanToCont matrix you entered should contain just 1 row with means (or NA's for those discrete stages from which no individuals move to the continuous class")
	if (sum(dimnames(meanToCont)[[2]]==stages[1:nDiscreteClasses])<nDiscreteClasses) stop("Error - the column names of the meanToCont matrix you entered should be in alphabetical order and match the column names of the discrete classes in discreteTrans (so without continuous)")
	#define the sd size of individuals coming from discrete stages to the continuous stage
	
if (distToCont=="negbin") { 
		if (class(thetaToCont)!="matrix") stop("Error - the thetaToCont matrix you entered should be a matrix")
		if (nrow(thetaToCont)!=1) stop("Error - the thetaToCont matrix you entered should contain just 1 row with means (or NA's for those discrete stages from which no individuals move to the continuous class")
		if (sum(dimnames(thetaToCont)[[2]]==stages[1:nDiscreteClasses])<nDiscreteClasses) stop("Error - the column names of the thetaToCont matrix you entered should be in alphabetical order and match the column names of the discrete classes in discreteTrans (so without continuous)")
		# make the regression to relate the probability of individuals moving to any of the discrete stages as a function of their size 
	}
	if (sum(discreteTrans[stages[1:nDiscreteClasses],"continuous"])==0) {
		moveToDiscrete <- glm(rep(0,21)~1, family = binomial)
	} else {
		subData <- subset(dataf, dataf$stage == "continuous" & dataf$surv == 1)
		subData$contToDiscrete <- 1
		subData$contToDiscrete[subData$stageNext == "continuous"] <- 0
		subData$size2 <- subData$size ^ 2
		subData$size3 <- subData$size ^ 3
		if (length(grep("expsize", as.character(continuousToDiscreteExplanatoryVariables))) > 0) subData$expsize <- exp(subData$size)
		if (length(grep("logsize", as.character(continuousToDiscreteExplanatoryVariables))) > 0) subData$logsize <- log(subData$size)
		moveToDiscrete <- glm(paste('contToDiscrete~',continuousToDiscreteExplanatoryVariables,sep=''), family = binomial, data = subData)
	}
	
	#define new object
	disTrans <- new("discreteTransInteger")
	disTrans@discreteTrans <- discreteTrans
	disTrans@meanToCont <- meanToCont
	disTrans@thetaToCont <- thetaToCont
	disTrans@moveToDiscrete <- moveToDiscrete
	disTrans@distToCont <- distToCont
	return(disTrans)
}

# =============================================================================
# =============================================================================
## FUNCTION FOR TURNING DATA INTO MATRIX DEFINING ENVIRONMENTAL TRANSITIONS ###
## data is vector of env level at t, and one timestep later, at t+1

makeEnvObj <- function(dataf){
  #turn into index starting at 1
  minval <-  min(c(dataf$covariate,dataf$covariateNext),na.rm=TRUE)
  startEnv <- dataf$covariate-minval+1
  nextEnv <- dataf$covariateNext-minval+1
  
  
  nEnvClass <- max(c(startEnv,nextEnv), na.rm=TRUE)
  desired.mat <- matrix(0,nEnvClass,nEnvClass) 
  mats<-table(startEnv,nextEnv)
  rx <- as.numeric(rownames(mats));#print(rx)
  cx <- as.numeric(colnames(mats))
  desired.mat[cbind(rep(rx,length(cx)),rep(cx,each=length(rx)))]=c(as.matrix(mats))
  
  rc <- new("envMatrix",
      nEnvClass = nEnvClass)
  
  rc@.Data <- t(t(desired.mat)/colSums(desired.mat))
  
  return(rc) 
  
}



