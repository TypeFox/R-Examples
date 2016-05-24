#' An example of the Mulstistrata (multi-state) model in which states are routes taken by migrating fish.
#' 
#' @name skagit
#' @author Megan Moore <megan.moore at noaa.gov> 
#' @examples
#'# There are just two states which correspond to route A and route B. There are 6 occasions
#'# which are the locations rather than times. After release at 1=A there is no movement 
#'# between states for the first segment, fish are migrating downriver together and all pass 2A. 
#'# Then after occasion 2, migrants go down the North Fork (3A) or the South Fork (3B), 
#'# which both empty into Skagit Bay. Once in saltwater, they can go north to Deception Pass (4A)
#'# or South to a receiver array exiting South Skagit Bay (4B). Fish in route A can then only go
#'# to the Strait of Juan de Fuca, while fish in route B must pass by Admiralty Inlet (5B). 
#'# Then both routes end with the array at the Strait of Juan de Fuca.
#'#
#'#       1A
#'#        |
#'#        2A
#'#      /     \
#'#    3A        3B
#'#   /  \      /  \ 
#'# 4A   4B    4A  4B       
#'#  |     \    /   |
#'#   5A    5B  5A   5B
#'#      \   \   /    /
#'#            6
#'# 
#'# from 3A and 3B they can branch to either 4A or 4B; branches merge at 6  
#'# 5A does not exist so p=0; only survival from 4A to 6 can be 
#'# estimated which is done by setting survival from 4A to 5A to 1 and
#'# estimating survival from 5A to 6 which is then total survival from 4A to 6.
#'\donttest{
#'pathtodata=paste(path.package("RMark"),"extdata",sep="/")
#'skagit=import.chdata(paste(pathtodata,"skagit.txt",sep="/"),field.types=c("f"),header=TRUE)
#'skagit.processed=process.data(skagit,model="Multistrata",groups=c("tag"))
#'skagit.ddl=make.design.data(skagit.processed)
#'#
#'# p
#'#
#'# Can't be seen at 5A or 2B,6B (the latter 2 don't exist)
#'skagit.ddl$p$fix=ifelse((skagit.ddl$p$stratum=="A"&skagit.ddl$p$time==5) | 
#'  (skagit.ddl$p$stratum=="B"&skagit.ddl$p$time%in%c(2,6)),0,NA)
#'# Estimated externally from current data to allow estimation of survival at last interval
#'skagit.ddl$p$fix[skagit.ddl$p$tag=="v7"&skagit.ddl$p$time==6&skagit.ddl$p$stratum=="A"]=0.687
#'skagit.ddl$p$fix[skagit.ddl$p$tag=="v9"&skagit.ddl$p$time==6&skagit.ddl$p$stratum=="A"]=0.975
#'#
#'# Psi
#'#
#'# only 3 possible transitions are A to B at time interval 2 to 3 and 
#'# for time interval 3 to 4 from A to B and from B to A
#'# rest are fixed values
#'skagit.ddl$Psi$fix=NA
#'# stay in A for intervals 1-2, 4-5 and 5-6
#'skagit.ddl$Psi$fix[skagit.ddl$Psi$stratum=="A"&
#'   skagit.ddl$Psi$tostratum=="B"&skagit.ddl$Psi$time%in%c(1,4,5)]=0
#'# stay in B for interval 4-5
#'skagit.ddl$Psi$fix[skagit.ddl$Psi$stratum=="B"&skagit.ddl$Psi$tostratum=="A"
#'   &skagit.ddl$Psi$time==4]=0
#'# leave B to go to A for interval 5-6
#'skagit.ddl$Psi$fix[skagit.ddl$Psi$stratum=="B"&skagit.ddl$Psi$tostratum=="A"&
#'  skagit.ddl$Psi$time==5]=1
#'# "stay" in B for interval 1-2 and 2-3 because none will be in B
#'skagit.ddl$Psi$fix[skagit.ddl$Psi$stratum=="B"&skagit.ddl$Psi$tostratum=="A"&
#'  skagit.ddl$Psi$time%in%1:2]=0
#'# 
#'# S
#'#
#'# None in B, so fixing S to 1
#'skagit.ddl$S$fix=ifelse(skagit.ddl$S$stratum=="B"&skagit.ddl$S$time%in%c(1,2),1,NA)
#'skagit.ddl$S$fix[skagit.ddl$S$stratum=="A"&skagit.ddl$S$time==4]=1
#'# fit model
#'p.timexstratum.tag=list(formula=~time:stratum+tag,remove.intercept=TRUE)
#'Psi.sxtime=list(formula=~-1+stratum:time)
#'S.stratumxtime=list(formula=~-1+stratum:time)
#'#
#'S.timexstratum.p.timexstratum.Psi.sxtime=mark(skagit.processed,skagit.ddl,
#'  model.parameters=list(S=S.stratumxtime,p= p.timexstratum.tag,Psi=Psi.sxtime))
#'# calculation of cummulative survival for entire route
#'Sest=plogis(coef(S.timexstratum.p.timexstratum.Psi.sxtime)$estimate)
#'# A
#'prod(Sest[c(1:3,6)])
#'#[1] 0.1644
#'# B
#'prod(Sest[c(1,2,4,5,7)])
#'#[1] 0.1154
#' }

NULL

#' #' Burnham Live-Dead Model
#' 
#' An example of the Burnham live-dead model using simulated data LD1.inp from Chapter 9 of Cooch and White
#' 
#' @name Burnham
#' @author Luke Eberhart-Phillips<luke.eberhart at gmail.com> 
#' @examples
#' \donttest{
#'###############################################################################
#'#### RMARK script for conducting the Burnham model tutorial in Chapter 9.3 ####
#'#################### the of the Cooch and White MARK book #####################
#'###############################################################################
#'##################### Code by: Luke Eberhart-Phillips #########################
#'###### Dept. Animal Behaviour, Bielefeld University, Bielefeld, Germany #######
#'##################### email: luke.eberhart at gmail.com ##########################
#'###############################################################################
#'
#'# import/convert the simulated "LD1.inp" MARK capture history into an RMARK 
#'# dataframe, while defining the two groups as "Y" for individuals marked as 
#'# young, and "A" for individuals marked as adults
#'# NOTE: the "LD1.inp" file is found in the zipped folder downloaded when you
#'# click on "Example data files" in the drop-down menu of the MARK book webpage
#'# (http://www.phidot.org/software/mark/docs/book/)
#' pathtodata=paste(path.package("RMark"),"extdata",sep="/")
#' LD=convert.inp(paste(pathtodata,"ld1",sep="/"),
#'            group.df=data.frame(age_marked=c("Y","A")))
#'# process the data by defining the model type as "Burnham" and the groups in
#'# the data.  In this case the only group is the age at which individuals were 
#'# marked
#'LD.proc=process.data(data = LD, 
#'		model = "Burnham",
#'		groups=c("age_marked"),
#'		age.var=1,
#'		initial.age=c(1,0))
#'
#'# make the design data from the process data above
#'LD.ddl=make.design.data(LD.proc)
#'
#'# add the correct binning to the design data so that individuals that were
#'# marked as young are adults in their second year of life, where as those
#'# marked as adults are adults for their entire life.
#'LD.ddl=add.design.data(data = LD.proc,
#'		ddl = LD.ddl,
#'		parameter="S", 
#'		type = "age",
#'		bins = c(0,1,8),
#'		right = FALSE,
#'		name = "age",
#'		replace = TRUE)
#'
#'# do the same to the F parameter
#'LD.ddl=add.design.data(data = LD.proc,
#'		ddl = LD.ddl,
#'		parameter="F", 
#'		type = "age",
#'		bins = c(0,1,8),
#'		right = FALSE,
#'		name = "age",
#'		replace = TRUE)
#'
#'# check parameter matrix to see if groups were binned correctly in the S matrix
#'PIMS(mark(data = LD.proc,
#'				ddl = LD.ddl,
#'				model.parameters=list(S=list(formula=~age)),
#'				output=FALSE,
#'				model = "Burnham"),
#'		"S")
#'
#'# Create the formulas that describe variation in the parameter we want to test.
#'# In this case we want to test for an age effect on survival and fidelity,
#'# while keeping recapture and recovery probabilities constant.
#'S.age=list(formula=~age) # S(age)
#'p.dot=list(formula=~1) # p(.)
#'F.age=list(formula=~age) # F(age)
#'r.dot=list(formula=~1) # r(.)
#'
#'# Run the model
#'LD.model.age.F.S=mark(data = LD.proc,
#'		ddl = LD.ddl, 
#'		model.parameters = list(S = S.age, p = p.dot, 
#'				F =F.age, r = r.dot), 
#'		invisible = FALSE, 
#'		model = "Burnham")
#'
#'# Check the paramter estimates, they should be the same as those generated
#'# when doing the tutorial in chapter 9.3 of the in MARK Book (table on pg 9-8)
#'LD.model.age.F.S$results$real
#'
#'# Clean your working directory
#'cleanup(ask=FALSE)
#'}

NULL

#' Lark Sparrow
#' 
#' An example of Multiple Scale Occupancy model for some lark sparrow data that was contributed by David Pavlacky 
#' at Rocky Mountain bird observatory. The study design was a GRTS selection of paired "Deferred" and "Grazed" pastures.  
#' The point count locations within each pasture were a random selection of systematic point count locations separated by 250 m. Each point count had
#' a radius of 125m. A removal design was used to the set the data to missing after the first detection.  
#' The point count data were set to missing when fewer than nine points were surveyed.     
#' 
#' @name larksparrow
#' @aliases LASP
#' @docType data
#' @format A data frame with 52 observations on the following 20 variables.
#' \describe{\item{ch}{a character vector containing the encounter history}
#' \item{ceap}{a factor with two levels "Deferred" and "Grazed" corresponding to a rest rotation grazing system with pastures either rested (Deferred) or grazed (Grazed) during the spring bird breeding season.}
#' \item{cwx}{a continuous covariate for primary occasion x, representing an ocular estimate of the proportion of area  covered by crested wheatgrass in a 50-m radius around the point count location.}
#' \item{tdx}{ a continuous covariate for primary occasion x, representing the starting time (h) of each 6-min point count survey measured on the ratio scale (1.5 h = 1 h 30 min).}
#' }
#' @keywords datasets
#' @examples
#' \donttest{
#' 
#' # This example is excluded from testing to reduce package check time
#'# Create dataframe
#'data(LASP)
#' mscale=LASP
#'
#'# Process data with MultScalOcc model and use group variables
#'
#'mscale.proc=process.data(mscale,model="MultScalOcc",groups=c("ceap"),begin.time=1,mixtures=3)
#'
#'# Create design data
#'
#'ddl=make.design.data(mscale.proc)
#'
#'# Create function to build models
#'
#'do.Species=function()
#'{
#'	p.1=list(formula=~1)   
#'	p.2=list(formula=~ceap)    
#'	p.3=list(formula=~td)
#'	
#'	Theta.1=list(formula=~1)    
#'	Theta.2=list(formula=~ceap)   
#'	Theta.3=list(formula=~cw)
#'	
#'	Psi.1=list(formula=~1)    
#'	Psi.2=list(formula=~ceap)    
#'	
#'	cml=create.model.list("MultScalOcc")
#'	return(mark.wrapper(cml,data=mscale.proc,ddl=ddl,adjust=FALSE,realvcv=TRUE))
#'}
#'
#'# Run function to get results
#'
#'Species.results=do.Species()
#'
#'# Output model table and estimates
#'
#'Species.results$model.table
#'
#'Species.results[[as.numeric(rownames(Species.results$model.table[1,]))]]$results$real
#'Species.results[[as.numeric(rownames(Species.results$model.table[1,]))]]$results$beta
#'
#'write.csv(Species.results$model.table,file="lasp_model_selection.csv",row.names=FALSE)
#'
#'write.csv(Species.results[[as.numeric(rownames(Species.results$model.table[1,]))]]$results$real,
#'  file="lasp_m1_real.csv")
#'write.csv(Species.results[[as.numeric(rownames(Species.results$model.table[1,]))]]$results$beta,
#'  file="lasp_m1_beta.csv")
#'
#'# Covariate prediction and model averaging
#'
#'# p(time of day)
#'
#'mintd <- min(mscale[,12:20])
#'maxtd <- max(mscale[,12:20])
#'td.values <- mintd+(0:100)*(maxtd-mintd)/100
#'
#'PIMS(Species.results[[1]],"p",simplified=FALSE)
#'
#'td <- covariate.predictions(Species.results,data=data.frame(td1=td.values),indices=c(21))
#'
#'write.table(td$estimates,file="lasp_cov_pred_p_td.csv",sep=",",col.names=TRUE,row.names=FALSE)
#'
#'# Theta(crested wheatgrass cover)
#'
#'mincw <- min(mscale[,3:11])
#'maxcw <- max(mscale[,3:11])
#'cw.values <- mincw+(0:100)*(maxcw-mincw)/100
#'
#'PIMS(Species.results[[1]],"Theta",simplified=FALSE)
#'
#'cw <- covariate.predictions(Species.results,data=data.frame(cw1=cw.values),indices=c(3))
#'
#'write.table(cw$estimates,file="lasp_cov_pred_theta_cw.csv",sep=",",col.names=TRUE,row.names=FALSE)
#'
#'# Psi(ceap grazing for wildlife practice)
#'
#'ceap.values <- as.data.frame(matrix(c(1,2),ncol=1))
#'names(ceap.values) <- c("index")
#'
#'PIMS(Species.results[[1]],"Psi",simplified=FALSE)
#'
#'ceap <- covariate.predictions(Species.results,data=data.frame(ceap=ceap.values))
#'
#'write.table(ceap$estimates,file="lasp_cov_pred_psi_ceap.csv",sep=",",col.names=TRUE,row.names=FALSE)
#'
#' }

NULL

#' San Luis Valley mallard data
#' 
#' A recovery data set for mallards in San Luis Valley Colorado
#' 
#' This is a data set that accompanies program MARK as an example for Brownie and
#' Seber recovery model.  In those input files it is in a summarized format. Here
#' it is in the LD encounter history format. The data can be stratified using ReleaseAge (Adult, Young) 
#' as a grouping variable.  
#'  
#' Note that in the MARK example the variable is named Age.  In the R code, the
#' fields "age" and "Age" have specific meanings in the design data related to
#' time since release.  These will override the use of a field with the same
#' name in the individual covariate data, so the names "time", "Time",
#' "cohort", "Cohort", "age", and "Age" should not be used in the individual
#' covariate data with possibly the exception of "cohort" which is not defined
#' for models with "Square" PIMS such as POPAN and other Jolly-Seber type
#' models.
#' 
#' @name brownie
#' @docType data
#' @format A data frame with 108 observations on the following 5 variables.
#' \describe{ \item{ch}{a character vector containing the encounter history of
#' each bird} \item{freq}{frequency of that encounter history} \item{ReleaseAge}{the age of the bird when it was
#' released}}
#' @keywords datasets
#' @examples
#' # brownie=import.chdata("brownie.inp",field.types=c("n","f"))
#' \donttest{
#' # This example is excluded from testing to reduce package check time
#'  data(brownie)
#' # default ordering of ReleaseAge is alphabetic so it is 
#' # Adult, Young which is why initial.ages=c(1,0)
#' # Seber Recovery
#' br=process.data(brownie,model="Recovery",groups="ReleaseAge",age.var=1,initial.ages=c(1,0))
#' br.ddl=make.design.data(br,parameters=list(S=list(age.bins=c(0,1,10)),
#'                                            r=list(age.bins=c(0,1,10))),right=FALSE)
#' mod=mark(br,br.ddl,model.parameters=list(S=list(formula=~-1+age:time,link="sin"),
#'                                            r=list(formula=~-1+age:time,link="sin")))
#' # Brownie Recovery
#' br=process.data(brownie,model="Brownie",groups="ReleaseAge",age.var=1,initial.ages=c(1,0))
#' br.ddl=make.design.data(br,parameters=list(S=list(age.bins=c(0,1,10)),
#'                                f=list(age.bins=c(0,1,10))),right=FALSE)
#' mod=mark(br,br.ddl,model.parameters=list(S=list(formula=~-1+age:time,link="sin"),
#'                                f=list(formula=~-1+age:time,link="sin")))
#' mod=mark(br,br.ddl,model.parameters=list(S=list(formula=~-1+age,link="sin"),
#'                                f=list(formula=~-1+age,link="sin")))
#' #Random effects Seber recovery
#' br=process.data(brownie,model="REDead",groups="ReleaseAge",age.var=1,initial.ages=c(1,0))
#' br.ddl=make.design.data(br,parameters=list(S=list(age.bins=c(0,1,10)),
#'                                        r=list(age.bins=c(0,1,10))),right=FALSE)
#' mod=mark(br,br.ddl,model.parameters=list(S=list(formula=~age),r=list(formula=~age)))
#' #Pledger Mixture Seber recovery
#' br=process.data(brownie,model="PMDead",groups="ReleaseAge",
#'                            mixtures=3,age.var=1,initial.ages=c(1,0))
#' br.ddl=make.design.data(br,parameters=list(S=list(age.bins=c(0,1,10)),
#'                             r=list(age.bins=c(0,1,10))),right=FALSE)
#' mod=mark(br,br.ddl,model.parameters=list(pi=list(formula=~mixture),
#'                      S=list(formula=~age+mixture),r=list(formula=~age)))
#' br=process.data(brownie,model="PMDead",groups="ReleaseAge",
#'                      mixtures=2,age.var=1,initial.ages=c(1,0))
#' br.ddl=make.design.data(br,parameters=list(S=list(age.bins=c(0,1,10)),
#'                       r=list(age.bins=c(0,1,10))),right=FALSE)
#' mod=mark(br,br.ddl,model.parameters=list(pi=list(formula=~age),
#'                       S=list(formula=~age+mixture),r=list(formula=~age)))
#' }

NULL

#' Black duck known fate data
#' 
#' A known fate data set on Black ducks that accompanies MARK as an example
#' analysis using the Known model.
#' 
#' This is a data set that accompanies program MARK as an example for Known
#' fate. The data can be stratified using BirdAge as a grouping variable.  The
#' function \code{run.Blackduck} defined below in the examples creates some of
#' the models used in the dbf file that accompanies MARK.
#' 
#' Note that in the MARK example the variable is named Age.  In the R code, the
#' fields "age" and "Age" have specific meanings in the design data related to
#' time since release.  These will override the use of a field with the same
#' name in the individual covariate data, so the names "time", "Time",
#' "cohort", "Cohort", "age", and "Age" should not be used in the individual
#' covariate data with possibly the exception of "cohort" which is not defined
#' for models with "Square" PIMS such as POPAN and other Jolly-Seber type
#' models.
#' 
#' @name Blackduck
#' @docType data
#' @format A data frame with 48 observations on the following 5 variables.
#' \describe{ \item{ch}{a character vector containing the encounter history of
#' each bird} \item{BirdAge}{the age of the bird: a factor with levels \code{0}
#' \code{1} for young and adult} \item{Weight}{the weight of the bird at time
#' of marking} \item{Wing_Len}{the wing-length of the bird at time of marking}
#' \item{condix}{the condition index of the bird at time of marking} }
#' @keywords datasets
#' @examples
#' 
#' data(Blackduck)
#' # Change BirdAge to numeric; starting with version 1.6.3 factor variables are
#' # no longer allowed.  They can work as in this example but they can be misleading
#' # and fail if the levels are non-numeric.  The real parameters will remain 
#' # unchanged but the betas will be different.
#' Blackduck$BirdAge=as.numeric(Blackduck$BirdAge)-1
#' run.Blackduck=function()
#' {
#' #
#' # Process data
#' #
#' bduck.processed=process.data(Blackduck,model="Known")
#' #
#' # Create default design data
#' #
#' bduck.ddl=make.design.data(bduck.processed)
#' #
#' #  Add occasion specific data min < 0; I have no idea what it is
#' #
#' bduck.ddl$S$min=c(4,6,7,7,7,6,5,5)
#' #
#' #  Define range of models for S
#' #
#' S.dot=list(formula=~1)
#' S.time=list(formula=~time)
#' S.min=list(formula=~min)
#' S.BirdAge=list(formula=~BirdAge)
#' #
#' # Note that in the following model in the MARK example, the covariates
#' # have been standardized.  That means that the beta parameters will be different
#' # for BirdAge, Weight and their interaction but the likelihood and real parameter
#' # estimates are the same.
#' #
#' S.BirdAgexWeight.min=list(formula=~min+BirdAge*Weight)
#' S.BirdAge.Weight=list(formula=~BirdAge+Weight)
#' #
#' # Create model list and run assortment of models
#' #
#' model.list=create.model.list("Known")
#' bduck.results=mark.wrapper(model.list,data=bduck.processed,ddl=bduck.ddl,
#'                invisible=FALSE,threads=1)
#' 
#' #
#' # Return model table and list of models
#' #
#' return(bduck.results)
#' }
#' bduck.results=run.Blackduck()
#' bduck.results
#' 
#' 
#' 
NULL

#' White-tailed deer double observer spotlight capture-recapture analysis
#' 
#' This data represents a set of independent double observer road-transect survey data of white-tailed deer
#' on Brosnan Forest, South Carolina surveyed in August, 2005-2009.  The primary reason for
#' this package is to provide a completely reproducible example of the analysis from Collier et al. (2012).  We used a
#' Huggins closed capture model implemented in MARK \url{http://www.phidot.org/software/mark/} via RMark 
#' both of which will need to be installed on the system to use this package.  The data have 2 time periods (primary observer (t1) was a thermal imager, secondary observer (t2) was
#' a spotlight observer in the same vehicle on the same side) with the primary objective of the study being to evaluate 
#' the detection (recapture) rates of white-tailed deer using spotlights as a survey method. 
#'
#' @details  In addition to detailing the analysis used by Collier et al. (2012), this example documents the 
#' use of the \code{share} argument in the RMark parameter specification because there is presently very little
#' documentation on the use of \code{share}. Parameters in MARK models rarely share columns of the design matrix. For
#' example while you might want to use the same covariate for survival and capture probability, you would never use the
#' same beta (same column of the design matrix) for each parameter.  However, there are exceptions when the parameters 
#' represent similar quantities and that is when the \code{share} argument is useful.  For example, in the closed capture models
#' p is initial capture probability and c is recapture probability. In this case, it would make perfect sense to use the same
#' column of the design matrix for both parameters. The most obvious case is to fit a model in which p=c.
#'
#' In RMark, certain pairs of parameters have been identified as similar and shareable.  These can be found in the file parameters.txt
#' which is in the RMark directory in your R library. With each pair that is shareable, the first one listed is the primary parameter.
#' When you want to share columns in the design matrix, share=TRUE is added to the specification of the primary parameter. A parameter
#' specification is not given for the other secondary parameter when they are shared.  When RMark, sees that the parameters are to be shared it
#' creates a pooled set of design data and adds a column with the name of the secondary parameter and its value is 0 for the rows
#' for the primary parameter and 1 for the rows for the secondary parameter. For example, with the closed capture model if share=TRUE is
#' added to the parameter specification for p, a model is not specified for c, and the pooled design data set contains a field called c.
#' The added field allows construction of models where there are restricted differences between the parameters. For example, p=list(formula=~time+c,share=TRUE)
#' will fit a model in which capture probability varies by time and recapture probability includes an additive difference on the link scale.
#' Because the design data are pooled when you share parameters, if you modify design data for one of the parameters, the other most be modified as
#' as well, so the columns of the design data for both parameters are the same or RMark will give an error.
#' 
#' The argument \code{share} is used in all the candidate models in the below example analysis.  As a simplified example of how 
#' \code{share} works, look at the candidate models in the \code{bfrun{}} function call named \code{mod.2} and \code{mod.2a} (note that
#' \code{mod.2a} was not included in the supplemental file available from the Journal of Wildlife Management and is only included in this package).  Both
#' of these models are conducting the exact same analysis, with the first \code{mod.2}, we used the formula \code{~time}  (if you don't 
#' know what this means go read the MARKBOOK at \url{http://www.phidot.org/software/mark/}.  Notice, however, we used the argument
#' \code{share} in \code{mod.2}, which tells RMark to share columns of the MARK design matrix.  For comparison, so you can evaluate how 
#' \code{share} works for yourself, \code{mod.2a} recreates the same analysis as \code{mod.2}, but uses the approach more typical to MARK 
#' analyses where each parameter is specified independently and uniquely.
#'
#' @name deer
#' @docType data
#' @format The format is a data frame with 4508 observations on the following 7 variables.
#' \describe{ \item{SL (spotlight)}{0/1 whether deer was missed/seen by the spotlight observer}
#' \item{TI (thermal imager)}{0/1 whether deer was missed/seen by the thermal imager observer}
#' \item{Group}{Factor with 79 levels representing each unique paired (TI-SL) survey conducted}
#' \item{Year}{Factor with 5 levels for year of survey}
#' \item{MaxCount}{Count of maximum number of deer seen for each survey, only needed
#' for bootstrapp analysis in MARK, not used in bfdeeR package}
#' \item{Cluster}{Value assigning each deer to a specific observation cluster, only needed
#' for bootstrapp analysis in MARK, not used in bfdeeR package}
#' \item{MgmtUnit}{Management unit identification} }
#' @references Collier, B. A., S. S. Ditchkoff, J. B. Raglin, and C. R. Ruth.  2012. Spotlight surveys 
#' for white-tailed deer: monitoring panacea or exercise in futility? Journal of Wildlife Management, In Press.  
#' @keywords datasets
#' @author Bret Collier
#' @examples 
#' 
#' \donttest{
#' # This example is excluded from testing to reduce package check time
#' data(deer)
#' x=data.frame(ch=paste(deer$TI, deer$SL, sep=""), Survey=factor(deer$Group), 
#'      Year=factor(deer$Year), Cluster=deer$Cluster, MgtUnit=factor(deer$MgmtUnit))
#' x$ch=as.character(x$ch)
#' bfrun=function(){
#' x.proc=process.data(x, model="Huggins", groups=c("Survey", "Year", "MgtUnit"))
#' x.ddl=make.design.data(x.proc)
#' 
#' #Silly Null model, constant p & c sharing 1 parameter (one detection estimate)
#' p.shared=list(formula=~1,share=TRUE)
#' mod.1=mark(x.proc, x.ddl, model.parameters=list(p=p.shared), invisible=FALSE)
#'  
#' #2 Parameter Null Model, constant p, constant c, different p and c (one estimate for each; p ne c)
#' #p(time), c(-), share=TRUE, detection is time dependent, with recapture parameter shared
#' p.sharetime=list(formula=~time, share=TRUE)
#' mod.2=mark(x.proc, x.ddl, model.parameters=list(p=p.sharetime), invisible=FALSE)
#' 
#' #2a Parameter Null Model, constant p, constant c,
#' # different p and c (one estimate for each; p ne c) not using share
#' mod.2a=mark(x.proc, x.ddl, model.parameters=list(p=list(formula=~1), c=list(formula=~1)))
#' 
#' #Fully parameterized model, different p and c for each survey transect replicate, 
#' # management unit, method (TI or SL) and any observers
#' p.survey=list(formula=~Survey*time, share=TRUE)
#' mod.3=mark(x.proc, x.ddl, model.parameters=list(p=p.survey), invisible=FALSE)
#' 
#' #p(MU), c(MU), initial detection and recapture differ and are management unit dependent
#' p.mu=list(formula=~MgtUnit*time, share=TRUE)
#' mod.4=mark(x.proc, x.ddl, model.parameters=list(p=p.mu), invisible=FALSE)
#' 
#' #p(MU) detection is management unit dependent
#' p.mu=list(formula=~MgtUnit, share=TRUE)
#' mod.5=mark(x.proc, x.ddl, model.parameters=list(p=p.mu), invisible=FALSE)
#' 
#' #p(Yr + MgtUnit),  detection is year + MgtUnit
#' p.yearMgtUnit=list(formula=~Year*time+MgtUnit, share=TRUE)
#' mod.6=mark(x.proc, x.ddl, model.parameters=list(p=p.yearMgtUnit), invisible=FALSE)
#' 
#' #p(Year), initial detection and recapture are year dependent
#' p.year=list(formula=~Year*time, share=TRUE)
#' mod.7=mark(x.proc, x.ddl, model.parameters=list(p=p.year), invisible=FALSE)
#' 
#' return(collect.models())
#' }
#' bf.out=bfrun()
#' bf.out
#' 
#' #export function to send dataset and covariates data to MARK for bootstrap analysis 
#' #(not run but here for completeness)
#' #export.MARK(x.proc, "BFdeer", mod.3, replace=TRUE, ind.covariates="all")
#' }
NULL




#' Example data for Closed Robust Design Multistrata
#' 
#' Data and Script to simulate the MSCRD example of 15.7.1 from the MARK book
#' Cooch and White
#' 
#' 
#' @name crdms
#' @docType data
#' @format A data frame with 557 observations on the following 2 variables.
#' \describe{ \item{ch}{a character vector of encounter histories}
#' \item{freq}{a numeric vector of frequencies of each history} }
#' @references For Cooch and White book see
#' \url{http://www.phidot.org/software/mark/}
#' @source This example was constructed by Andrew Paul who is with Fish and
#' Wildlife Division of the Alberta Provincial Government
#' @keywords datasets
#' @examples
#' \donttest{
#' # This example is excluded from testing to reduce package check time
#' #Script to simulate the MSCRD 
#' #example of 15.7.1 from the MARK
#' #book
#' #created by AJP 21 Dec 2010
#' 
#' 
#' #convert .inp data - only needed to create crdms
#' #ch.data<-convert.inp("rd_simple1.inp")
#' data(crdms)
#' #set time intervals
#' #4 primary periods each with 3 secondary occasions
#' t.int<-c(rep(c(0,0,1),3),c(0,0))
#' 
#' #process data for RMark
#' crdms.data<-process.data(crdms,model="CRDMS",time.interval=t.int,
#' 			strata.labels=c("1","U"))
#' #change Psi parameters that are obtained by subtraction
#' crdms.ddl<-make.design.data(crdms.data,
#' 		parameters=list(Psi=list(subtract.stratum=c("1","1"))))
#' 
#' #create grouping index for unobserved p and c (i.e., always zero)
#' up=as.numeric(row.names(crdms.ddl$p[crdms.ddl$p$stratum=="U",]))
#' 
#' #create grouping index to fix Psi for unobs to unbos at time 1
#' #this isn't necessary but it allows this Psi to be fixed to a value
#' #that can be flagged and not erroneously interpreted
#' Psiuu1=as.numeric(row.names(crdms.ddl$Psi[crdms.ddl$Psi$stratum=="U"&
#' 			crdms.ddl$Psi$time==1,]))
#' 
#' #create dummy variable for constraining last Psi in Markovian model
#' #variable is called ctime for constrained time
#' crdms.ddl$Psi$ctime=crdms.ddl$Psi$time
#' crdms.ddl$Psi$ctime[crdms.ddl$Psi$time==3]=2
#' do_example=function()
#' {
#' #Initial assumptions
#' S.dot=list(formula=~1)  #S equal for both states and constant over time
#' p.session=list(formula=~session, share=TRUE,  #p=c varies with session 
#' 		fixed=list(index=up,value=0)) #p set to zero for unobs
#' 
#' #Model 1 - Markovian movement
#' Psi.markov=list(formula=~ctime+stratum,
#' 		fixed=list(index=Psiuu1,value=9e-99)) #9e-99 is a flag
#' model.1=mark(crdms.data,ddl=crdms.ddl,
#' 	model.parameters=list(S=S.dot,
#' 		p=p.session,
#' 		Psi=Psi.markov),threads=2)
#' 
#' #Model 2 - Random movement
#' Psi.rand=list(formula=~time)
#' model.2=mark(crdms.data,ddl=crdms.ddl,
#' 	model.parameters=list(S=S.dot,
#' 		p=p.session,
#' 		Psi=Psi.rand),threads=2)
#' 
#' #Model 3 - No movement
#' Psi.fix=list(formula=~1,fixed=0)
#' model.3=mark(crdms.data,ddl=crdms.ddl,
#' 	model.parameters=list(S=S.dot,
#' 		p=p.session,
#' 		Psi=Psi.fix),threads=2)
#' 		
#' #collect and store models
#' crdms.res<-collect.models()
#' 
#' print(crdms.res)
#' invisible()
#' }
#' do_example()
#' }
NULL


#' Dipper capture-recapture data
#' 
#' A capture-recapture data set on European dippers from France that
#' accompanies MARK as an example analysis using the CJS and POPAN models.  The
#' dipper data set was orginally described as an example by Lebreton et al
#' (1992).
#' 
#' This is a data set that accompanies program MARK as an example for CJS and
#' POPAN analyses.  The data can be stratified using sex as a grouping
#' variable.  The functions \code{run.dipper}, \code{run.dipper.alternate},
#' \code{run.dipper.popan} defined below in the examples mimic the models used
#' in the dbf file that accompanies MARK. Note that the models used in the MARK
#' example use PIM coding with the sin link function which is often better at
#' identifying the number of estimable parameters.  The approach used in the R
#' code uses design matrices and cannot use the sin link and is less capable at
#' counting parameters.  These differences are illustrated by comparing the
#' results of \code{run.dipper} and \code{run.dipper.alternate} which fit the
#' same set of "CJS" models.  The latter fits the models with constraints on
#' some parameters to achieve identifiability and the former does not. Although
#' it does not influence the selection of the best model it does infleunce
#' parameter counts and AIC ordering of some of the less competitive models. In
#' using design matrices it is best to constrain parameters that are confounded
#' (e.g., last occasion parameters in Phi(t)p(t) CJS model) when possible to
#' achieve more reliable counts of the number of estimable parameters.  See
#' \code{\link{adjust.parameter.count}} for more dicussion on this point.
#' 
#' Note that the covariate "sex" defined in dipper has values "Male" and
#' "Female".  It cannot be used directly in a formula for MARK without using it
#' do define groups because MARK.EXE will be unable to read in a covariate with
#' non-numeric values.  By using \code{groups="sex"} in the call the
#' \code{\link{process.data}} a factor "sex" field is created that can be used
#' in the formula.  Alternatively, a new covariate could be defined in the data
#' with say values 0 for Female and 1 for Male and this could be used without
#' defining groups because it is numeric.  This can be done easily by
#' translating the values of the coded variables to a numeric variable.  Factor
#' variables are numbered 1..k for k levels in alphabetic order.  Since Female
#' < Male in alphabetic order then it is level 1 and Male is level 2.  So the
#' following will create a numeric sex covariate.
#' 
#' \preformatted{ dipper$numeric.sex=as.numeric(dipper$sex)-1 }
#' 
#' See \code{\link{export.chdata}} for an example that creates a .inp file for
#' MARK with sex being used to describe groups and a numeric sex covariate.
#' 
#' @name dipper
#' @docType data
#' @format A data frame with 294 observations on the following 2 variables.
#' \describe{ \item{ch}{a character vector containing the encounter history of
#' each bird} \item{sex}{the sex of the bird: a factor with levels
#' \code{Female} \code{Male}} }
#' @source Lebreton, J.-D., K. P. Burnham, J. Clobert, and D. R. Anderson.
#' 1992. Modeling survival and testing biological hypotheses using marked
#' animals: case studies and recent advances. Ecol. Monogr. 62:67-118.
#' @keywords datasets
#' @examples
#' \donttest{
#' # This example is excluded from testing to reduce package check time
#' data(dipper)
#' dipper.model=mark(dipper)
#' run.dipper=function()
#' {
#' #
#' # Process data
#' #
#' dipper.processed=process.data(dipper,groups=("sex"))
#' #
#' # Create default design data
#' #
#' dipper.ddl=make.design.data(dipper.processed)
#' #
#' # Add Flood covariates for Phi and p that have different values
#' #
#' dipper.ddl$Phi$Flood=0
#' dipper.ddl$Phi$Flood[dipper.ddl$Phi$time==2 | dipper.ddl$Phi$time==3]=1
#' dipper.ddl$p$Flood=0
#' dipper.ddl$p$Flood[dipper.ddl$p$time==3]=1
#' #
#' #  Define range of models for Phi
#' #
#' Phidot=list(formula=~1)
#' Phitime=list(formula=~time)
#' Phisex=list(formula=~sex)
#' Phisextime=list(formula=~sex+time)
#' Phisex.time=list(formula=~sex*time)
#' PhiFlood=list(formula=~Flood)
#' #
#' #  Define range of models for p
#' #
#' pdot=list(formula=~1)
#' ptime=list(formula=~time)
#' psex=list(formula=~sex)
#' psextime=list(formula=~sex+time)
#' psex.time=list(formula=~sex*time)
#' pFlood=list(formula=~Flood)
#' #
#' # Run assortment of models
#' #
#' dipper.phidot.pdot          =mark(dipper.processed,dipper.ddl,
#'                  model.parameters=list(Phi=Phidot,p=pdot))
#' dipper.phidot.pFlood      	=mark(dipper.processed,dipper.ddl,
#'                  model.parameters=list(Phi=Phidot,p=pFlood))
#' dipper.phidot.psex        	=mark(dipper.processed,dipper.ddl,
#'                  model.parameters=list(Phi=Phidot,p=psex))
#' dipper.phidot.ptime       	=mark(dipper.processed,dipper.ddl,
#'                  model.parameters=list(Phi=Phidot,p=ptime))
#' dipper.phidot.psex.time		=mark(dipper.processed,dipper.ddl,
#'                  model.parameters=list(Phi=Phidot,p=psex.time))
#' dipper.phitime.ptime      	=mark(dipper.processed,dipper.ddl,
#'                  model.parameters=list(Phi=Phitime, p=ptime))
#' dipper.phitime.pdot       	=mark(dipper.processed,dipper.ddl,
#'                  model.parameters=list(Phi=Phitime,p=pdot))
#' dipper.phitime.psex		=mark(dipper.processed,dipper.ddl,
#'                  model.parameters=list(Phi=Phitime,p=psex))
#' dipper.phitime.psex.time	=mark(dipper.processed,dipper.ddl,
#'                  model.parameters=list(Phi=Phitime,p=psex.time))
#' dipper.phiFlood.pFlood    	=mark(dipper.processed,dipper.ddl,
#'                  model.parameters=list(Phi=PhiFlood, p=pFlood))
#' dipper.phisex.pdot        	=mark(dipper.processed,dipper.ddl,
#'                  model.parameters=list(Phi=Phisex,p=pdot))
#' dipper.phisex.psex        	=mark(dipper.processed,dipper.ddl,
#'                  model.parameters=list(Phi=Phisex,p=psex))
#' dipper.phisex.psex.time        	=mark(dipper.processed,dipper.ddl,
#'                  model.parameters=list(Phi=Phisex,p=psex.time))
#' dipper.phisex.ptime       	=mark(dipper.processed,dipper.ddl,
#'                  model.parameters=list(Phi=Phisex,p=ptime))
#' dipper.phisextime.psextime	=mark(dipper.processed,dipper.ddl,
#'                  model.parameters=list(Phi=Phisextime,p=psextime))
#' dipper.phisex.time.psex.time	=mark(dipper.processed,dipper.ddl,
#'                  model.parameters=list(Phi=Phisex.time,p=psex.time))
#' dipper.phisex.time.psex 	=mark(dipper.processed,dipper.ddl,
#'                  model.parameters=list(Phi=Phisex.time,p=psex))
#' dipper.phisex.time.pdot		=mark(dipper.processed,dipper.ddl,
#'                  model.parameters=list(Phi=Phisex.time,p=pdot))
#' dipper.phisex.time.ptime	=mark(dipper.processed,dipper.ddl,
#'                  model.parameters=list(Phi=Phisex.time,p=ptime))
#' #
#' # Return model table and list of models
#' #
#' return(collect.models() )
#' }
#' 
#' dipper.results=run.dipper()
#' 
#' run.dipper.alternate=function()
#' {
#' #
#' # Process data
#' #
#' dipper.processed=process.data(dipper,groups=("sex"))
#' #
#' # Create default design data
#' #
#' dipper.ddl=make.design.data(dipper.processed)
#' #
#' # Add Flood covariates for Phi and p that have different values
#' #
#' dipper.ddl$Phi$Flood=0
#' dipper.ddl$Phi$Flood[dipper.ddl$Phi$time==2 | dipper.ddl$Phi$time==3]=1
#' dipper.ddl$p$Flood=0
#' dipper.ddl$p$Flood[dipper.ddl$p$time==3]=1
#' #
#' #  Define range of models for Phi
#' #
#' Phidot=list(formula=~1)
#' Phitime=list(formula=~time)
#' Phitimec=list(formula=~time,fixed=list(time=6,value=1))
#' Phisex=list(formula=~sex)
#' Phisextime=list(formula=~sex+time)
#' Phisex.time=list(formula=~sex*time)
#' PhiFlood=list(formula=~Flood)
#' #
#' #  Define range of models for p
#' #
#' pdot=list(formula=~1)
#' ptime=list(formula=~time)
#' ptimec=list(formula=~time,fixed=list(time=7,value=1))
#' psex=list(formula=~sex)
#' psextime=list(formula=~sex+time)
#' psex.time=list(formula=~sex*time)
#' psex.timec=list(formula=~sex*time,fixed=list(time=7,value=1))
#' pFlood=list(formula=~Flood)
#' #
#' # Run assortment of models
#' #
#' dipper.phidot.pdot          =mark(dipper.processed,dipper.ddl,
#'                   model.parameters=list(Phi=Phidot,p=pdot))
#' dipper.phidot.pFlood      	=mark(dipper.processed,dipper.ddl,
#'                   model.parameters=list(Phi=Phidot,p=pFlood))
#' dipper.phidot.psex        	=mark(dipper.processed,dipper.ddl,
#'                   model.parameters=list(Phi=Phidot,p=psex))
#' dipper.phidot.ptime       	=mark(dipper.processed,dipper.ddl,
#'                   model.parameters=list(Phi=Phidot,p=ptime))
#' dipper.phidot.psex.time		=mark(dipper.processed,dipper.ddl,
#'                   model.parameters=list(Phi=Phidot,p=psex.time))
#' dipper.phitime.ptimec      	=mark(dipper.processed,dipper.ddl,
#'                   model.parameters=list(Phi=Phitime, p=ptimec))
#' dipper.phitime.pdot       	=mark(dipper.processed,dipper.ddl,
#'                   model.parameters=list(Phi=Phitime,p=pdot))
#' dipper.phitime.psex		=mark(dipper.processed,dipper.ddl,
#'                   model.parameters=list(Phi=Phitime,p=psex))
#' dipper.phitimec.psex.time	=mark(dipper.processed,dipper.ddl,
#'                   model.parameters=list(Phi=Phitimec,p=psex.time))
#' dipper.phiFlood.pFlood    	=mark(dipper.processed,dipper.ddl,
#'                   model.parameters=list(Phi=PhiFlood, p=pFlood))
#' dipper.phisex.pdot        	=mark(dipper.processed,dipper.ddl,
#'                   model.parameters=list(Phi=Phisex,p=pdot))
#' dipper.phisex.psex        	=mark(dipper.processed,dipper.ddl,
#'                   model.parameters=list(Phi=Phisex,p=psex))
#' dipper.phisex.psex.time        	=mark(dipper.processed,dipper.ddl,
#'                   model.parameters=list(Phi=Phisex,p=psex.time))
#' dipper.phisex.ptime       	=mark(dipper.processed,dipper.ddl,
#'                   model.parameters=list(Phi=Phisex,p=ptime))
#' dipper.phisextime.psextime	=mark(dipper.processed,dipper.ddl,
#'                   model.parameters=list(Phi=Phisextime,p=psextime),adjust=FALSE)
#' dipper.phisex.time.psex.timec	=mark(dipper.processed,dipper.ddl,
#'                   model.parameters=list(Phi=Phisex.time,p=psex.timec))
#' dipper.phisex.time.psex 	=mark(dipper.processed,dipper.ddl,
#'                   model.parameters=list(Phi=Phisex.time,p=psex))
#' dipper.phisex.time.pdot		=mark(dipper.processed,dipper.ddl,
#'                   model.parameters=list(Phi=Phisex.time,p=pdot))
#' dipper.phisex.time.ptimec	=mark(dipper.processed,dipper.ddl,
#'                   model.parameters=list(Phi=Phisex.time,p=ptimec))
#' #
#' # Return model table and list of models
#' #
#' return(collect.models() )
#' }
#' dipper.results.alternate=run.dipper.alternate()
#' #
#' # Merge two sets of models into a single model list and include the 
#' # initial model as a demo for merge.mark
#' #
#' dipper.cjs=merge.mark(dipper.results,dipper.results.alternate,dipper.model)
#' dipper.cjs
#' #
#' # next delete some of the models to show how this is done with remove.mark
#' #
#' dipper.cjs=remove.mark(dipper.cjs,c(2,4,9))
#' dipper.cjs
#' 
#' run.dipper.popan=function()
#' {
#' #
#' # Process data
#' #
#' dipper.processed=process.data(dipper,model="POPAN",group="sex")
#' #
#' # Create default design data
#' #
#' dipper.ddl=make.design.data(dipper.processed)
#' #
#' # Add Flood covariates for Phi and p that have different values
#' #
#' dipper.ddl$Phi$Flood=0
#' dipper.ddl$Phi$Flood[dipper.ddl$Phi$time==2 | dipper.ddl$Phi$time==3]=1
#' dipper.ddl$p$Flood=0
#' dipper.ddl$p$Flood[dipper.ddl$p$time==3]=1
#' #
#' #  Define range of models for Phi
#' #
#' Phidot=list(formula=~1)
#' Phitime=list(formula=~time)
#' Phisex=list(formula=~sex)
#' Phisextime=list(formula=~sex+time)
#' Phisex.time=list(formula=~sex*time)
#' PhiFlood=list(formula=~Flood)
#' #
#' #  Define range of models for p
#' #
#' pdot=list(formula=~1)
#' ptime=list(formula=~time)
#' psex=list(formula=~sex)
#' psextime=list(formula=~sex+time)
#' psex.time=list(formula=~sex*time)
#' pFlood=list(formula=~Flood)
#' #
#' #  Define range of models for pent
#' #
#' pentsex.time=list(formula=~sex*time)
#' #
#' #  Define range of models for N
#' #
#' Nsex=list(formula=~sex)
#' #
#' # Run assortment of models
#' #
#' dipper.phisex.time.psex.time.pentsex.time=mark(dipper.processed,dipper.ddl,
#' model.parameters=list(Phi=Phisex.time,p=psex.time,pent=pentsex.time,N=Nsex),
#' invisible=FALSE,adjust=FALSE)
#' dipper.phisex.time.psex.pentsex.time=mark(dipper.processed,dipper.ddl,
#' model.parameters=list(Phi=Phisex.time,p=psex,pent=pentsex.time,N=Nsex),
#' invisible=FALSE,adjust=FALSE)
#' #
#' # Return model table and list of models
#' #
#' return(collect.models() )
#' }
#' 
#' dipper.popan.results=run.dipper.popan()
#'
#'# *****************************************************************
#'# Here is an example of user specified links for each real parameter
#'  data(dipper)
#'  dipper.proc=process.data(dipper)
#'  dipper.ddl=make.design.data(dipper.proc)
#'# dummy run of make.mark.model to get links and design data. 
#'# parm.specific set to TRUE so it will create a link for 
#'# each parameter because for this model they are all the
#'# same (logit) and if this was not specified you'ld get a vector with one element
#'  dummy=make.mark.model(dipper.proc,dipper.ddl,simplify=FALSE,parm.specific=TRUE)
#'  input.links=dummy$links
#'# get model indices for p where time=4
#'  log.indices=dipper.ddl$p$model.index[dipper.ddl$p$time==4]
#'# assign those links to log
#'  input.links[log.indices]="Log"
#'# Now these can be used with any call to mark
#'  mymodel=mark(dipper.proc,dipper.ddl,input.links=input.links)
#'  summary(mymodel)
#' }
NULL



#' Exercise 7 example data
#' 
#' An example occupancy data set used as exercise 7 in the occupancy website
#' developed by Donovan and Hines.
#' 
#' This is a data set from exercise 7 of Donovan and Hines occupancy web site
#' (\url{http://www.uvm.edu/rsenr/vtcfwru/spreadsheets/?Page=occupancy/occupancy.htm}).
#' 
#' @name Donovan.7
#' @docType data
#' @format A data frame with 20 observations (sites) on the following 2
#' variables.  \describe{ \item{ch}{a character vector containing the presence
#' (1) and absence (0) for each of 5 visits to the site} \item{freq}{frequency
#' of sites (always 1)} }
#' @keywords datasets
#' @examples
#' 
#' # Donovan.7 can be created with
#' # Donovan.7=convert.inp("Donovan.7.inp")
#' 
#' do.exercise.7=function()
#' {
#'   data(Donovan.7)
#' # Estimates from following agree with estimates on website but the
#' # log-likelihood values do not agree.  Maybe a difference in whether the
#' # constant binomial coefficients are included.
#'   Donovan.7.poisson=mark(Donovan.7,model="OccupRNPoisson",invisible=FALSE,threads=1)
#' # THe following model was not in exercise 7.
#'   Donovan.7.negbin=mark(Donovan.7,model="OccupRNNegBin",invisible=FALSE,threads=1)
#'   return(collect.models())
#' }
#' exercise.7=do.exercise.7()
#' # Remove # to see output
#' # print(exercise.7)
#' 
#' 
#' 
NULL





#' Exercise 8 example data
#' 
#' An example occupancy data set used as exercise 8 in the occupancy website
#' developed by Donovan and Hines.
#' 
#' This is a data set from exercise 8 of Donovan and Hines occupancy web site
#' (\url{http://www.uvm.edu/rsenr/vtcfwru/spreadsheets/?Page=occupancy/occupancy.htm}).
#' In MARK, it uses 2 digits to allow a count of 0 to 99 at each site, so the
#' history has 10 digits for 5 visits (occasions).
#' 
#' @name Donovan.8
#' @docType data
#' @format A data frame with 20 observations (sites) on the following 2
#' variables.  \describe{ \item{ch}{a character vector containing the counts
#' for each of 5 visits to the site} \item{freq}{frequency of sites (always 1)}
#' }
#' @keywords datasets
#' @examples
#' \donttest{
#' # This example is excluded from testing to reduce package check time
#' # Donovan.8 can be created with
#' # Donovan.8=convert.inp("Donovan.8.inp")
#' do.exercise.8=function()
#' {
#'   data(Donovan.8)
#' # Results agree with the values on the website.
#'   Donovan.8.poisson=mark(Donovan.8,model="OccupRPoisson",invisible=FALSE,threads=2)
#' # The following model was not in exercise 8. The NegBin model does 
#' # better if it is initialized with the r and lambda from the poisson.
#'   Donovan.8.negbin=mark(Donovan.8,model="OccupRNegBin",
#'     initial=Donovan.8.poisson,invisible=FALSE,threads=2)
#'   return(collect.models())
#' }
#' exercise.8=do.exercise.8()
#' # Remove # to see output
#' # print(exercise.8)
#' }
#' 
NULL





#' Rabbit capture-recapture data
#' 
#' A capture-recapture data set on rabbits derived from Edwards and Eberhardt
#' (1967) that accompanies MARK as an example analysis using the closed
#' population models.
#' 
#' This data set is used in MARK to illustrate the various closed population
#' models including "Closed", "HetClosed", "FullHet","Huggins","HugHet", and
#' "FullHugHet".  The first 3 include N in the likelihood whereas the last 3
#' are based on the Huggins approach which does not use N in the likelihood.
#' The Het... and FullHet... models are based on the Pledger mixture model
#' approach. Some of the examples demonstrate the use of the \code{share}
#' argument in the \code{model.parameters} list for parameter \code{p} which
#' allows sharing common values for \code{p} and \code{c}.
#' 
#' @name edwards.eberhardt
#' @docType data
#' @format A data frame with 76 observations on the following variable.
#' \describe{ \item{ch}{a character vector} }
#' @source Edwards, W.R. and L.L. Eberhardt 1967.  Estimating cottontail
#' abundance from live trapping data. J. Wildl. Manage. 31:87-96.
#' @keywords datasets
#' @examples
#' \donttest{
#' # This example is excluded from testing to reduce package check time
#' #
#' # get data
#' #
#' data(edwards.eberhardt)
#' #
#' # create function that defines and runs the analyses as defined in 
#' # MARK example dbf file
#' #
#' run.edwards.eberhardt=function()
#' {
#' #
#' #  Define parameter models
#' #
#' pdotshared=list(formula=~1,share=TRUE)
#' ptimeshared=list(formula=~time,share=TRUE)
#' ptime.c=list(formula=~time+c,share=TRUE)
#' ptimemixtureshared=list(formula=~time+mixture,share=TRUE)
#' pmixture=list(formula=~mixture)
#' #
#' # Run assortment of models
#' #
#' #
#' #   Capture Closed models
#' #
#' #  constant p=c
#' ee.closed.m0=mark(edwards.eberhardt,model="Closed",
#'                    model.parameters=list(p=pdotshared))
#' #  constant p and constant c but different
#' ee.closed.m0c=mark(edwards.eberhardt,model="Closed")
#' #  time varying p=c
#' ee.closed.mt=mark(edwards.eberhardt,model="Closed",
#'                    model.parameters=list(p=ptimeshared))
#' #
#' #  Closed heterogeneity models
#' #
#' #  2 mixtures Mh2
#' ee.closed.Mh2=mark(edwards.eberhardt,model="HetClosed",
#'                    model.parameters=list(p=pmixture))
#' #  Closed Mth2 - p different for time; mixture additive
#' ee.closed.Mth2.additive=mark(edwards.eberhardt,model="FullHet",
#'                    model.parameters=list(p=ptimemixtureshared),adjust=TRUE)
#' #
#' #    Huggins models
#' #
#' # p=c constant over time
#' ee.huggins.m0=mark(edwards.eberhardt,model="Huggins",
#'                    model.parameters=list(p=pdotshared))
#' # p constant c constant but different; this is default model for Huggins
#' ee.huggins.m0.c=mark(edwards.eberhardt,model="Huggins")
#' # Huggins Mt
#' ee.huggins.Mt=mark(edwards.eberhardt,model="Huggins",
#'                    model.parameters=list(p=ptimeshared),adjust=TRUE)
#' #
#' #    Huggins heterogeneity models
#' #
#' #  Mh2 - p different for mixture
#' ee.huggins.Mh2=mark(edwards.eberhardt,model="HugHet",
#'                    model.parameters=list(p=pmixture))
#' #  Huggins Mth2 - p different for time; mixture additive
#' ee.huggins.Mth2.additive=mark(edwards.eberhardt,model="HugFullHet",
#'                    model.parameters=list(p=ptimemixtureshared),adjust=TRUE)
#' #
#' # Return model table and list of models
#' #
#' return(collect.models() )
#' }
#' #
#' # fit models in mark by calling function created above
#' #
#' ee.results=run.edwards.eberhardt()
#' }
NULL





#' Simulated data from Cormack-Jolly-Seber model
#' 
#' A simulated data set from CJS model to demonstrate use of grouping variables
#' and individual covariates in an analysis of a \code{mark} model. The true
#' model for the simulated data is S(age+year)p(year).
#' 
#' The \code{weight}, \code{age} and \code{region} are static variables that
#' are defined based on the values when the animal was released. \code{age} is
#' a factor variable representing an age level.  The actual ages (at time of
#' release) are 0,1,2 for the 3 levels respectively.
#' 
#' @name example.data
#' @docType data
#' @format A data frame with 6000 observations on the following 5 variables.
#' \describe{ \item{ch}{a character vector} \item{weight}{a numeric vector}
#' \item{age}{a factor with levels \code{1} \code{2} \code{3}} \item{sex}{a
#' factor with levels \code{F} \code{M}} \item{region}{a factor with levels
#' \code{1} \code{2} \code{3} \code{4}} }
#' @keywords datasets
#' @examples
#' \donttest{
#' # This example is excluded from testing to reduce package check time
#' data(example.data)
#' run.example=function()
#' {
#' PhiTime=list(formula=~time)
#' pTimec=list(formula=~time,fixed=list(time=7,value=1))
#' pTime=list(formula=~time)
#' PhiAge=list(formula=~age)
#' Phidot=list(formula=~1)
#' PhiweightTime=list(formula=~weight+time)
#' PhiTimeAge=list(formula=~time+age)
#' mod1=mark(example.data,groups=c("sex","age","region"),
#'                            initial.ages=c(0,1,2))
#' mod2=mark(example.data,model.parameters=list(p=pTimec,Phi=PhiTime),
#'           groups=c("sex","age","region"),initial.ages=c(0,1,2))
#' mod3=mark(example.data,model.parameters=list(Phi=Phidot,p=pTime),
#'           groups=c("sex","age","region"),initial.ages=c(0,1,2))
#' mod4=mark(example.data,model.parameters=list(Phi=PhiTime),
#'           groups=c("sex","age","region"),initial.ages=c(0,1,2))
#' mod5=mark(example.data,model.parameters=list(Phi=PhiTimeAge),
#'           groups=c("sex","age","region"),initial.ages=c(0,1,2))
#' mod6=mark(example.data,model.parameters=list(Phi=PhiAge,p=pTime),
#'           groups=c("sex","age","region"),initial.ages=c(0,1,2))
#' mod7=mark(example.data,model.parameters=list(p=pTime,Phi=PhiweightTime),
#'           groups=c("sex","age","region"),initial.ages=c(0,1,2))
#' mod8=mark(example.data,model.parameters=list(Phi=PhiTimeAge,p=pTime),
#'           groups=c("sex","age","region"),initial.ages=c(0,1,2))
#' return(collect.models())
#' }
#' example.results=run.example()
#' }
NULL





#' Example of Immigration-Emigration LogitNormal Mark-Resight model
#' 
#' Data and example illustrating Immigration-Emigration LogitNormal
#' Mark-Resight model
#' 
#' 
#' @name IELogitNormalMR
#' @docType data
#' @format A data frame with 34 observations on the following variable.
#' \describe{ \item{ch}{a character vector} }
#' @keywords datasets
#' @examples
#' \donttest{
#' # This example is excluded from testing to reduce package check time
#' data(IELogitNormalMR)
#' IElogitNor.proc=process.data(IELogitNormalMR,model="IELogitNormalMR",
#' 	counts=list("Marked Superpopulation"=c(28, 29, 30, 30, 30, 33, 33, 33, 33, 34, 34, 34),
#' 	"Unmarked Seen"=c(264, 161, 152, 217, 217, 160, 195, 159, 166, 152, 175, 190),
#' 	"Marked Unidentified"=c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)),
#' 		         time.intervals=c(0,0,0,1,0,0,0,1,0,0,0))
#' IElogitNor.ddl=make.design.data(IElogitNor.proc)
#' mod1=mark(IElogitNor.proc,IElogitNor.ddl,
#' 		model.parameters=list(p=list(formula=~-1+session),
#' 				              sigma=list(formula=~session),
#' 							  alpha=list(formula=~-1+session:time),
#' 							  Nstar=list(formula=~session),
#' 							  Nbar=list(formula=~session)))
#' summary(mod1)
#' # You can use the initial value to get a better estimate.
#'  mod2=mark(IElogitNor.proc,IElogitNor.ddl,
#' 		  model.parameters=list(p=list(formula=~-1+session:time),
#' 							  sigma=list(formula=~session),
#' 							  alpha=list(formula=~-1+session:time),
#' 							  Nstar=list(formula=~session),
#' 							  Nbar=list(formula=~session)),
#' 							  initial=mod1)
#'  summary(mod2)			  
#' }
NULL





#' Killdeer nest survival example data
#' 
#' A data set on killdeer that accompanies MARK as an example analysis for the
#' nest survival model.
#' 
#' This is a data set that accompanies program MARK as an example for nest
#' survival. The data structure for the nest survival model is completely
#' different from the capture history structure used for most MARK models.  To
#' cope with these data you must import them into a dataframe using R commands
#' and assign the specific variable names shown above. The id and Freq fields
#' are optional.  Freq is assumed to be 1 if not given.  You cannot import the
#' MARK .inp file structure directly into R without some manipulation.  Also
#' note that \code{\link{import.chdata}} and \code{\link{convert.inp}} do NOT
#' work for nest survival data. In the examples section below, the first
#' section of code provides an example of converting the killdeer.inp file into
#' a dataframe for RMark.
#' 
#' If your dataframe contains a variable AgeDay1, which is the age of the nest
#' on the first occasion then you can use a variable called NestAge which will
#' create a set of time-dependent covariates named NestAge1,NestAge2
#' ...NestAge(nocc-1) which will provide a way to incorporate the age of the
#' nest in the model.  This was added because the age covariate in the design
#' data for S assumes all nests are the same age and is not particularly
#' useful. This effect could be incorporated by using the add() function in the
#' design matrix but RMark does not have any capability for doing that and it
#' is easier to create a time-dependent covariate to do the same thing.
#' 
#' @name killdeer
#' @docType data
#' @format A data frame with 18 observations on the following 6 variables.
#' \describe{ \item{id}{a MARK comment field with a nest id}
#' \item{FirstFound}{the day the nest was first found} \item{LastPresent}{the
#' last day that chicks were present} \item{LastChecked}{the last day the nest
#' was checked} \item{Fate}{the fate of the nest; 0=hatch and 1 depredated}
#' \item{Freq}{the frequency of nests with this data; usually 1} }
#' @keywords datasets
#' @examples
#' \donttest{
#' # This example is excluded from testing to reduce package check time
#' # EXAMPLE CODE FOR CONVERSION OF .INP TO NECESSARY DATA STRUCTURE
#' # read in killdeer.inp file
#' #killdeer=scan("killdeer.inp",what="character",sep="\n")
#' # strip out ; and write out all but first 2 lines which contain comments
#' #write(sub(";","",killdeer[3:20]),"killdeer.txt")
#' # read in as a dataframe and assign names
#' #killdeer=read.table("killdeer.txt")
#' #names(killdeer)=c("id","FirstFound","LastPresent","LastChecked","Fate","Freq")
#' #
#' # EXAMPLE CODE TO RUN MODELS CONTAINED IN THE MARK KILLDEER.DBF
#' data(killdeer)
#' # produce summary
#' summary(killdeer)
#' # Define function to run models that are in killdeer.dbf
#' # You must specify either the number of occasions (nocc) or the time.intervals 
#' # between the occasions.
#' run.killdeer=function()
#' {
#'    Sdot=mark(killdeer,model="Nest",nocc=40)
#'    STime=mark(killdeer,model="Nest",
#'        model.parameters=list(S=list(formula=~I(Time+1))),nocc=40,threads=2)
#'    STimesq=mark(killdeer,model="Nest",
#'        model.parameters=list(S=list(formula=~I(Time+1)+I((Time+1)^2))),nocc=40,threads=2)
#'    STime3=mark(killdeer,model="Nest",
#'       model.parameters=list(S=list(formula=~I(Time+1)+I((Time+1)^2)+I((Time+1)^3))),
#'                    nocc=40,threads=2)
#'    return(collect.models())
#' }
#' # run defined models
#' killdeer.results=run.killdeer()
#' }
NULL





#' Example of LogitNormal Mark-Resight model
#' 
#' Data and example illustrating LogitNormal Mark-Resight model.
#' 
#' 
#' @name LogitNormalMR
#' @docType data
#' @format A data frame with 35 observations on the following variable.
#' \describe{ \item{ch}{a character vector} }
#' @keywords datasets
#' @examples
#' 
#' data(LogitNormalMR)
#' logitNor.proc=process.data(LogitNormalMR,model="LogitNormalMR",
#' 		counts=list("Unmarked seen"=c(96,68,59),
#' 				    "Marked Unidentified"=c(0,0,0,0,1,1,1,0,0,3,0,1)),
#' 			         time.intervals=c(0,0,0,1,0,0,0,1,0,0,0))
#' logitNor.ddl=make.design.data(logitNor.proc)
#' mod=mark(logitNor.proc,logitNor.ddl,threads=1)
#' summary(mod)
#' 
NULL





#' Mallard nest survival example
#' 
#' A nest survival data set on mallards.  The data set and analysis is
#' described by Rotella et al.(2004).
#' 
#' The first 5 fields must be named as they are shown for nest survival models.
#' \code{Freq} and the remaining fields are optional.  See
#' \code{\link{killdeer}} for more description of the nest survival data
#' structure and the use of the special field \code{AgeDay1}. The habitat
#' variables \code{Native},\code{Planted},\code{Wetland},\code{Roadside} were
#' originally coded as 0/1 dummy variables to enable easier modelling with
#' MARK.  A better alternative in RMark is to create a single variable
#' \code{habitat} with values of \code{"Native"},\code{"Planted"},
#' \code{"Wetland"}, and \code{"Roadside"}. Then the Hab model in the example
#' below would become:
#' 
#' \preformatted{ Hab=mark(mallard,nocc=90,model="Nest",
#' model.parameters=list(S=list(formula=~habitat)), groups="habitat") } For
#' this example, that doesn't make a big difference but if you have more than
#' one factor and you want to construct an interaction or you have a factor
#' with many levels, then it is more efficient to use factor variables rather
#' than dummy variables.
#' 
#' @name mallard
#' @docType data
#' @format A data frame with 565 observations on the following 13 variables.
#' \describe{ \item{FirstFound}{the day the nest was first found}
#' \item{LastPresent}{the last day that chicks were present}
#' \item{LastChecked}{the last day the nest was checked} \item{Fate}{the fate
#' of the nest; 0=hatch and 1 depredated} \item{Freq}{the frequency of nests
#' with this data; always 1 in this example} \item{Robel}{Robel reading of
#' vegetation thickness} \item{PpnGrass}{proportion grass in vicinity of nest}
#' \item{Native}{dummy 0/1 variable; 1 if native vegetation}
#' \item{Planted}{dummy 0/1 variable; 1 if planted vegetation}
#' \item{Wetland}{dummy 0/1 variable; 1 if wetland vegetation}
#' \item{Roadside}{dummy 0/1 variable; 1 if roadside vegetation}
#' \item{AgeFound}{age of nest in days the day the nest was found}
#' \item{AgeDay1}{age of nest at beginning of study} }
#' @author Jay Rotella
#' @source Rotella, J.J., S. J. Dinsmore, T.L. Shaffer.  2004.  Modeling
#' nest-survival data: a comparison of recently developed methods that can be
#' implemented in MARK and SAS.  Animal Biodiversity and Conservation
#' 27:187-204.
#' @keywords datasets
#' @examples
#' 
#'# Last updated June 2, 2011
#'# Read in data, which are in a simple text file that
#'# looks like a MARK input file but (1) with no comments or semicolons and
#'# (2) with a 1st row that contains column labels
#'# mallard=read.table("mallard.txt",header=TRUE)
#'
#'# The mallard data set is also incuded with RMark and can be retrieved with
#'# data(mallard)

#'
#'#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#'# Example of use of RMark for modeling nest survival data -   #
#'# Mallard nests example                                       #
#'# The example runs the 9 models that are used in the Nest     #
#'# Survival chapter of the Gentle Introduction to MARK and that# 
#'# appear in Table 3 (page 198) of                             #
#'# Rotella, J.J., S. J. Dinsmore, T.L. Shaffer.  2004.         #
#'# Modeling nest-survival data: a comparison of recently       #
#'# developed methods that can be implemented in MARK and SAS.  #
#'#   Animal Biodiversity and Conservation 27:187-204.          #
#'#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#' \donttest{
#' # This example is excluded from testing to reduce package check time
#' require(RMark)
#' require(plotrix)
#' 
#'# Retrieve data
#'data(mallard)
#'
#'# Treat dummy variables for habitat types as factors
#'mallard$Native=factor(mallard$Native)
#'mallard$Planted=factor(mallard$Planted)
#'mallard$Wetland=factor(mallard$Wetland)
#'mallard$Roadside=factor(mallard$Roadside)
#'
#'# Examine a summary of the dataset
#'summary(mallard)
#'
#'# Write a function for evaluating a set of competing models
#'run.mallard=function()
#'{
#'# 1. A model of constant daily survival rate (DSR)
#'	Dot=mark(mallard,nocc=90,model="Nest",
#'			model.parameters=list(S=list(formula=~1)))
#'	
#'# 2. DSR varies by habitat type - treats habitats as factors
#'#  and the output provides S-hats for each habitat type
#'	Hab=mark(mallard,nocc=90,model="Nest",
#'			model.parameters=list(S=list(formula=~Native+Planted+Wetland)),
#'			groups=c("Native","Planted","Wetland"))
#'	
#'# 3. DSR varies with vegetation thickness (Robel reading)
#'# Note: coefficients are estimated using the actual covariate
#'# values. They are not based on standardized covariate values.
#'	Robel=mark(mallard,nocc=90,model="Nest",
#'			model.parameters=list(S=list(formula=~Robel)))
#'	
#'# 4. DSR varies with the amount of native vegetation in the surrounding area
#'# Note: coefficients are estimated using the actual covariate
#'# values. They are not based on standardized covariate values.
#'	PpnGr=mark(mallard,nocc=90,model="Nest",
#'			model.parameters=list(S=list(formula=~PpnGrass)))
#'	
#'# 5. DSR follows a trend through time
#'	TimeTrend=mark(mallard,nocc=90,model="Nest",
#'			model.parameters=list(S=list(formula=~Time)))
#'	
#'# 6. DSR varies with nest age
#'	Age=mark(mallard,nocc=90,model="Nest",
#'			model.parameters=list(S=list(formula=~NestAge)))
#'	
#'# 7. DSR varies with nest age & habitat type
#'	AgeHab=mark(mallard,nocc=90,model="Nest",
#'			model.parameters=list(S=list(formula=~NestAge+Native+Planted+Wetland)),
#'			groups=c("Native","Planted","Wetland"))
#'	
#'# 8. DSR varies with nest age & vegetation thickness
#'	AgeRobel=mark(mallard,nocc=90,model="Nest",
#'			model.parameters=list(S=list(formula=~NestAge+Robel)))
#'	
#'# 9. DSR varies with nest age & amount of native vegetation in surrounding area
#'	AgePpnGrass=mark(mallard,nocc=90,model="Nest",
#'			model.parameters=list(S=list(formula=~NestAge+PpnGrass)))
#'	
#'#
#'# Return model table and list of models
#'#
#'	return(collect.models() )
#'}
#'
#'# The next line runs the 9 models above and takes a minute or 2
#'mallard.results=run.mallard()  
#'#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#'# Examine table of model-selection results #
#'#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#'
#'export.MARK(mallard.results$Age$data,"MallDSR",mallard.results,replace=TRUE,ind.covariates="all")
#'mallard.results                        # print model-selection table to screen
#'options(width=100)                     # set page width to 100 characters
#'sink("results.table.txt")              # capture screen output to file
#'print(mallard.results)                 # send output to file
#'sink()                                 # return output to screen
#'
#'# remove "#" on next line to see output in notepad                                            
#'# system("notepad results.table.txt",invisible=FALSE,wait=FALSE) 
#'
#'#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#'# Examine output for constant DSR model #
#'#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#'# Remove "#" on next line to see output
#'# mallard.results$Dot                  # print MARK output to designated text editor
#'mallard.results$Dot$results$beta       # view estimated beta for model in R
#'mallard.results$Dot$results$real       # view estimated DSR estimate in R
#'
#'#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#'# Examine output for 'DSR by habitat' model #
#'#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#'# Remove "#" on next line to see output
#'# mallard.results$Hab                  # print MARK output to designated text editor
#'mallard.results$Hab$design.matrix      # view the design matrix that was used
#'mallard.results$Hab$results$beta       # view estimated beta for model in R
#'mallard.results$Hab$results$beta.vcv   # view variance-covariance matrix for beta's
#'mallard.results$Hab$results$real       # view the estimates of Daily Survival Rate
#'
#'#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#'# Examine output for best model #
#'#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#'# Remove "#" on next line to see output
#'# mallard.results$AgePpnGrass            # print MARK output to designated text editor
#'mallard.results$AgePpnGrass$results$beta # view estimated beta's in R
#'mallard.results$AgePpnGrass$results$beta.vcv # view estimated var-cov matrix in R
#'
#'# To obtain estimates of DSR for various values of 'NestAge' and 'PpnGrass'
#'#   some work additional work is needed.
#'
#'# Store model results in object with simpler name
#'AgePpnGrass=mallard.results$AgePpnGrass
#'# Build design matrix with ages and ppn grass values of interest
#'# Relevant ages are age 1 to 35 for mallards
#'# For ppngrass, use a value of 0.5
#'fc=find.covariates(AgePpnGrass,mallard)
#'fc$value[1:35]=1:35                    # assign 1:35 to 1st 35 nest ages
#'fc$value[fc$var=="PpnGrass"]=0.1       # assign new value to PpnGrass
#'design=fill.covariates(AgePpnGrass,fc) # fill design matrix with values
#'# extract 1st 35 rows of output
#'AgePpn.survival=compute.real(AgePpnGrass,design=design)[1:35,] 
#'# insert covariate columns
#'AgePpn.survival=cbind(design[1:35,c(2:3)],AgePpn.survival)     
#'colnames(AgePpn.survival)=c("Age","PpnGrass","DSR","seDSR","lclDSR","uclDSR")
#'AgePpn.survival       # view estimates of DSR for each age and PpnGrass combo
#'
#'# Plot results
#'with(data=AgePpn.survival,plot(Age,DSR,'l',ylim=c(0.88,1)))
#'grid()
#'axis.break(axis=2,breakpos=0.879,style='slash')
#'with(data=AgePpn.survival,points(Age,lclDSR,'l',lty=3))
#'with(data=AgePpn.survival,points(Age,uclDSR,'l',lty=3))
#'
#'# assign 17 to 1st 50 nest ages
#'fc$value[1:89]=17                     
#'# assign range of values to PpnGrass
#'fc$value[fc$var=="PpnGrass"]=seq(0.01,0.99,length=89) 
#'design=fill.covariates(AgePpnGrass,fc) # fill design matrix with values
#'AgePpn.survival=compute.real(AgePpnGrass,design=design)
#'# insert covariate columns
#'AgePpn.survival=cbind(design[,c(2:3)],AgePpn.survival)     
#'colnames(AgePpn.survival)=
#'		c("Age","PpnGrass","DSR","seDSR","lclDSR","uclDSR")
#'# view estimates of DSR for each age and PpnGrass combo    
#'AgePpn.survival   
#'
#'# Plot results
#'# open new graphics window for new plot
#'dev.new()                          
#'with(data=AgePpn.survival,plot(PpnGrass,DSR,'l',ylim=c(0.88,1)))
#'grid()
#'axis.break(axis=2,breakpos=0.879,style='slash')
#'with(data=AgePpn.survival,points(PpnGrass,lclDSR,'l',lty=3))
#'with(data=AgePpn.survival,points(PpnGrass,uclDSR,'l',lty=3))
#'
#'# The "rm" command can be used to remove all objects from the .Rdata file.
#'# Cleaning up objects as shown in this script is good programming 
#'# practice and a good idea as long as the computations are not time consuming.  
#'# However, if your analysis is taking several hours or days in MARK then 
#'# clearly you'll want to hang onto the resulting objects in R and you 
#'# won't want to use the following command. It has been commented out for
#'# this example; the "#" on the next line needs to be removed to do the clean up.
#'# rm(list=ls(all=TRUE))
#'
#'# The next line deletes orphaned output files from MARK.  
#'# ?cleanup will give a more complete description of this function.
#'cleanup(ask=FALSE)
#' }
NULL



#' Multistrata example data
#' 
#' An example data set which appears to be simulated data that accompanies MARK
#' as an example analysis using the Multistrata model.
#' 
#' This is a data set that accompanies program MARK as an example for the
#' Multistrata model. The models created by RMark are all "Parm-specific"
#' models by default. The sin link is not allowed because all models are
#' specified via the design matrix. Although you can set links for the
#' parameters, usually the default values are preferable.  See
#' \code{\link{make.mark.model}} for additional help building formula for Psi
#' using the remove.intercept argument.
#' 
#' @name mstrata
#' @docType data
#' @format A data frame with 255 observations on the following 2 variables.
#' \describe{ \item{ch}{a character vector containing the encounter history of
#' each bird with strata} \item{freq}{the number of birds with that capture
#' history} }
#' @keywords datasets
#' @examples
#'  \donttest{
#' # This example is excluded from testing to reduce package check time
#' data(mstrata)
#' run.mstrata=function()
#' {
#' #
#' # Process data
#' #
#' mstrata.processed=process.data(mstrata,model="Multistrata")
#' #
#' # Create default design data
#' #
#' mstrata.ddl=make.design.data(mstrata.processed)
#' #
#' #  Define range of models for S; note that the betas will differ from the output
#' #  in MARK for the ~stratum = S(s) because the design matrix is defined using
#' #  treatment contrasts for factors so the intercept is stratum A and the other
#' #  two estimates represent the amount that survival for B abd C differ from A.
#' #  You can use force the approach used in MARK with the formula ~-1+stratum which
#' #  creates 3 separate Betas - one for A,B and C.
#' #
#' S.stratum=list(formula=~stratum)
#' S.stratumxtime=list(formula=~stratum*time)
#' #
#' #  Define range of models for p
#' #
#' p.stratum=list(formula=~stratum)
#' #
#' #  Define range of models for Psi; what is denoted as s for Psi
#' #  in the Mark example for Psi is accomplished by -1+stratum:tostratum which
#' #  nests tostratum within stratum.  Likewise, to get s*t as noted in MARK you
#' #  want ~-1+stratum:tostratum:time with time nested in tostratum nested in
#' #  stratum.
#' #
#' Psi.s=list(formula=~-1+stratum:tostratum)
#' #
#' # Create model list and run assortment of models
#' #
#' model.list=create.model.list("Multistrata")
#' #
#' # Add on specific model that is paired with fixed p's to remove confounding
#' #
#' p.stratumxtime=list(formula=~stratum*time)
#' p.stratumxtime.fixed=list(formula=~stratum*time,fixed=list(time=4,value=1))
#' model.list=rbind(model.list,c(S="S.stratumxtime",p="p.stratumxtime.fixed",
#'   Psi="Psi.s"))
#' #
#' # Run the list of models
#' #
#' mstrata.results=mark.wrapper(model.list,data=mstrata.processed,ddl=mstrata.ddl,threads=2)
#' #
#' # Return model table and list of models
#' #
#' return(mstrata.results)
#' }
#' mstrata.results=run.mstrata()
#' mstrata.results
#'
#' # Example of reverse Multistratum model
#' data(mstrata)
#' mod=mark(mstrata,model="Multistrata")
#' mod.rev=mark(mstrata,model="Multistrata",reverse=TRUE)
#' Psilist=get.real(mod,"Psi",vcv=TRUE)
#' Psilist.rev=get.real(mod.rev,"Psi",vcv=TRUE)
#' Psivalues=Psilist$estimates
#' Psivalues.rev=Psilist.rev$estimates
#' TransitionMatrix(Psivalues[Psivalues$time==1,])
#' TransitionMatrix(Psivalues.rev[Psivalues.rev$occ==1,])
#' }
NULL





#' Multi-state occupancy example data
#' 
#' An occupancy data set for modelling multi-state data (0,1,2).
#' 
#' This is a data set from Nichols et al (2007).
#' 
#' @name NicholsMSOccupancy
#' @docType data
#' @format A data frame with 40 records for 54 observations (sites) on the
#' following 2 variables.  \describe{ \item{ch}{a character vector containing
#' the presence (state 1), presence (state 2), and absence (0) for each visit
#' to the site, and a "." if the site was not visited} \item{freq}{frequency of
#' sites with that history} }
#' @references Nichols, J. D., J. E. Hines, D. I. MacKenzie, M. E. Seamans, and
#' R. J. Gutierrez.  2007. Occupancy estimation and modeling with multiple
#' states and state uncertainty.  Ecology 88:1395-1400.
#' @keywords datasets
#' @examples
#' \donttest{
#' # This example is excluded from testing to reduce package check time
#' # To create the data file use:
#' # NicholsMSOccupancy=convert.inp("NicholsMSOccupancy.inp")
#' #
#' # Create a function to fit the 12 models in Nichols et al (2007).
#' do.MSOccupancy=function()
#' {
#' #  Get the data
#'    data(NicholsMSOccupancy)
#' # Define the models; default of Psi1=~1 and Psi2=~1 is assumed
#'    # p varies by time but p1t=p2t
#'    p1.p2equal.by.time=list(formula=~time,share=TRUE)  
#'    # time-invariant p p1t=p2t=p1=p2
#'    p1.p2equal.dot=list(formula=~1,share=TRUE)    
#'    #time-invariant p1 not = p2
#'    p1.p2.different.dot=list(p1=list(formula=~1,share=FALSE),p2=list(formula=~1))  
#'    # time-varying p1t and p2t
#'    p1.p2.different.time=list(p1=list(formula=~time,share=FALSE),p2=list(formula=~time)) 
#'    # delta2 model with one rate for times 1-2 and another for times 3-5; 
#'    #delta2 defined below
#'    Delta.delta2=list(formula=~delta2) 
#'    Delta.dot=list(formula=~1)  # constant delta
#'    Delta.time=list(formula=~time) # time-varying delta
#' # Process the data for the MSOccupancy model
#'    NicholsMS.proc=process.data(NicholsMSOccupancy,model="MSOccupancy")
#' # Create the default design data
#'    NicholsMS.ddl=make.design.data(NicholsMS.proc)
#' # Add a field for the Delta design data called delta2.  It is a factor variable
#' # with 2 levels: times 1-2, and times 3-5.
#'    NicholsMS.ddl=add.design.data(NicholsMS.proc,NicholsMS.ddl,"Delta",
#'      type="time",bins=c(0,2,5),name="delta2")
#' # Create a list using the 4 p modls and 3 delta models (12 models total)
#'    cml=create.model.list("MSOccupancy")
#' # Fit each model in the list and return the results
#'    return(mark.wrapper(cml,data=NicholsMS.proc,ddl=NicholsMS.ddl))
#' }
#' # Call the function to fit the models and store it in MSOccupancy.results
#' MSOccupancy.results=do.MSOccupancy()
#' # Print the model table for the results
#' print(MSOccupancy.results)
#' # Adjust model selection by setting chat=1.74
#' MSOccupancy.results=adjust.chat(chat=1.74,MSOccupancy.results)
#' # Print the adjusted model selection results table
#' print(MSOccupancy.results)
#' #
#' # To fit an additive model whereby p1 and p2 differ by time and p2 differs from
#' # p1 a constant amount on the logit scale, use
#' #
#' # p varies by time logit(p1t)=logit(p2t)+constant
#' p1.plust.p2.by.time=list(formula=~time+p2,share=TRUE) 
#' }
NULL





#' Example of Poisson Mark-Resight model
#' 
#' Data and example illustrating Poisson Mark-Resight model with 2 groups and
#' one occasion.
#' 
#' 
#' @name Poisson_twoMR
#' @docType data
#' @format A data frame with 93 observations on the following 2 variables.
#' \describe{ \item{ch}{a character vector} \item{pg}{a factor
#' with levels \code{group1} \code{group2}} }
#' @keywords datasets
#' @examples
#'  \donttest{
#' # This example is excluded from testing to reduce package check time
#' data(Poisson_twoMR)
#' pois.proc=process.data(Poisson_twoMR,model="PoissonMR",groups="pg",
#' 		counts=list("Unmarked Seen"=matrix(c(1237,588),nrow=2,ncol=1),
#' 				    "Marked Unidentified"=matrix(c(10,5),nrow=2,ncol=1),
#' 					"Known Marks"=matrix(c(60,0),nrow=2,ncol=1)))
#' pois.ddl=make.design.data(pois.proc)
#' mod=mark(pois.proc,pois.ddl,
#' 		model.parameters=list(alpha=list(formula=~1),
#' 				              U=list(formula=~-1+group),
#' 				              sigma=list(formula=~1,fixed=0)),
#' 	 	                      initial=c(0.9741405 ,0.0000000 ,6., 5.),threads=2)
#' summary(mod)	 	                 
#' }
NULL





#' Example of Poisson Mark-Resight model
#' 
#' Data and example illustrating Poisson Mark-Resight model.
#' 
#' 
#' @name PoissonMR
#' @docType data
#' @format A data frame with 68 observations on the following 1 variables.
#' \describe{ \item{ch}{a character vector} }
#' @keywords datasets
#' @examples
#' \donttest{
#' # This example is excluded from testing to reduce package check time
#' data(PoissonMR)
#' pois.proc=process.data(PoissonMR,model="PoissonMR",
#' 		counts=list("Unmarked Seen"=c(1380, 1120, 1041, 948),
#' 				    "Marked Unidentified"=c(8,10,9,11),
#' 					"Known Marks"=c(45,67,0,0)))
#' pois.ddl=make.design.data(pois.proc)
#' mod=mark(pois.proc,pois.ddl,
#' 		model.parameters=list(Phi=list(formula=~1,link="sin"),
#' 			      GammaDoublePrime=list(formula=~1,share=TRUE,link="sin"),
#' 			      alpha=list(formula=~-1+time,link="log"),
#' 				  U=list(formula=~-1+time,link="log"),
#' 				  sigma=list(formula=~-1+time,link="log")),
#' 		          initial=c(1,1,1,1,-1.4,-.8,-.9,-.6,6,6,6,6,2,-1),threads=2)
#' summary(mod)
#' }
NULL





#' Robust Design occupancy example data
#'
#' A simulated data set on a breeding bird as an example of robust design
#' occupancy modeling.
#'
#' These are simulated data for an imaginary situation with 35 independent
#' 'sites' on which presence/absence of a breeding bird is recorded 3 times
#' annually for 3 years. Potential variables influencing site occupancy are
#' the size of the site in hectares (samplearea) and canopy cover percentage
#' (cover). The timing of the surveys within the year is thought to influence
#' the detection of occupancy, so the week the survey was conducted is included
#' in 9 variables that are named as occps where p is the primary session (year)
#' number and s is the secondary session (visit) number. Using
#' \code{data(RDOccupancy)} will retrieve the completed dataframe and using
#' \code{example(RDOccpancy)} will run the example code. However, in this
#' example we also show how to import the raw data and how they were modified
#' to construct the \code{RDOccupancy} dataframe.
#'
#' For this example, the raw data are shown below and the code below assumes
#' the file is named \code{RD_example.txt}.
#'
#' \preformatted{ ch samplearea cover occ11 occ12 occ13 occ21 occ22 occ23 occ31
#' occ32 occ33 11011.100 12 0.99 1 5 6 2 4 . 1 5 8 000110100 9 0.64 4 5 8 1 2 7
#' 2 5 9 10.100110 9 0.21 1 2 . 1 5 8 2 3 6 110000100 8 0.54 2 5 9 5 8 11 2 5 8
#' 111101100 15 0.37 1 3 5 6 8 9 5 7 12 11..11100 10 0.04 1 2 . . 2 3 5 8 14
#' 100000100 17 0.58 2 3 8 5 6 7 2 . 9 100110000 9 0.38 5 8 14 1 2 8 5 8 16
#' 1001.0100 6 0.25 4 6 8 1 . 3 1 5 6 1.110000. 17 0.34 1 . 4 3 5 9 4 5 .
#' 111100000 3 0.23 1 2 3 4 5 6 7 8 9 000000000 15 0.87 1 2 8 2 5 6 3 7 11
#' 1111.0010 8 0.18 1 2 4 1 . 3 2 3 . 10011011 . 7 0.72 2 4 5 2 6 7 1 2 .
#' 110001010 14 0.49 2 5 6 4 8 9 11 12 13 101.10100 13 0.31 1 2 3 . 2 5 1 4 6
#' 100000010 10 0.6 1 5 7 8 9 10 5 8 9 010100010 12 0.67 1 4 5 2 6 8 3 4 7
#' 110.01110 11 0.71 1 2 3 . 4 6 1 2 7 10.11.100 10 0.26 1 2 . 1 2 . 1 5 6
#' 110100.10 9 0.56 1 4 7 2 3 4 . 2 7 010000000 10 0.16 1 5 7 8 9 11 6 7 8
#' 000000.00 10 0.46 1 2 5 2 5 8 . 3 4 1.0000100 12 0.69 2 . 4 5 7 9 1 2 4
#' 100010000 11 0.42 1 2 3 4 5 6 7 8 9 000000000 12 0.42 2 5 6 5 8 9 1 3 4
#' 0.1100110 8 0.72 1 . 5 2 5 8 1 5 7 11.100100 11 0.51 1 5 . 1 2 4 4 5 6
#' 000000000 11 0.37 1 2 3 4 5 6 7 8 9 001100111 12 0.54 1 2 3 1 2 3 1 2 3
#' 10.1.1100 9 0.37 1 2 . 3 . 5 1 6 8 000000000 7 0.38 1 5 7 6 8 11 1 9 14
#' 1011.0100 8 0.35 1 5 7 2 . 5 1 3 4 100110000 9 0.86 1 2 4 2 3 6 1 2 4
#' 11.100111 8 0.57 1 5 . 2 6 7 1 3 5 }
#'
#' The data could be read into a dataframe with code as follows:
#' \preformatted{RDOccupancy<-read.table("RD_example.txt",
#' colClasses=c("character", rep("numeric",2), rep("character", 9)),
#' header=TRUE)}
#'
#' Note that if the file was not in the same working directory as your
#' workspace (.RData) then you can set the working directory to the directory
#' containing the file by using the following command before the
#' \code{read.table}.
#'
#' \code{setwd(your working directory location here)}
#'
#' In the data file "." represents a site that was not visited on an occasion.
#' Those "." values are read in fine because \code{ch} is read in as a
#' character string. However, "." has also been used in the file in place of
#' numeric values of the \code{occ} variable. Because "." is not numeric, R
#' will coerce the input value to an NA value for each "." and will treat
#' the column they are in as a factor.  Thus, the "NA" will not be a
#' valid numeric value for MARK, so we need to change it to a number. To avoid
#' the coercion, the \code{occ} values were read in as characters and the
#' following code changes all "." to "0" and then coverts the fields to numeric
#' values:
#'
#' \preformatted{ for (i in 4:12) { RDOccupancy[RDOccupancy[,i]==".",i]="0"
#' RDOccupancy[,i]=as.numeric(RDOccupancy[,i]) } }
#'
#' It is fine to use zero (or any numeric value) in place of missing values for
#' session-dependent covariates as the "0's" provide no information for
#' modeling as they are tied to un-sampled occasions. However, all values of a
#' site-specific covariate (e.g., cover) are used, so there cannot be any
#' missing values.  Note, however that use of "0's" in the time-dependent covariates
#' will influence predictions output by MARK for that parameter, as they will be 
#' biased low due to the zero's being included in estimating the mean for that 
#' parameter.
#'
#' The code below and associated comments provide a self contained example for
#' importing, setting up, and evaluating the any of the general robust design
#' type models (RDOccupEG, RDOccupPE, RDOccupPG) using RMARK. Unlike standard
#' occupancy designs, robust designs require the user to designate primary and
#' secondary occasions using the argument \code{time.intervals}. For this
#' example, we have 3 primary occasions (year) with 3 secondary sampling
#' occasions within each year, thus, we would set our \code{time.intervals} as
#' follows to represent 0 interval between secondary occasions and interval of
#' 1 (years in this case) between primary occasions:
#'
#' \preformatted{ time.intervals=c(0,0,1,0,0,1,0,0) }
#'
#' The first 0 designates the interval between the first and second sampling
#' occasion in year 1, the second 0 designates the interval between the second
#' and third sampling occasion in year 1, and the 1 indicated the change from
#' primary period 1 to primary period 2. See \code{\link{process.data}} for
#' more information on the use of \code{time.intervals}.
#'
#' @name RDOccupancy
#' @docType data
#' @format A data frame with 35 observations on the following 12 variables
#' \describe{ \item{ch}{A character vector containing the presence (1) and
#' absence (0) or (.) not visited for each of 3 visits (secondary occasions)
#' over 3 years (primary occasions)} \item{cover}{percentage canopy cover at
#' each sampled habitat} \item{occ11}{one of 9 session-dependent variables
#' occ11 to occ33 containing the week the survey was conducted; p is the
#' primary session number and s is the secondary session number}
#' \item{occ12}{one of 9 session-dependent variables occ11 to occ33 containing
#' the week the survey was conducted; p is the primary session number and s is
#' the secondary session number} \item{occ13}{one of 9 session-dependent
#' variables occ11 to occ33 containing the week the survey was conducted; p is
#' the primary session number and s is the secondary session number}
#' \item{occ21}{one of 9 session-dependent variables occ11 to occ33 containing
#' the week the survey was conducted; p is the primary session number and s is
#' the secondary session number} \item{occ22}{one of 9 session-dependent
#' variables occ11 to occ33 containing the week the survey was conducted; p is
#' the primary session number and s is the secondary session number}
#' \item{occ23}{one of 9 session-dependent variables occ11 to occ33 containing
#' the week the survey was conducted; p is the primary session number and s is
#' the secondary session number} \item{occ31}{one of 9 session-dependent
#' variables occ11 to occ33 containing the week the survey was conducted; p is
#' the primary session number and s is the secondary session number}
#' \item{occ32}{one of 9 session-dependent variables occ11 to occ33 containing
#' the week the survey was conducted; p is the primary session number and s is
#' the secondary session number} \item{occ33}{one of 9 session-dependent
#' variables occ11 to occ33 containing the week the survey was conducted; p is
#' the primary session number and s is the secondary session number}
#' \item{samplearea}{continuous variable indicating area size (ha) of the
#' sampled habitat} }
#' @author Bret Collier
#' @keywords datasets
#' @examples
#'
#' \donttest{
#' # This example is excluded from testing to reduce package check time
#' data(RDOccupancy)
#' #
#' # Example of epsilon=1-gamma
#' test_proc=process.data(RDOccupancy,model="RDOccupEG",time.intervals=c(0,0,1,0,0,1,0,0))
#' test_ddl=make.design.data(test_proc)
#' test_ddl$Epsilon$eps=-1
#' test_ddl$Gamma$eps=1
#' p.dot=list(formula=~1)
#' Epsilon.random.shared=list(formula=~-1+eps, share=TRUE)
#' model=mark(test_proc,test_ddl,model.parameters=list(Epsilon=Epsilon.random.shared, p=p.dot))
#' #
#' # A self-contained function for evaluating a set of user-defined candidate models
#' run.RDExample=function()
#' {
#' # Creating list of potential predictor variables for Psi
#'
#' Psi.area=list(formula=~samplearea)
#' Psi.cover=list(formula=~cover)
#' Psi.areabycover=list(formula=~samplearea*cover)
#' Psi.dot=list(formula=~1)
#' Psi.time=list(formula=~time)
#'
#' # Creating list of potential predictor variables for p
#' # When coding formula with session-dependent (primary or secondary)
#' # covariates, you do NOT have to include the session identifiers (
#' # the ps of occps) in the model formula. You only need to specify ~occ.
#' # The variable suffix can be primary occasion numbers or
#' # primary and secondary occasion numbers.
#'
#' p.dot=list(formula=~1)
#' p.occ=list(formula=~occ)
#' p.area=list(formula=~sample.area)
#' p.coverbyocc=list(formula=~occ*cover)
#'
#' # Creating list of potential predictor variables for Gamma
#' # and/or Epsilon (depending on which RDOccupXX Parameterization is used)
#'
#' gam.area=list(formula=~samplearea)
#' epsilon.area=list(formula=~samplearea)
#' gam.dot=list(formula=~1)
#' epsilon.dot=list(formula=~1)
#'
#' # setting time intervals for 3 primary sessions with
#' # secondary session length of 3,3,3
#'
#' time_intervals=c(0,0,1,0,0,1,0,0)
#'
#' # Initial data processing for RMARK RDOccupPG
#' # (see RMARK appendix C-3 for list of RDOccupXX model paramterizations)
#'
#' RD_process=process.data(RDOccupancy, model="RDOccupPG",
#' time.intervals=time_intervals)
#' RD_ddl=make.design.data(RD_process)
#' # Candidate model list
#' # 1. Occupancy, detection, and colonization are constant
#'
#' model.p.dot.Psi.dot.gam.dot<-mark(RD_process, RD_ddl,
#' model.parameters=list(p=p.dot, Psi=Psi.dot, Gamma=gam.dot),
#' invisible=TRUE)
#'
#' # 2. Occupancy varies by time, detection is constant,
#' # colonization is constant
#'
#' model.p.dot.Psi.time.gam.dot<-mark(RD_process, RD_ddl,
#' model.parameters=list(p=p.dot, Psi=Psi.time, Gamma=gam.dot),
#' invisible=TRUE)
#'
#' # 3. Occupancy varies by area, detection is constant,
#' # colonization varies by area
#'
#' model.p.dot.Psi.area.gam.area<-mark(RD_process,
#' RD_ddl, model.parameters=list(p=p.dot, Psi=Psi.area,
#' Gamma=gam.area), invisible=TRUE)
#'
#' # 4. Occupancy varies by cover, detection is constant,
#' # colonization varies by area
#'
#' model.p.dot.Psi.cover.gam.area<-mark(RD_process, RD_ddl,
#' model.parameters=list(p=p.dot, Psi=Psi.cover, Gamma=gam.area),
#' invisible=TRUE)
#'
#' # 5. Occupancy is constant, detection is session dependent,
#' # colonization is constant
#'
#' model.p.occ.Psi.dot.gam.dot<-mark(RD_process, RD_ddl,
#' model.parameters=list(p=p.occ, Psi=Psi.dot, Gamma=gam.dot),
#' invisible=TRUE)
#'
#' # 6. Occupancy varied by area, detection is session
#' # dependent, colonization is constant
#' model.p.occ.Psi.area.gam.dot<-mark(RD_process, RD_ddl,
#' model.parameters=list(p=p.occ, Psi=Psi.area, Gamma=gam.dot),
#' invisible=TRUE)
#' #
#' # Return model table and list of models
#' #
#' return(collect.models())
#' }
#' # This runs the 6 models above-Note that if you use
#' # invisible=FALSE in the above model calls
#' # then the mark.exe prompt screen will show as each model is run.
#'
#' robustexample<-run.RDExample() #This runs the 6 models above
#'
#' # Outputting model selection results
#' robustexample 	# This will print selection results
#' options(width=150)	# Sets page width to 100 characters
#' sink("results.table.txt") # Captures screen output to file
#'
#' # Remove comment to see output
#' #print(robustexample) # Sends output to file
#' sink() # Returns output to screen
#' #
#' # Allows you to view results in notepad;remove # to see output
#' # system("notepad results.table.txt", invisible=FALSE, wait=FALSE)
#'
#' # Examine the output for Model 1: Psi(.), p(.), Gamma(.)
#' # Opens MARK results file in text editor
#' #robustexample$model.p.dot.Psi.dot.gam.dot
#'
#' # View beta estimates for specified model in R
#' robustexample$model.p.dot.Psi.dot.gam.dot$results$beta
#'
#' # View real estimates for specified model in R
#' robustexample$model.p.dot.Psi.dot.gam.dot$results$real
#'
#' # Examine the best fitting model which has a time-dependent
#' # effect on detection
#' # (Model 5: Psi(.), p(occ), Gamma(.))
#'
#' # View beta estimates for specified model in R
#' robustexample$model.p.occ.Psi.dot.gam.dot$results$beta
#'
#' # View real estimates for specified model in R
#' robustexample$model.p.occ.Psi.dot.gam.dot$results$real
#'
#' # View estimated variance/covariance matrix in R
#' robustexample$model.p.occ.Psi.dot.gam.dot$results$beta.vcv
#'
#' 
#' # View model averages estimates for session-dependent
#' # detection probabilities
#' model.average(robustexample, "p", vcv=TRUE)
#'
#' # View model averaged estimate for Psi (Occupancy)
#' model.average(robustexample, "Psi", vcv=TRUE)
#'
#' # View model averaged estimate for Gamma (Colonization)
#' model.average(robustexample, "Gamma", vcv=TRUE)
#'
#' #
#' # Compute real estimates across the range of covariates
#' # for a specific model parameter using Model 6
#' #
#' # Identify indices we are interested in predicting
#' # see covariate.predictions for information on
#' # index relationship to real parameters
#'
#' summary.mark(robustexample$model.p.occ.Psi.area.gam.dot, se=TRUE)
#' # Define data frame of covariates to be used for analysis
#'
#' ha<-sort(RDOccupancy$samplearea)
#'
#' # Predict parameter of interest (Psi) across the
#' # range of covariate data of interest
#'
#' Psi.by.Area<-covariate.predictions(robustexample,
#' data=data.frame(samplearea=ha), indices=c(1))
#'
#' # View dataframe of real parameter estimates without var-cov
#' # matrix printing (use str(Psi.by.Area) to evaluate structure))
#'
#' Psi.by.Area[1]
#'
#' #Create a simple plot using plot() and lines()
#'
#' plot(Psi.by.Area$estimates$covdata, Psi.by.Area$estimates$estimate,
#' type="l", xlab="Patch Area", ylab="Occupancy", ylim=c(0,1))
#' lines(Psi.by.Area$estimates$covdata, Psi.by.Area$estimates$lcl, lty=2)
#' lines(Psi.by.Area$estimates$covdata, Psi.by.Area$estimates$ucl, lty=2)
#'
#' # For porting graphics directly to file, see pdf() or png(),
#'
#' } 
NULL




#' Robust design salamander occupancy data
#' 
#' A robust design occupancy data set for modelling presence/absence data for
#' salamanders.
#' 
#' This is a data set that I got from Gary White which is suppose to be
#' salamander data collected with a robust design.
#' 
#' @name RDSalamander
#' @docType data
#' @format A data frame with 40 observations (sites) on the following 2
#' variables.  \describe{ \item{ch}{a character vector containing the presence
#' (1) and absence (0) with 2 primary occasions with 48 and 31 visits to the
#' site} \item{freq}{frequency of sites (always 1)} }
#' @keywords datasets
#' @examples
#' \donttest{
#' # This example is excluded from testing to reduce package check time
#' fit.RDOccupancy=function()
#' {
#'    data(RDSalamander)
#'    occ.p.time.eg=mark(RDSalamander,model="RDOccupEG",
#'       time.intervals=c(rep(0,47),1,rep(0,30)),
#'       model.parameters=list(p=list(formula=~session)),threads=2)
#'    occ.p.time.pg=mark(RDSalamander,model="RDOccupPG",
#'       time.intervals=c(rep(0,47),1,rep(0,30)),
#'       model.parameters=list(Psi=list(formula=~time),
#'       p=list(formula=~session)),threads=2)
#'    occ.p.time.pe=mark(RDSalamander,model="RDOccupPE",
#'       time.intervals=c(rep(0,47),1,rep(0,30)),
#'       model.parameters=list(Psi=list(formula=~time),
#'       p=list(formula=~session)),threads=2)
#' return(collect.models())
#' }
#' RDOcc=fit.RDOccupancy()
#' print(RDOcc)
#' }
#' 
NULL





#' Robust design example data
#' 
#' A robust design example data set that accompanies MARK as an example
#' analysis using the various models for the robust design.
#' 
#' This is a data set that accompanies program MARK as an example for robust
#' models. The data are entered with the summary format using the variable
#' \code{freq} which represents the number of critters with that capture
#' (encounter) history.  The data set represents a robust design with 5 primary
#' occasions and within each primary occasion the number of secondary occasions
#' is 2,2,4,5,2 respectively.  This is represented with the
#' \code{time.intervals} argument of \code{\link{process.data}} which are
#' 0,1,0,1,0,0,0,1,0,0,0,0,1,0. The 0 time intervals represent the secondary
#' sessions in which the population is assumed to be closed. The non-zero
#' values are the time intervals between the primary occasions.  They are all 1
#' in this example but they can have different non-zero values.  The code
#' determines the structure of the robust design based on the time intervals.
#' The intervals must begin and end with at least one 0 and there must be at
#' least one 0 between any 2 non-zero elements. The number of occasions in a
#' secondary session is one plus the number of contiguous zeros.
#' 
#' @name robust
#' @docType data
#' @format A data frame with 668 observations on the following 2 variables.
#' \describe{ \item{ch}{a character vector containing the encounter history }
#' \item{freq}{the number of critters with that capture history} }
#' @keywords datasets
#' @examples
#' \donttest{
#' # This example is excluded from testing to reduce package check time
#' data(robust)
#' run.robust=function()
#' {
#' #
#' # data from Robust.dbf with MARK
#' # 5 primary sessions with secondary sessions of length 2,2,4,5,2
#' #
#' time.intervals=c(0, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0)
#' #
#' # Random emigration, p=c varies by time and session, S by time
#' #
#' S.time=list(formula=~time)
#' p.time.session=list(formula=~-1+session:time,share=TRUE)
#' GammaDoublePrime.random=list(formula=~time,share=TRUE)
#' model.1=mark(data = robust, model = "Robust", 
#'             time.intervals=time.intervals,
#'             model.parameters=list(S=S.time,
#'             GammaDoublePrime=GammaDoublePrime.random,p=p.time.session),threads=2)
#' #
#' # Random emigration, p varies by session, uses Mh but pi fixed to 1, 
#' # S by time.This model is in the example Robust with MARK but it is
#' # a silly example because it uses the heterogeneity model but then fixes 
#' # pi=1 which means there is no heterogeneity.Probably the data were 
#' # not generated under Mh.  See results of model.2.b
#' #
#' pi.fixed=list(formula=~1,fixed=1)
#' p.session=list(formula=~-1+session,share=TRUE)
#' model.2.a=mark(data = robust, model = "RDHet", 
#'             time.intervals=time.intervals,
#'             model.parameters=list(S=S.time,
#'             GammaDoublePrime=GammaDoublePrime.random,
#'             p=p.session,pi=pi.fixed),threads=2)
#' #
#' # Random emigration, p varies by session, uses Mh and in this 
#' # case pi varies and so does p across
#' # mixtures with an additive session effect.
#' #
#' pi.dot=list(formula=~1)
#' p.session.mixture=list(formula=~session+mixture,share=TRUE)
#' model.2.b=mark(data = robust, model = "RDHet", 
#'             time.intervals=time.intervals,
#'             model.parameters=list(S=S.time,
#'             GammaDoublePrime=GammaDoublePrime.random,
#'             p=p.session.mixture,pi=pi.dot),threads=2)
#' #
#' # Markov constant emigration rates, pi varies by session, 
#' # p=c varies by session, S constant
#' # This model is in the example Robust with MARK 
#' # but it is a silly example because it
#' # uses the heterogeneity model but then fixes pi=1 
#' # which means there is no heterogeneity.
#' # Probably the data were not generated under Mh.  
#' # See results of model.3.b
#' #
#' S.dot=list(formula=~1)
#' pi.session=list(formula=~session)
#' p.session=list(formula=~-1+session,share=TRUE)
#' GammaDoublePrime.dot=list(formula=~1)
#' GammaPrime.dot=list(formula=~1)
#' model.3.a=mark(data = robust, model = "RDHet", 
#'             time.intervals=time.intervals,
#'             model.parameters=list(S=S.dot,
#'             GammaPrime=GammaPrime.dot,
#'             GammaDoublePrime=GammaDoublePrime.dot,
#'             p=p.session,pi=pi.session),threads=2)
#' #
#' # Markov constant emigration rates, pi varies by session, 
#' # p=c varies by session+mixture, S constant. This is model.3.a 
#' # but allows pi into the model by varying p/c by mixture.
#' #
#' S.dot=list(formula=~1)
#' pi.session=list(formula=~session)
#' GammaDoublePrime.dot=list(formula=~1)
#' GammaPrime.dot=list(formula=~1)
#' model.3.b=mark(data = robust, model = "RDHet", 
#'             time.intervals=time.intervals,
#'             model.parameters=list(S=S.dot,
#'             GammaPrime=GammaPrime.dot,
#'             GammaDoublePrime=GammaDoublePrime.dot,
#'             p=p.session.mixture,pi=pi.session),threads=2)
#' #
#' # Huggins Random emigration, p=c varies by time and session, 
#' # S by time
#' # Beware that this model is not quite the same 
#' # as the others above that say random emigration because
#' # the rates have been fixed for the last 2 occasions.  
#' # That was done with PIMS in the MARK example and
#' # here it is done by binning the times so that times 3 and 4 
#' # are in the same bin, so the time model
#' # has 3 levels (1,2, and 3-4).  By doing so the parameters 
#' # become identifiable but this may not be
#' # reasonable depending on the particulars of the data.  
#' # Note that the same time binning must be done both for
#' # GammaPrime and GammaDoublePrime because the parameters are 
#' # the same in the random emigration model.  If you
#' # forget to bin one of the parameters across time it will fit 
#' # a model but it won't be what you expect as it will
#' # not share parameters.  Note the use of the argument "right".  
#' # This controls whether binning is inclusive on the right (right=TRUE) 
#' # or on the left (right=FALSE).  Using "right" nested in the list
#' # of design parameters is equivalent to using it as a calling 
#' # argument to make.design.data or add.design.data.
#' #
#' S.time=list(formula=~time)
#' p.time.session=list(formula=~-1+session:time,share=TRUE)
#' GammaDoublePrime.random=list(formula=~time,share=TRUE)
#' model.4=mark(data = robust, model = "RDHuggins", 
#'         time.intervals=time.intervals,design.parameters=
#'         list(GammaDoublePrime=list(time.bins=c(1,2,5))),
#'         right=FALSE, model.parameters=
#'         list(S=S.time,GammaDoublePrime=GammaDoublePrime.random,
#'         p=p.time.session),threads=2)
#' 
#' return(collect.models())
#' }
#' robust.results=run.robust()
#' #
#' #  You will receive a warning message that the model list 
#' #  includes models of different types which are not compatible 
#' #  for comparisons of AIC.  That is because
#' #  the runs include closed models which include N 
#' #  in the likelihood and Huggins models which don't include 
#' #  N in the likelihood.  That can be avoided by running
#' #  the two types of models in different sets.
#' #
#' robust.results
#' }
NULL





#' Salamander occupancy data
#' 
#' An occupancy data set for modelling presence/absence data for salamanders.
#' 
#' This is a data set that accompanies program PRESENCE and is explained on
#' page 99 of MacKenzie et al. (2006).
#' 
#' @name salamander
#' @docType data
#' @format A data frame with 39 observations (sites) on the following 2
#' variables.  \describe{ \item{ch}{a character vector containing the presence
#' (1) and absence (0) for each visit to the site} \item{freq}{frequency of
#' sites (always 1)} }
#' @references MacKenzie, D.I., Nichols, J. D., Royle, J.A., Pollock, K.H.,
#' Bailey, L.L., and Hines, J.E.  2006. Occupancy Estimation and Modeling:
#' Inferring Patterns and Dynamics of Species Occurence. Elsevier, Inc. 324p.
#' @keywords datasets
#' @examples
#' \donttest{
#' # This example is excluded from testing to reduce package check time
#' do.salamander=function()
#' {
#'    data(salamander)
#'    occ.p.dot=mark(salamander,model="Occupancy")
#'    occ.p.time=mark(salamander,model="Occupancy",
#'          model.parameters=list(p=list(formula=~time)))
#'    occ.p.mixture=mark(salamander,model="OccupHet",
#'          model.parameters=list(p=list(formula=~mixture)))
#'    return(collect.models())
#' }
#' salamander.results=do.salamander()
#' print(salamander.results)
#' }
NULL


#' Occupancy data for Mahoenui Giant Weta
#' 
#' An occupancy data set for modelling presence/absence data for salamanders.
#' 
#' This is a data set that accompanies program PRESENCE and is explained on
#' pages 116-122 of MacKenzie et al. (2006).
#' 
#' @name weta
#' @docType data
#' @format A data frame with 72 observations (sites) on the following 7
#' variables.  \describe{ \item{ch}{a character vector containing the presence
#' (1) and absence (0), or (.) not visited for each of 5 visits to the site}
#' \item{Browse}{0/1 dummy variable to indicate browsing} \item{Obs1}{observer
#' number for visit 1; . used when site not visited} \item{Obs2}{observer
#' number for visit 2; . used when site not visited} \item{Obs3}{observer
#' number for visit 3; . used when site not visited} \item{Obs4}{observer
#' number for visit 4; . used when site not visited} \item{Obs5}{observer
#' number for visit 5; . used when site not visited} }
#' @references MacKenzie, D.I., Nichols, J. D., Royle, J.A., Pollock, K.H.,
#' Bailey, L.L., and Hines, J.E.  2006. Occupancy Estimation and Modeling:
#' Inferring Patterns and Dynamics of Species Occurence. Elsevier, Inc. 324p.
#' @keywords datasets
#' @examples
#' 
#' #  The data can be imported with the following command using the
#' #  tab-delimited weta.txt file in the data subdirectory.
#' #   weta=import.chdata("weta.txt",field.types=c(rep("f",6)))
#' #  Below is the first few lines of the data file that was constructed
#' #  from the .xls file that accompanies PRESENCE.
#' #ch	Browse	Obs1	Obs2	Obs3	Obs4	Obs5
#' #0000.	1	1	3	2	3	.
#' #0000.	1	1	3	2	3	.
#' #0001.	1	1	3	2	3	.
#' #0000.	0	1	3	2	3	.
#' #0000.	1	1	3	2	3	.
#' #0000.	0	1	3	2	3	.
#' #
#' \donttest{
#' # This example is excluded from testing to reduce package check time
#' # retrieve weta data
#' data(weta)
#' # Create function to fit the 18 models in the book
#' fit.weta.models=function()
#' {
#' #  use make.time.factor to create time-varying dummy variables Obs1 and Obs2
#' #  observer 3 is used as the intercept
#'    weta=make.time.factor(weta,"Obs",1:5,intercept=3)
#' #  Process data and use Browse covariate to group sites; it could have also
#' #  been used an individual covariate because it is a 0/1 variable.
#'    weta.process=process.data(weta,model="Occupancy",groups="Browse")
#'    weta.ddl=make.design.data(weta.process)
#' #  time factor variable copied to Day to match names used in book
#'    weta.ddl$p$Day=weta.ddl$p$time
#' # Define p models
#'    p.dot=list(formula=~1)
#'    p.day=list(formula=~Day)
#'    p.obs=list(formula=~Obs1+Obs2)
#'    p.browse=list(formula=~Browse)
#'    p.day.obs=list(formula=~Day+Obs1+Obs2)
#'    p.day.browse=list(formula=~Day+Browse)
#'    p.obs.browse=list(formula=~Obs1+Obs2+Browse)
#'    p.day.obs.browse=list(formula=~Day+Obs1+Obs2+Browse)
#' # Define Psi models
#'    Psi.dot=list(formula=~1)
#'    Psi.browse=list(formula=~Browse)
#' # Create model list
#'    cml=create.model.list("Occupancy")
#' # Run and return marklist of models
#'    return(mark.wrapper(cml,data=weta.process,ddl=weta.ddl))
#' }
#' weta.models=fit.weta.models()
#' # Modify the model table to show -2lnl and use AIC rather than AICc
#' weta.models$model.table=model.table(weta.models,use.AIC=TRUE,use.lnl=TRUE)
#' # Show new model table which duplicates the results except they have
#' # some type of error with the model Psi(.)P(Obs+Browse) which should have
#' # 5 parameters rather than 4 and the -2lnl also doesn't agree with the results here
#' weta.models
#' #
#' # display beta vcv matrix of the Psi parameters (intercept + browse=1)
#' # matches what is shown on pg 122 of Occupancy book
#' weta.models[[7]]$result$beta.vcv[8:9,8:9]
#' # compute variance-covariance matrix of Psi0(6; unbrowsed) ,Psi1(7; browsed)
#' vcv.psi=get.real(weta.models[[7]],"Psi",vcv=TRUE)$vcv.real
#' vcv.psi
#' # Compute proportion unbrowsed and browsed
#' prop.browse=c(37,35)/72
#' prop.browse
#' # compute std error of overall estimate as shown on pg 121-122
#' sqrt(sum(prop.browse^2*diag(vcv.psi)))
#' # compute std error and correctly include covariance between Psi0 and Psi1
#' sqrt( t(prop.browse) %*% vcv.psi %*% prop.browse )
#' # show missing part of variance 2 times cross-product of prop.browse * covariance
#' 2*prod(prop.browse)*vcv.psi[1,2]
#' sqrt(sum(prop.browse^2*diag(vcv.psi))+2*prod(prop.browse)*vcv.psi[1,2])
#' 
#' }
NULL

#' White-winged dove Jolly-Seber POPAN Analysis Details
#' 
#' This dataset represents 2 years of capture-mark-recapture data collected on uniquely identifiable leg-banded (size 4) white-winged dove
#' captured in Alice, Texas, USA (Latitude 27.25, Longitude -98.07) between mid-February and mid-September during 2009 and 2010.  The package was developed such that others
#' could recreate the analysis developed by Collier, B. A., S. R. Kremer, C. D. Mason, J. Stone, K. W. Calhoun, and M. J. Peterson.  2012. Immigration and 
#' recruitment in an urban white-winged dove breeding colony.  Journal of Fish and Wildlife Management, In review., and see how the data and results 
#' were used to estimate population level recruitment (number juveniles in population over number adults in population).
#' 
#' @details White-winged doves were captured (aged:  AHY=after-hatch year, HY=hatch year) continuously in baited walk-in dove traps beginning in Februrary and ending in September in each year (2009 and 2010).  
#' During 2009 5,101 white-winged doves were captured (2,894 AHY, 2,207 HY) while in 2010 3,502 white-winged doves were captured (3,106 AHY, 486 HY). We used approximately 2 week
#' date windows to categorized our encounter histories for analysis in MARK \url{http://www.phidot.org/software/mark/} via RMark  
#' using these dates:  27 Feb; 13 March; 27 March; 10 April; 24 April; 8 May; 22 May; 5 June; 19 June;	3 July;	17 July; 31 July; 14 August-End; giving us 13 encounter occasions.
#' 
#' I wanted to force b0=0 for the first time frame, as none could be there when we started as they had not arrived yet in any real number.  
#' If you take the first column out, the numbers get ridiculously screwy for the super population size and the entry parameters, 
#' because the initial population has individuals in it, thus the JSPOPAN estimates something like 40% of the birds were already in 
#' the population pre-trapping, which is biologically impossible.
#' 
#' Initially, because of the parameter structure in a JS-POPAN model and the fact that the initial entry probability is 1 minus the sum of the resultant entry probabilities 
#' for subsequent sampling occasions, and because occasionally a couple of doves were captured during the initial time frame, we were getting entry values for the initial time period
#' representing >40% of the total doves captured. This was not plausible because these doves migrate in from Mexico/S. America and very few birds have arrived prior to trapping.  
#' To obtain reasonable results a 'fake' encounter occasion (time -1) with no captures ("0") was appended to beginning of each capture history to force b0=0.
#' As such, the data for \code{wwdo.09} and \code{wwdo.10} will have 14, not 13 encounter histories.
#' 
#' Important to note, in case you don't read Collier et al. (2012), is the fact that the 2010 dataset is kind of screwy relative to the estimation of the entry parameters for HY wwdo.  Basically,
#' what happened was we caught a bunch of AHY birds, but when we captured HY birds, we only caught a few (~400 in 2010) and of those we captured, we rarely, if ever had any recaptures.  Without getting
#' into a bunch of speculation on what happened, as we don't really know, we suspect it had something to do with the fact that 2009 was a extreme drought in Texas, as to where 2010 was extremely wet, 
#' so mast based food sources (mulberry tree's are everywhere) and such were readily available in the urban environment, as such, lower trapping success.  So, when you fit the 2010 dataset using the
#' same candidate model set at the 2009, you get pretty nonsensical answers for the entry \code{b} parameters.  As such, we did not use a group specific entry model for 2010.  If you want more 
#' detail, see Collier et al. (2012).
#' 
#' Within the model code below, you can see that we fixed the both the \code{Phi} and \code{p} parameters for those periods in the analysis when HY birds could not be
#' in the population (e.g., when all birds migrated into the breeding colony and no HY birds had been produced yet).  There is probably a little slack in the range, as 
#' it is possible that there were some HY white-winged doves in the population in the last period we fixed, but we did not catch any there so we opted to fix it as well. 
#' 
#' Note that the R function \code{wwdo.popan} is set up in RMark speak, and will run either the \code{wwdo.09} or \code{wwdo.10} datasets as long as you specify one
#' in the lines where I have listed \code{data(wwdo.09)} and \code{wwdo=wwdo.09}.  If you want to see the 2010 analysis, just change those lines to \code{wwdo.10}.
#' 
#' @name wwdo.popan
#' @aliases wwdo.09 wwdo.10
#' @docType data
#' @format The format is 2 data frames, one for 2009 and one for 2010. 2009 has 5101 unique captures, 2010 has 3502 unique captures.
#' \describe{ \item{Prefix}{Unique band prefix identifier (usually 0914 or 0984)}
#' \item{Suffix}{Unique band suffix (numeric value)}
#' \item{E0-E13}{0/1 representing whether a dove was captured (1) or not captured (0) during that (E) sampling occasion}
#' \item{Age}{Age class with AHY=after-hatch year and HY=hatch year} }
#' @references Collier, B. A., S. R. Kremer, C. D. Mason, J. Stone, K. W. Calhoun, and M. J. Peterson.  2012. Immigration and 
#' recruitment in an urban white-winged dove breeding colony.  Journal of Fish and Wildlife Management, In review. 
#' @keywords datasets
#' @examples
#' 
#' \donttest{
#' # This example is excluded from testing to reduce package check time
#' data(wwdo.09)
#' wwdo=wwdo.09
#' wwdo.popan=function(){
#' wwdo.proc=process.data(wwdo, model="POPAN", groups="Age")
#' wwdo.ddl=make.design.data(wwdo.proc)
#' 
#' #Fixing Phi Parameters for sampling periods where HY WWDO were not available in population
#' 	hy.phi1=as.numeric(row.names(wwdo.ddl$Phi[wwdo.ddl$Phi$group=="HY" & 
#'                       wwdo.ddl$Phi$time==1,]))
#' 	hy.phi2=as.numeric(row.names(wwdo.ddl$Phi[wwdo.ddl$Phi$group=="HY" &
#'                       wwdo.ddl$Phi$time==2,]))
#' 	hy.phi3=as.numeric(row.names(wwdo.ddl$Phi[wwdo.ddl$Phi$group=="HY" & 
#'                       wwdo.ddl$Phi$time==3,]))
#' 	hy.phi4=as.numeric(row.names(wwdo.ddl$Phi[wwdo.ddl$Phi$group=="HY" & 
#'                       wwdo.ddl$Phi$time==4,]))
#' 	hy.phi5=as.numeric(row.names(wwdo.ddl$Phi[wwdo.ddl$Phi$group=="HY" & 
#'                       wwdo.ddl$Phi$time==5,]))
#' 	hy.phi6=as.numeric(row.names(wwdo.ddl$Phi[wwdo.ddl$Phi$group=="HY" & 
#'                       wwdo.ddl$Phi$time==6,]))
#' 	hy.phi7=as.numeric(row.names(wwdo.ddl$Phi[wwdo.ddl$Phi$group=="HY" & 
#'                       wwdo.ddl$Phi$time==7,]))
#' 	hy.phi.fix=c(hy.phi1, hy.phi2, hy.phi3, hy.phi4, hy.phi5, hy.phi6, hy.phi7)
#' 	
#' #Fixing PENT Parameters for sampling period where HY WWDO were not available in population
#' 	hy.pent2=as.numeric(row.names(wwdo.ddl$pent[wwdo.ddl$pent$group=="HY" & 
#'                        wwdo.ddl$pent$time==2,]))
#' 	hy.pent3=as.numeric(row.names(wwdo.ddl$pent[wwdo.ddl$pent$group=="HY" & 
#'                        wwdo.ddl$pent$time==3,]))
#' 	hy.pent4=as.numeric(row.names(wwdo.ddl$pent[wwdo.ddl$pent$group=="HY" & 
#'                        wwdo.ddl$pent$time==4,]))
#' 	hy.pent5=as.numeric(row.names(wwdo.ddl$pent[wwdo.ddl$pent$group=="HY" & 
#'                        wwdo.ddl$pent$time==5,]))
#' 	hy.pent6=as.numeric(row.names(wwdo.ddl$pent[wwdo.ddl$pent$group=="HY" & 
#'                        wwdo.ddl$pent$time==6,]))
#' 	hy.pent7=as.numeric(row.names(wwdo.ddl$pent[wwdo.ddl$pent$group=="HY" & 
#'                        wwdo.ddl$pent$time==7,]))
#' 	hy.pent.fix=c(hy.pent2, hy.pent3, hy.pent4, hy.pent5, hy.pent6, hy.pent7)
#' 	
#' #####
#' #Real Parameter Definitions
#' #####
#' #Detection process
#' 	p.dot=list(formula=~1)
#' 	p.time=list(formula=~time)
#' 	p.group=list(formula=~group)
#' 	p.g.time=list(formula=~group:time)
#' 	
#' #Survival process
#' 	Phi.dot.fix=list(formula=~1, fixed=list(index=hy.phi.fix, value=c(0,0,0,0,0,0,0)))
#' 	Phi.time.fix=list(formula=~time, fixed=list(index=hy.phi.fix, value=c(0,0,0,0,0,0,0)))
#' 	Phi.age.fix=list(formula=~group, fixed=list(index=hy.phi.fix, value=c(0,0,0,0,0,0,0)))
#' 	Phi.timeage.fix=list(formula=~time:group, 
#'      fixed=list(index=hy.phi.fix, value=c(0,0,0,0,0,0,0)))
#' 	
#' #Entry Process-always time dependent, otherwise makes no sense in my situation
#' 	pent.time.fix=list(formula=~time, fixed=list(index=hy.pent.fix, value=c(0,0,0,0,0,0)))
#' 	
#' 	Model.1=mark(wwdo.proc, wwdo.ddl, 
#'    model.parameters=list(Phi=Phi.dot.fix, p=p.dot, pent=pent.time.fix, N=list(formula=~group)),
#'    invisible=FALSE)
#' 	Model.2=mark(wwdo.proc, wwdo.ddl, 
#'   model.parameters=list(Phi=Phi.time.fix, p=p.dot, pent=pent.time.fix, N=list(formula=~group)),
#'    invisible=FALSE)
#' 	Model.3=mark(wwdo.proc, wwdo.ddl, 
#'   model.parameters=list(Phi=Phi.age.fix, p=p.dot, pent=pent.time.fix, N=list(formula=~group)), 
#'     invisible=FALSE)
#' 	Model.4=mark(wwdo.proc, wwdo.ddl, 
#'   model.parameters=list(Phi=Phi.timeage.fix, p=p.dot, pent=pent.time.fix, N=list(formula=~group)), 
#'     invisible=FALSE)
#' 	Model.5=mark(wwdo.proc, wwdo.ddl,  
#'   model.parameters=list(Phi=Phi.timeage.fix, p=p.time, pent=pent.time.fix, N=list(formula=~group)), 
#'    invisible=FALSE)
#' 	Model.6=mark(wwdo.proc, wwdo.ddl, 
#'    model.parameters=list(Phi=Phi.timeage.fix,p=p.g.time, pent=pent.time.fix,N=list(formula=~group)), 
#'    invisible=FALSE)
#' 	collect.models()
#' }
#' wwdo.out=wwdo.popan()
#' wwdo.out
#' }
NULL

#' Summary of changes by version
#' 
#' A good place to look for changes.  Often I'll add changes here but don't
#' always get to it in the documentation for awhile.  They are ordered from
#' newest to oldest.
#' 
#' Version 2.0.9 (1 Dec 2011) \itemize{ \item Patch was made to  \code{\link{make.mark.model}} to fix bug in
#' PIM creation for a multi-session model and there was just 1 session. Thanks to Erin Roche for helping to 
#' identify this bug.
#'  \item Patch was made to  \code{\link{make.mark.model}} to fix bug in
#' handling of mlogits for pi and Omega parameters with more than one group. These parameters were introduced
#' with the RDMSMisClass and other new models that were recently added.
#' \item Additional changes were made to \code{\link{export.MARK}} to re-fix changes for robust and nest survival
#' model export to MARK.
#' \item A function \code{\link{mark.wrapper.parallel}} written by Eldar Rakhimberdiev provides a parallel processing
#' version of mark.wrapper.  See the example in the help for the function.  The parallel version is functionally the same
#' and can be used in place of \code{mark.wrapper} to run sequentially or in parallel. It does not include the run argument
#' however.
#' \item A bug was fixed in \code{\link{get.real}} which caused an R error for a triangular PIM that had only a single entry.
#' Thanks to Amanda Goldberg for discovering and reporting this error. 
#' }
#' Version 2.0.8 (7 Oct 2011) \itemize{ \item Both \code{\link{setup.model}}
#' and \code{\link{setup.parameters}} were re-written to use data files
#' models.txt and parameters.txt to define models and parameters which should
#' make it easier to add new models.  The latter function is now much simpler
#' and smaller. \item Model RDOccupEG now allows sharing Epsilon and Gamma
#' parameters.  Epsilon is the dominant parameter which gets the share=TRUE
#' argument. See \code{\link{RDOccupancy}} for an example. Thanks to Jake Ivan
#' for an example of what was needed. \item To avoid confusion, the arguments
#' \code{component} and \code{component.name} were removed from the parameter
#' specification because these have not been required since v1.3 when full
#' support for individual covariates were included.  Likewise the argument
#' \code{covariates} in \code{\link{mark}} and \code{\link{make.mark.model}}
#' was removed because it was only needed to support the component approach.
#' \item  \code{\link{export.MARK}} was modified so that if all individual covariates
#' are output, it excludes factor covariates. \item Many more models were added to those
#' supported in RMark. Now 92 of the 137 models in MARK are supported. See MarkModels.pdf in the RMark directory of your R library
#' to see which models are supported (in red). Most of the remaining unsupported models
#' are versions with mis-identification error and they are not shown in MarkModels.pdf.  
#' \item Previously RMark stored the input file in a temporary file Markxxx.tmp. Using a common
#' filename caused problems when more than one model were spawned to different CPUs, so
#' now it uses a random temporary file name. It also no longer uses common file markxxx.vcv.
#' Thanks to Glenn Stauffer for testing these changes.
#' \item Sessions are now labelled using value of session time rather than numerically from 1
#' to largest session number in robust designs.  Thanks to Tommy Garrison for this suggestion.
#' \item A bug was fixed in \code{\link{get.real}} which caused incorrect assignments of fixed
#' parameter values in the unusual case where a fixed parameter had a non-zero design matrix row.
#' \item \code{\link{mark.wrapper}} was modified so it returns a list of the models that were
#' constructed if run==FALSE.  Thanks to Eldar Rakhimberdiev for the suggestion and code.
#' \item Code was added to \code{\link{extract.mark.output}} to extract deviance degrees of freedom.  Thanks 
#'  again to Eldar for contributing this code.
#' \item A bug was fixed in \code{\link{make.mark.model}} that prevented use of sin link on within session parameters
#' in a robust design model. Thanks to Tommy Garrison for reporting this bug.
#' \item A bug was fixed in \code{\link{make.mark.model}} which prevented the use of time varying covariates
#' with shared parameters. Thanks to Andre Breton for reporting this bug.
#' \item simplify argument was removed from functions because I have not found a reason not to simplify and
#' I have not been testing code with simplify=FALSE.
#' }
#' Version 2.0.7 (25 August 2011) \itemize{ \item Change to
#' \code{\link{make.mark.model}} to fix bug in which mlogits were incorrectly
#' assigned in ORDMS model when both Psi and pent used mlogit links. Thanks to
#' Glenn Stauffer for identifying and tracking down this error.  \item Fixed
#' \code{\link{export.chdata}} and \code{\link{export.MARK}} so Nest survival
#' models can be exported.  Also changed default of argument ind.covariates to
#' "all" which will use all individual covariates in the data in the file sent
#' to MARK.  Thanks to Jay Rotella for his help. \item Made change to
#' \code{\link{process.data}} and \code{\link{mark}} to include a new argument
#' reverse, which if set to TRUE with model="Multistratum" will reverse the
#' timing of transition and survival. See \code{\link{mstrata}} for an example
#' of a reverse multistratum model. \item Made change to
#' \code{\link{make.design.data}} to allow for zero time intervals in
#' non-robust design model.  This was needed to allow use of the reverse time
#' structure in multi-state models. In addition, for the reverse multistate
#' model the function now adds an occasion (occ) field to the design data
#' because the time field will be constant when with a 0 time interval.  The
#' row names of the design matrix and real parameters was extended for this
#' model to include occasion (o) and occasion cohort (oc) to create unique
#' labels because cohort and time are not unique with 0 time intervals.  }
#' Version 2.0.6 (1 July 2011) \itemize{ \item Change to MR resight examples so
#' the output does not appear in notepad which was causing problem with check
#' on CRAN submission } Version 2.0.5 (29 June 2011) \itemize{ \item Order of
#' arguments for \code{\link{model.average.marklist}} were switched incorrectly
#' such that ... was the second argument.  This has been fixed to the original
#' format where ... is at the end.  This resulted in the value of any arguments
#' other than the first to be ignored unless specifically named e.g.,
#' \code{parameter="Phi"}. Thanks to Rod Towell for reporting this error. \item
#' Additional changes were made to .First.lib to 1) examine any MarkPath
#' setting, 2) look in C:/program files/mark or c:/program files (x86)/mark, 3)
#' or to search the path.  If MARK executable cannot be found in any of those
#' ways then a warning is issued that the user needs to set MarkPath to the
#' path specification for the location of the MARK executable.  Thanks to Bryan
#' Wright for help with tracking down a bug.  \item A bug in
#' \code{\link{add.design.data}} was fixed where the function gave incorrect
#' results in pim.type was anything other than "all".  Thanks to Jeff Hostetler
#' for reporting this error. \item The mark-resight models PoissonMR,
#' LogitNormalMR, and IELogitNormalMR were added.  This required some changes
#' to \code{\link{make.mark.model}},\code{\link{make.design.data}} and
#' \code{\link{compute.design.data}} and the addition of an argument
#' \code{counts} to \code{\link{process.data}} to provide mark-resight count
#' data that are not in the capture history format. The format for
#' \code{counts} is a named list of vectors (one group) or matrices (each group
#' is a row) where the names of the list elements are specific to the model.
#' Currently there is no checking to make sure these are named correctly. Some
#' of these models are very sensitive to starting values; thus the use of
#' initial values in the examples.  Thanks to Brett McClintock for his help
#' incorporating these models. } Version 2.0.4 (1 June 2011) \itemize{ \item
#' Change made to .First.lib to check for availability of mark software that
#' depends on operating system. It now provides a warning rather than stopping
#' package attachment.  This allows the user to set MarkPath to some location
#' outside of default location Program Files or Program Files (x86) without
#' changing path which requires administrator privilege on Windows. \item
#' Change made to \code{\link{run.mark.model}} to use shQuote because link to
#' mark.exe was not working in some cases.  } Version 2.0.3 (17 May 2011)
#' \itemize{ \item Now requires R 2.13 to use path.package function \item
#' Change made to .First.lib to check for availability of mark software that
#' depends on operating system. } Version 2.0.2 (16 May 2011) \itemize{ \item
#' MarkPath no longer needs to be set if mark.exe is no longer in the default
#' location (c:/Program Files/mark) but is specified in the path. The code now
#' uses the R function Sys.which to find the correct location for mark.exe, if
#' it is contained in the Path. \item Code for crm models was removed from
#' RMark and moved to a different R package called marked that is
#' under-development.  This removes the FORTRAN code and accompanying dll and
#' some functions and help files that were extraneous to RMark capabilities.
#' There is still some code in some functions for crm that could be removed at
#' some point.  None of this matters to those who use RMark for its original
#' purpose as an interface to mark.exe. \item Many superficial changes were
#' made to code so it could be posted on CRAN.  Three changes that may be
#' noticed by users involved renaming deriv.inverse.link, summary.ch and
#' merge.design.covariates to deriv_inverse.link, summary_ch, and
#' merge_design.covariates.  These names conflicted with the generic functions
#' deriv, summary and merge.  It is only the last 2 that you may have in your
#' scripts and you will have to rename them. The deprecated function
#' merge.occasion.data was removed. Sorry for any inconvenience. \item The
#' dependency on Hmisc for the examples was removed by replacing errbar with
#' plotCI. } Version 2.0.1 (21 Feb 2011) \itemize{ \item Made a change to
#' \code{\link{run.mark.model}} to handle output filenames that exceeded
#' mark9999. \item Added an argument \code{prefix} in \code{\link{mark}},
#' \code{\link{run.mark.model}} and \code{\link{cleanup}}. Like other
#' parameters for \code{\link{mark}}, it can also be used in
#' \code{\link{mark.wrapper}} as one of the ... arguments. Previously the mark
#' files have always been named "marknnn.*". By specifying \code{prefix} you
#' can now create sets of models with different prefixes. For example,
#' \code{prefix="cu"} would result in cu001.*(cu001.out, cu001.inp, etc),
#' cu002.* etc.  This provides the ability to name files to do things like
#' naming them based on the species being analyzed.  In general, there is no
#' need to work with these files directly because the filename.* is stored with
#' each \code{mark} object and the various R functions use that link to provide
#' the information from the files. If you use prefixes other than "mark",
#' you'll need to call call \code{\link{cleanup}} with each prefix to remove
#' unused files.  See \code{\link{run.mark.model}} for an example that shows
#' use of the prefix argument to split the dipper data into separate analyses
#' for each sex.  Note that use of prefix was not mandatory here to separate
#' the analyses but it provided a useful example.  \item Additional usefulness
#' has been coded for argument \code{initial} for assigning initial values to
#' beta parameters.  Previously, the options were either a vector of the same
#' length as the new model to be run or a previously run model in a \code{mark}
#' object from which equivalent betas are extracted based on their names in the
#' design matrix.  Now if the vector contains names for the elements they will
#' be matched with the new model like with the model option for initial.  If
#' any betas in the new model are not matched, they are assigned 0 as their
#' initial value, so the length of the \code{initial} vector no longer needs to
#' match the number of parameters in the new model as long as the elements are
#' named. The names can be retrieved either from the column names of the design
#' matrix or from \code{rownames(x$results$beta)} where \code{x} is the name of
#' the \code{mark} object.  \item Using the feature above, I added a new
#' argument \code{use.initial} to \code{\link{mark.wrapper}}.  If
#' \code{use.initial=TRUE}, prior to running a model it looks for the first
#' model that has already been run (if any) for each parameter formula and
#' constructs an \code{initial} vector from that previous run. For example, if
#' you provided 5 models for p and 3 for Phi in a CJS model, as soon as the
#' first model for p is run, in the subsequent 2 models with different Phi
#' models, the initial values for p are assigned based on the run with the
#' first Phi model.  At the outset this seemed like a good idea to speed up
#' execution times, but from the one set of examples I ran where several
#' parameters were at boundaries, the results were discouraging because the
#' models converged to a sub-optimal likelihood value than the runs using the
#' default initial values.  I've left this option in but set its default value
#' to FALSE. \item A possibly more useful argument and feature was added to
#' \code{\link{mark.wrapper}} in the argument \code{initial}.  Previously, you
#' could use \code{initial=model} and it would use the estimates from that
#' model to assign initial values for any model in the set defined in
#' \code{\link{mark.wrapper}}. Now I've defined \code{initial} as a specific
#' argument and it can be used as above but you can also use it to specify a
#' \code{marklist} of previously run models. When you do that, the code will
#' lookup each new model to be run in the set of models specified by
#' \code{initial} and if it finds one with the matching name then it will use
#' the estimates for any matching parameters as initial values in the same way
#' as \code{initial=model} does.  The model name is based on concatenating the
#' names of each of the parameter specification objects.  To make this useful,
#' you'll want to adapt to an approach that I've started to use of naming the
#' objects something like p.1,p.2 etc rather than naming them something like
#' p.dot, p.time as done in many of the examples.  I've found that using
#' numeric approach is much less typing and cumbersome rather than trying to
#' reflect the formula in the name. By default, the formula is shown in the
#' model selection results table, so it was a bit redundant.  Now where I see
#' this being the most benefit.  Individual covariate models tend to run rather
#' slowly. So one approach is to run the sequence of models (eg results stored
#' in initial_marklist), including the set of formulas with all of the
#' variables other than individual covariates.  Then run another set with the
#' same numbering scheme, but adding the individual covariates to the formula
#' and using \code{initial=initial_marklist} That will work if each parameter
#' specification has the same name (eg., p.1=list(formula=~time) and then
#' p.1=list(formula=~time+an_indiv_covariate)).  All of the initial values will
#' be assigned for the previous run except for any added parameters (eg.
#' an_indiv_covariate) which will start with a 0 initial value. \item I added a
#' new function \code{\link{search.output.files}} to the set of utility
#' functions. This can be useful to search all of the output files in a
#' \code{marklist} for a specific string like "numerical convergence suspect"
#' or just "WARNING".  The function returns the model numbers in the
#' \code{marklist} that contain that string in the output file. \item Further
#' changes were needed to \code{\link{popan.derived}} to handle data that are
#' summarized (i.e. frequency of capture history >1). Thanks to Carl Schwarz
#' for reporting and finding the change needed. \item A bug was fixed in
#' \code{create.dm} was fixed so it would return a matrix instead of a vector
#' with the formula ~1 } Version 2.0.0 (14 Jan 2011) \itemize{ \item The
#' packages msm, Hmisc, nlme, plotrix are now explicitly required to install
#' RMark.  These were used in examples or for specialty functions, but to avoid
#' problems these must be installed as well. \item Added example for "CRDMS"
#' model that was created by Andrew Paul. See \code{\link{crdms}}. \item Change
#' was made for gamma link in v1.9.3 was only made to Pradel model and not
#' Pradsen like it stated.  Both now correctly use the logit link as the
#' default for gamma.  Thanks to Gina Barton for bringing this to my attention.
#' \item Change was made in \code{\link{make.mark.model}} that prevented use of
#' groups with Nest survival models.  Thanks to Jeff Warren for bringing this
#' to my attention.  \item Change was made in \code{\link{make.design.data}}
#' because re-ordering of parameters caused issues with the CRDMS model because
#' the \code{subtract.stratum} were not being set for \code{Psi}.  } Version
#' 1.9.9 (2 Nov 2010) \itemize{ \item This version was built with R 2.12 and
#' will not work with earlier versions of R.  It contains both 32 and 64 bit
#' versions and R will automatically ascertain which to use. \item Parameter
#' ordering for some models (RDHet, RDFullHet, RDHHet, RDHFHet, OccupHet,
#' RDOccupHetPE, RDOccupHetPG, RDOccupHetEG, MSOccupancy, ORDMS, and CRDMS) had
#' to be changed such that the models could be imported into the MARK
#' interface. This change can influence any code you have written for those
#' models if you specified parameter indices because the ordering of the
#' parameters were changed. For example, see change in example for
#' \code{\link{RDOccupancy}} to use indices=c(1) from c(10).  Thanks to Gary
#' White for helping me work this out.  \item Modified code in
#' \code{\link{extract.mark.output}} to handle cases with more than 9999 real
#' parameters because MARK outputs **** when it exceeds 9999. } Version 1.9.8
#' (15 Sept 2010) \itemize{ \item Added \code{model="CRDMS"} with parameters S,
#' Psi, N, p, and c for closed robust design multi-state models \item Patched
#' \code{\link{popan.derived}} which produced incorrect abundance estimates
#' with unequal time intervals. Thanks to Andy Paul for finding this error and
#' testing for me. \item Patched \code{\link{export.MARK}} which failed for
#' robust design models.  Thanks to Dave Hewitt and Gary White for discovering
#' and isolating the problem. } Version 1.9.7 (14 April 2010) \itemize{ \item
#' Added \code{model="Brownie"} with parameters S and f which is the Brownie et
#' al. parameterization of the recovery model. "Recovery only" (in RMark)
#' \code{model="Recovery"} which is also encounter type "dead" in MARK uses the
#' Seber parameterization with parameters S and r which is also used in the
#' models for live and dead encounter models. \item Added
#' \code{model="MSLiveDead"} with parameters S, r, Psi and p.  It is the
#' multistate version of the Burnham model in which F=1. \item Added
#' \code{\link{compute.Sn}} to utility functions for computation of natural
#' survival from total survival when all harvest is reported. Patched
#' \code{\link{nat.surv}} which was incorrectly rejecting based on model type.
#' \item Added \code{\link{var.components.reml}} to provide an alternate
#' variance components estimation using REML or maximum likelihood. It allows a
#' random component that is not iid which is all that
#' \code{\link{var.components}} can do. \item Replaced all T/F values with
#' TRUE/FALSE to avoid conflicts with objects named TRUE or FALSE. } Version
#' 1.9.6 (1 February 2010) \itemize{ \item Writing the
#' \code{\link{popan.derived}} function has led me down all sorts of paths.  I
#' had to make one small change to this function to handle externally saved
#' models.  However, various changes listed below were brought on by using this
#' function with a relatively large POPAN model. \item Most importantly I
#' discovered an error in computations of real parameters with an mlogit link
#' in which some of the real parameters involved in the mlogit link were fixed.
#' This is NOT a problem if you simply used the real parameter values extracted
#' from MARK; however, if you were using either \code{\link{compute.real}}
#' (model.average uses this function) or \code{\link{covariate.predictions}} to
#' compute those real parameters (with an mlogit link) then they may be
#' incorrect.  This would have been apparent because their value would have
#' changed relative to the original values extracted from the MARK output.
#' Correcting this error involved changes in \code{\link{compute.real}},
#' \code{\link{convert.link.to.real}} and \code{\link{covariate.predictions}}
#' to correct the real parameter estimates and their standard errors.  I've not
#' found it in the MARK documentation but deduced that if you use the mlogit
#' link and fix real parameters it uses those fixed real parameter values in
#' the calculation. A simple example will make it clear.  Consider pent for 5
#' occasions where the first is computed by subtraction and you then have 4
#' real parameters.  Let's assume that the 3rd and fourth parameters were fixed
#' to 0.  Then the real parameters are calculated as follows:
#' pent2=exp(beta2)/(1+ exp(beta2)+exp(beta3)+exp(0)+exp(0)),
#' pent3=exp(beta3)/(1+ exp(beta2)+exp(beta3)+exp(0)+exp(0)), pent4=0, pent5=0
#' and pent1=1-pent2-pent3-0-0 (in this case). Obviously you would not want to
#' fix any real parameters to be >1, <0 or to have the sum to be <0 or >1.
#' This structure also had implications on how the standard error was
#' calculated. \item In addition the coding was made more efficient in
#' \code{\link{covariate.predictions}} for the case where
#' \code{data(index=somevector)} is used without any data entries for covariate
#' values in the design matrix.  While the task performed with that use of the
#' function could be done with \code{\link{compute.real}}, it is useful to have
#' the capability in \code{\link{covariate.predictions}} as well because it
#' will then model average over the listed set of parameters. The previous
#' approach to coding was inefficient and led to very large matrices that were
#' unnecessary and could cause failure with insufficient memory for large
#' analyses.  \item Calculation of NGross was added to
#' \code{\link{popan.derived}} and logical arguments \code{N} and \code{NGross}
#' were added to control what was computed in the call. In addition, argument
#' \code{drop} was added which is passed to \code{\link{covariate.predictions}}
#' to control whether models are dropped when variance of betas are not all
#' positive.  \item \code{\link{var.components}} was modified to use qr matrix
#' inversion. The returned value for beta is now a dataframe that includes the
#' std errors which are extracted from the vcv matrix. Also, if the design
#' matrix only uses a portion of the vcv matrix, the appropriate rows and
#' columns are now extracted. Prior to this change, the standard errors would
#' have been unreliable if the design matrix didn't use the entire set of
#' thetas and vcv matrix. } Version 1.9.5 (4 December 2009) \itemize{ \item A
#' bug in \code{\link{covariate.predictions}} was fixed that would assign fixed
#' values incorrectly if the parameter indices were specified in anything but
#' ascending order.  The error would have been obvious to anyone that may have
#' encountered it because estimated parameters would likely have been assigned
#' a fixed value.  In most cases indices would be passed in order if they were
#' selected from the design data unless indices were chosen from more than one
#' parameter type.  I discovered it using the \code{\link{popan.derived}}
#' function I added in v1.9.4 because it requests indices for multiple
#' parameters in a single function call.  } Version 1.9.4 (6 November 2009)
#' \itemize{ \item Note that this version was built with R 2.10 which no longer
#' supports compiled help files (chtml). If you were using compiled help files
#' to get help with RMark, you'll need to switch to regular html files by using
#' options(help_type="html") in R.  You can put this command in your
#' RProfile.site file so it is set up that way each time you start R. The
#' functionality is the same but it is not as pretty.  You can find the index
#' (what used to be in a window on the left) as a link at the bottom of each
#' help page. \item Changed \code{\link{export.MARK}} so it will not allow
#' selection of a project name that would over-write an existing .inp file.
#' \item Added the function \code{\link{popan.derived}} which for POPAN models
#' computes derived abundance estimates by group and occasion and sum of group
#' abundances for each occasion. For some reason RMark is unable to extract all
#' of the derived parameters from the MARK binary file for POPAN models.  This
#' function provides the derived abundance estimates and their var-cov matrix
#' and adds the abundance estimate sum across groups for each occasion which is
#' not provided by MARK. Note that by default confidence intervals are based on
#' a normal distribution to match the output of MARK, but if you want
#' log-normal intervals use \code{normal=FALSE}. \item The
#' \code{\link{model.average.list}} and \code{\link{model.average.marklist}}
#' functions were modified to use revised estimator for the unconditional
#' standard error (eq 6.12 of Burnham and Anderson (2002)) which is now the
#' default in MARK.  To use eq 4.9 (the prior formula) set the argument
#' \code{revised=FALSE}. \item Fixed bug in \code{\link{model.average.list}}
#' which in some cases failed when the list of var-cov matrices were specified.
#' \item Code in \code{\link{model.average.marklist}} was changed to set
#' standard error to 0 if the variance is negative.  The same is done in the
#' var-cov matrix for the variance and any corresponding covariances.  The
#' results from RMark will now match the model average results from MARK for
#' this case.  It is not entirely clear that this is the best approach when
#' ill-fitted models are included.  \item Changed code in \code{\link{cleanup}}
#' to handle case in which a model did not run. \item Changed use of
#' \code{grep} and \code{regexpr} in \code{\link{convert.inp}} and
#' \code{\link{extract.mark.output}} to accomodate change in R.2.10.  The code
#' should work in earlier versions of R but if not update to R2.10. \item
#' Created a function \code{\link{adjust.value}} and kept special case of
#' \code{adjust.chat}.  For any field other than \code{chat} it will adjust the
#' value in \code{model$results}.  As an example, to adjust the effective
#' sample size (ESS) use \code{model.list=adjust.value("n",value,model.list)}
#' where value is replaced with the ESS you want to use. As part of the change
#' \code{\link{model.table}} was changed to recompute AICc.  } Version 1.9.3
#' (24 September 2009) \itemize{ \item Default link function for Gamma in the
#' Pradel seniority model had been incorrectly set to log and has now been
#' changed to logit to restrict it to be a probability. \item A bug was fixed
#' in \code{\link{compute.real}} and \code{\link{covariate.predictions}} in
#' which confidence intervals were being incorrectly scaled by c (chat
#' adjustment) instead of sqrt(c).  The reported standard errors were correctly
#' using sqrt(c) and only the confidence intervals for the real predictions
#' were too large.  Simply re-running the prediction computations for a model
#' will provide the correct results. There is no reason to re-run the models.
#' Also this in no way affects model selection. \item A related bug was fixed
#' in \code{\link{compute.real}} and \code{\link{covariate.predictions}} which
#' created invalid confidence intervals for real parameters if a probability
#' link other than logit was used and a single type of link was used for all
#' real parameters (e.g., sin).  } Version 1.9.2 (10 August 2009) \itemize{
#' \item Added the function \code{\link{export.MARK}} which creates a .Rinp,
#' .inp and optionally renames one or more output files for import to MARK.
#' The July 2009 version of MARK now contains a File/RMARK Import menu item
#' which will automatically create the MARK project using the information in
#' these files.  This prevents problems that have been encountered in creating
#' MARK projects with RMark output because the data/group structures are setup
#' exactly in MARK as they were in RMark. See \code{\link{export.MARK}} for an
#' example and instructions.  \item Fixed a problem in
#' \code{\link{make.design.data}} which prevented use of \code{remove.unused}
#' with unequal time intervals and more than one group. \item At least one
#' person has encountered a problem with a very large number of parameters in
#' which RMark created the input file with PIMs written in exponential notation
#' for the larger indices.  MARK will not accept that format and it will fail.
#' The solution to this is to set the R option scipen to a positive number.
#' Start with options(scipen=1) and increase if necessary. } Version 1.9.1 (2
#' June 2009) \itemize{ \item Fixed a problem \code{\link{make.design.data}}
#' which was not using \code{begin.time} to label the session values \item Made
#' a change in \code{\link{export.chdata}} like the change in
#' \code{\link{make.mark.model}} to accomodate change with release of version
#' R2.9.0. \item Made a change in \code{\link{process.data}} so that
#' \code{strata.labels} can be specified for Multistrata designs like with
#' ORDMS so an unobserved strata can be included. \item A warning was added to
#' the help file for \code{\link{export.chdata}} and \code{\link{export.model}}
#' so it is clear that the MARK database must be created correctly and with the
#' .inp file created by \code{\link{export.chdata}} from the processed data
#' that was used to create the models that are being exported.  This is to
#' ensure that the group structure is setup such that the assumed model
#' structure for groups matches the model structure setup in the .inp file. }
#' Version 1.9.0 (30 April 2009) \itemize{ \item Fixed a bug in
#' \code{\link{summary.mark}} which occasionally produced erroneous results
#' with \code{showall=FALSE}. \item Made a change in
#' \code{\link{make.mark.model}} to accomodate change with release of version
#' R2.9.0. \item RMark now requires R version 2.8.1 or higher. } Version 1.8.9
#' (9 March 2009) \itemize{ \item Changed \code{\link{model.average.marklist}}
#' and \code{\link{covariate.predictions}} to set NaN or Inf results in v-c
#' matrix to 0 to cope with poorly determined models.  Also, for each function
#' the dropping of models is now restricted to cases in which there are
#' negative variances for the betas being used in the averaged parameter
#' estimates. Unused betas are ignored.  For example, if
#' \code{\link{model.average}} is called with \code{parameter="Phi"}, then the
#' model will only be dropped if there is a negative variance for one of the
#' betas associated with "Phi". \item In \code{\link{var.components}} the
#' tolerance value (\code{tol}) in the call to \code{uniroot} was reduced to
#' 1e-15 which should provide better estimates of the process variance when it
#' is small.  Previously a process variance less than 1e-5 would be treated as
#' 0. \item Made changes to \code{\link{cleanup}}, \code{\link{coef.mark}}, and
#' \code{\link{make.mark.model}} to accomodate externally saved model objects
#' (\code{external=TRUE}). } Version 1.8.8 (5 December 2008) \itemize{ \item An
#' error was fixed in \code{\link{make.time.factor}} which created incorrect
#' assignments when only some of the time dependent variables contained a "."
#' for occasions with no data. \item Patched \code{\link{compute.design.data}}
#' which was not creating the design data in the same order as the PIM
#' construction for the newly added \code{ORDMS} model. \item Generalized
#' section of code in \code{\link{make.mark.model}} to handle \code{mlogit}
#' structure for \code{ORDMS} model. \item Fixed a bug in
#' \code{\link{process.data}} in which the initial ages were not correctly
#' assigned in some situations with multiple grouping variables.  Note that it
#' is always a good idea to examine the design data after it is created to make
#' sure it is structured properly because it relates the data and model
#' structure via the grouping variables and the pre-defined variables (ie age,
#' time etc), While I've done a lot of testing, I have certainly not tried
#' every possible example and there is always the potential for an error to
#' occur in a circumstance that I've not encountered. } Version 1.8.7 (13
#' November 2008) \itemize{ \item An argument \code{common.zero} was added to
#' function \code{\link{make.design.data}} and
#' \code{\link{compute.design.data}}.  It can be set to TRUE to make the
#' \code{Time} variable have a common time origin of \code{begin.time} which is
#' useful for shared parameters like \code{p} and \code{c} in closed capture
#' and similar models.  \item The function \code{\link{read.mark.binary}} was
#' patched to work with the newer versions of MARK.EXE since 1 Oct 2008. \item
#' The model type \code{ORDMS} for open robust design multi-state models was
#' added. An example data set will be added at a later date after further
#' testing has been completed. \item Some patches were made to fix some aspects
#' of profile intervals and to fix adjustment by chat in
#' \code{\link{summary.mark}} when \code{showall=F}. The notation for profile
#' intervals is now included in the field \code{model$results$real$note} where
#' \code{model} is the name of a mark model.  Previously an incomplete notation
#' was kept in \code{model$results$real$fixed} but that field is now used
#' exclusively to denote fixed parameters.  It is important to realize that
#' profile intervals computed by MARK are only found in
#' \code{model$results$real$note} and are not changed by a \code{chat}
#' adjustment unless the model is re-run. None of the intervals computed by
#' \code{RMark} and displayed by \code{\link{summary.mark}} are profile
#' intervals. } Version 1.8.6 (28 October 2008) \itemize{ \item A bug in an
#' error message for \code{initial.ages} in \code{\link{process.data}} was
#' fixed. \item A new function \code{\link{var.components}} was added to
#' provide variance components capability as in the MARK interface except that
#' shrinkage estimators are not computed currently.  \item Fixed parameter
#' values are now being reported correctly by
#' \code{\link{covariate.predictions}}. Also over-dispersion (c>1) was not
#' being included in the variances for parameters except those using the mlogit
#' link.  \item Some utility functions were added including
#' \code{\link{pop.est}},\code{\link{nat.surv}}, and
#' \code{\link{extract.indices}}. \item The function
#' \code{\link{model.average}} has been changed to a generic function.
#' Currently it supports 2 classes: 1) list, and 2) marklist.  The latter was
#' the original \code{model.average} which has been renamed
#' \code{\link{model.average.marklist}} and the first argument has been renamed
#' \code{x} instead of \code{model.list} to match the standard generic function
#' approach.  The previous syntax \code{model.average(...)} will work as long
#' as the usage does not name the first agument as in the example
#' \code{model.average(model.list=dipper.results,...)}.  The list formulation
#' (\code{\link{model.average.list}}) was created to enable a generic model
#' averaging of estimates instead of just real parameter estimates from a
#' \code{mark} model.  It could be used with any set of estimates, model
#' weights and estimates of precision. \item A change is needed to
#' \code{\link{read.mark.binary}} to accomodate the change to mark.exe with the
#' version dated 1 Oct 2008. Some data types (notably Nest survival) may not
#' work with the new version of mark.exe.  Working with Gary to make the patch.
#' If you need an older version of mark.exe contact me. } Version 1.8.5 (8
#' October 2008) \itemize{ \item A bug in \code{\link{process.data}} was fixed
#' that prevented use of a dataframe contained in a list while using the
#' \code{groups} argument. \item Profile intervals on the real parameters can
#' now be obtained from MARK using the arguments \code{profile.int} and
#' optionally \code{chat} in \code{\link{mark}}. The argument
#' \code{profile.int} can be set to \code{TRUE} and a profile interval will be
#' constructed for all real parameters, or a vector of parameter indices can be
#' specified to restrict the profiling to certain parameters.  The value
#' specified by \code{chat} is passed to MARK for over-dispersion. \item
#' References to cjs, js etc have been removed from here because this code was
#' removed 5/11/11. \item Yet another fix to \code{\link{summary.ch}} which
#' gave incorrect results for the number recaptured at least once when
#' \code{marray=F} and the data contained non-unity values for \code{freq}. }
#' Version 1.8.4 (29 August 2008) \itemize{ \item A generic function
#' \code{\link{coef.mark}} was added to extract the table of betas from the
#' model with the expression \code{coef(model)} where \code{model} is a
#' \code{mark} model that has been run and contains output.  The table includes
#' standard errors and confidence intervals. \item An argument \code{brief} was
#' added to \code{\link{summary.mark}}. If \code{brief=TRUE} the real
#' parameters are not included in the summary. \item References to cjs, js etc
#' have been removed from here because this code was removed 5/11/11. \item A
#' bug in \code{\link{summary.ch}} was fixed.  It would produce erroneous
#' results when the data contained a non-constant \code{freq} field. Results
#' with the default of \code{freq=1} were fine. \item The function
#' \code{\link{adjust.chat}} and its help file were changed such that it was
#' clear that a \code{model.list} argument was needed. } Version 1.8.3 (25 July
#' 2008) \itemize{ \item For robust design models, an error trap was added to
#' \code{\link{process.data}} to make sure that the capture history length
#' matches the specification for the \code{time.intervals}.  This error was
#' already trapped for non-robust models. \item Fixed an error in
#' \code{\link{make.mark.model}} that prevented interaction model of
#' session/time-specific individual covariates in a robust design model. \item
#' Fixed an error in \code{\link{process.data}} so that the field \code{freq}
#' is optional for nest survival data sets. \item \code{\link{print.mark}} was
#' modified to add an argument \code{input} which if set to \code{input=TRUE}
#' will have the MARK input file be displayed rather than the output file.
#' Also, \code{wait=FALSE} was set in the system command which means the viewer
#' window will be opened and you can carry on with R.  Before you had to close
#' the viewer window before proceeding with R. \item An example
#' \code{\link{RDOccupancy}} provided by Bret Collier was added for the Robust
#' Occupancy model which shows the use of session and time-varying individual
#' covariates in a robust design model. } Version 1.8.2 (26 June 2008)
#' \itemize{ \item \code{\link{summary.ch}} was modified to allow missing
#' cohorts (no captures/recaptures) for an occasion and to fix a bug in which
#' \code{bygroup=FALSE} did not work when groups were defined. \item To avoid
#' running out of memory, an argument \code{external} has been added to
#' \code{\link{collect.models}}, \code{\link{mark}}, \code{\link{rerun.mark}},
#' and \code{\link{run.mark.model}}.  As with all arguments of
#' \code{\link{mark}}, \code{external} can also be set in
#' \code{\link{mark.wrapper}}.  Likewise, \code{external} can also be set in
#' \code{\link{run.models}} and it is passed to \code{\link{run.mark.model}}.
#' The default is \code{external=FALSE} but if it is set to \code{TRUE} then
#' the mark model object is saved in an external file with an extension
#' \code{.rda} and the same base filename as its matching MARK output files.
#' The mark object in the workspace is a character string which is the name of
#' the file with the saved image (e.g., "mark001.rda").  If
#' \code{external=TRUE} with \code{\link{mark.wrapper}} then the resulting
#' marklist contains a list entry for each mark model which is only the
#' filename and then the last entry is the \code{model.table}.  All of the
#' functions recognize the dual nature of the mark object (i.e., filename or
#' mark object) in the workspace.  So even if the mark object only contains the
#' filename, functions like \code{\link{print.mark}} or
#' \code{\link{summary.mark}} will work.  However, if you have used
#' \code{external=TRUE} and you want to look at part of a mark object without
#' using one of the functions, then use the function \code{load.model}.
#' Whereas, before you may have typed \code{mymark$results}, if you use
#' \code{external=TRUE}, you would replace the above with
#' \code{load.model(mymark)$results}. \item Functions \code{\link{store}} and
#' \code{\link{restore}} were created to \code{store} externally and
#' \code{restore} models from external storage into the R workspace.  They work
#' on a \code{marklist} and are only needed to \code{store} externally existing
#' marklist models or ones originally created with \code{external=FALSE} or to
#' \code{restore} if you change your mind and decide to keep them in the R
#' workspace.  \item Error in setup for robust design occupancy models with
#' more than 2 primary sessions was fixed.  The error resulted in mark.exe
#' crashing.  \item The concept of time-varying individual covariates has been
#' expanded to include robust design models which have both primary (session)
#' and secondary (time) occasion-specific data.  For a robust design, a
#' time-varying individual covariate can be either session-dependent or
#' session-time dependent.  As an example, if there are 3 primary sessions and
#' each has 4 secondary occasions, then the individual covariates can be named
#' x1,x2,x3 to be primary session-dependent or named
#' x11,x12,x13,x14,x21,x22,x23,x24,x31,x32,x33,x34. The value of x can be any
#' name for the covariate.  In the formula only the base name is used (e.g.,
#' \code{~x}) and RMark fills in the individual covariate names that it finds
#' that match either the session or session-time individual covariates. }
#' Version 1.8.1 (19 May 2008) \itemize{ \item Added function
#' \code{\link{summary.ch}} to provide summaries of the capture history data
#' (resighting matrices and m-arrays). It will not work with all types of
#' models at present.  It will work with CJS and Jolly-type models. \item Added
#' argument \code{model.name} to \code{\link{model.table}} to be able to use
#' alternate names in the model table. It can use either the model name with
#' each mark object which uses a formula notation (the current approach) or it
#' can use the name of the R object containing the mark model
#' (\code{model.name=FALSE}). See \code{\link{model.table}} for an example.
#' Also, the help file for \code{model.table} was updated to reflect the code
#' changes implemented in version 1.7.3. \item Code in
#' \code{\link{mark.wrapper}} was modified to output the number of columns and
#' column names of the design matrix for each model if \code{run=FALSE}.  This
#' allows a check of each of the columns included in the model.  By reviewing
#' these you can assess whether the model was constucted as you intended.  If
#' there is any question you can either use \code{\link{model.matrix}} or
#' \code{\link{make.mark.model}} to examine the design matrix more thoroughly.
#' \item A bug in the new function \code{\link{merge.design.covariates}} was
#' fixed in which \code{merge} was sorting the design data which does obvious
#' bad things.  Adding \code{sort=FALSE} does not appear to mean that the data
#' frame is left in its original order. To prevent this, the dataframe is
#' forced to remain in its original order by adding a sequence field for
#' re-sorting after the merge. } Version 1.8.0 (8 May 2008) \itemize{ \item
#' Fixed a bug in \code{\link{model.average}} that caused it to fail and issue
#' an error when any of the models included a time dependent covariate in the
#' parameter being averaged. \item Added \code{\link{merge.design.covariates}}
#' which is meant to replace \code{merge.occasion.data}. This new function
#' allows covariates to be assigned by \code{time}, \code{time} and
#' \code{group}, or just \code{group}. It also uses a simplified list of
#' arguments and works with individual design dataframes rather than the entire
#' ddl. It uses the R function \code{\link{merge}} which can be used on its own
#' to merge design covariates into the design data.  You can use
#' \code{\link{merge}} directly as this function only checks for some common
#' mistakes before it calls \code{merge} and it handles reassignemnt of row
#' names in the case were design data have been deleted.  An example, where you
#' might want to use \code{merge} instead of this function would be situations
#' where the design data are not just group, time or group-time specific.  For
#' example, if groups were specified by two different factor variables say
#' initial age and region and the design covariates were only region-specific.
#' It would be more efficient to use \code{merge} directly rather than this
#' function which would require an entry for each group which would be each
#' pairing of inital age and region.  If you use \code{\link{merge}} and you
#' deleted design data prior to merging, save the row.names, merge and then
#' reassign the row.names. \item An argument \code{run} was added to
#' \code{\link{mark.wrapper}}.  If set to FALSE, then it will run through each
#' set of models in \code{model.list} and try to build each model but does not
#' attempt to run it.  This is useful to check for and fix any errors in the
#' formula before setting off a large run. If you use \code{run=FALSE} do not
#' include arguments that are meant to be passed to
#' \code{\link{run.mark.model}} like \code{adjust}. } Version 1.7.9 (7 April
#' 2008) \itemize{ \item \code{\link{make.design.data}} was fixed so that
#' \code{remove.unused=T} will work properly when different \code{begin.time}
#' values are specified for each group. } Version 1.7.8 (12 March 2008)
#' \itemize{ \item Changed the default link for N to log in the
#' \code{\link{setup.parameters}} for the \code{HetClosed} and \code{FullHet}.
#' It was incorrectly set to logit which created incorrect estimates to be
#' computed in \code{\link{model.average}} because MARK forces the log link for
#' N regardless of what is set in the input file. } Version 1.7.7 (6 March
#' 2008) \itemize{ \item Supressed warning message that occurred with code to
#' check the validity of the sin link in \code{\link{make.mark.model}}. \item
#' Fixed a couple of bugs in \code{\link{covariate.predictions}} that prevented
#' it from working for some cases after including code for the sin link. \item
#' Added function \code{\link{release.gof}} to construct the RELEASE goodness
#' of fit test and extract the TEST2 and TEST3 final chi-square, df and
#' P-values. } Version 1.7.6 (26 Feb 2008) \itemize{ \item
#' \code{\link{make.mark.model}} was modified to change the capitalization of
#' the link functions and to remove all spaces after "=" in the input file for
#' mark.exe. These differences were preventing the MARK interface from fully
#' importing the model. Although the model would be imported and could be run
#' inside the MARK interface, median c-hat would not run and would give an
#' error stating "Invalid Link" for any model imported from RMark.  Now
#' transfering a model from RMark to the MARK interface is fully functional (I
#' hope).  If you want to import an output file that was created with a prior
#' version of RMark without re-running it, use a text editor on the output file
#' and remove any spaces before and after an = sign.  Then change the
#' capitalization of the links to "Logit", "MLogit", "Log", "LogLog",
#' "CLogLog", "Identity". \item The sin link is now supported if the formula
#' for the parameter generates an identity matrix for the parameter. For
#' example, if you use ~-1+time instead of ~time then the resulting design
#' matrix will be an identity for time.  Likewise, for interactions use
#' ~-1+group:time instead of ~group*time.  If you select the sin link and the
#' resulting design matrix is not an identity for the parameter, an error will
#' be given and the run will stop. \item To match the output from MARK, the
#' confidence intervals for real parameters using any 0-1 link including
#' loglog,cloglog,logit and sin are now computed using the logit
#' transformation.  For previous versions this will only affect any results
#' that were using loglog and cloglog. Previously, it was using the chosen link
#' to compute the se and the interval endpoints. The latter is still used for
#' the log and identity links which are not bounded in 0-1. \item The model
#' "Jolly" was added to the supported list of models. Parameters include
#' Phi,p,Lambda,N.  It is not a particularly numerically stable model and often
#' will not converge. Use of options="SIMANNEAL" in call to \code{\link{mark}}
#' is recommended for better convergence.  It will take much longer to converge
#' but is mroe reliable. } Version 1.7.5 (24 Jan 2008) \itemize{ \item
#' \code{\link{model.average}} was modified to ignore any models that did not
#' run and either had no attached output file or no results. \item
#' \code{\link{read.mark.binary}} and \code{\link{extract.mark.output}} were
#' modified to extract and store the real.vcv matrix (var-cov matrix of the
#' simplified real parameters) in the mark object if realvcv=TRUE. The default
#' is realvcv=FALSE. This argument has been added to functions
#' \code{\link{mark}}, \code{\link{run.mark.model}} and
#' \code{\link{rerun.mark}}. \item An argument delete has also been added to
#' \code{\link{mark}} and \code{\link{run.mark.model}}.  The default value is
#' FALSE but if set to TRUE it deletes all output files created by MARK after
#' extracting the results.  This is most useful for simulations that could
#' easily create thousands of output files and after extracting the results the
#' model objects are no longer needed. This is just a convenience to replace
#' the need to call \code{\link{cleanup}}. } Version 1.7.4 (10 Jan 2008)
#' \itemize{ \item A bug in \code{\link{make.mark.model}} was fixed.  It was
#' preventing creation of individual (site) covariate models for parameters
#' with only a single parameter (single index) in certain circumstances like
#' Psi1 in the MSOccupancy model. \item The fix to \code{merge.occasion.data}
#' in version 1.7.1 did not work when design data had been deleted.  That has
#' been remedied. \item Various functions with some operating specific calls
#' have been modified so they will work on either Windows or Linux.  Thus, the
#' there is a single file for all source/help for both operating systems in
#' RMarkSource.zip. It can be downloaded to either Windows or Linux to build
#' the package. You need to build the package for Linux but not for Windows.
#' For Windows, you only need RMark.zip which contains the pre-built package
#' which only needs to be installed. Currently, with Linux the variable
#' MarkPath is ignored and mark.exe is assumed to be in the path.  Also, for
#' Linux the default for MarkViewer is "pico" (an editor on some Linux
#' machines). This can be modified in \code{\link{print.mark}} or by setting
#' MarkViewer to a different value.  The one Linux specific function is
#' \code{read.mark.binary.linux}.  The function
#' \code{\link{extract.mark.output}} calls either \code{read.mark.binary.linux}
#' or \code{\link{read.mark.binary}} depending on the operating system. A Linux
#' version of mark.exe (32 or 64 bit) can be obtained from Evan. } Version
#' 1.7.3 (4 Jan 2008) \itemize{ \item In working with the occupancy models, it
#' became apparent that it would be useful to have a new function called
#' \code{\link{make.time.factor}} which creates time-varying dummy variables
#' from a time-varying factor variable.  An example is given using observer
#' with the occupancy dataset \code{\link{weta}} from the MacKenzie et al
#' Occupancy modelling book. \item To match the results in the book, I added
#' arguments \code{use.AIC} and \code{use.lnl} to function
#' \code{\link{model.table}} to construct a results table with AIC rather than
#' AICc and -2LnL values. The latter is more useful with a mix of models some
#' using individual covariates and others not. \item A modification was made to
#' \code{\link{make.mark.model}} with the \code{MSOccupancy} model to fix the
#' name of the added data for parameter \code{p1} when \code{share=TRUE} to be
#' \code{p2}. For an example which uses p2 to construct an additive model, see
#' \code{\link{NicholsMSOccupancy}}. } Version 1.7.2 (20 Dec 2007) \itemize{
#' \item In changing code for the occupancy models, a brace was misplaced which
#' prevented the nest survival models from working. This has been fixed.  Also,
#' the example code for \code{\link{mallard}} and \code{\link{killdeer}} was
#' modified to exclude the calls to process the input file.  This enables use
#' of the function \code{example()} to run the example code (e.g.
#' \code{example(mallard)}). From now on as I add examples they are being
#' included in my test set to avoid this type of problem in the future. }
#' Version 1.7.1 (14 Dec 2007) \itemize{ \item If you update with this version
#' of RMark make sure to update MARK also, so you get the fixes for some of the
#' occupancy models. \item A minor bug was fixed in function
#' \code{merge.occasion.data} that created duplicate row names and prevented
#' the design data from being used in a model. \item Thirteen different
#' occupancy models were added. Models in the following list use the
#' designation from MARK: \code{Occupancy, OccupHet, RDOccupEG, RDOccupPE,
#' RDOccupPG, RDOccupHetEG, RDOccupHetPE, RDOccupHetPG,OccupRNPoisson},
#' \code{OccupRNNegBin,OccupRPoisson,OccupRNegBin,MSOccupancy}. \code{Het}
#' means it uses the Pledger mixture and those with \code{RD} are the robust
#' design models. The 2 letter designations for the RD models are shorthand for
#' the parameters that are estimated.  For \code{EG}, Psi, Epsilon, and Gamma
#' are estimated, for \code{PE} gamma is dropped and for \code{PG}, Epsilon is
#' dropped.  For the latter 2 models, Psi can be estimated for each primary
#' occasion. The last 5 models include the Royle/Nichols count (RPoisson) and
#' presence (RNPoisson) models and the multi-state occupancy model.  See
#' \code{\link{salamander}} for an example of \code{Occupancy, OccupHet},
#' \code{\link{Donovan.7}} for an example of
#' \code{OccupRNPoisson,OccupRNNegBin}, \code{\link{Donovan.8}} for an example
#' of \code{OccupRPoisson,OccupRNegBin}, see \code{\link{RDSalamander}} for an
#' example of the robust design models and \code{\link{NicholsMSOccupancy}} for
#' an example of \code{MSOccupancy}. \code{\link{salamander}} data. \item The
#' functions \code{\link{create.model.list}} and \code{\link{mark.wrapper}}
#' were modified so that a list of parameters can be used to loop.  This is
#' useful in the situation with shared parameters such as \code{p1} and
#' \code{p2} in the \code{MSOccupancy} model, \code{closed} models etc. See
#' \code{p1.p2.different.dot} in \code{\link{NicholsMSOccupancy}} for an
#' example.  It can also be useful if the model definitions are linked
#' conceptually (e.g., when one parameter is time dependent, the other should
#' also be time dependent). \item The "." value in an encounter history is now
#' acceptable to RMark and gets passed to MARK for interpretation as a missing
#' value. \item \code{\link{print.marklist}} was fixed to show the model table
#' properly after a c-hat adjustment was made. The change in the code in
#' version 1.6.5 to add parameter specific values to the model table had the
#' side-effect of dropping the model name if c-hat was adjusted. } Version
#' 1.7.0 (7 Nov 2007) \itemize{ \item A function
#' \code{\link{deltamethod.special}} for computation of delta method variances
#' of some special functions was added.  It uses the function
#' \code{deltamethod} from the package \code{msm}.  You need to install the
#' package \code{msm} from CRAN to use it. \item A more complete example
#' (\code{\link{mallard}}) created by Jay Rotella was added for the nest
#' survival model. His script provides a nice tutorial for RMark and the
#' utility of R to provide a wide-open capability to calculate/plot etc with
#' the results. It also demonstrates the advantages of scripting in R to
#' document your analysis and enable it to be repeated.  Before you use his
#' tutorial you need to install the package plotrix from CRAN. At a later date,
#' Jay has said he will add some additional examples to demonstrate use of the
#' \code{deltamethod} function to create variances for functions of the results
#' from MARK. \item Various changes were made to help files.  A more complete
#' description of \code{\link{cleanup}} was given to tie into
#' \code{\link{mallard}} example. } Version 1.6.9 (10 Oct 2007) \itemize{ \item
#' Nest survival model was added to list of MARK models supported by RMark. See
#' \code{\link{killdeer}} for an example.  Note that the data structure for
#' nest models is completely different from the standard capture history so the
#' functions \code{import.chdata}, \code{export.chdata} and \code{convert.inp}
#' do not work with nest data structure. \item Slight change was made to
#' \code{\link{run.mark.model}} and \code{\link{print.mark}}to accomodate
#' change in R 2.6.0. } Version 1.6.8 (2 Oct 2007) \itemize{ \item Changes were
#' made to \code{merge.occasion.data} to enable group and time-specific design
#' covariates to be added to the design data. \item Change was made to
#' \code{\link{setup.parameters}} to use a log-link for N in the closed-capture
#' models.  MARK forces that link for N but the change was needed for
#' \code{\link{model.average}} which does the inverse-link computation.  Note
#' that the reported N in \code{\link{model.average}} is actually f0 (number
#' not seen).  To get the correct values for N simply add M_t+1 (unique number
#' captured) to f0.  That is the way MARK computes N.  The std error and
#' confidence interval is on f0 such that the lower ci on N will never be less
#' than M_t+1. \item An error was fixed in the output of
#' \code{\link{model.average}}.  When you selected a specific parameter, it was
#' giving a UCL which was a copy from one of the models and not the UCL from
#' the model averaging.  If you didn't specify \code{vcv=T} it only showed the
#' errant UCL and if you did specify \code{vcv=T} then it showed the correct
#' LCL and UCL but then added the errant UCL in a column. This occurred because
#' it was adding covariate data for the specific parameter and was shifted a
#' column because of a change in 1.6.1.  } Version 1.6.7 (7 Aug 2007) \itemize{
#' \item Changes were made in \code{\link{print.mark}},
#' \code{\link{print.summary.mark}} and \code{\link{compute.design.data}} to
#' acommodate changes in V2.5.1 of R. When upgrading versions of R problems may
#' occur if RMark was built with an earlier version of R.  The version of R
#' that was used to build RMark is listed on the screen each time it is loaded
#' with library(RMark) \code{This is RMark 1.6.7 Built: R 2.5.1;
#' i386-pc-mingw32; 2007-08-07 09:00:33; windows} \item The help file for
#' \code{\link{import.chdata}} was expanded to clarify the differences between
#' it and \code{\link{convert.inp}} and the use of the \code{freq} field. }
#' Version 1.6.6 (14 May 2007) \itemize{ \item Function
#' \code{\link{make.mark.model}} was fixed so that the real label indices were
#' properly written when \code{simplify=FALSE} is used. \item Function
#' \code{\link{make.mark.model}} was also changed to remove the parameter
#' simplification for \code{mlogit} parameters that was added in v1.4.5.  I
#' mistakenly assumed that the \code{mlogit} parameters were setup such that
#' the normalization to sum to 1 was done will all the real parameters in the
#' set (i.e., all PSI for a single stratum).  In fact, the \code{mlogit} values
#' are only specified for the unique real parameters so if there is any
#' simplification and the sum of the probabilities is close to 1 (excluding
#' subtraction value) the values will not be properly constrained.  For
#' example, with the \code{\link{mstrata}} data if the problem was constrained
#' such that PSI from AtoB was equal to AtoC, it is still necessary to have
#' these as separate real parameters and constrain them with the design matrix.
#' As it turns out, with the \code{\link{mstrata}} example it does not matter
#' because the problem is such that the sum of Psi for AtoB and AtoC is not
#' close to 1 (same for other strata) and any link will work. This change will
#' only be noticeable in situations in which the constraint matters (i.e., the
#' probability for the subtraction parameter is near 0).  The change back to
#' non-simplification for \code{mlogit} parameters may increase execution times
#' because the design matrix size has been increased.  Previous users of the
#' Multistrata design will see very little difference in there results if they
#' only used models containing stratum:tostratum because that will create an
#' all-different PIM within each \code{mlogit} set. When I ran the
#' \code{\link{mstrata}} examples with this version and compared them to v1.6.5
#' the results were different but they were differences in the 5th or smaller
#' decimal point due to differences in numerical optimization. } Version 1.6.5
#' ( 3 May 2007) \itemize{ \item Function \code{\link{model.table}} was
#' modified to include parameter formula fields in the \code{model.table}
#' dataframe of a \code{marklist}. Previously only the \code{model.name} was
#' included which is a concatenation of the individual parameter formulas. The
#' additional fields allows extracting the model table results based on one of
#' the parameter formulas or to create a matrix of model AICc or other values
#' with rows as one parameter and colums as the other.  See
#' \code{\link{model.table}} for an example. \item Function
#' \code{\link{process.data}} was modified such that factor variables used for
#' grouping retain the ordering of the factor levels in the data file.
#' Previously they would revert back to default ordering and the re-leveling
#' would also have to be repeated on the design data also. \item An argument
#' \code{brief} was added to \code{\link{mark}} to control amount of summary
#' output. \item Fixed a bug in \code{\link{get.real}} that prevented
#' computation for models without the stored covariate values. \item Added code
#' to \code{\link{make.mark.model}} that prevents constructing models with
#' empty rows in the design matrix unless the parameter is fixed.  For example,
#' if you were to try ~-1+Time for the dipper data, it will fail now because
#' there is no value for the intercept (Time=0). \item Function
#' \code{\link{mark.wrapper}} outputs the model name to the screen before
#' running the model which helps associate any error messages to the model if
#' \code{output=F}. } Version 1.6.4 (7 March 2007) \itemize{ \item A new
#' function \code{\link{export.model}} was created to copy the output files
#' into the naming convention needed to append them into a MARK .dbf/.fpt
#' database so they can be used with the MARK interface features. This is
#' useful to be able to use some of the features not contained in RMark such as
#' median c-hat and variance component estimation. To create a MARK database,
#' first use \code{\link{export.chdata}} to create a .inp file to pass the data
#' into MARK.  Start MARK and use File/New to create a new database.  Select
#' the appropriate Data Type (model in RMark) and fill in the appropriate
#' values for encounter occasions etc.  For the Encounter Histories File Name,
#' select the file you created with \code{\link{export.chdata}}. Once you have
#' created the database in the Program MARK interface, click on the Browse menu
#' item and then Output/Append and select the output file(s) (i.e those with a
#' Y.tmp) that you exported with \code{\link{export.model}}.  Note that this
#' will not work with output files run with versions of RMark prior to this one
#' because the MARK interface will give a parse error for the design matrix.
#' To get around that you can edit the output file and remove the spaces in the
#' line with the design matrix header.  For example, it should look as follows
#' \code{design matrix constraints=7 covariates=7} without spaces around the =
#' sign. \item The minor change described above was made in the input file with
#' spacing on the design matrix line to enable proper appending of the output
#' into a MARK .dbf/.fpt database. \item The function \code{\link{cleanup}} was
#' modified to delete all mark*.tmp files.  Do not use \code{cleanup} until you
#' have appended any exported models. \item An argument \code{use.comments} was
#' added to \code{\link{import.chdata}} to enable comment fields to be used as
#' row names in the data frame.  A comment is indicated as in MARK with /*
#' comment */.  They can be anywhere in the record but they must be unique and
#' they can not have a column header (field name). \item Function
#' \code{\link{create.model.list}} was modified such that it only includes
#' lists with a \code{formula} element.  This prevents collecting other objects
#' that are named similarly but are not model defintions. } Version 1.6.3 (5
#' March 2007) \itemize{ \item A minor change in \code{\link{make.mark.model}}
#' and \code{\link{find.covariates}} was made to accommodate use of the same
#' covariate in different formulas (e.g. Phi and p). Previous code worked
#' except any call to \code{\link{get.real}} would fail. Previously a duplicate
#' of the covariate was entered in the data file to MARK.  Now only a single
#' copy is passed. \item An argument \code{default} has been added to the model
#' definition (\code{parameters} in \code{\link{make.mark.model}} and
#' \code{model.parameters} in \code{\link{mark}}). The argument sets the
#' default value for parameters represented by design data that have been
#' deleted. \item Checks were added in \code{\link{make.mark.model}} to fail if
#' any of the individual covariates used are either factor variables or contain
#' NAs. Both could fail in the MARK.EXE run but the error message would be less
#' obvious.  Factor variables can work as an individual covariate, if the
#' levels are numeric.  But it was easier to exclude all factor variables from
#' being individual covariates.  They can easily be converted to a continuous
#' version (e.g. Blackduck$BirdAge=as.numeric(Blackduck$BirdAge)-1). The code
#' for the \code{\link{Blackduck}} was changed to make BirdAge a continuous
#' rather than factor variable.  Factor variables can still be used to define
#' groups and then used in the formula.  They just can't be used as individual
#' covariates. This change was made because a factor variable was in the data
#' but not defined in \code{groups} and when it was used in the formula it
#' would create a float error in MARK.EXE and that would be confusing and hard
#' to track down. } Version 1.6.2 (28 Feb 2007) \itemize{ \item The fix in
#' 1.6.1 to avoid the incorrect design matrix was not sufficiently general and
#' created a parse error in R if you attempted to use any design data
#' covariates that were created with a cut function to create factor variables
#' by binning a variable.  This has been corrected in this version. \item The
#' code in \code{\link{read.mark.binary}} has been changed to skip over the v-c
#' matrix for the derived parameters if it is not found in the file.  This was
#' causing an error with the PRADREC model type. } Version 1.6.1 (17 Jan 2007)
#' \itemize{ \item An important bug was fixed in \code{\link{make.mark.model}}
#' in which an incorrect design matrix would be created if you used two
#' individual covariates in the same formula whereby one of the covariate names
#' was contained within the other.  For example, if you used ~mass+mass2 where
#' mass2=mass^2, it would actually create a design matrix with columns mass
#' product(mass,mass2) which would be the model mass+mass^3.  This happened due
#' to the way the code identified columns where it needed to replace dummy
#' values with individual covariate names. Since mass was contained in mass2 it
#' added mass to the column as a product. The code now does exact matching so
#' the error can no longer occur. \item An argument \code{indices} was added to
#' the function \code{\link{model.average}} which enables restricting the model
#' averaging to a specific set of parameters as identified by the all-different
#' parameter indices.  This is most useful in large models with many different
#' indices such that memory limitations are encountered in constructing the
#' variance-covariance matrix of the real parameters.  For example, with a CJS
#' analysis of data with 18 groups and 26 years of data, the number of
#' parameter indices exceeds 22,000.  Even by restricting the parameters to
#' either Phi or p with the \code{parameter} argument there are still 11,000
#' which would attempt to create a matrix containing 11,000 x 11,000 elements
#' which can exceed the memory limit. In most cases, there are far fewer unique
#' parameters and this argument allows you to select which parameters to
#' average. \item Time-varying covariates are no longer needed for all times if
#' the formula is correctly written to exclude them in the resulting design
#' matrix.  \code{\link{make.mark.model}} still reports missing time-varying
#' covariates but will continue to try and fit the model but if the missing
#' variables are used in the design matrix the model will fail.  As an example
#' consider a time varying covariate x for recapture times 1990 to 1995.  The
#' code expects to find variables x1990, x1991, x1992, x1993, x1994, x1995.
#' However, lets say that the values are only known for 1993-1995.  If you
#' define a variable I'll call recap in the design data which has a value 1 for
#' 1993-1995 and a value 0 for 1990-1992 then if you use the formula ~recap:x
#' the resulting design matrix will only use the known variables for 1993-1995
#' but you will still be warned that the other values (x1990 - x1992) are
#' missing. \item A bug was fixed in \code{\link{extract.mark.output}} which
#' prevented it from obtaining more than the last mean covariate value from the
#' MARK output. \item \code{\link{fill.covariates}} was modified such that only
#' a partial list of covariate values need to be specified with \code{data} and
#' the remainder are filled in with default values depending on argument
#' \code{usemean}. \item The output from \code{\link{summary.mark}} was
#' modified for real parameters when \code{se=T} to include
#' \code{all.diff.index} to provide the indices of each real parameter in the
#' all-different PIM structure.  They are useful to restrict
#' \code{\link{covariate.predictions}} and \code{\link{model.average}} to a
#' specific set of real parameters. \item A new function
#' \code{\link{covariate.predictions}} was created to compute real parameter
#' values for multiple covariate values and their variance-covariance matrix.
#' It will also model average those values if a marklist is passed to the
#' function. Two examples from chapter 12 of Cooch and White are provided to
#' give examples of models with individual covariates and the use of this
#' function. \item The default value of \code{vcv} in
#' \code{\link{model.average}} has been changed to \code{FALSE}. } Version
#' 1.6.0 (27 Nov 2006) \itemize{ \item A bug was fixed in \code{\link{PIMS}}
#' which prevented it from working with Multistrata models. \item Bugs were
#' fixed in \code{\link{make.design.data}} which prevented use of argument
#' \code{remove.unused=T} with Multistrata models and also for any type of
#' model when there were no grouping variables. \item Bugs were fixed in
#' \code{\link{process.data}} which gave incorrect ordering of intial ages if
#' the factor variable for the age group was numeric and more than two digits.
#' Also, the number of groups in the data was not correct if the number of loss
#' on capture records exceeded the number without loss on capture within a
#' group. \item Bugs were fixed in \code{\link{setup.parameters}} and
#' \code{\link{setup.model}} that prevented use of the Barker model and that
#' reported an erroneous list of model names when an incorrect type of model
#' was selected. } Version 1.5.9 (26 June 2006) \itemize{ \item A bug was fixed
#' in \code{\link{convert.inp}} which prevented the code from working with
#' groups and two or more covariates.  Note that there are limitations to this
#' function which may require some minor editing of the file.  The limitations
#' have been added to the help file (\code{\link{convert.inp}}). } Version
#' 1.5.8 (22 June 2006) \itemize{ \item Argument \code{options} was added to
#' \code{\link{mark}} and \code{\link{make.mark.model}} with a default NULL
#' value. It is simply a character string that is tacked onto the \code{Proc
#' Estimate} statement for the MARK .inp file.  It can be used to request
#' options such as NoStandDM (to not standardize the design matrix) or
#' SIMANNEAL (to request use of the simulated annealing optimization method) or
#' any existing or new options that can be set on the estimate proc. \item A
#' bug in \code{\link{model.table}} was fixed so it would accomodate the change
#' from v1.3 to a marklist in which the model.table was switched to the last
#' entry in the list. \item A bug in \code{\link{summary.mark}} was fixed so it
#' would properly display QAICc when chat > 1. \item Function
#' \code{\link{adjust.chat}} was modified such that it returns a marklist with
#' each model having a new chat value and the model.table is adjusted for the
#' new chat value. \item Function \code{\link{adjust.parameter.count}} was
#' modfied so it returns the mark model object rather than using eval to modify
#' the object in place.  The latter does not work with models in a marklist and
#' calls made within functions. } Version 1.5.7 (8 June 2006) \itemize{ \item
#' Argument \code{data} was added to function \code{\link{model.average}} to
#' enable model averaging parameters at specific covariate values rather than
#' the mean value of the observed data.  An example is given in the help file.
#' \item Argument \code{parameter} of function \code{\link{model.average}} now
#' has a default of NULL and if it is not specified then all of the real
#' parameters are model averaged rather than those for a particular type of
#' parameter (eg p or Psi). \item A bug was fixed in function
#' \code{\link{compute.real}} that caused the function to fail for computations
#' of \code{Psi}. } Version 1.5.6 (6 June 2006) \itemize{ \item
#' \code{print.summary.mark} was modified so fixed parameters are noted. \item
#' Argument \code{show.fixed} was added to \code{\link{summary.mark}} to
#' control whether fixed parameters are shown as NA (FALSE) or as the value at
#' which they were fixed.  If \code{se=T} the default is \code{show.fixed=T}
#' otherwise \code{show.fixed=F}.  The latter is most useful in displaying
#' values in PIM format (without std errors), so fixed values are displayed as
#' blanks instead of NA. \item Argument \code{links} was added to
#' \code{\link{convert.link.to.real}} and the default value for argument
#' \code{model} is now NULL.  One or the other must be given. If the value for
#' \code{links} is given then they are used in place of the links specified in
#' the \code{model} object.  This provides for additional flexibility in
#' changing link values for computation (eg use of log with mlogit). \item
#' Argument \code{drop} was added to \code{\link{model.average}}.  If
#' \code{drop=TRUE} (the default), then any model with one or more non-positive
#' (0 or negative) variances is not used in the model averaging. \item An error
#' in computation of the v-c matrix of mlogit link values in
#' \code{compute.links.from.reals} was fixed. This did not affect confidence
#' intervals for real parameters (eg Psi) in \code{model.average} because it
#' uses the logit transformation for confidence intervals on real parameters
#' that use mlogit link (eg Psi). \item \code{\link{get.real}} was unable to
#' extract a single parameter value(eg constant Phi model).  This was fixed.
#' \item The argument \code{parm.indices} was removed from the functions
#' \code{\link{compute.real}} and \code{\link{convert.link.to.real}} because
#' the subsetting can be done easily with the complete results returned by the
#' functions.  This changed the examples in \code{\link{fill.covariates}}.
#' \item \code{\link{compute.real}} and subsequently \code{\link{get.real}}
#' return a field \code{fixed} when se=TRUE that denotes whether a real
#' parameter is a fixed parameter or an estimated parameter at a boundary which
#' is identified by having a standard error=0. } Version 1.5.5 (1 June 2006)
#' \itemize{ \item \code{model} has been deleted from the arguments in
#' TransitionMatrix.  It was only being used to ascertain whether the model was
#' a Multistrata model.  This is now determined more accurately by looking for
#' the presence of \code{tostratum} in the argument \code{x} which is a
#' dataframe created for \code{Psi} from the function \code{\link{get.real}}.
#' The function also works with the estimates dataframe generated from
#' \code{\link{model.average}}. See help for \code{\link{TransitionMatrix}} for
#' an example. \item An argument \code{vcv} was added to function
#' \code{\link{model.average}}.  If the argument is TRUE (the default value)
#' then the var-cov matrix of the model averaged real parameters is computed
#' and returned and the confidence intervals for the model averaged parameters
#' are constructed. Models with non-positive variances for betas are reported
#' and dropped from model averaging and the weights are renormalized for the
#' remaining models. \item A new function
#' \code{\link{compute.links.from.reals}} was added to the library to transform
#' real parameters to its link space.  It has 2 functions both related to model
#' averaged estimates. Firstly, it is used to transform model averaged
#' estimates so the normal confidence interval can be constructed on the link
#' values and then back-transformed to real space.  The second function is to
#' enable parametric bootstrapping in which the error distbution is assumed to
#' be multivariate normal for the link values. From a single model, the link
#' values are easily constructed from the betas and design matrix so this
#' function is not needed.  But for model averaging there is no equivalent
#' because the real parameters are averaged over a variety of models with the
#' same real parameter structure but differing design structures.  This
#' function allows for link values and their var-cov matrix to be created from
#' the model averaged real estimates. } Version 1.5.4 (30 May 2006) \itemize{
#' \item In function \code{\link{mark}} an argument \code{retry} was added to
#' enable the analysis to be re-run up to the number of times specified.  An
#' analysis is only re-run if there are "singular" beta parameters which means
#' that they are either non-estimable (confounded) or they are at a boundary.
#' Beginning with this version, \code{\link{extract.mark.output}} was modified
#' such that the singular parameters identified by MARK are extracted from the
#' output (if any) and the indices for the beta parameters are stored in the
#' list element \code{model$results$singular}. The default value for
#' \code{retry} is 0 which means it will not retry.  When the model is re-run
#' the initial values are set to the values at the completion of the last run
#' except for the "singular" parameters which are set to 0.  Using \code{retry}
#' will not help if the parameters are non-estimable.  However, if the
#' parameters are at a boundary because the optimization "converged" to a
#' sub-optimal set of parameters, then setting \code{retry} to 1 or a suitably
#' small value will often help it find the MLEs by moving away from the
#' boundary. If the parameters are estimable and setting \code{retry} does not
#' work, then it may be better to set new initial parameters by either
#' specifying their values or using a model with similar parameters that did
#' converge. \item A new function \code{\link{rerun.mark}} was created to
#' simplify the process of refitting models with new starting values when the
#' models were initially created with \code{\link{mark.wrapper}} which runs a
#' list of models by using all combinations of the formulas defined for the
#' various parameters in the model.  Thus, individual calls to \code{mark} are
#' not constructed by the user and re-running an analysis from the resulting
#' list would require constructing those calls. The argument
#' \code{model.parameters} is now stored in the model object and it is used by
#' this new function to avoid constructing calls to rerun the analysis.  With
#' this new function you only need to specify the model list element to be
#' refitted, the processed dataframe, the design data and the model list
#' element (or different model) to be used for initial values. See
#' \code{\link{rerun.mark}} for an example. \item To make
#' \code{\link{rerun.mark}} a viable approach for all circumstances, the
#' functions \code{\link{mark.wrapper}} and \code{\link{model.table}} were
#' modified such that models that fail to converge at the outset (i.e., does
#' not provide estimates in the output file) are stored in the model list
#' created by the former function and they are reported as models that did not
#' run and are skipped in the model.table by the latter function.  This enables
#' a failed model to be reanalyzed with \code{\link{rerun.mark}} using another
#' model that converged for starting values. } Version 1.5.3 (25 May 2006)
#' \itemize{ \item In function \code{\link{get.real}} a fix was made to
#' accommodate constant pims and a warning is given if the v-c matrix for the
#' betas has non-positive variances. \item In function
#' \code{\link{make.mark.model}}, the argument \code{initial} can now be a
#' single value which is then assigned as the initial value for all betas.  I
#' have found this useful for POPAN models.  For some models I have run, the
#' models fail to converge in MARK with the default initial values it uses (I
#' believe it uses \code{initial=0}). I have had better luck using
#' \code{initial=1}.  By allowing the use of a single value you can use the
#' same generic starting value for each model without figuring out the number
#' of betas in each model. Also note that you can specify another model that
#' has already been run to use as initial values for a new model and it will
#' match parameter values. \item A bad bug was fixed in \code{\link{cleanup}}
#' which was unfortunately deleting files containing "out", "inp", "res" or
#' "vcv" rather than those having these as extensions.  This happened without
#' your knowledge if you chose ask=FALSE.  Good thing I had a backup.  Anyhow,
#' I have now restricted it to files that are named by RMark with markxxxx.inp
#' etc where xxxx is a numeric value.  Thus if you assign your own basefile
#' name for output files you'll have to delete them manually.  Better safe than
#' sorry. } Version 1.5.2 (18 May 2006) \itemize{ \item Two new functions were
#' added in this version. \code{\link{convert.inp}} converts a MARK encounter
#' history input file to an RMark dataframe.  This will be particularly useful
#' for those folks who have already been using MARK. Instead of converting and
#' importing their data with \code{\link{import.chdata}} they can use the
#' \code{\link{convert.inp}} to import their .inp file directly.  It can also
#' be used to directly import any of the example .inp files that accompany MARK
#' and the MARK electronic book
#' (\url{http://www.phidot.org/software/mark/docs/book/}).  The second new
#' function is only useful for tutorials and for first time users trying to
#' understand the way RMark works.  The function \code{\link{PIMS}} displays
#' the full PIM structure or the simplified PIM structure for a parameter in a
#' model.  The user does not directly manipulate PIMS in RMark and they are
#' essentially transparent to the user but for those with MARK experience being
#' able to look at the PIMS may help with the transition. } Version 1.5.1 (11
#' May 2006) \itemize{ \item Functions \code{\link{compute.link}} and
#' \code{\link{get.link}} were added to compute link values rather than the
#' parameter estimates. \item A function \code{\link{convert.link.to.real}} was
#' added to convert link estimates to real parameter estimates. Previously a
#' similar internal function was used within \code{compute.real} but to provide
#' more flexibility it was put into a separate function. \item An argument
#' \code{beta} was added to \code{\link{get.real}} to enable it to be changed
#' in the computation of the real parameters rather than always using the
#' values in \code{results$beta}. \item A function
#' \code{\link{TransitionMatrix}} was added to create a transition matrix for
#' the Psi values.  It is provided for all strata including the
#' \code{subtract.stratum}.  Standard errors and confidence intervals can also
#' be returned. \item \code{\link{make.mark.model}} was modified to include
#' \code{time.intervals} as an element in the mark object. } Version 1.5.0 (9
#' May 2006) \itemize{ \item If output file already exists user is given option
#' to create mark model from existing files.  Only really useful if a bug
#' occurs (which occurred to me from 1.4.9 changes) and once fixed any models
#' already run can be brought into R by running the same model over and
#' specifying the existing base \code{filename}.  Base \code{filename} values
#' are no longer prefixed with MRK to enable this change. \item On occasion
#' MARK will complete the analysis but fail to create the v-c matrix and v-c
#' file.  The code has been modified to skip over the file if it is missing and
#' output a warning. \item Two new functions have been added to ease handling
#' of marklist objects. \code{\link{merge.mark}} merges an unspecified number
#' of marklist and mark model objects into a new marklist with an optional
#' model.table. \code{\link{remove.mark}} can be used to remove mark models
#' from a marklist.  See \code{\link{dipper}} for examples of each function.
#' \item Various changes were made to functions that compute real parameter
#' estimates, their standard errors, confidence intervals and
#' variance-covariance matrix.  The functions that were changed include
#' \code{\link{compute.real}},\code{\link{find.covariates}},\code{\link{get.real}},\code{\link{fill.covariates}}.
#' For examples, see help for latter two functions. } Version 1.4.9 (3 May
#' 2006) \itemize{ \item Argument \code{initial} of
#' \code{\link{make.mark.model}} was not working after model simplification was
#' added in v1.2.  This was modified to select initial values from the model
#' based on names of design matrix columns rather than column contents which
#' have different numbers of rows depending on the simplification. \item
#' \code{\link{extract.mark.output}} was fixed to extract the correct -2LnL
#' from the output file in situations in which initial values were specified. }
#' Version 1.4.8 (25 April 2006) \itemize{ \item Argument \code{silent} was
#' added to \code{\link{mark}} and \code{\link{mark.wrapper}} with a default
#' value of \code{FALSE}.  This overcomes the problem described above in 1.4.7.
#' \item Code was added to \code{\link{collect.model.names}} to prevent it from
#' tripping up when files contain an asterisk which R uses for special names.
#' \item Use of T and F was properly changed to TRUE and FALSE in various
#' functions to prevent errors when T or F are R objects. \item Code for naming
#' files was modified to avoid problems when more than 999 analyses were run in
#' the same directory. \item Bug in setting fixed parameters with argument
#' \code{fixed=list(index=,value=)} was corrected. \item Argument
#' \code{remove.intercept} was added to parameter definition to force removal
#' of intercept in designs with nested factor interactions with additional
#' factor variables (e.g.,
#' \code{Psi=list(formula=~sex+stratum:tostratum,remove.intercept=TRUE)}). }
#' Version 1.4.7 (10 April 2006) \itemize{ \item An error was fixed in the
#' \code{Psi} simplification code. Note that with the fix in 1.4.2 to trap
#' errors, a side effect is that non-trapped errors that occur in the R code
#' will now fail without any error messages.  If the error occurs in making the
#' model, then the model will not be run, but you will not receive a message
#' that the model failed.  I may have to make the error trapping a
#' user-settable option to provide better error tracking. } Version 1.4.6 (7
#' April 2006) \itemize{ \item Assurance code was added to test that the
#' mlogits were properly assigned.  An error message will be given if there has
#' been any unforeseen problem created by the simplification.  This eliminates
#' any need for the user to check them as described under 1.4.5 above. }
#' Version 1.4.5 (6 April 2006) \itemize{ \item For multistrata models, the
#' code for creating the mlogit links for Psi was not working properly if there
#' was more than one group.  This was fixed in this version. \item
#' Simplification of the PIMS has now been extended to include mlogit
#' parameters.  That was not a trivial exercise and while I feel confident it
#' is correct, double check the assignment of mlogit links for complex models,
#' as I have not checked many examples at present.  Within a stratum, the
#' corresponding elements for Psi for each of the tostratum (movement from
#' stratum to each of the other strata excluding the \code{subtract.stratum})
#' should have the same mlogit(xxx) value such that it can properly compute the
#' value for \code{subtract.stratum} by subtraction. } Version 1.4.4 (4 April
#' 2006) \itemize{ \item By including the test on model failure, errors that
#' would stop program were not being displayed.  This has been fixed in this
#' version. \item An error was fixed in using time-varying covariates when some
#' of the design data had been deleted. } Version 1.4.3 (30 March 2006)
#' \itemize{ \item Problem with pop up window has been fixed.  It will no
#' longer appear if the model does not converge but the model will show as
#' having failed. \item An error was fixed in extracting output from the MARK
#' output file when for some circumstances the label for beta parameters
#' included spaces.  This now works properly. } Version 1.4.2 (14 March 2006)
#' \itemize{ \item Errors in the FORTRAN code were preventing completion of
#' large batch jobs.  Now these errors are caught and models that fail are
#' reported and skipped over.  Unfortunately, it does require user intervention
#' to close the popup window.  Make sure you select Yes to close the window
#' especially if you use the default \code{invisible=FALSE} such that the
#' window does not appear.  If you select No, you will not able to close the
#' window and R will hang. \item A new list element was added to
#' \code{parameters} in \code{\link{make.design.data}} for parameters such as
#' \code{Psi} to set the value of \code{tostratum} that is computed by
#' subtraction. The default is to compute the probabilitity of remaining in the
#' stratum.  The following is an example with strata A to D and setting A to be
#' computed by subtraction for each stratum: \cr
#' ddl=make.design.data(data.processed, \cr
#' parameters=list(Psi=list(pim.type="constant",subtract.stratum=c("A","A","A","A")),
#' \cr p=list(pim.type="constant"),S=list(pim.type="constant"))) } Version
#' 1.4.1 (11 March 2006) \itemize{ \item A value "constant" was added for the
#' argument \code{pim.type}. Note that \code{pim.type} is only used for
#' triangular PIMS. See \code{\link{make.design.data}} \item Some code changes
#' were made to \code{\link{make.mark.model}} which lessen time to create the
#' MARK input file for large models. \item Function
#' \code{\link{add.design.data}} was modified to accomodate robust design and
#' deletion of design data; this was missed in v1.4 changes. \item
#' \code{model.name} argument in \code{\link{mark}} and
#' \code{\link{make.mark.model}} was not working.  This was fixed. } Version
#' 1.4 (9 March 2006) \itemize{ \item Robust design models added. See
#' \code{\link{robust}} for an example. \item Function \code{\link{cleanup}}
#' was modified so warning messages/errors do not occur if no models/files are
#' found. \item Parameters in the design matrix are now ordered in the same
#' consistent arrangement.  In prior versions they were arranged based on their
#' order in the argument call. \item Argument \code{right} was added to
#' \code{\link{make.design.data}}, \code{\link{add.design.data}} and in
#' \code{design.parameters} of \code{\link{make.mark.model}} to control whether
#' bins are inclusive on the right (default).  The \code{\link{robust}} example
#' uses this argument in a call to \code{\link{mark}}. } Version 1.3 (22 Feb
#' 2006) \itemize{ \item Time varying covariates can now be included in the
#' model formula. See \code{\link{make.mark.model}} for details. \item New
#' model types for Known (Known-fate) and Multistrata (CJS with different
#' strata) were added. See \code{\link{Blackduck}} and \code{\link{mstrata}}
#' for examples. \item Specific rows of the design data can now be removed for
#' parameters that should not be estimated.  Default fixed values can be
#' assigned.  The function \code{\link{make.design.data}} now accepts an
#' argument \code{remove.unused} which can be used to automatically remove
#' unused design data for nested models. It's behavior is also determined by
#' the new argument \code{default.fixed} in \code{\link{make.mark.model}}.
#' \item \code{\link{summary.mark}} now produces a summary object and
#' \code{\link{print.summary.mark}} prints the summary object. Changes were
#' made to output when \code{se=T}. \item A new function
#' \code{merge.occasion.data} was created to add occasion specific covariates
#' to the design data. \item New functions \code{\link{mark.wrapper}} and
#' \code{\link{create.model.list}} were created to automate running models from
#' a set of model specifications for each model parameter. \item The argument
#' \code{begin.time} in \code{\link{process.data}} can now be a vector to
#' enable a different beginning time for each group. \item An argument
#' \code{pim.type} was added to parameter specification to enable using pims
#' with time structure for data sets with a single release cohort for CJS. See
#' \code{\link{make.design.data}} \item Model lists created with
#' \code{\link{collect.models}} are now given the class "marklist" which is
#' used with \code{\link{cleanup}} and \code{\link{print.marklist}} (see
#' \code{\link{print.mark}}). \item The function \code{\link{collect.models}}
#' now places the model.table at the end of the returned list such that each
#' model number in model.table is now the element number in the returned list.
#' Previously it was 1+ that number. \item Input, output, v-c and residual
#' results files from MARK are now stored in the directory containing the
#' \code{.Rdata} workspace.  They are numbered consecutively and the field
#' \code{output} contains the base filename.  The function
#' \code{\link{cleanup}} was created to delete files that are no longer linked
#' to \code{mark} or \code{marklist} objects. \item Model averaged estimates
#' and standard errors of real parameters can be obtained with the function
#' \code{\link{model.average}}. } 
#' Version 1.2 (4 Oct 2005) \itemize{ \item By
#' default the PIM structure is simplified to use the fewest number of unique
#' parameters.  This reduces the size of the design matrix and should reduce
#' run times. \item The above change was made in some versions still numbered
#' 1.1, but it contained an error that caused the links command for MARK to be
#' constructed incorrectly. \item \code{adjust} argument has been added to
#' \code{\link{collect.models}} to enable control of number of parameters and
#' resulting AIC values. \item \code{model.list} in \code{\link{model.table}}
#' and \code{\link{adjust.chat}} can now also be a list of models created by
#' \code{\link{collect.models}} which allows operating on sets of models. }
#' 
#' @name Whatsnew
#' @author Jeff Laake
NULL


#' A beginners introduction and guide to RMark
#' 
#' The RMark package is a collection of R functions that can be used as an
#' interface to MARK for analysis of capture-recapture data.
#' 
#' The library contains various functions that import/export capture data,
#' build capture-recapture models, run the FORTRAN program MARK.EXE, and
#' extract and display output.  Program MARK has its own user interface;
#' however, model development can be rather tedious and error-prone because the
#' parameter structure and design matrix are created by hand. This interface in
#' R was created to use the formula and design matrix functions in R to ease
#' model development and reduce errors. This R interface has the following
#' advantages: 1) Uses model notation to create design matrices rather than
#' designing them by hand in MARK or in EXCEL, which makes model development
#' faster and more reliable. All-different PIMS are automatically created for
#' each group (if any). 2) Allows models based on group (factor variables) and
#' individual covariates with groups created on the fly. Age, cohort, group and
#' time variables are pre-defined for use in formulas. 3) Both real and beta
#' labels are automatically added for easy output interpretation.  4) Input,
#' output and specific results (eg parameter estimates, AICc etc) are stored in
#' an R object where they can be manipulated as deemed useful (eg plotting,
#' further calculations, simulation etc). 5) Parameter estimates can be
#' displayed in triangular PIM format (if appropriate) for ease of
#' interpretation. 6) Easy setup of batch jobs and the calls to the R functions
#' document the model specifications and allow models to be easily reproduced
#' or re-run if data are changed. 7) Covariate-specific estimates of real
#' parameters can be computed within R without re-running the analysis.
#' 
#' The MARK capture-recapture models that are currently
#' supported are provided in MarkModels.pdf which is installed in the RMark directory of your R library.
#' You can also find a list in MARK under Help/Data Types.
#' ' 
#' There is one limitation of this interface. All models
#' in this interface are developed via a design matrix approach rather than
#' coding the model structure via parameter index matrices (PIMS). In most
#' cases, a logit or other link is used by default which has implications for
#' ability of MARK to count the number of identifiable parameters (see
#' \code{\link{dipper}} for an example).  However, beginning with v1.7.6 the
#' sin link is now supported if the formula specifies an identity design matrix
#' for the parameter.
#' 
#' Before you begin, you must have installed MARK
#' (\url{http://www.cnr.colostate.edu/~gwhite/mark/mark.htm}) on your computer
#' or at least have a current copy of MARK.EXE. As long as you selected the
#' default location for the MARK install (c:/Program Files/Mark), the
#' \code{RMark} library will be able to find it.  If for some reason, you choose
#' to install it in a different location, see the note section in
#' \code{\link{mark}} for instructions on setting the variable MarkPath to
#' specify the path.  In addition to installing MARK, you must have installed
#' the \code{RMark} library into the R library directory.  Once done with those
#' tasks, run R and enter library(RMark) (or put it in your .First function) to
#' attach the library of functions.
#' 
#' The following is a categorical listing of the functions in the library with
#' a link to the help for each function. To start, read the help for functions
#' \code{\link{import.chdata}} and \code{\link{mark}} to learn how to import
#' your data and fit a simple model.  The text files for the examples shown in
#' \code{import.chdata} are in the subdirectory data within the R Library
#' directory in RMark. Next look at the example data sets and analyses
#' \code{\link{dipper}}, \code{\link{edwards.eberhardt}}, and
#' \code{\link{example.data}}. After you see the structure of the examples and
#' the use of functions to fit a series of analyses, explore the remaining
#' functions under Model Fitting, Batch Analyses, Model Selection and Summary
#' and Display.  If your data and models contain individual covariates, read
#' the section on Real Parameter Computation to learn how to compute estimates
#' of real parameters at various covariate values.
#' 
#' Input/Output data & results
#' 
#' \code{\link{import.chdata}},\code{\link{read.mark.binary}},
#' \code{\link{extract.mark.output}}
#' 
#' Exporting Models to MARK interface
#' 
#' \code{\link{export.chdata}}, \code{\link{export.model}}
#' 
#' Model Fitting
#' 
#' \code{\link{mark}}, \code{\link{process.data}},
#' \code{\link{make.design.data}}, \code{\link{add.design.data}},
#' \code{\link{make.mark.model}}, \code{\link{run.mark.model}}
#' \code{\link{merge_design.covariates}}
#' 
#' Batch analyses with functions
#' 
#' \code{\link{run.models}}, \code{\link{collect.models}},
#' \code{\link{create.model.list}}, \code{\link{mark.wrapper}}
#' 
#' Summary and display
#' 
#' \code{\link{summary.mark}}, \code{\link{print.mark}},
#' \code{\link{print.marklist}}, \code{\link{get.real}},
#' \code{\link{compute.real}}, \code{\link{print.summary.mark}}
#' 
#' Model Selection/Goodness of fit
#' 
#' \code{\link{adjust.chat}}, \code{\link{adjust.parameter.count}},
#' \code{\link{model.table}} , \code{\link{release.gof}},
#' \code{\link{model.average}}
#' 
#' Real Parameter computation
#' 
#' \code{\link{find.covariates}}, \code{\link{fill.covariates}},
#' \code{\link{compute.real}} , \code{\link{covariate.predictions}}
#' 
#' Utility and internal functions
#' 
#' \code{\link{collect.model.names}}, \code{\link{compute.design.data}},
#' \code{\link{extract.mark.output}}, \code{\link{inverse.link}},
#' \code{\link{deriv.inverse.link}}, \code{\link{setup.model}},
#' \code{\link{setup.parameters}}, \code{\link{valid.parameters}},
#' \code{\link{cleanup}}
#' 
#' For examples, see \code{\link{dipper}} for CJS and POPAN, see
#' \code{\link{example.data}} for CJS with multiple grouping variables, see
#' \code{\link{edwards.eberhardt}} for various closed-capture models, see
#' \code{\link{mstrata}} for Multistrata, and see \code{\link{Blackduck}} for
#' known fate. The latter two are examples of the use of
#' \code{\link{mark.wrapper}} for a shortcut approach to creating a series of
#' models. Other examples have been added for the various other models. In MarkModels.pdf it also
#' lists the name of examples that are provided for each model.
#' 
#' @name ABeginnersGuide
#' @aliases ABeginnersGuide
#' @author Jeff Laake
#' @references MARK: Dr. Gary White, Department of Fishery and Wildlife
#' Biology, Colorado State University, Fort Collins, Colorado, USA
#' \url{http://www.cnr.colostate.edu/~gwhite/mark/mark.htm}
#' @keywords utility
NULL

