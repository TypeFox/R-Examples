  
######### Test script for DivE ###########
# This script runs the DiveMaster function from the DivE package on the example dataset 'Bact1' (using only 2 models)
# Required package: DivE

###########################################
#### Step 1. Load package ###
###########################################
require(DivE)

###########################################
#### Step 2. Load data files ####
###########################################
data(Bact1) 		# Example sample dataset
data(ModelSet) 		# 58 models to fit
data(ParamRanges) 	# Parameter ranges for the 58 models
data(ParamSeeds) 	# 58 sets of candidate initial parameters


###########################################
#### Step 3. Truncate 58 models to 2 models (for quick demonstration purposes only) ####
###########################################
testmodels 	<- list()
testmeta 	<- list()
paramranges <- list()
Mods 		<- c(1,2)		

for (i in 1:length(Mods)) 
{
	testmodels 			<- c(testmodels 	,	ModelSet[  Mods[i] ]		)
	testmeta	[[i]]	<- ParamSeeds	[[ Mods[i] ]]
	paramranges	[[i]] 	<- ParamRanges	[[ Mods[i] ]]
}


###########################################
#### Step 4 (simple version). Run the DiveMaster function ####
###########################################


################
# SIMPLEST METHOD, USING DiveMaster 
################

# With two samples (Main sample + one subset)
result <- DiveMaster(models=testmodels, init.params=testmeta, param.ranges = paramranges,
		main.samp=Bact1, subsizes=2, NResamples=50, nrf=10, fitloop=1, numit=50) # default parameters are: nrf=1, NResamples=1e3, numit=1e5. Values here chosen to speed example up

################
# LONGER METHOD, using component functions 
################

#### Step 4 (component function version). ####
# Again, with two samples (Main sample + one subset)

###########################################
### 4.1 create rarefaction data (divsubsamples object)
###########################################

Bact1length = sum(Bact1$Count)
dss_1 	<- divsubsamples(Bact1, nrf=2, minrarefac=1, maxrarefac=0.5*Bact1length	, NResamples=30)
dss_2 	<- divsubsamples(Bact1, nrf=2, minrarefac=1, maxrarefac=Bact1length		, NResamples=30)
dss 	<- list(dss_2, dss_1)

###########################################
### 4.2 fit models (create fitsingleMod object)
###########################################
fmm 		<- list()	## list of fitted model data
for (i in 1:length(Mods)) 
{
	fsm.temp <- fitsinglemod(model.list = testmodels[i]	, init.param = testmeta[[i]],  
			
			param.range = paramranges[[i]], 
			numit 		= 10^2, 
			varleft 	= 1e-8, 
			fitloops 	= 1,
			minplaus 	= 10 ,
			tot.pop 	= 100 * Bact1length,			
			nrf   		= 2, 
			minrarefac	= 1, 
			NResamples 	= 30, 
			main.samp 	= Bact1, 	
      
      # Option 1:	provide previously calculated rareafaction data - keeps rarefaction data consistent across models
			subsizes   	= c(Bact1length, 0.5*Bact1length), 
      		data.default = FALSE,
			dssamps   	= dss 
      
      ## Option 2:	calculate rarefaction data separately for each model - not recommended.
			#subsizes 	= 2, 
			#data.default = TRUE
	)  
	fmm[[fsm.temp$modelname]] <- fsm.temp
}

###########################################
### 4.3 score models (create list of class scoresinglemod)
###########################################

num.mod 		<- length(Mods)
ssm				<- matrix(rep(NA, num.mod*4), nrow=num.mod, ncol=4)	## list of model scores
colnames(ssm) 	<- c("fit", "accuracy", "similarity", "plausibility")

mod.rownames	<- matrix(rep(NA, num.mod), nrow = num.mod, ncol = 1)
mod.score 		<- matrix(rep(NA, num.mod), nrow = num.mod, ncol = 1)
TopX			<- 5		## 

for (i in 1:length(Mods)) 
{
	ssm.temp <- scoresinglemod(fsm = fmm[[i]], precision.lv = c(1e-04, 
					0.005, 0.005), plaus.pen = 500 )
	mod.score.temp <- combine.criteria(ssm = ssm.temp, crit.wts = c(1, 1, 1, 1))
	ssm[i, 1] <- ssm.temp$fit
	ssm[i, 2] <- ssm.temp$accuracy
	ssm[i, 3] <- ssm.temp$similarity
	ssm[i, 4] <- ssm.temp$plausibility
	
	mod.score[i, ] 	<- mod.score.temp
	mod.rownames[i] <- fmm[[i]]$modelname
}

rownames(mod.score) <- mod.rownames
colnames(mod.score) <- "Combined score"
rownames(ssm) 		<- mod.rownames

###########################################
### 4.3 Produce Estimates
###########################################

m <- min(TopX, num.mod)
topX_scores <- sort(unique(mod.score))[1:m]		
lenX <- length(topX_scores)
topX_index <- which(mod.score %in% topX_scores)
predX.vector <- c()
for (i in 1:lenX)
{
	tmp <- ((fmm[[topX_index[i]]])$global)[1, ]
	predX.vector <- c(predX.vector, tmp)
}
PointEstimate 	<- geo.mean(predX.vector)
UpperBound 		<- max(predX.vector)
LowerBound 		<- min(predX.vector)



###########################################
#### Step 5. View outputs ####
###########################################
result # Comparison of combined scores
summary(result) # Summary of combined score
result$ssm # Detailed comparison of scores
result$fmm$logistic # Fit details (model 1 - logistic model)
result$fmm$negexp # Fit details (model 2 - negative exponential model)
summary(result$fmm$logistic) # Fit summary (model 1)
result$fmm$logistic$param # Individual fit details ($param example)
plot(result$fmm$logistic) # Local plot of model 1 fit
plot(result$fmm$logistic, range="global") # Global plot of model 1 fit
plot(result$fmm$negexp) # Local plot of model 2 fit
plot(result$fmm$negexp, range="global") # Global plot of model 2 fit

popdiversity(result, 10^6)	## calculate diversity at another population size. 

###########################################
#### Step 6. Miscellaneous ####
###########################################
?DiveMaster # View main help file
?DivE # Package summary

# Component functions
?divsubsamples
?fitsinglemod
?scoresinglemod







