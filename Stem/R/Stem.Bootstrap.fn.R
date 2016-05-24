`Stem.Bootstrap.fn` <-
function(x, StemModel, seed.list.out, output.kriging=NULL){

###########
#1) fix a seed x
#2) simulate a new data set z
#3) estimate the parameter vector phi
#4) predict in new spatial locations
##########
	##STEP 1: FIX THE SEED 
	#x = seed
	set.seed(x)
	position   = which(unlist(seed.list.out)== x)
	cat(paste("=========================== BOOTSTRAP ITERATION N.",position,"=========================== "),"\n")
	#write(position , file = "position.out", append=T)


        ##STEP 2: SIMULATIOM
        StemModel$skeleton$phi = StemModel$estimates$phi.hat   #it follows that the initial values are given by the ML Estimates  
	simulated.z	= Stem.Simulation(StemModel = StemModel) 
	StemModel$data$z = simulated.z

	###STEP 3: PARAMETER ESTIMATION
	MLE 	= Stem.Estimation(StemModel = StemModel)
		
        ##STEP4: SPATIAL PREDICTION
	#data.newlocations	= output.kriging$data.newlocations
	#time.point 			= output.kriging$time.point
	#~ spat.predictions  =  Stem.Kriging(StemModel 		= StemModel,
						#~ output.estimation		= MLE,
						#~ coord.newlocations 	= data.newlocations$coord.newlocations,
						#~ covariates.newlocations= data.newlocations$covariates.newlocations,
						#~ K.newlocations 		= data.newlocations$K.newlocations,
						#~ time.point 			= time.point,
						#~ regular.grid 		= output.kriging$regular.grid)

     
	    

	return(MLE = MLE) #, spat.predictions=spat.predictions))
}

