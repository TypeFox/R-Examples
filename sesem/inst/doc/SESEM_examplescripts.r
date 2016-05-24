#Example Script for Spatial SEM

# Updated Jan 22, 2014

# Citation: Lamb, E.G., K. Mengersen, K.J. Stewart, U. Attanayake, S.D. Siciliano. 2014. Spatially explicit 
# structural equation modeling. Ecology xx:xxx-xxx.

# Functions here were tested using R 3.0.2
# libraries "lavaan", "gplots", and "mgcv" are required


# The following examples will walk you through the example SEM analyses provided in the text and appendices to Lamb et al. 2014. Ecology

#Step 1 Load library and data

library(sesem)

truelove#Truelove lowland datafile
plantcomp #plant competition datafile
alexfiord #alexandra fiord datafile

# Note that if you are using the script files provided in the appendices to Lamb et al 2014
# you will need to download the following files to your working directory
#	truelove <- read.csv("rawdata_truelove.csv") # read Truelove lowland datafile
# 	alexfiord <- read.csv("rawdata_Alexfiord.csv")# read Alexandra Fiord datafile
#	plantcomp <- read.csv("rawdata_plantcomp.csv")# read Plant Competition datafile

#Step 2 calculate distance matrix and specify lag distance bin sizes

	# function calc.dist(datafile) calculates and stores the distance matrix for all X-Y pairs
	#	requires a dataset where the first two columns are x and y coordinates

truelovedist<-calc.dist(truelove)
alexfiorddist<-calc.dist(alexfiord)
plantcompdist<-calc.dist(plantcomp)

	# Need to specify 2 vectors for bin sizes:
	#	binsize is a vector of lag distances starting at 0. The following values are the upper limits of each distance bin
	#		binsize should have n+1 elements where n is the number of lag distance bins desired
	#	bin.name is a vector of unique names for each lag distance bin. there should be n elements, one for each bin

	# Bin sizes and names can be specified manually or using function make.bin

#make.bin(dist.mat,type="ALL",p.dist=50,n.bins=10,s.size=100)

	# generates cut off values for lag distance bins and corresponding bin names
	# The function has three default parameter values available, if user does not want to specify:
		#Inference distance as a percentage(p.dist) = 50%
		#Number of bins (n.bins) = 10
		#Sample size (s.size) = 100
	#dist.mat is a distance matrix produced by calc.dist
	#Can use type="ALL","n.bins" OR "s.size" to control parameter values.
	#The function produces a list object containing (1.)binsize and (2.)binname
	#These two vectors (binsize and binname) will be used by make.covar to calculate variance covariance matrices for each lag distance bin
	
	#Special note:
	#User specified number of lag distance bins OR sample size 
	#will be used to calculate initial cutoff value of each lag distance bin.
	#However, if the cutoff point is in between a lag distance bin, 
	#real cutoff will apply at the upper margin of the particular bin. 
	#Therefore, resulting number of bins are less than or equal AND 
	#resulting sample sizes are greater than or equal to the value specified by the user.
	
# function make.bin to generate all possible bins the Truelove lowland dataset to a maximum inference distance of 50%

#	note the need to extract sizes and names from the bin size object as shown below

Truelove_bins<-make.bin(truelovedist,type="ALL") #inference distance=50% AND number of bins=ALL
truebinsize<-Truelove_bins[1][[1]] #truelove lowland bin sizes
truebinname<-Truelove_bins[2][[1]] #truelove lowland bin names

# plotbin(distancematrix,binsize) provides a histogram of all pairs in the distance matrix and a plot of sample 
# 	size for the selected bin size.
	
plotbin(truelovedist,truebinsize)

	
# manual bin sizes for the Alexandra Fiord dataset
#	note that the bin name vector must have one less element than bin size
	
alexbinsize = c(0,1,2,2.2,4,5,8,16,32,64,96,128,160) #alex fiord bin sizes
alexbinname=c("Bin1","Bin2","Bin3","Bin4","Bin5","Bin6","Bin7","Bin8","Bin9","Bin10","Bin11","Bin12")#alex fiord Bin names

plotbin(alexfiorddist,alexbinsize)

# function make.bin to generate bins with 200 samples each to a maximum inference distance of 50% 
#	for the Plant Competition dataset

#	note the need to extract sizes and names from the bin size object as shown below

Plant_bins<-make.bin(plantcompdist,type="s.size",s.size=200) #inference distance=50% AND sample size=200
plantbinsize<-Plant_bins[1][[1]] #plant competition bin sizes
plantbinname<-Plant_bins[2][[1]] #plant competition bin names

plotbin(plantcompdist,plantbinsize)

	#User can use type="ALL","n.bins" OR "s.size"
	#make.covar uses binsize and binname, only if attibutes=NULL and is.vector=TRUE.


#Step 3 Calculate covariance matrices
	
	# function make.covar(datafile,distancematrix,binsize,binname) calculates variance covariance matrices
	#	for each lag distance bin and for a flat (non-spatial) bin
	#	Values in the bin summary include the low and high cut points for each bin and the count "aa.lagcnt" of 
	#	samples in each bin
	
	# this step can be computationally time consuming; the resulting object can be saved for later sesem analysis
	
truelove_covar<-make.covar(truelove,truelovedist,truebinsize,truebinname)
alex_covar<-make.covar(alexfiord,alexfiorddist,alexbinsize,alexbinname)
plant_covar<-make.covar(plantcomp,plantcompdist,plantbinsize,plantbinname)


# Step 4 Run Models

# function runModels()  
# 	lavaan library must be installed. 
#   need to specifify path model using lavaan syntax. 
#	See http://cran.r-project.org/web/packages/lavaan/lavaan.pdf for details
#	list object created by runModels() used in subsequent steps

#Truelove  and Alexandra Fiord Path Models
spatial_model<-'
	Gram ~ Moisture
	N_Fix ~ Bryoph + Lich + SoilCrust
	SoilCrust ~ Bryoph + Lich + Gram + Shrubs + Forbs	
	Bryoph ~ Gram + Shrubs + Forbs + Moisture
	Lich ~ Moisture + Forbs + Gram + Shrubs + Bryoph
	Forbs ~ Moisture
	Gram ~~ Forbs
	Shrubs ~ Moisture	
	Gram ~~ Shrubs
	Shrubs ~~ Forbs
	'

#Plant Competition Path Model
spatial_model_plant<-'
	SoilMoist ~ Topog.Pos
	ShootBio ~ TotalN + Topog.Pos + SoilMoist
	RootBio ~ TotalN + Topog.Pos + SoilMoist
	ShootBio ~~ RootBio
	LightInt ~ ShootBio + Topog.Pos
	Comp.Intensity ~ LightInt + RootBio + TotalN + SoilMoist
	SpRich ~ Comp.Intensity + SoilMoist + RootBio + ShootBio + Topog.Pos
	'

spatial_model_results_true<-runModels(spatial_model,truelove_covar)
spatial_model_results_alex<-runModels(spatial_model,alex_covar)
spatial_model_results_plant<-runModels(spatial_model_plant,plant_covar)

# Step 5 Extract analysis results

	# modelsummary() extracts basic model summary information from the object
	# created by runModels in a readable format
	
modelsummary(spatial_model_results_true)
modelsummary(spatial_model_results_alex)
modelsummary(spatial_model_results_plant)

# bin.results(spatial_model_results,bin="binflat") 
#	extracts path coefficients, standard errors, and standardized coefficients for a particular bin in a readable format

bin.results(spatial_model_results_true,bin="Bin2")# prints results for bin2 
bin.results(spatial_model_results_alex,bin="binflat") # prints results for flat (nonspatial) bin


# bin.rsquare(spatial_model_results,bin="binflat") 
#	extracts r-square values for each dependent variable for a particular bin

bin.rsquare(spatial_model_results_plant,bin="binflat")

# Step 6a plot model fit results

#	plotmodelfit(spatial_model_results,plots="all",add.line="none",rmsea_err=T,pch=16,lwd=2,lty=1)
#	model fit plotting with optional lines, and formatting

#	spatial_model_results = object output by runModels
#	plots = options for selecting fit indices to plot
#		"all" indicats all of chi square, cfi, rmsea, and srmr to be plotted 
#		"chi", "cfi", "rmsea", "srmr" select individual metrics for plotting
#	add.lines. options for plotting a fit line 
#		"none" indicates no line
#		"step" plots straight line segments between points
#		"smooth" plots a smoothed curve fit using function lowess. Smoothed lines do not include the flat model (lag distance zero)
#	rmsea_err should confidence limits for rmsea be plotted?
#	pch, lwd, lty options for formatting points and fit lines.

plotmodelfit(spatial_model_results_true) # plots all metrics without lines and with rmsea errors
											# note that warnings will arise from the plotting of rmsea errors with zero length
											#these can be ignored
plotmodelfit(spatial_model_results_alex,add.line="step",rmsea_err=F,lwd=2) # plots stepped fit line
plotmodelfit(spatial_model_results_plant,add.line="smooth",rmsea_err=F,pch=16,lty=1) # plots smoothed fit line

plotmodelfit(spatial_model_results_true,plot="cfi",add.line="smooth",pch=16,lty=1,lwd=4,cex.lab=1.5,cex=2,cex.axis=1.5) #plots only cfi

# Step 6b plot path coefficent changes with lag distances

#	plotpath(spatial_model_results,path.type="directed",add.line="none",rmsea_err=T,pch=16,lwd=2,lty=1)
#   
#	spatial_model_results = object output by runModels
#	path.type= options for selecting which paths to plot
#		"directed" = only directed paths plotted
#		"undirected" = only undirected correlations plotted
#		"both" = all paths plotted
#		"user" = allows user to specify particular paths and a particular order for plotting
#			argument selectpath must also be provided with path.type="user"
#	selectpath = required when path.type="user"
#		usage is selectpath==c(5,18,16,23,29) where values refer to path numbers
#		path numbers can be obtained using spatial_model_results[2]
#	add.lines = options for plotting a fit line 
#		"none" indicates no line
#		"step" plots straight line segments between points
#		"smooth" plots a smoothed curve fit using function lowess
#	add.error = should standard error bars be added for each path coefficient
#	pcut = p-value cutoff above which points with non significant p-values are shaded grey
#		set pcut=1 to have all points black
#	pch, lwd, lty cex, cex.lab, cex.axis, cex.main options for formatting points and fit lines.

plotpath(spatial_model_results_true,pch=11) # plots 
plotpath(spatial_model_results_plant,pch=16,path.type="directed",add.line="smooth")
plotpath(spatial_model_results_alex,path.type="both",add.error=F,add.line="step",pcut=0.01,lwd=1)
plotpath(spatial_model_results_true,path.type="user",selectpath=c(5,18,16,23,29),add.error=F,add.line="smooth",pcut=0.01,lwd=1)

#
#Step 6c plot changes in path coefficents using generalized additive models

#	gam.path(spatial_model_results,path.type="directed",selectpath="none selected",
#		plot.points=T,se.plot=T,lwd.pred=2,lty.pred=1,lwd.se=2,lty.se=3,cex=1,cex.axis=1,cex.lab=1
#		,xlab="Lag Distance",ylab="Unst. Path Coeff.",yaxt="s",xaxt="s")
#
#	library "mgcv" is required
   
#	This function prints some of the gam results for each path in the model as well as producing figures
   
#	spatial_model_results = object output by runModels
#	path.type= options for selecting which paths to plot
#		"directed" = only directed paths plotted
#		"undirected" = only undirected correlations plotted
#		"both" = all paths plotted
#		"user" = allows user to specify particular paths and a particular order for plotting
#			argument selectpath must also be provided with path.type="user"
#	selectpath = required when path.type="user"
#		usage is selectpath==c(5,18,16,23,29) where values refer to path numbers
#		path numbers can be obtained using spatial_model_results[2]
#	plot.points = should points be plotted on figure
#	se.plot = should standard errors for the prediction lines be printed
#	additional arguments formatting figure

#displays gam relationships for all directed paths in model
par(mar=c(2.5,2.5,1.5,1))
gam.path(spatial_model_results_true,plot.points=F,se.plot=T)

#displays gam relationships for a selection of paths in model
gam.path(spatial_model_results_plant,path.type="user",selectpath=c(5,18,16,14,7),plot.points=T,se.plot=T)

gam.path(spatial_model_results_alex,path.type="user",selectpath=c(5,18,16,14,7),plot.points=F,se.plot=F,lwd.pred=2,ylab="",xlab="",yaxt="n",xaxt="n")


#Step 7. Optimize model

# avg.modindices() extracts modification indices for all models in the object produced by runModels() 
# 		and summarizes the mod indices by taking the mean for each possible additional path
#		across all lag distance bins. The flat model is not included in these calculations
#	 modcut eliminates printing of average MI values below the cutoff. The default is 4
# 
# Based on the modification indices create a new lavaan path model description and return to step 4.


avg.modindices(spatial_model_results_true,modcut=3)# display modindices
avg.modindices(spatial_model_results_alex)# display modindices
avg.modindices(spatial_model_results_plant,modcut=1)# display modindices

