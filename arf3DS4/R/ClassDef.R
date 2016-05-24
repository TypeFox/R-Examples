#############################################
# arf3DS4 S4 CLASS DEFINITIONS				#
# Wouter D. Weeda							#
# University of Amsterdam					#
#############################################

#[CONTAINS]
#version			[user]
#settings			[user]
#registration		[user]
#functional			[user]
#experiment			[user]
#options			[user]
#nifti.fileinfo		[user]
#nifti.header		[user]
#fmri.data			[user]
#wald				[user]
#data				[user]
#model				[user]
#sequence			[user]

## arf version class (version is set here)
setClass(
	Class='version',
	representation=representation(
		version='numeric',
		build='numeric',
		update='numeric',
		svnrev='numeric'
	),
	prototype=prototype(
		version=2,
		build=5,
		update=10,
		svnrev=219
	),
	package='arf3DS4'
)

## arf settings class (containing all standard directory names)
setClass(
	Class='settings',
	representation=representation(
		expRda='character',				#name of experimentobject
		optionsRda='character',			#name of optionsobject
		startRda='character',			#name of startvec file
		dataRda='character',			#name of dataobject in each datadir
		modelRda='character',			#name of modeldatafile
		statsRda='character',			#name of statsfile
		regRda='character',
		funcRda='character',
		subjectPrefix='character', 		#subjects prefix
		conditionPrefix='character', 	#condition prefix
		modelPrefix='character',		#model prefix
		subjectDir='character', 		#subjects prefix
		conditionDir='character', 		#condition prefix
		dataDir='character',			#data directory name
		weightsDir='character',			#weights directory name
		avgDir='character',				#avg directory name
		regDir='character',
		funcDir='character',
		betaDir='character',			#beta directory
		modelDir='character',			#model directory
		statsDir='character',			#stats directory
		modeldatDir='character', 		#Residual/Derivatives/Weights directory
		avgdatFile='character',			#averageDataFileName
		avgWFile='character',			#averageWeightFileName
		avgtstatFile='character',		#averagetstatFile
		modelDataFile='character',		#name of modelDataiFile
		modelnamesRda='character',		#modelNames file
		residualFile='character',		#name of residual binary
		derivativeFile='character',		#name of derivative binary
		weightFile='character',
		lowresFile='character',
		lowresAvg='character',
		logFile='character',
		version='ANY'
	),
	prototype=prototype(
		expRda='experiment.Rda',		#name of experimentobject
		optionsRda='options.Rda',		#name of optionsobject
		startRda='start.Rda',			#startvec
		dataRda='data.Rda',				#name of dataobject in each data
		modelRda='model.Rda',			#name of modelfile
		statsRda='stats.Rda',			#name of statsfile
		regRda='registration.Rda',
		funcRda='functional.Rda',
		subjectPrefix='', 				#subjects prefix
		conditionPrefix='', 			#condition prefix
		subjectDir='subjects', 			#subjects prefix
		conditionDir='conditions', 		#condition prefix
		modelPrefix='model',			#model prefix
		dataDir='data',					#data directory name
		weightsDir='weights',			#weights directory name
		avgDir='avg',					#avg directory name
		regDir='reg',
		funcDir='funcs',
		betaDir='beta',					#beta directory
		modelDir='models',				#model directory
		statsDir='stats',				#stats directory
		modeldatDir='data', 			#Residual/Derivatives/Weights directory	
		avgdatFile='avgdata',			#averageDataFileName
		avgWFile='avgweight',			#averageWeightFileName
		avgtstatFile='avgtstat',		#average tstat filename
		modelDataFile='avgmodel',		#name of modelNiftiFile
		modelnamesRda='modelnames.Rda', #modelNames file
		residualFile='residuals.bin',	#name of residual binary
		derivativeFile='derivs.bin',	#name of derivative binary
		weightFile='weights.bin',
		lowresFile='lowres_struct',
		lowresAvg='avglowres_struct',
		logFile='arfprocess.log',
		version=new('version')
		),
		package='arf3DS4'
)

## arf/FSL registration class
setClass(
	Class='registration',
	representation=representation(
		fullpath='character',		#fullpath of regfile
		filename='character',		#filename of regfile
		linkedfile='character',		#link to betafile
		examp2high='character',		#examp2high ASCII mat (affine transform)
		high2stand='character',		#high2standard ASCII mat (affine transform)
		examp2stand='character',	#examp2high ASCII mat (affine transform)
		example='character',		#example.nii.gz (corr to betafile)
		highres='character',		#highres.nii.gz (subject anatomical scan)
		standard='character',		#standard.nii.gz (FSL MNI152)
		Dex='matrix',
		Dhi='matrix',
		Dst='matrix',
		SXhi ='matrix',
		Aex2hi='matrix', 
		Ahi2st ='matrix',
		OXst ='matrix',
		version='ANY'
	),
	prototype=prototype(
		version=new('version')
	),
	package='arf3DS4'
	
)

## arf/FSL functional_data class
setClass(
	Class='functional',
	representation=representation(
		fullpath='character',			#fullpath of functionalvolume
		functionaldata='character',		#functionalvolumefilename
		filename='character',			#full filename of functional.Rda
		linkedfiles='character',		#link to betafilename
		timings='numeric',				#stimulus timings for condition/contrast
		version='ANY'
	),
	prototype=prototype(
		version=new('version')
	),
	package='arf3DS4'
)


## arf experiment class (containing, rootdirectories, condition numbers etc. and preferred sequences)
setClass(
	Class='experiment',
	contains='settings',
	representation=representation(
		path='character', 			#rootdirectory of experiment
		name='character',			#experiment name
		subject.num='numeric',		#number of subjects
		subject.names='character',	#names of subjects
		condition.num='numeric',	#number of conditions
		condition.names='character' #names of condition directories
	),
	package='arf3DS4'
)

## arf (general) settings class (containing general settings on NLM, starting values and range checks)
setClass(
	Class='options',
	representation=representation(
		nlm.gradtol='numeric',		#NLM gradient tolerance
		nlm.steptol='numeric',		#NLM stepsize tolerance
		opt.method='character',		#optim method
		opt.lower='numeric',		#optim lowerbound
		opt.upper='numeric',		#optim upperbound
		min.analyticalgrad='logical', #use analytical gradient
		min.iterlim='numeric',		#minimization iteration limit
		min.routine='character',	#which routine is used
		min.boundlim='numeric',		#number of iterations before exiting with persistent bound error	
		start.method='character',	#which method of determining starting values is used ('rect','load','none')
		start.maxfac='numeric',		#fallOff factor in the determination of region width
		start.vector='numeric',		#vector containing startingvalues (if !start.method=='fixed' only t5 and t6 are used) vector is recycled for regions.
		chk.method='character',		#which method is used to check the range of parameter values
		chk.range='numeric',		#vector containing ranges for each parameter (vector is recycled for regions)		
		sw.type='character',		#method to use with Residuals ('diag','full')
		output.mode='character',	#output mode
		version='ANY'


	),
	prototype=prototype(
		nlm.gradtol=1e-6,		#NLM gradient tolerance
		nlm.steptol=1e-3,		#NLM stepsize tolerance
		opt.method='L-BFGS-B',	#optim method
		opt.lower=c(rep(0,3),rep(0.1,3),rep(-.90,3),-Inf),	#L-BFGS-U lower bound
		opt.upper=c(rep(Inf,3),rep(1,3),rep(.90,3),Inf),	#L-BFGS-U upper bound
		min.analyticalgrad=T, 	#use analytical gradient
		min.iterlim=2000,		#iteration limit
		min.boundlim=50,		#boundary iterations limit
		min.routine=c('optim','vpv'),	#which routine
		start.method='use',	
		start.maxfac=1,
		start.vector=c(0,0,0,0,0,0,.05,.05,.05,10),
		chk.method='imagedim',
		chk.range=c(0,0,0,0,0,0,-.95,-.95,-.95,-1e+64,0,0,0,0,0,0,.95,.95,.95,1e+64),
		sw.type='diag',
		output.mode=c('none'),
		version=new('version')
	),
	package='arf3DS4'
)

## make nifti.fileinfo class
setClass(
	Class='nifti.fileinfo',
	representation=representation(
		fullpath='character',		#Full path of nifti file
		filename='character',		#name of nifti file
		filetype='character',		#type (nifti or analyze)
		extension='character',		#extension (.nii, .hdr/.img)
		gzipped='logical',			#gzipped or not
		endian='character',			#endian
		version='ANY'
	),
	prototype=prototype(
		gzipped=TRUE,
		extension='nii',
		filetype='nifti',
		endian=.Platform$endian,
		version=new('version')
	),
	package='arf3DS4'
)

## nifti.header class (inherits nifti.fileinfo)
setClass(
	Class='nifti.header',
	contains='nifti.fileinfo',
	representation=representation(
		sizeof_hdr = 'numeric',			#size of header (must be 348))
		data_type =  'character',		
		db_name =  'character',			
		extents = 'numeric',			
		session_error = 'numeric',
		regular = 'character',
		dim_info = 'character',
		dims = 'numeric',				#dimensions (num of dim, dimx, dimy,...)
		intent_p1 = 'numeric',
		intent_p2 = 'numeric',
		intent_p3 = 'numeric',
		intent_code = 'numeric',
		datatype = 'numeric',			#storage data type
		bitpix = 'numeric',				#bits per pixel
		slice_start = 'numeric',		
		pixdim = 'numeric',				#voxel dimensions
		vox_offset = 'numeric',			#offset of data in .nii file
		scl_slope = 'numeric',
		scl_inter = 'numeric',
		slice_end = 'numeric',
		slice_code = 'numeric',
		xyzt_units = 'character',
		cal_max = 'numeric',
		cal_min = 'numeric',
		slice_duration = 'numeric',
		toffset = 'numeric',
		glmax = 'numeric',
		glmin = 'numeric',
		descrip = 'character',
		aux_file = 'character',
		qform_code = 'numeric',
		sform_code = 'numeric',
		quatern_b = 'numeric',
		quatern_c = 'numeric',
		quatern_d = 'numeric',
		qoffset_x = 'numeric',
		qoffset_y = 'numeric',
		qoffset_z = 'numeric',
		srow_x = 'numeric',
		srow_y = 'numeric',
		srow_z = 'numeric',
		intent_name = 'character',		#meaning of data
		magic = 'character',			#magicstring
		data.type='character',			#type of data
		data.signed='logical'			#signed data
	),
	prototype=prototype(
		sizeof_hdr = 348,			#size of header (must be 348))
		data_type =  '',		
		db_name =  '',			
		extents = 0,			
		session_error = 0,
		regular = 'r',
		dim_info = '',
		dims = c(3,0,0,0,0,0,0,0),	#dimensions (num of dim, dimx, dimy,...)
		intent_p1 = 0,
		intent_p2 = 0,
		intent_p3 = 0,
		intent_code = 5,
		datatype = 16,			#storage data type
		bitpix = 32,			#bits per pixel
		slice_start = 0,		
		pixdim = c(1,1,1,1,1,0,0,0),
		vox_offset = 0,			#offset of data in .nii file
		scl_slope = 1,
		scl_inter = 0,
		slice_end = 0,
		slice_code = 0,
		xyzt_units = '\n',
		cal_max = 0,
		cal_min = 0,
		slice_duration = 0,
		toffset = 0,
		glmax = 0,
		glmin = 0,
		descrip = 'ARF3DS4',
		aux_file = '',
		qform_code = 1,
		sform_code = 1,
		quatern_b = 0,
		quatern_c = 0,
		quatern_d = 0,
		qoffset_x = 0,
		qoffset_y = 0,
		qoffset_z = 0,
		srow_x = c(1,0,0,0),
		srow_y = c(0,1,0,0),
		srow_z = c(0,0,1,0),
		intent_name = 'estimates',	#meaning of data
		magic = 'n+1',					#magicstring
		data.type = 'double',			#type of data
		data.signed = TRUE				#signed data
	),
	package='arf3DS4'
)

## fmri data class (inherits nifti.header)
setClass(
	Class='fmri.data',
	contains='nifti.header',
	representation=representation(
		datavec='numeric'
	),
	package='arf3DS4'
)

## wald statistics class
setClass(
	Class='wald',
	representation=representation(
		consts='matrix',	#constants for tests
		stats='matrix',		#statistic
		df1='numeric',		#df1
		df2='numeric',		#df2
		pvalues='matrix'	#pvalues
	),
	package='arf3DS4'
)

## arf data class (containing info on the locations of the data and weightfiles, the dimensions, number of runs, and relevant nifti parameters.)
setClass(
	Class='data',
	representation=representation(
		name='character',			#indicator of data (subject,condition etc.)
		fullpath='character',		#fullpath of files
		betafiles='character',		#vector of char containing run datafiles
		weightfiles='character',	#vector of char containing weight datafiles
		avgdatfile='character',		#filename of average data
		avgWfile='character',		#filename of average weights
		avgtstatFile='character',	#filename of avgtstatFile
		n='numeric',				#number of voxels of the images
		mask='numeric',				#vector of length n defining a mask
		ss='numeric',				#sums of squares of the data
		regDir='character',			#directory of registration dirs
		regRda='character',			#filename of registration file 
		funcDir='character',		#directory of functional volume
		funcRda='character',		#filename of functional file
		runs='numeric',				#number of runs
		dataHeader='ANY',
		version='ANY'
	),
	prototype=prototype(
		version=new('version'),
		avgdatfile='',
		avgWfile='',
		avgtstatFile='',
		dataHeader=new('nifti.header')
	),
	package='arf3DS4'
)


## arf model class (containing information on the fitted arf model, it extends the data class)
setClass(
	Class='model',
	contains='data',
	representation=representation(
		modelname='character',		#name of the model (default is region_n)
		modelpath='character',		#path to the model directory
		modeldatapath='character',  #path to the modeldata dir
		residualFile='character',	#Residual Filename
		derivativeFile='character',	#Derivative Filename
		weightFile='character',		#weightfilename
		modelDataFile='character',	#niftifilename
		fullmodelDataFile='character',#full niftifilename
		modelFile='character',		#modelFilename
		optionsFile='character',	#optionsFilename
		startFile='character',		#startvalueFilename
		logFile='character',		#logFileName
		convergence='character', 	#convergence information
		iterates='numeric',			#number of iterations
		minimum='numeric',			#minimum of objective function
		estimates='numeric',		#vector of parameter estimates (t1r1..t6r1,t1r2..t6r2,t1rR..t6rR)
		gradient='numeric',			#gradient of solution
		hessian='matrix',			#hessian matrix
		params='numeric',			#number of parameters in the model
		modeltype='character',		#type of model (currently simple or gauss)
		sandwichmethod='character',	#sandwichmethod
		varcov='matrix',			#variance covariance matrix (full form)
		warnings='character',		#warnings (pos def var/covar etc.)
		fit='matrix',				#fit value 
		wald='ANY',					#object of class 'wald'
		regions='numeric',			#number of fitted regions
		startval='numeric',			#vector of starting values
		proctime='matrix',			#processing time
		valid='logical'				#is model valid (converged and no warnings)
	),
	prototype=prototype(
			valid=FALSE,
			proctime=matrix(0,1,2,dimnames=list(c(''),c('mintime','swtime'))),
			fit=matrix(0,1,2,dimnames=list(c(''),c('BIC','RMSEA'))),
			modeltype='gauss',
			params=10,
			wald=new('wald')
	),
	package='arf3DS4'
)

## arf sequence class (containing info (fit, valid) on a sequence of models)
setClass(
		Class='sequence',
		representation=representation(
				current='numeric',			#current minimum
				regions='numeric',			#vector of regions to fit (can be sequential or any other combination)
				mnames='character',			#vector of names of models
				fit='numeric',				#vector of fit measures (to evaluate best fit)
				minimum='numeric',			#vector of minima for the fit functions 
				best='logical',				#which model has the best fit
				valid='logical'				#vector of validity of solutions (all estimates and variances ok)
		),
		package='arf3DS4'
)

## arf modelnames class, contains modelnames for a subject/condition
setClass(
		Class='mnames',
		representation=representation(
				experiment='ANY',
				subject='character',
				condition='character',
				mnames='character'
		),
		package='arf3DS4'
)

## arf correlation class, contains correlations and partial correlations
setClass(
		Class='arfcorrelation',
		representation=representation(
				timebyreg='matrix',
				corr='matrix',
				corr.pval='matrix',
				pacorr='matrix',
				pacorr.pval='matrix',
				num.corr='numeric'
				
		),
		package='arf3DS4'
)
