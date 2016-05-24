#############################################
# arf3DS4 S4 DISPLAY FUNCTIONS				#
# Wouter D. Weeda							#
# University of Amsterdam					#
#############################################


makeExample <- function(path) {
	
	if(missing(path)) path=system.file('data',package='arf3DS4')
	
	#Definitions
	expname = 'example-experiment'
	subject = 'wouter'
	condition = 'A'
	p = .Platform$file.sep
	
	#load example datafiles
	example = get('.example',pos=-1)
		
	#makeExperimentDirs
	makeExpDirs(paste(path,p,sep=''),expname,subject,condition)
	
	#write Nifti objects to path
	fullpath = paste(path,p,expname,p,'subjects',p,subject,p,'conditions',p,condition,p,'data',p,'beta',sep='')
	example$tstat1@fullpath = example$tstat2@fullpath = fullpath
	writeData(example$tstat1)
	writeData(example$tstat2)
	
	fullpath = paste(path,p,expname,p,'subjects',p,subject,p,'conditions',p,condition,p,'data',p,'funcs',sep='')
	example$sevents@fullpath = fullpath
	writeData(example$sevents)
	
	#load and set experiment
	ex = loadExp(paste(path,p,expname,sep=''),'set')
	createAllAverages(ex)
	
	#assign .experiment to the .arfInternal  
	assign('.experiment',ex,envir=.arfInternal)
	
	return(invisible(ex))
	
}

