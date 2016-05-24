#############################################
# arf3DS4 S4 EXPERIMENT FUNCTIONS			#
# Wouter D. Weeda							#
# University of Amsterdam					#
#############################################

#[CONTAINS]
#checkExp
#makeExpDirs		[user]
#chngRootExp
#setExp			
#getExp				[user]
#loadExp			[user]
#checkVersion
#getFSLdata


checkExp <- 
function(experiment) 
# checkExp determines if the experiment structure is valid
{
	
	#set separator and 'all is well' flag
	sp <- .Platform$file.sep
	allIsWell=TRUE
	
	#check if input is of class 'experiment'
	if(class(experiment)!='experiment') stop('Input must be of class \'experiment\'')
	
	subd <- paste(.experiment.path(experiment),sp,.experiment.subjectDir(experiment),sep='')
	
	for(sdirs in 1:.experiment.subject.num(experiment)) {
		
		sn <- paste(subd,sp,.experiment.subjectPrefix(experiment),.experiment.subject.names(experiment)[sdirs],sep='')
		subc <- paste(sn,sp,.experiment.conditionDir(experiment),sep='')
		
		for(cdirs in 1:.experiment.condition.num(experiment )) {
			
			cn <- paste(subc,sp,.experiment.conditionPrefix(experiment),.experiment.condition.names(experiment)[cdirs],sep='')
			
			#check if data,weights,avg and models directories exist.
			if(!file.exists(paste(cn,sp,.experiment.modelDir(experiment),sep=''))) {warning(paste('In:',cn,' modelDir does not exits or does not match settings',sep=''));allIsWell=F}
			if(!file.exists(paste(cn,sp,.experiment.dataDir(experiment),sep=''))) {warning(paste('In:',cn,' dataDir does not exits or does not match settings',sep=''));allIsWell=F}
			if(!file.exists(paste(cn,sp,.experiment.dataDir(experiment),sp,.experiment.betaDir(experiment),sep=''))) {warning(paste('In:',cn,sp,.experiment.dataDir(experiment),' betaDir does not exits or does not match settings',sep=''));allIsWell=F}
			if(!file.exists(paste(cn,sp,.experiment.dataDir(experiment),sp,.experiment.weightsDir(experiment),sep=''))) {warning(paste('In:',cn,sp,.experiment.dataDir(experiment),' weightsDir does not exits or does not match settings',sep=''));allIsWell=F}
			if(!file.exists(paste(cn,sp,.experiment.dataDir(experiment),sp,.experiment.avgDir(experiment),sep=''))) {warning(paste('In:',cn,sp,.experiment.dataDir(experiment),' avgDir does not exits or does not match settings',sep=''));allIsWell=F}
			
		}
	}
	
	#return FALSE if not correct, return TRUE if all is well
	return(invisible(allIsWell))
}


makeExpDirs <- 
function(path=getwd(),name='default_experiment',subjectind=1,conditionind=1,settings=new('settings'))
# makeExpDirs creates a directory structure given the path, number of subjects, conditions and a settings object
# by default creates a structure in the working directory with 1 subject 1 condition and default settings
{
	#set separator
	sp <- .Platform$file.sep
	
	#check if path exists
	if(path=='') path=paste(getwd(),sp,name,sep='') else path=paste(path,sp,name,sep='')
	
	if(file.exists(path)) {
		
		valAns=F
		while(!valAns) {
			answ = readline('Directory already exists, overwrite? [y/n] ')
			if(tolower(answ)=='y' | tolower(answ)=='n') valAns=T
		}
		
		if(tolower(answ)=='y') file.remove(list.files(path,'.Rda',full.names=T))	
		if(tolower(answ)=='n') stop('User stop: not overwriting directory.')
				
	}
	
	if(!file.exists(path)) dir.create(path,recursive=T)
	
	#set path separator at end
	path <- file.path(dirname(path),basename(path),'')
	
	#create new experimentclass
	experiment <- new('experiment',settings)
	
	#set the experiment values
	.experiment.path(experiment)=path
	.experiment.name(experiment)=name
	.experiment.subject.num(experiment)=length(subjectind)
	.experiment.subject.names(experiment)=as.character(subjectind)
	.experiment.condition.num(experiment)=length(conditionind)
	.experiment.condition.names(experiment)=as.character(conditionind)
		
	#save experiment Rda
	save(experiment,file=paste(.experiment.path(experiment),sp,.settings.expRda(settings),sep=''))
	
	#create subjects directory
	subd <- paste(path,sp,.settings.subjectDir(settings),sep='')
	dir.create(subd,showWarnings=F)
	
	for(sdirs in 1:.experiment.subject.num(experiment)) {
		
		#create individual subject dirs
		sn <- paste(subd,sp,.settings.subjectPrefix(settings),.experiment.subject.names(experiment)[sdirs],sep='')
		dir.create(sn,showWarnings=F)
		
		#create conditions dir
		subc <- paste(sn,sp,.settings.conditionDir(settings),sep='')
		dir.create(subc,showWarnings=F)		

		#create funcDir
		dir.create(paste(sn,sp,.settings.funcDir(settings),sep=''),showWarnings=F)
		
		for(cdirs in 1:.experiment.condition.num(experiment)) {
			
			#create individual condition dirs
			cn <- paste(subc,sp,.settings.conditionPrefix(settings),.experiment.condition.names(experiment)[cdirs],sep='')
			dir.create(cn,showWarnings=F)
			
			#models dir and modelnamesfile
			dir.create(paste(cn,sp,.settings.modelDir(settings),sep=''),showWarnings=F)
			mnames=''
			save(mnames,file=paste(cn,sp,.settings.modelDir(settings),sp,.settings.modelnamesRda(settings),sep=''),showWarnings=T)
			
			#data,weights,avg
			dn <- paste(cn,sp,.settings.dataDir(settings),sep='')
			dir.create(dn,showWarnings=F)
						
			dir.create(paste(dn,sp,.settings.betaDir(settings),sep=''),showWarnings=F)
			dir.create(paste(dn,sp,.settings.weightsDir(settings),sep=''),showWarnings=F)
			dir.create(paste(dn,sp,.settings.avgDir(settings),sep=''),showWarnings=F)
			dir.create(paste(dn,sp,.settings.regDir(settings),sep=''),showWarnings=F)
			dir.create(paste(dn,sp,.settings.funcDir(settings),sep=''),showWarnings=F)
			
		
		}
	}
	
	#save experiment
	save(experiment,file=paste(.experiment.path(experiment),sp,.experiment.expRda(experiment),sep=''))

	#final check
	if(!checkExp(experiment)) {
		warning('Experiment did not pass checks! Check warnings!')
	} else return(invisible(experiment))
	
}


chngRootExp <- 
function(path=getwd(),quiet=F) 
# resets and checks the experiment path
{
	
	#set separator
	sp <- .Platform$file.sep
	
	#check if one file is in the experiment root
	if(length(list.files(path,'.Rda',full.names=T))!=1) stop('No experiment rda file found or multiple rda files found.')
	experiment <- loadRda(list.files(path,'.Rda',full.names=T))
	
	#set the correct path
	path <- path.expand(path)
	path <- file.path(dirname(path),basename(path),'')
	
	.experiment.path(experiment) <- path
	
	#check the experiment dirs, if good save and exit. if not good stop
	if(checkExp(experiment)) {
		save(experiment,file=paste(.experiment.path(experiment),sp,.experiment.expRda(experiment),sep=''))
		if(!quiet) cat('Experiment path changed. Experiment saved to',paste(.experiment.path(experiment),sp,.experiment.expRda(experiment),sep=''),'\n')
	} else {
		stop('Experiment structure not valid,check warnings.')
	}
	
	return(invisible(experiment))
}


setExp <- 
function(path=getwd(),tempsub=1,tempcond=1,auto=TRUE,createWeights=TRUE,overwrite=F,load=T) 
#setExp creates an experiment-class object based on existing directories
{
	#set separator
	sp <- .Platform$file.sep
	
	#check if only dir else pre-append working directory
	if(length(grep(sp,path))==0) path <- paste(getwd(),sp,path,sep='')
	
	#set path separator at end
	path <- file.path(dirname(path),basename(path),'')
	
	#set expname
	expname <- strsplit(path,.Platform$file.sep)[[1]]
	expname <- expname[length(expname)]
	
	#if experiment file exists open it, else create a new one
	exp <-  list.files(path,pattern='Rda',full.names=T)
	if(length(exp)==1 & overwrite==FALSE) settings <- loadRda(exp) else settings <- new('settings') 
		
	#create new experimentclass
	experiment <- new('experiment',settings)
	.experiment.name(experiment) <- expname
	.experiment.path(experiment) <- path
	
	#cat('[',toupper(expname),']\n')
	#cat(' Experiment root:',path,'\n')
	
	#determine directories subjects
	whichdirs <- file.info(list.files(.experiment.path(experiment),full.names=T))$isdir
	fileList <- list.files(.experiment.path(experiment),full.names=F)[whichdirs]
	if(length(fileList)!=1) stop('Multiple subject directories found in rootdir') 
		
	#set subjectsdirname
	.experiment.subjectDir(experiment) <- fileList
	#cat('  \\Subjects directory:',fileList,'\n')
	
	#get subjects data
	subd <- paste(.experiment.path(experiment),sp,.experiment.subjectDir(experiment),sep='')
	whichdirs <- file.info(list.files(subd,full.names=T))$isdir
	fileList <- list.files(subd,full.names=F)[whichdirs]
	.experiment.subject.num(experiment) <- length(fileList)
	.experiment.subject.names(experiment) <- fileList
	
	#go to the first subject		
	sn <- paste(subd,sp,.experiment.subjectPrefix(experiment),.experiment.subject.names(experiment)[tempsub],sep='')
	whichdirs <- file.info(list.files(sn,full.names=T))$isdir
	fileList <- list.files(sn,full.names=F)[whichdirs]
	
	#get conditions dirname
	if(length(fileList)!=1) fileList <- .experiment.conditionDir(experiment) 
	.experiment.conditionDir(experiment) <- fileList
	#cat('   \\Condition directory:',fileList,'\n')
			
	#get condition data
	subc <- paste(sn,sp,.experiment.conditionDir(experiment),sep='')
	whichdirs <- which(file.info(list.files(subc,full.names=T))$isdir)
	fileList <- list.files(subc,full.names=F)[whichdirs]
	.experiment.condition.num(experiment) <- length(fileList)
	.experiment.condition.names(experiment) <- fileList
	
	#go to first condition 	
	cn <- paste(subc,sp,.experiment.conditionPrefix(experiment),.experiment.condition.names(experiment)[tempcond],sep='')
	whichdirs <- which(file.info(list.files(cn,full.names=T))$isdir)
	fileList <- list.files(cn,full.names=F)[whichdirs]
		
	#make uniform weights when no weights exist
	if(createWeights) makeWeights(experiment)
	
	#setAllObjects
	setAllObjects(experiment,overwrite)

	#create functional files
	setFuncFiles(func_data='',experiment=experiment)
	
	#check the experiment dirs, if good save and exit. if not good stop
	if(checkExp(experiment)) {
		save(experiment,file=paste(.experiment.path(experiment),.settings.expRda(settings),sep=''))
		cat('Experiment correctly set. Experiment saved to',paste(.experiment.path(experiment),.settings.expRda(settings),sep=''),'\n')
		if(load) loadExp(paste(.experiment.path(experiment),sep=''))
				
		assign('.experiment',experiment,envir=.arfInternal)
		save(experiment,file=paste(.experiment.path(experiment),sp,.experiment.expRda(experiment),sep=''))

	} else {
		stop('Experiment structure not valid,check warnings.')
	}
	
	return(invisible(experiment))
	
}

getExp <-
function()
#return the experiment-object
{
	experiment <- try(get('.experiment',envir=.arfInternal),silent=T)
	if(attr(experiment,'class')=='try-error') stop('Experiment not loaded. Run loadExp first.')
	return(experiment)
}


loadExp <- 
function(path=getwd(),method=c('fast','set','rda'))
#loadExp loads an experiment 
{
	
	method = match.arg(method)
	
	if(method=='fast') fast=TRUE
	if(method=='set') {
		fast=FALSE
		overwrite=TRUE
		set=TRUE
	}
	if(method=='rda') {
		fast=FALSE
		overwrite=TRUE
		set=FALSE
	}
	#set separator
	sp <- .Platform$file.sep
	
	allIsWell = TRUE
	
	#check if only dir else pre-append working directory
	if(length(grep(sp,path))==0) path <- paste(getwd(),sp,path,sep='')
	
	#set filename and load experiment
	filename <- list.files(path,'.Rda',full.names=T)
	if(length(filename)!=1) stop(paste('No experiment rda file found or multiple rda files found in ',path,sep=''))
	.experiment <- experiment <- loadRda(filename)
	
	if(!fast) {
		#remove filename to obtain root-path
		fn <- strsplit(filename,.Platform$file.sep)[[1]]
		path <- sub(fn[length(fn)],'',filename)
		path <- path.expand(path)
		path <- file.path(dirname(path),basename(path),'')
		
		if(!set)  {
			#change root for the experiment-file to the current root
			.experiment <- experiment <- chngRootExp(path,quiet=T)
			#set and check all objects based on subjects/condition info and settings
			allIsWell <- setAllObjects(experiment,overwrite=overwrite) 
		} else .experiment <- experiment <- setExp(path,1,1,TRUE,TRUE,TRUE,FALSE)
		
	}  
	
	#save experiments
	assign('.experiment',.experiment,envir=.arfInternal)
	save(experiment,file=paste(.experiment.path(experiment),sp,.experiment.expRda(experiment),sep=''))
		
	#return loaded info 
	if(allIsWell) cat('Loaded experiment ',.experiment.name(experiment),' (version ',.version.version(.experiment.version(experiment)),'.',.version.build(.experiment.version(experiment)),'-',.version.update(.experiment.version(experiment)),')\n',sep='') else cat('Loaded experiment ',.experiment.name(experiment),' (version',.version.version(.experiment.version(experiment)),'.',.version.update(.experiment.version(experiment)),'-',.version.build(.experiment.version(experiment)),') with warnings!\n',sep='')
	
	return(invisible(experiment))
		
}

checkVersion <- 
function(curversion,version=1,build=0,update=0) 
#checks the version of the ARF code (check if curversion larger than given version)
{
	allIsWell=FALSE
		
	if(version < .version.version(curversion)) {
		allIsWell=TRUE
	} else {
		if(version == .version.version(curversion)) {
			if(build < .version.build(curversion)) {
				allIsWell=TRUE
			} else {
				if(build == .version.build(curversion)) {
					if(update < .version.update(curversion)) {
						allIsWell=TRUE
					}						
				}
			}	
		}
	}
	
	return(allIsWell)
}


getFSLdata <- function(fsldir=getwd(),subjectlist,contrastnums,expname='FSLtoARF',expdir=getwd(),settings=new('settings'),quiet=F,doReg=T,doFunc=F) 
{
	#set separator
	sp = .Platform$file.sep
	
	cont = TRUE
	
	cat('getFSLdata featdir',fsldir,'\n')
	cat('getFSL for subjects:\n')
	for(i in 1:length(subjectlist)) cat(' >',subjectlist[i],'\n')
	
	#standard FSL files DO NOT CHANGE!!!!
	functional = c('filtered_func_data.nii.gz')
	registration = c('example_func2standard.mat','example_func2highres.mat','highres2standard.mat','example_func.nii.gz','highres.nii.gz','standard.nii.gz')
	cope = 'cope'
	varcope = 'varcope'
	
	subjects <- length(subjectlist)
	subinf <- vector(subjects,mode='list')
	
	for(i in 1:subjects) if(!file.exists(paste(fsldir,sp,subjectlist[i],'.feat',sep=''))) {warning(paste('subject',subjectlist[i],'does not exist'));cont=FALSE}
	
	if(!cont) stop('Cannot run getFSLdata, reason: one or more subjects does not exist!')
	
	cat('fetching fsl filenames...')
	#get all subject files and info for copying
	for(subj in 1:subjects) {
			
		dirs = list.files(fsldir,subjectlist[subj],full.names=T)
		runs = length(dirs)
				
		func_file=reg_file=cope_file=vector(runs,mode='list')
				
		for(tri in 1:runs) {
		
			#get functional data
			func_file[[tri]] = list.files(dirs[tri],functional,full.names=T)
			
			#get registration data
			rdir = paste(dirs[tri],sp,'reg',sep='')
			reg_file[[tri]] = list(ex2stand=list.files(rdir,registration[1],full.names=T),ex2hi=list.files(rdir,registration[2],full.names=T),hi2st=list.files(rdir,registration[3],full.names=T),ex=list.files(rdir,registration[4],full.names=T),hi=list.files(rdir,registration[5],full.names=T),st=list.files(rdir,registration[6],full.names=T))
			
			reg_file[[tri]]$hi = reg_file[[tri]]$hi[-1]
			reg_file[[tri]]$st = reg_file[[tri]]$st[-c(1,2)]
			
			
			#get copes and varcopes
			sdir = paste(dirs[tri],sp,'stats',sep='')
			all_copes = list.files(sdir,cope,full.names=T)
			
			use_cope=use_varcope=numeric(0)
			
			for(i in contrastnums){
				use_cope = c(use_cope,grep(paste(cope,i,sep=''),all_copes))
				use_varcope = c(use_varcope,grep(paste(varcope,i,sep=''),all_copes))
			}
			
			use = unique(c(use_cope,use_varcope))
			if(length(use)>0) all_copes = all_copes[use]
			vcs = grep(varcope,all_copes)
			
			#get varcopes out of it and put all in list 
			cope_file[[tri]] = list(copes=all_copes[-vcs],varcopes=all_copes[vcs])
		}
		
		subinf[[subj]] = list(func_file=func_file,reg_file=reg_file,cope_file=cope_file,runs=runs,contrasts=length(cope_file[[1]]$copes),name=subjectlist[subj])
	
	}
	cat('ok\n')
	
	if(length(contrastnums)==0) contrastnums = 1:subinf[[1]]$contrasts
	
	
	cat('making experiment dirs...')
	#make new experiment
	experiment <- makeExpDirs(path=expdir,name=expname,subjectind=subjectlist,conditionind=paste('contrast',contrastnums,sep=''),settings=settings)
	subd <- paste(.experiment.path(experiment),sp,.settings.subjectDir(settings),sep='')
	cat('ok\n')
	
	cat('copying images [status]:\n')
	
	#fill the experiment with files
	for(sdirs in 1:length(subinf)) {
		
		#get to subjectsdir
		sn <- paste(subd,sp,.settings.subjectPrefix(settings),subinf[[sdirs]]$name,sep='')
		subc <- paste(sn,sp,.settings.conditionDir(settings),sep='')
			
		for(cdirs in 1:.experiment.condition.num(experiment)) {
			
			#get to conditiondir and datadir
			cn <- paste(subc,sp,.settings.conditionPrefix(settings),.experiment.condition.names(experiment)[cdirs],sep='')
			dn <- paste(cn,sp,.settings.dataDir(settings),sep='')
			
			#get the runs for the conditions 
			for(i in 1:subinf[[sdirs]]$runs) {
				copefile = avgfile = readData(subinf[[sdirs]]$cope_file[[i]]$copes[cdirs])
				varcopefile = readData(subinf[[sdirs]]$cope_file[[i]]$varcopes[cdirs])
				tstat = .fmri.data.datavec(copefile)/sqrt(.fmri.data.datavec(varcopefile))
				tstat[is.nan(tstat)]=0
				.fmri.data.fullpath(avgfile) = paste(dn,sp,.settings.betaDir(settings),sep='')
				.fmri.data.filename(avgfile) = paste('arf_tstat',i,sep='')
				.fmri.data.descrip(avgfile) = 'arf t_stat image'
				writeData(avgfile,tstat)
				if(!quiet) cat('',subinf[[sdirs]]$cope_file[[i]]$copes[cdirs],'>',paste(dn,sp,.settings.betaDir(settings),sp,'arf_tstat',i,'\n',sep=''))
			}
			
			cropVolumeAuto(paste(dn,sp,.settings.betaDir(settings),sep=''),quiet=T)
		}
	}
	cat('> ready\n')
	
	
	#set experiment and create weights
	experiment <- setExp(.experiment.path(experiment),1,1,T,T,T,F)
	
	if(doReg) {
		cat('setting registration files [status]:\n')
		#create and fill regs and funcs
		for(sdirs in 1:length(subinf)) {
			
			#get to subjectsdir
			sn <- paste(subd,sp,.settings.subjectPrefix(settings),subinf[[sdirs]]$name,sep='')
			subc <- paste(sn,sp,.settings.conditionDir(settings),sep='')
			
			for(cdirs in 1:.experiment.condition.num(experiment)) {
				
				#get to conditiondir and datadir
				cn <- paste(subc,sp,.settings.conditionPrefix(settings),.experiment.condition.names(experiment)[cdirs],sep='')
				dn <- paste(cn,sp,.settings.dataDir(settings),sep='')
				
				arfdata <- loadRda(paste(dn,sp,.settings.dataRda(settings),sep=''))
				createRegs(arfdata)
			
				
				filelist_reg = list.files(paste(dn,sp,.settings.regDir(settings),sep=''))
					
				for(i in 1:subinf[[sdirs]]$runs) {
					
					registration <- loadRda(paste(dn,sp,.settings.regDir(settings),sp,filelist_reg[i],sp,.settings.regRda(settings),sep=''))
								
					for(f in 1:length(subinf[[sdirs]]$reg_file[[i]])) {
						fn = as.character(subinf[[sdirs]]$reg_file[[i]][f])
						fns = strsplit(fn,sp)[[1]]
						file.copy(fn,paste(.registration.fullpath(registration),sp,fns[length(fns)],sep=''))
						if(!quiet) cat('',fn,'>',paste(.registration.fullpath(registration),sp,fns[length(fns)],sep=''),'\n')
						
					}
					
					registration <- setRegFiles(registration)
					registration <- setRegParams(registration)
							
				}
			}
		}
		cat('> ready\n')
	} else cat('registration files not requested.\n')
	
	if(doFunc) {
		
		cat('setting functional files [status]:\n')
		#create and fill regs and funcs
		for(sdirs in 1:length(subinf)) {
			
			#get to subjectsdir
			sn <- paste(subd,sp,.settings.subjectPrefix(settings),subinf[[sdirs]]$name,sep='')
			subc <- paste(sn,sp,.settings.conditionDir(settings),sep='')
			
			for(cdirs in 1:.experiment.condition.num(experiment)) {
				cn <- paste(subc,sp,.settings.conditionPrefix(settings),.experiment.condition.names(experiment)[cdirs],sep='')
				dn <- paste(cn,sp,.settings.dataDir(settings),sep='')
				arfdata <- loadRda(paste(dn,sp,.settings.dataRda(settings),sep=''))
				createFuncs(arfdata)
			}
			
			filelist_func = list.files(paste(dn,sp,.settings.funcDir(settings),sep=''))
			
			for(i in 1:subinf[[sdirs]]$runs) {
				functional <- loadRda(paste(dn,sp,.settings.funcDir(settings),sp,filelist_func[i],sp,.settings.funcRda(settings),sep=''))
				fn = as.character(subinf[[sdirs]]$func_file[[i]][1])
				fns = strsplit(fn,sp)[[1]]
				funcdatapath <- paste(sn,sp,.settings.funcDir(settings),sp,filelist_func[i],sep='')
				if(!file.exists(funcdatapath)) dir.create(funcdatapath)
				file.copy(fn,paste(funcdatapath,sp,fns[length(fns)],sep=''))
				if(!quiet) cat('',fn,'>',paste(funcdatapath,sp,fns[length(fns)],sep=''),'\n')

				
			}
		}
				
		setFuncFiles(experiment)
		cat('> ready\n')
		
		
	} else cat('functional files not requested.\n')
	
	#set experiment and create weights
	experiment <- setExp(.experiment.path(experiment),1,1,T,T,T,T)
	
	return(invisible(experiment))
}
