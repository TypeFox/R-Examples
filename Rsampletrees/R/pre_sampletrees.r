## This file contains the functions that setup the parameter files
## for a sampletrees run

# Create a default object containing the arguments/settings to pass to sampletrees
defaultArgs=function()
{
	args=list(RunName="Run",  Seed=NA, DataType="h", DataFile=NA, LocationFile=NA,WeightFile=NA,
			FocalPoint=NA, ChainLength=1000, BurnIn=0, Thinning=1,
			InitialTheta=1,MinTheta=0.0001,MaxTheta=10,
			InitialRho=0.0004,ScaleRho=0.1,ShapeRho=1,
			InitialTreeFile=NA, 
			RandomTree=FALSE, HaploFreqFile=NA, InitialHaploFile=NA,  
			HaploListFile=NA, clean=FALSE)
	
	class(args)="pars"
  	return(args)	
}


# Set up a new object containing the arguments/settings and optionally 
# pass some of the settings for initialization
newArgs=function(...)
{
	args=defaultArgs()
	otherargs=list(...)

	if (length(otherargs)!=0){
		args=changeArgs(args,...)
	}
	
	return(args)	
	
}

# This function is used to change any of the options in the
# settings object. Because the settings object is a list, 
# one could just change the values of the list. But it is 
# safer to use this function because it will test that
# the option that is being changed is actually valid. 
changeArgs=function(object, ...) 
{
	UseMethod("changeArgs", object)
}

changeArgs.default=function(object, ...)
{
	return(paste("changeArgs method for class",class(object),"is not defined"))
}

changeArgs.pars=function(object,...)
{
	args=object
	otherargs=list(...)

	if (length(otherargs)!=0){
		
		inboth=intersect(names(args),names(otherargs))
		args[inboth]=otherargs[inboth]
		notfound <- setdiff(names(otherargs),names(args))
		
		if (length(notfound)!=0) {
			warning(paste("The following are not valid settings: ", notfound))
		}
		
	}
	
	return(args)	
}


# This function is used to set the weights for each of the proposal
# types in sampletrees. These are written to a file and the 
# name of the file is set in the settings object
setWeights=function(args, WeightFile=NULL, weightmat=NULL){
	
	if (is.null(WeightFile)){
		# Use a default filename
		WeightFile="weights"
	}
	args$WeightFile=WeightFile
	
	if (is.null(weightmat)){
		# Use default weightings
		
		if (args$DataType=="h"){
			weightmat=cbind(c(0.1,0.1,0.4,0.2,0.2),1:5)
		} else if (args$DataType=="g"){
			weightmat=cbind(c(0.1,0.1,0.35,0.1,0.1,0.125,0.125),1:7)
		} else {
			stop("Datatype must be 'g' or 'h' to set default weight values\n")
		}
	} else {
		# Set weights using weightmat, but check that they have been set to
		# valid values
		
		# Weights must sum to 1
		if (sum(weightmat[,1])!=1){
			stop("Weights must sum to 1\n")
		} 
		
		# The valid update types are 1:7
		validnums=1:7
		if ( length(setdiff(weightmat[,2],validnums))>0 ){
			stop("Update indices must be between 1 and 7\n")
		}
		
		
		# If type is genotype then one of the update types must be 6 or 7
		if ( (args$DataType=="g")&&
				(sum(is.element( c(6,7), weightmat[,2] ))==0 ) ){
			stop("If data is genotype either update type 6 or 7 must be done\n")
		}
	}
	
	write.table(weightmat,WeightFile,row.names=F,col.names=F,quote=F,na="")
	return(args)
	
}




# Nicely print the parameter object with the tags so it is 
# obvious how to access elements of the list
print.pars=function(x, ...)
{
	args=x
	outmat=matrix(nrow=length(args),ncol=2)
	outmat=data.frame(outmat)
	colnames(outmat)=c("Name","Value")
	
	outmat[1:length(args),"Name"]=names(args)
	outmat[1:length(args),"Value"]=unlist(args)	

	cat("Setting name (Name) and its current value (Value):\n")
	print(outmat)
}



	

# Write the settings object to a file. This file can then be used
# in a sampletrees run
writeArgs=function(args, outfile=NULL)
{
	if (is.null(outfile)){
		# Write to default file name
		outfile = "run1.par"
		cat("Options written to file: run1.par\n")
	} else {
		cat(paste("Options written to file: ",outfile,"\n",sep=""))
	}
	
	if (args$clean==FALSE){
	 	warning("Settings haven't been checked and/or Settings contain errors.\n  Sampletrees may not run without errors. Run checkArgs().\n")
 	}
	
	cat(paste("RunName ",args$RunName,"\n",sep=""),file=outfile,append=F)
	if (!is.na(args$Seed)){
		cat(paste("Seed ",args$Seed,"\n",sep=""),file=outfile,append=T)
	}
	
	for (i in 3:16){
		cat(paste(names(args)[i]," ",format(args[[i]], scientific=FALSE),"\n",sep=""),file=outfile,append=T)
	}
	
	if (!is.na(args$InitialTreeFile)){
		cat(paste("InitialTree ",1,"\n",sep=""),file=outfile,append=T)
		cat(paste("InitialTreeFile ",args$InitialTreeFile,"\n",sep=""),file=outfile,append=T)
	}
	
	if (args$RandomTree==TRUE){
		cat(paste("RandomTree ",1,"\n",sep=""),file=outfile,append=T)
	}
	
	if (args$DataType=="g"){
		if (is.na(args$HaploFreqFile)){
			warning("WARNING: Default value for haplofreqfile
					\t\t\t Sampler will not run with default values for this option\n")
		} 
		cat(paste("HaploFreqFile ",args$HaploFreqFile,"\n",sep=""),file=outfile,append=T)
		
		if (!is.na(args$InitialHaploFile)){
			cat(paste("InitialHaplos ",1,"\n",sep=""),file=outfile,append=T)
			cat(paste("InitialHaploFile ",args$InitialHaploFile,"\n",sep=""),file=outfile,append=T)
			
		}
		
		if (!is.na(args$HaploListFile)){
			cat(paste("HaploList ",1,"\n",sep=""),file=outfile,append=T)
			cat(paste("HaploListFile ",args$HaploListFile,"\n",sep=""),file=outfile,append=T)
		}
		
	}

}


# This function is used to start a new job where a previous one finished
restartRun=function(newrunname, oldargs=NULL, argfile=NULL, extrait=NULL, totalsamples=NULL)
{
	if ( is.null(extrait) && is.null(totalsamples) ) {
		stop("Must specify one of extrait (additional iterations to perform) or totalit (number of iterations desired)")
	} else if ( (!is.null(extrait)) && (!is.null(totalsamples)) ){
		stop("Only one of extrait (additional iterations to perform) or totalit (number of iterations desired) can be specified")
	} else if ( is.null(argfile)&&is.null(oldargs)){
		stop("Must specify one of args (an object of class pars ) or argfile (a file that can be read in to an object of class pars)")
	} else if ( (!is.null(argfile))&&(!is.null(oldargs)) ){
		stop("Only one of args (an object of class pars ) or argfile (a file that can be read in to an object of class pars) can be specified")
	} else {
		
		if (!is.null(argfile)){
			if (!file.exists(argfile)){
				stop("Parameter file to read in does not exist\n")
			} else {
				args=readArgs(argfile)
			}
		} else {
			args=oldargs
		}
		oldrunname=args$RunName
		
		
		# Set up new parameter values based on previous run
		args$InitialTreeFile=paste(oldrunname,"_lasttree.out",sep="")
		
		samplefile=paste(oldrunname,"_samples.out",sep="")
		nlines=length(count.fields(samplefile))
		oldres=scan(samplefile,skip=(nlines-1))
		
		if (!is.null(extrait)){
			args$ChainLength=extrait
		} else {
			args$ChainLength=(totalsamples-(nlines-1)*args$Thinning)*args$Thinning
		}
		
		args$BurnIn=0
		args$InitialTheta=oldres[2]
		args$InitialRho=oldres[3]
		
		
		if (args$DataType=="g"){
			
			# Initial haplotypes are in the lasttree results
			args$InitialHaploFile=NA
			
			# Will this file for sure exist?
			fname=paste(oldrunname,"_haplolist.out",sep="")
			if (file.exists(fname)){
				args$HaploListFile=fname
			}
		}
		
		args$Seed=NA
		args$RunName=newrunname	
		
		return(args)
	}
	
}


# This function reads in parameter values from a file
# Must make sure that type of input is correct!
readArgs=function(filename, check=TRUE)
{
	args=newArgs()
	
	if ( !file.exists(filename) ){
		stop("Filename does not exist\n")
	} else {
		
		tryval=try(read.table(filename),silent=TRUE)
		
		if (class(tryval)=="try-error"){
			
			tryval=try(readLines(filename),silent=TRUE)
			if (class(tryval)=="try-error"){
				stop(paste("Check format of input file:",filename))
			} else {
				vals=readLines(filename)
				vals=strsplit(vals,split=" ")
				pad=sapply(vals,length)==1
				
				for (i in which(pad==TRUE)){
					temp=vals[i][[1]]
					temp=c(temp,NA)
					vals[i][[1]]=temp
				}
				
				vals=matrix(unlist(vals),ncol=2,byrow=T)
			}
		} else {
			vals=read.table(filename, colClasses="character")
		}
		
		
		# Vals is a matrix with the first column the tags 
		# and the second column the value. 
		for (i in 1:nrow(vals)){
			tag=names(args)==vals[i,1]
			if (sum(tag)!=0){
				index=which(tag)
				args[[index]]=vals[i,2]
			}
			index=0
		}

		numericargs=c(2,7:16)
		
		for (i in numericargs){
			args[[i]]=as.numeric(args[[i]])
		}
		
	}

	
	# Checking the options that were read in
	if (check==TRUE){
		cat("Checking the validity of the options specified in the parameter file:\n\n")
		args=checkArgs(args)
	}
	
	return(args)
		
}



## This function does a check of all the options specified in args
## to ensure that there are no obvious errors that will cause sampletrees
## to not work.
checkArgs=function(args)
{
	args$clean=TRUE
	
	cat("Error messages:\n")
	
	# Check that the values with no defaults have been specified
	# If they have, ensure the files actually exist
	if (is.na(args$DataFile)||is.na(args$LocationFile)||is.na(args$WeightFile)||
			is.na(args$FocalPoint)){
		cat("\tDataFile, LocationFile, WeightFile  and/or FocalPoint are set to default values\n")
		args$clean=FALSE
	} else {
		if (!file.exists(args$DataFile)){
			cat("\tDataFile does not exist\n")
			args$clean=FALSE
		}
		if (!file.exists(args$LocationFile)){
			cat("\tLocationFile does not exist\n")
			args$clean=FALSE
		}
		if (!file.exists(args$WeightFile)){
			cat("\tWeightFile does not exist\n")
			args$clean=FALSE
		}
	}
	
	if (args$DataType=="g"){
		if (is.na(args$HaploFreqFile)){
			cat("\tHaploFreqFile set to default value\n")
			args$clean=FALSE
		}
	}
	
	
	# Check that the optional files have been specified
	if ( !is.na(args$InitialTreeFile) && !file.exists(args$InitialTreeFile) ){
		cat("\tInitialTreeFile does not exist\n")
		args$clean=FALSE
	}
	
	if ( (args$DataType=="g")&&(!is.na(args$InitialHaploFile)) ){
		if (!file.exists(args$InitialHaploFile)){
			cat("\tInitialHaploFile does not exist\n")
			args$clean=FALSE
		}	
	}
	
	if ( (args$DataType=="g")&&(!is.na(args$HaploList)) ){
		if (!file.exists(args$HaploListFile)){
			cat("\tHaploListFile does not exist\n")
			args$clean=FALSE
		}	
	}
	
	
	# Check that the arguments that are supposed to be numeric actually are numeric
	numericargs=c(args$FocalPoint, args$ChainLength, args$BurnIn,args$Thinning,
			args$InitialTheta, args$MinTheta, args$MaxTheta, args$InitialRho,
			args$ScaleRho, args$ShapeRho,args$Start)
	
	for (i in numericargs){
		if (!is.numeric(i)){
			cat(paste("Non-numeric value for:",i,"\n"))
			args$clean=FALSE
		}
	}
	
	
	# Check that the genotype/sequence file is in the right format
	nchars=-1
	if (!is.na(args$DataFile)&&file.exists(args$DataFile)){
		
		nfields=unique(count.fields(args$DataFile))
		
		if (length(nfields)!=1) {
			cat("All rows of DataFile must have only 1 column\n")
			args$clean=FALSE
		} else {
			
			if (args$DataType=="g"){
				dat=read.table(args$DataFile)
				not012=(dat!=0)&(dat!=1)&(dat!=2)
				if (sum(not012)>0){
					cat("\tAt least one row has a genotype that is not 0, 1 or 2\n")
					args$clean=FALSE
				}
				nchars=nfields
				
			} else if (args$DataType=="h"){
				
				if (nfields!=1){
					cat("\tDataFile must have only 1 column\n")
					args$clean=FALSE 
				} else {
					dat=read.table(args$DataFile,colClasses="character")
					nchars=unique(apply(dat,1,nchar))
					if (length(nchars)!=1){
						cat("\tVariable sequences of length:",nchars,"\n")
						args$clean=FALSE
					}
					temp=unlist(strsplit(dat[,1],split=""))
					not01=(temp!=0)&(temp!=1)
					if (sum(not01)>0){
						cat("\tAt least one sequence is not composed of 0 and 1\n")
						args$clean=FALSE
					}
				}
				
			} else {
				cat("DataType must be 'g' or 'h'\n")
				args$clean=FALSE
			}
		}
		
	}
	
	
	# Check the locations file 
	if (!is.na(args$LocationFile)&&file.exists(args$LocationFile)){
		locations=scan(args$LocationFile,quiet=T)
		if (is.numeric(locations)!=TRUE){
			cat("\tSNP locations must be numeric\n")
			args$clean=FALSE
		} else if (sum(order(locations)!=(1:length(locations)))>0){
			cat("\tSNP locations must be in increasing order\n")
			args$clean=FALSE
		} else if ( (args$FocalPoint<locations[1])||
				(args$FocalPoint>locations[length(locations)]) ){
			cat(paste("\tFocal point must be >=",locations[1],"or","<=",locations[length(locations)],"\n"))
			args$clean=FALSE 
		} 
		if ( (nchars!=-1)&&(length(locations)!=nchars) ){
			cat("\tNumber of loci in data file does not match number of loci in locations file\n")
			args$clean=FALSE
		} 
	}
	
	
	# Check the weights file
	if (!is.na(args$WeightFile)&&(file.exists(args$WeightFile))){
		weights=strsplit(readLines(args$WeightFile),split=" ")
		sumweights=sum(as.numeric(sapply(weights,"[[",1)))
		if (sumweights!=1){
			cat("\tWeights do not sum to 1\n")
			args$clean=FALSE
		} 
		veclengths=sapply(weights,length)
		if (args$DataType=="g"){
			patt="[1-7]"
		} else {
			patt="[1-5]"
		}
		all=NULL
		for (i in 1:length(weights)){
			temp=weights[[i]][2:veclengths[i]]
			if (length(grep(patt,temp,invert=TRUE))>0){
				cat("\tUpdate types must be 1-5 for haplotype or 1-7 for genotype\n")
				args$clean=FALSE	
			}
			all=c(all,temp)
		}
		if (args$DataType=="g"){
			if (length(grep("[6-7]",all))==0){
				cat("\tUpdate type must include 6 or 7 for genotype\n")
				args$clean=FALSE
			}
		} else if (args$DataType=="h"){
			if (length(grep("[6-7]",all))!=0){
				cat("\tUpdate type can not includ 6 or 7 for haplotype\n")
				args$clean=FALSE
			}
			
		}
	}
	
	
	# Check the intial haplotype file if option is used
	if ( (args$DataType=="g")&&(!is.na(args$InitialHaploFile))&&
			file.exists(args$InitialHaploFile) ) {
		
		ncols=unique(count.fields(args$InitialHaploFile))
		if (ncols!=1){
			cat("\tInitialHaploFile must have only 1 column\n")
			args$clean=FALSE 
		} else {
			initdat=read.table(args$InitialHaploFile,colClasses="character")
			ncols=unique(apply(initdat,1,nchar))
			if (length(ncols)!=1){
				cat("\tVariable sequence lengths in InitialHaploFile\n")
				args$clean=FALSE
			} else if (ncols!=nchars){
				cat("\tInitialHaploFile must have same number of loci as DataFile\n")
				args$clean=FALSE
			} 
			temp=matrix(as.numeric(unlist(strsplit(initdat[,1],split=""))),ncol=nchars,byrow=T)
			not01=(temp!=0)&(temp!=1)
			if (sum(not01)>0){
				cat("\tAt least one sequence in InitialHaploFile is not composed of 0 and 1\n")
				args$clean=FALSE
			}
			genos=temp[c(TRUE,FALSE),]+temp[c(FALSE,TRUE),]
			if (sum(genos!=dat)>0){
				cat("\tInitial haplotypes in InitialHaploFile do not match genotypes in DataFile\n")
				args$clean=FALSE
			}
			
		}
		
	}
	
	# Check the haplolist file
	if ( (args$DataType=="g")&&(!is.na(args$HaploListFile))&&
			file.exists(args$HaploListFile) ) {
		
		ncols=unique(count.fields(args$HaploListFile))
		if (ncols!=1){
			cat("\tHaploListFile must have only 1 column\n")
			args$clean=FALSE 
		} else {
			
			initdat=read.table(args$HaploListFile)
			if (sum(apply(initdat,1,is.numeric))!=nrow(initdat)){
				cat("\tHaploListFile contains non-numeric value\n")
				args$clean=FALSE
			} else {
				
				initdat=read.table(args$HaploListFile,colClasses="character")
				ncols=unique(apply(initdat,1,nchar))
				if (length(ncols)!=1){
					cat("\tVariable sequence lengths in HaploListFile\n")
					args$clean=FALSE
				} else if (ncols!=nchars){
					cat("\tHaploListFile must have same number of loci as DataFile\n")
					args$clean=FALSE
				}  else {
					temp=matrix(as.numeric(unlist(strsplit(initdat[,1],split=""))),ncol=nchars,byrow=T)
					not01=(temp!=0)&(temp!=1)
					if (sum(not01)>0){
						cat("\tAt least one sequence in HaploListFile is not composed of 0 and 1\n")
						args$clean=FALSE
					}
				}
			}
		}
		
	}
	
	
	# Check haplotype frequency file
	if ( (args$DataType=="g")&&(!is.na(args$HaploFreqFile))&&
			file.exists(args$HaploFreqFile) ) {
		freqdat=read.table(args$HaploFreqFile)
		if (nrow(freqdat)!=(nchars-1)){
			cat("\tRows in HaploFreqFile don't match number of loci\n")
			args$clean=FALSE
		} else if (!isTRUE(all.equal(rowSums(freqdat),rep(1,nrow(freqdat))))){
			cat("\tRows of HaploFreqFile must sum to 1\n")
			args$clean=FALSE
		}
	}
	
	
	
	
	if (args$clean==TRUE){
		cat("\tNo error messages. Sampler will run.\n")
	} else {
		cat("\tMust fix errors before running sampler.\n")
	}
	
	cat("Warning messages:\n")
	if (args$RunName=="Run"){
		cat("\tDefault run name used. May overwrite previous results\n")
	} else if (file.exists(paste(args$Runname,"_samples.out",sep=""))){
		cat("\tOutput file already exists. May overwrite previous results\n")
	} else{ 
		cat("\tNo warning messages.\n")
	}	
	
	return(args)	
	
}