## This file contains the functions for processing the sampletrees 
## output


## This function is used to add the trees from a sampletrees run
## to the sampletrees output object
addTrees=function(output, all=TRUE, lines=NULL, start=1, end=NULL, nlines=NULL){
	
	#require(ape)
	
	
	# Check that trees aren't already read in
	if (class(output$rawdata$Trees)=="multiPhylo"){
		stop("Trees already saved in the tree output object.\n")
	} else {
		output$rawdata$Trees=readTrees(output=output, all=all,start=start,lines=lines,end=end,nlines=nlines)
	}
	return(output)
}


## This function is used to read in trees from a run of
## sampletrees. It is expected to mostly be used as a helper function by
## addTrees and treeapply
readTrees=function(output=NULL, filenames=NULL, all=TRUE,lines=NULL,start=1, end=NULL, nlines=NULL)
{
	#require(ape)

	if (class(output$rawdata$Trees)=="multiPhylo"){
		cat("Trees already read in to the tree output object")
		return(output$rawdata$Trees)
	} else {
	
		if (is.null(filenames)&&is.null(output)){
			stop("At least one of filenames or output must not be null\n")
		} 

		if (!is.null(output)){
			filenames=output$rawdata$Trees	
		}
	
		# Check that file(s) exists
		if ( sum(!file.exists(filenames))>0 ){
			stop(paste("Tree file(s) ",filenames," don't exist\n"))
		} 
	
		# Get data on the files
		numfiles=length(filenames)
		numlines=vector(length=length(numfiles))
		lineindex=vector(length=length(numfiles))
		total=0
		found=FALSE
		filestart=start
	
		for (i in 1:numfiles){
			numlines[i]=length(count.fields(filenames[i]))
			lineindex[i]=total
			total=total+numlines[i]
			if ((start<total)&&(found==FALSE)){
				startfileindex=i
				if (startfileindex>1){
					filestart=start-sum(numlines[1:(i-1)])
				}
				found=TRUE
			}
		}
	
		if ((start<0)||(start>sum(numlines))){
			stop("The first line to be read in is either too low or too high\n")
		}

	
		if (all==TRUE){
			# Read in all the trees

			for (i in 1:length(filenames)){
				
				if (i==1){
					trees=read.tree(filenames[i])
				} else {
					trees=c(trees, read.tree(filenames[i]))
				}
	
			}
			names(trees)=as.character(output$rawdata$i)
				
		} else if (!is.null(lines)){
		
			# The lines to be read in are stored in the vector called lines
			if (!is.numeric(sum(lines))){
				# If a subset is to be read in, ensure that the number
				# of lines has to read in has been specified
				stop("Line values must be numeric\n")
			} else {
			
				toread=vector(length=length(lines))
				treedat=vector(length=length(lines))
				treenames=vector(length=length(lines))
			
				lower=1
				for (i in 1:length(numlines)){
					upper=lower+numlines[i]
					inrange=(lines>=lower)&(lines<upper)
					toread[inrange]=i
					lower=upper
				}
				if (sum(toread==0)>0){
					stop("More lines to be read in than are in tree files\n ")
				}
			

				for (i in 1:length(lines)){
					temp=scan(file=filenames[toread[i]], skip=(lines[i]-lineindex[toread[i]]-1), 
						nlines=1, what="character",quiet=TRUE)
					treedat[i]=temp[2]
				}
				trees=read.tree(text=treedat, tree.names=output$rawdata$i[lines])
			}
		
		} else if ((!is.null(nlines))||(!is.null(end))) {
		
		
			# The lines to be read in start at 'start', and end when 'nlines' has been
			# read in
			if (!is.null(end)){
				nlines=end-start+1
				if ( (end<start) || (end>sum(numlines)) ){
					stop("Last line to be read in must be greater than first and less than max number of trees in the files\n")
				}
			} else {
				end=start+nlines-1
				if ( (nlines>sum(numlines)) || (end>sum(numlines))){
					stop("More lines to be read in than are in tree files\n ")
				}
			}
		
			numtoread=nlines
			fileindex=startfileindex
			treedat=NULL
		
			while (numtoread>0){
				numfromfile=min(numlines[fileindex]-filestart+1,numtoread)
				temp=scan(file=filenames[fileindex], skip=filestart-1, 
					nlines=2*numfromfile, what="character",quiet=TRUE)
				treedat=c(treedat,temp[2*(1:numfromfile)])
				fileindex=fileindex+1
				numtoread=numtoread-numfromfile
				filestart=1
			}
		
			trees=read.tree(text=treedat, tree.names=output$rawdata$i[start:end])
		
		} else {
			stop("If all=FALSE must specify values for lines OR nlines OR stop")
		}
	
		return(trees)
	}

}



# This function is used to read in sampletrees output. Note that it does not
# read in the trees as these are processed separately.
readOutput=function(argobj=NULL, argfile=NULL, addtrees=FALSE)
{

	if (is.null(argobj)&&is.null(argfile)){
		stop("One of 'myargs' or 'argfile' must not be NULL")
	} else if (is.null(argobj)){
		argobj=readArgs(argfile, check=FALSE)
	}
	
	argobj$clean=TRUE
	
	# Read in the raw results and create a list of the raw results
	dat=read.table(paste(argobj$RunName,"_samples.out",sep=""), header=T)
	rawResults=list(i=dat[,"i"])
	rawResults$Theta=dat[,"theta"]
	rawResults$Rho=dat[,"rho"]
	rawResults$Trees=paste(argobj$RunName,"_trees.out",sep="")
	
	
	# Read in acceptance rates and save it as a list.
	# This list will be appended to by other functions
	accept=read.table(paste(argobj$RunName,"_accept.out",sep=""))
	colnames(accept)=c("Type","Total","Accepted")
	procResults=list(Accept=accept)
	
	
	# Check if there are processed results and read them in 
	fname=paste(argobj$RunName,"_proctrees.out",sep="")
	if (file.exists(fname)){
		procResults$TreeStats=read.table(fname,header=T)
	}
	
	
	# Store all the results as a list of three elements. Each element of the list
	# is in turn a list
	# List consists of (1) metadata giving settings
	#                  (2) list of raw data from the run
	#                  (3) list of processed data from the run
	results=list(runinfo=argobj, rawdata=rawResults, procdata=procResults)
	class(results)="treeoutput"
	if (addtrees==TRUE){
		results=addTrees(results)
	}

	return(results)
}


# This function is used to change any of the options in the
# settings object. Because the settings object is a list, 
# one could just change the values of the list. But it is 
# safer to use this function because it will test that
# the option that is being changed is actually valid. 
changeArgs.treeoutput=function(object,...)
{
	results=object
	args=results$runinfo
	otherargs=list(...)
	if (length(otherargs)!=0){
		
		inboth=intersect(names(args),names(otherargs))
		args[inboth]=otherargs[inboth]
		notfound <- setdiff(names(otherargs),names(args))
		
		if (length(notfound)!=0) {
			warning(paste("The following are not valid settings: ", notfound))
		}
		
	}
	
	results$runinfo=args
	
	return(results)	
}


print.treeoutput=function(x, ...)
{
	output=x
	
	# Print the run settings
	info=output$runinfo
	class(info)="pars"
	print(info)	
	
	# Print the run statistics
	cat("\nSampletrees run information:\n")
	cat(paste("\tTotal number samples read in: ",length(output$rawdata$i),"\n"))
	cat(paste("\tFirst sample index: ",output$rawdata$i[1],"\n"))
	cat(paste("\tLast sample index: ",output$rawdata$i[length(output$rawdata$i)],"\n"))
	
	# Print information about processed results
		if (length(output$procdata)>1){
		cat(paste("\tNumber summary statistics available: ",ncol(output$procdata$TreeStats)-1,"\n"))
	}
	
}




summary.treeoutput=function(object, ...)
{
	
	output=object
	
	outmat=matrix(nrow=6,ncol=2)
	outmat=data.frame(outmat)
	colnames(outmat)=c("Name","Value")
	
	outmat[1:3,"Name"]=names(output$runinfo)[8:10]
	outmat[1:3,"Value"]=unlist(output$runinfo)[8:10]
	outmat[4,]=c("First Observation",output$rawdata$i[1])
	outmat[5,]=c("Last Observation",output$rawdata$i[length(output$rawdata$i)])
	outmat[6,]=c("Number Observations", length(output$rawdata$i))
	
	print(outmat)
	
	cat("\nAcceptance proportions:\n")
	temptable=output$procdata$Accept
	temptable[output$procdata$Accept[,"Type"]==3,"Total"]=1
	propns=temptable[,"Accepted"]/temptable[,"Total"]
	names(propns)=output$procdata$Accept[,"Type"]
	print(round(propns,5))
	
	cat("\nMCMC output summaries:\n")
	cat("Theta:\n")
	print(summary(output$rawdata$Theta))
	cat("\nRho:\n")
	print(summary(output$rawdata$Rho))
	
	if (!is.null(output$procdata$TreeStats)){
		for (i in 2:ncol(output$procdata$TreeStats)){
			cat(paste("\n",names(output$procdata$TreeStats)[i],":\n",sep=""))
			print(summary(output$procdata$TreeStats[,i]))
		}
		
	}
		
}



treeapply=function(output, myfunc, funcname=NULL, maxlines=1000, treerange=NULL)
{
	#require(ape)

	# Compute the statistic on the trees
	if (class(output$rawdata$Trees)=="multiPhylo"){
		
		# The trees have been read in so use those
		treestat=unlist(lapply(output$rawdata$Trees, myfunc))
		treenames=as.numeric(names(output$rawdata$Trees))
		rnames=1:length(treestat)
		
	} else{

		# Begin reading in trees. If there is a range, use the range; otherwise, apply 
		# to all trees in the file
		if (!is.null(treerange)){
			
			testrange=treerange[1]+0:(length(treerange)-1)
			if (sum(testrange-treerange)==0){
				
				# Range consists of consecutive numbers so faster use of scan:
				trees=readTrees(output=output,all=FALSE, lines=NULL, 
						        start=treerange[1],end=treerange[length(treerange)])
						
			} else {
				
				# Range consists of non-consecutive subset
				trees=readTrees(output=output,all=FALSE, lines=treerange)
				
			}
			
			treestat=unlist(lapply(trees, myfunc))
			treenames=output$rawdata$i[treerange]
			rnames=treerange
			
		} else {
			
			treestat=NULL
			
			# No special range used
			for (i in 1:length(output$rawdata$Trees)){
				
				con=file(output$rawdata$Trees[i],open="rt")
				
				while (length(input <- readLines(con, n=maxlines)) > 0){ 
					trees=read.tree(text=input)
					substat=unlist(lapply(trees, myfunc))
					treestat=c(treestat,substat)
				} 
				close(con)
			}
			
			treenames=output$rawdata$i
			rnames=1:length(treestat)
		}
			
	}


	df=data.frame(treenames,treestat)
	rownames(df)=rnames
	if (!is.null(funcname)){
		colnames(df)=c("index",funcname)
	} else {
		colnames(df)=c("index","stat")
	}
	
	return(df)
}


addTreeStat=function(output, myfunc, funcname=NULL, maxlines=1000, treerange=NULL){
	
	if (!is.data.frame(output$procdata$TreeStats)){
			
		# This is the first tree statistic to be saved
		treestat=treeapply(output, myfunc, funcname=funcname, maxlines=maxlines,treerange=treerange)
		ncols=ncol(treestat)
		output$procdata$TreeStats=treestat
			
	} else {
			
		# Tree stats have already been read in so use the same tree indices
		treeindex=as.numeric(rownames(output$procdata$TreeStats))
		nrows=nrow(output$procdata$TreeStats)
			
		if ( (!is.null(treerange))&&(!identical(treerange,treeindex)) ){
			warning("Provided treerange is different from those of the TreeStats data.frame. 
					 The tree indices will be read from TreeStats\n")
			treerange=treeindex
		} else if (!identical(1:nrows,treeindex)){
			treerange=treeindex
		}
		treestat=treeapply(output, myfunc, funcname=funcname, maxlines=maxlines,treerange=treerange)
		
		
		# Add the statistic to the output object
		newdf=merge(output$procdata$TreeStats,treestat, by="index")
		rownames(newdf)=rownames(treestat)
		ncols=ncol(newdf)
		output$procdata$TreeStats=newdf
	}
	
	if (!is.null(funcname)){
		colnames(output$procdata$TreeStats)[ncols]=funcname
	} else {
		colnames(output$procdata$TreeStats)[ncols]=paste("Stat.",ncols-1,sep="")
	}
			

	return(output)
	
}

plot.treeoutput=function(x,oneperpage=FALSE,
		           asktoplot=FALSE, layoutmat=NULL,statnames=NULL, ...)
{	
	
	output=x
	
	# To modify settings to hit enter before next plot. This
	# must be re-set to initial value at end of function
	oldpar=par()$ask
	par(ask=asktoplot)
	
	
	# Set up the matrices for plotting
	if (!is.null(layoutmat)){
		if (is.matrix(layoutmat)){
			m=layoutmat
		} else{
			stop("Variable layoutmat must be of type matrix\n")
		}
	} else {
		if (is.data.frame(output$procdata$TreeStats)){
			m=matrix(c(1,2,3,4), byrow=T, nrow = 2, ncol = 2)
		} else {
			m=matrix(c(1,1,2,2,0,3,3,0), byrow=T, nrow = 2, ncol = 4)
		}
	}
	
	if (oneperpage==TRUE){
		plotdims=c(1,1)
		m=matrix(c(1),nrow=1,ncol=1)
	}
	layout(m)
	
	# Plot the acceptance rate
	newcounts=output$procdata$Accept[,2]
	newcounts[output$procdata$Accept[,1]==3]=1
	
	barplot(output$procdata$Accept[,3]/newcounts,ylim=c(0,1),
			ylab="acceptance proportion",xlab="update label",
			names=output$procdata$Accept[,1],main="Update acceptance proportions")
	
	# Plot theta and rho
	xs=output$rawdata$i
	
	plot(xs,output$rawdata$Theta, main="Traceplot of Theta", xlab="Iterations",ylab="Theta", type="l")
	lines(lowess(xs,output$rawdata$Theta),col="red")
	plot(xs,output$rawdata$Rho, main="Traceplot of Rho", xlab="Iterations",ylab="Rho",type="l")
	lines(lowess(xs,output$rawdata$Rho),col="red")


	# Plot the tree stats
	if (is.data.frame(output$procdata$TreeStats)){
		if (is.null(statnames)){
			plottitles=colnames(output$procdata$TreeStats)[2:ncol(output$procdata$TreeStats)]
		} else {
			if (length(statnames)!=(ncol(output$procdata$TreeStats)-1)){
				stop("Length of vector of tree statistics names doesn't match number of tree statistics\n ")
			} else {
				plottitles=statnames
			}
		}
		xs=output$procdata$TreeStats[,1]
	
		for (i in 2:ncol(output$procdata$TreeStats)){
			plot(xs,output$procdata$TreeStats[,i], main=paste("Traceplot of",plottitles[i-1]), 
					ylab=plottitles[i-1], xlab="Iterations",type="l")
			lines(lowess(xs,output$procdata$TreeStats[,i]),col="red")
		}
	}

	
	par(ask=oldpar)

}


merge.treeoutput=function(x,y,runname=NULL,...)
{


	output1=x
	output2=y
	
	## For a merged job, it is assumed that output2 started where
	## output2 finished and so many of the options should be the 
	## same between them
	
	## Error checking: Make sure that options that should be the same
	## are the same
	checkindex=c(3:7,12:13,15:16)
	for (i in checkindex){
		if (output1$runinfo[[i]]!=output2$runinfo[[i]]){
			stop(paste(names(output1$runinfo[[i]])," must have same value in both settings objects\n"))
		}
	}
	
	if ((output1$runinfo$DataType=='g')&&(output1$runinfo$HaploFreqFile!=output2$runinfo$HaploFreqFile)){
			stop(paste(names(output1$runinfo$HaploFreqFile)," must have same value in both settings objects\n"))
	}
	
	## This is not an error, but it will mess up calculating the first/last values 
	## and plotting
	if (output1$runinfo$Thinning!=output2$runinfo$Thinning){
		warning("Two different values used for thinning interval\n")
		newoutput1$runinfo$Thinning=NA
	}
	
	if (is.null(runname)){
		runname="Merge"
	}
	
	
	## The initial values are all those from output1. 
	newoutput=output1
	newoutput$runinfo$ChainLength=output1$runinfo$ChainLength+output2$runinfo$ChainLength
	newoutput$runinfo$RunName=runname
	
	
	## Merge the acceptance tables
	newtable=output1$procdata$Accept+output2$procdata$Accept
	newtable[,"Type"]=output1$procdata$Accept[,"Type"]
	
	if (sum(newtable[,"Type"]==3)==1){
		index=which(output1$procdata$Accept[,"Type"]==3)
		newtotal=output1$procdata$Accept[index,"Total"]+output2$procdata$Accept[index,"Total"]
		newpropn=1/newtotal*(output1$procdata$Accept[index,"Total"]*output1$procdata$Accept[index,"Accepted"]+
					output2$procdata$Accept[index,"Total"]*output2$procdata$Accept[index,"Accepted"])
	 	newtable[index,"Total"]=newtotal
		newtable[index,"Accepted"]=newpropn
	}
	newoutput$procdata$Accept=newtable
	
	
	## Merge the i, Theta and Rho values
	newoutput$rawdata$Theta=c(output1$rawdata$Theta,output2$rawdata$Theta)
	newoutput$rawdata$Rho=c(output1$rawdata$Rho,output2$rawdata$Rho)
	newoutput$rawdata$i=c(output1$rawdata$i,output1$rawdata$i[length(output1$rawdata$i)]+output2$rawdata$i)
	
	
	## Merge the Tree values
	newoutput$rawdata$Trees=c(output1$rawdata$Trees,output2$rawdata$Trees)
	
	
	## Merge the tree statistics
	treestatflag=c(is.data.frame(output1$procdata$TreeStats),is.data.frame(output2$procdata$TreeStats))
	if (sum(treestatflag)!=0){
		if (sum(treestatflag)==1){
			stop("Error: Tree statistics must be computed on both runs or on neither run")	
		} else {
			if (ncol(output1$procdata$TreeStats)!=ncol(output2$procdata$TreeStats)){
				stop("Error: The same sets of tree statistics must be computed for both runs")	
			} else {
				newoutput$procdata$TreeStats=rbind(output1$procdata$TreeStats,output2$procdata$TreeStats)
			}
		}
	}
	
	return(newoutput)
	
}


writeTreeoutput=function(output, argfile=NULL, ask=TRUE)
{
	
	#require("ape")
	runname=output$runinfo$RunName
	
	# Write out the settings
	newpars=output$runinfo
	newpars$clean=TRUE
	class(newpars)="pars"
	
	if (is.null(argfile)){
		argfile=paste(runname,"_Settings.par",sep="")
	}
	
	if ((ask==FALSE)||checkfile(argfile)){	
		writeArgs(newpars, outfile=argfile)
	} else {
		cat(paste("No output written to file ",argfile,"\n",sep=""))
	}
	
	
	# Write out theta and rho to file 

	samples=data.frame(output$rawdata$i, output$rawdata$Theta, output$rawdata$Rho)
	colnames(samples)=c("i","Theta","Rho")
	resfilename=paste(runname,"_samples.out", sep="")
	
	if ((ask==FALSE)||checkfile(resfilename)){
		write.table(samples,resfilename,quote=F,row.names=F)
	} else {
		cat(paste("No output written to ",resfilename,"\n",sep=""))
	}
	
	
	# Write out trees to file
	resfilename=paste(runname,"_trees.out", sep="")
	
	if ((ask==FALSE)||checkfile(resfilename)){
	
		if (class(output$rawdata$Trees)=="multiPhylo"){
			write.tree(output$rawdata$Trees,resfilename)
		} else if (length(output$rawdata$Trees)==1){
			file.copy(from=output$rawdata$Trees, to=resfilename, overwrite = FALSE)
		} else if (length(output$rawdata$Trees)>1){
			for (i in 1:length(output$rawdata$Trees)){
				file.append(resfilename, output$rawdata$Trees[i])
			}	
		}
		
	} else {
		cat(paste("No output written to ",resfilename,"\n",sep=""))
	}
	
	
	# Write out acceptances to file
	resfilename=paste(runname,"_accept.out", sep="")
	if ((ask==FALSE)||checkfile(resfilename)){
		write.table(output$procdata$Accept,resfilename,quote=F,row.names=F,col.names=F)
	} else {
		cat(paste("No output written to ",resfilename,"\n",sep=""))
	}
	
	
	# Write out any processed results to file. This assumes they are
	# stored as a dataframe
	if (!is.null(output$procdata$TreeStats)){
		
		resfilename=paste(runname,"_proctrees.out",sep="")
		
		if ((ask==FALSE)||checkfile(resfilename)){
			write.table(output$procdata[[2]],resfilename,quote=F,row.names=F)
		} else {
			cat(paste("No output written to ",resfilename,"\n",sep=""))
		}
		
	}
	
	
}

# Helper function
checkfile=function(filename){

	check=TRUE
	
	if (file.exists(filename)){
		check=FALSE
		response='x'
		while( !identical(response, 'y')&&!identical(response, 'n') ){
			cat(paste("Overwrite ",filename," (y or n)?\n",sep=""))
			response=scan(n=1, what="character",quiet=TRUE)
			if (identical(response, 'y')){
				check=TRUE
			}	 
		}
	}
	
	return(check)
}	