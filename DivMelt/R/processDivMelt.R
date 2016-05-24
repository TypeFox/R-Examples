processDivMelt <- function(dir, singleFile="",selectList="", 
	sampleMask="",
	theta1=50., theta2=30., 
	format="pdf", showStats=TRUE, 
	excludeEarly=TRUE,earlyCutoff=3.0, 
	excludeLate=FALSE,lateCutoff=3.0, 
	includeShoulders=TRUE,shoulderCutoff=0.1, 
	showFlorSmooth=FALSE,
	showDFlorSmooth=TRUE,
	useLasso=TRUE,
	t1Cutoff=1.0,
	t2Cutoff=1.0,
	t1SlopeWindow=1.0,
	t2SlopeWindow=1.0,
	scoreFile="",
	dataFile="",
	modelFile="", 
	rejectsFile="rejected.csv",
	statsFile="",
	xrangeOverride=FALSE,
	xrangemin=75,
	xrangemax=100,
	progressFunc=NULL){

     if(missing(dir)) stop("Need dir argument to know where the data reside")

    # --------------------------------------------------------------------------- #
    # set some parameters
    # --------------------------------------------------------------------------- #
    # filesystem constants
    data.path <- dir

    # how tightly do we want to smooth the curve?
    loess.span <- 0.05

    data.pathparts <- strsplit(data.path,'/')
    data.dir <- data.pathparts[[1]][length(data.pathparts[[1]])]

    if (selectList != "")
	selected <- matrix (scan (file=selectList, skip=0, what='character', sep="\t"), 
                    ncol=1, byrow=TRUE)

    # read in all DivMelt scores so that we can plot them on the same graphs
    #file.DivMelt <- "DivMelt_scores.csv"

    file.DivMelt <- scoreFile

    #"Amplicon","P","R","C","PTID","Date","Region","NegFilter","Width","Tm","StartTm","EndTm"
    #"gag-pol",1,"A",1,"10001322","8/18/1995","GAG1","Positive",4.89,84.47,81.67,86.56

    asc <- function(x) { strtoi(charToRaw(x),16L) }

    chr <- function(n) { rawToChar(as.raw(n)) }

    colFromWell <- function (wellNum) {
	return (floor((as.integer(wellNum)-1) / 8) + 1)
    }

    rowFromWell <- function (wellNum) {
	return ((as.integer(wellNum)-1) %% 8)
    }

    load.DivMelt <- function (file) {
	#Plate,Sample,Row,Col,T1,T2,Peak,Score
	#015-men-batch-1-wit-091008,1_050205675_GAG2,A,1,88.77,92.58,90.683,3.82
	data <- matrix (scan (file=file, skip=1, what='character',
		    sep=","), ncol=8, byrow=TRUE)
	diffs.names <- c ("Plate","Sample","Row","Col","T1","T2","Peak","Score")
	names (data) <- diffs.names
	return (data)
    }

    if (file.DivMelt != "")
    {
	print(paste("Loading DivMelt file ",file.DivMelt,sep=""))
	data.DivMelt <- load.DivMelt (file.DivMelt)

	newdat<-matrix(nrow=nrow(data.DivMelt),ncol=7)

	# set up our working data
	for(i in 1: nrow(data.DivMelt))     
	{
	    newdat[i,1]=data.DivMelt[i,2]
	    newdat[i,2]=data.DivMelt[i,5]
	    newdat[i,3]=data.DivMelt[i,6]
	    newdat[i,4]=data.DivMelt[i,7]
	    newdat[i,5]=data.DivMelt[i,8]
	    newdat[i,6]=data.DivMelt[i,3]
	    newdat[i,7]=data.DivMelt[i,4]
	    print(paste("Sample:",newdat[i,1]))
	    print(paste("Row:",newdat[i,6]))
	    print(paste("Col:",newdat[i,7]))
	}
    }

    #specify aspect ratio of final plot (drives interpretation of theta)
    plot.asp <- 2.0
    plot.bins <- 7
    plot.minTemp <- 71.

    scan.uptrough_limit <- earlyCutoff
    scan.dntrough_limit <- lateCutoff

    scan.uptrough = 0
    scan.dntrough = 0

    # what normal angles should indicate the location of our inflection points?
    # (units = degrees, for human convenience)
    theta.low <- theta1
    theta.high <- -theta2

    # data points are 0.1 degrees apart -> convert from deg to num samples
    t1SlopeWS = t1SlopeWindow*10.0
    t2SlopeWS = t2SlopeWindow*10.0

    degtorad=(22.0/7.0)/180.0
    # (let's convert the above into radians)
    theta.low.rad <- degtorad *(theta.low)
    theta.high.rad <- degtorad *(theta.high)

    if (dataFile=="")
    {
	data.out  <- paste(data.dir,".csv",sep="")
    } else {
	data.out  <- dataFile
    }

    basemeasures=27
    if (file.DivMelt != "")
    {
	extraplaces=11
    } else {
	extraplaces=0
    }
    if (useLasso)
    {
	basemeasures=basemeasures+1
    }
    out<- vector(length=basemeasures+extraplaces)

    outmeasure=1
    out[outmeasure]="Plate"
    outmeasure=outmeasure+1
    out[outmeasure]="Sample"
    outmeasure=outmeasure+1
    out[outmeasure]="Row"
    outmeasure=outmeasure+1
    out[outmeasure]="Col"
    outmeasure=outmeasure+1
    out[outmeasure]="T1"
    outmeasure=outmeasure+1
    out[outmeasure]="T2"
    outmeasure=outmeasure+1
    out[outmeasure]="Peak"
    outmeasure=outmeasure+1
    out[outmeasure]="Score"
    outmeasure=outmeasure+1
    out[outmeasure]="theta1"
    outmeasure=outmeasure+1
    out[outmeasure]="theta2"
    outmeasure=outmeasure+1
    out[outmeasure]="shouldersOn"
    outmeasure=outmeasure+1
    out[outmeasure]="shoulderCutoff"
    outmeasure=outmeasure+1
    out[outmeasure]="earlyOn"
    outmeasure=outmeasure+1
    out[outmeasure]="earlyCutoff"
    outmeasure=outmeasure+1
    out[outmeasure]="lateOn"
    outmeasure=outmeasure+1
    out[outmeasure]="lateCutoff"
    outmeasure=outmeasure+1
    out[outmeasure]="t1HldInterval"
    outmeasure=outmeasure+1
    out[outmeasure]="t2HldInterval"
    outmeasure=outmeasure+1
    out[outmeasure]="t1SlpWindow"
    outmeasure=outmeasure+1
    out[outmeasure]="t2SlpWindow"
    outmeasure=outmeasure+1
    out[outmeasure]="altT1s[1]"
    outmeasure=outmeasure+1
    out[outmeasure]="altT1s[2]"
    outmeasure=outmeasure+1
    out[outmeasure]="altT1s[3]"
    outmeasure=outmeasure+1
    out[outmeasure]="altT2s[1]"
    outmeasure=outmeasure+1
    out[outmeasure]="altT2s[2]"
    outmeasure=outmeasure+1
    out[outmeasure]="altT2s[3]"
    outmeasure=outmeasure+1
    out[outmeasure]="errorCode"
    outmeasure=outmeasure+1

    if (file.DivMelt != "")
    {
	out[outmeasure]="File Sample"
	outmeasure=outmeasure+1
	out[outmeasure]="File T1"
	outmeasure=outmeasure+1
	out[outmeasure]="File T2"
	outmeasure=outmeasure+1
	out[outmeasure]="File Peak"
	outmeasure=outmeasure+1
	out[outmeasure]="File Score"
	outmeasure=outmeasure+1
	out[outmeasure]="File Row"
	outmeasure=outmeasure+1
	out[outmeasure]="File Col"
	outmeasure=outmeasure+1
	out[outmeasure]="dT1"
	outmeasure=outmeasure+1
	out[outmeasure]="dT2"
	outmeasure=outmeasure+1
	out[outmeasure]="dPeak"
	outmeasure=outmeasure+1
	out[outmeasure]="dScore"
	outmeasure=outmeasure+1
    }
    if (useLasso)
    {
	out[outmeasure]="Accepted"
	outmeasure=outmeasure+1
    }
    write (out,file=data.out,ncolumns=outmeasure,sep=",",append=FALSE)

    stats<- vector(length=12)
    statsnum=1
    stats[statsnum]="Sample"
    statsnum=statsnum+1
    stats[statsnum]="Row"
    statsnum=statsnum+1
    stats[statsnum]="Col"
    statsnum=statsnum+1
    stats[statsnum]="min"
    statsnum=statsnum+1
    stats[statsnum]="max"
    statsnum=statsnum+1
    stats[statsnum]="mean"
    statsnum=statsnum+1
    stats[statsnum]="stddev"
    statsnum=statsnum+1
    stats[statsnum]="sigma"
    statsnum=statsnum+1
    stats[statsnum]="dmin"
    statsnum=statsnum+1
    stats[statsnum]="dmax"
    statsnum=statsnum+1
    stats[statsnum]="dmean"
    statsnum=statsnum+1
    stats[statsnum]="dstddev"
    statsnum=statsnum+1

    if (statsFile != "")
    {
	stats.out  <- statsFile
	write (stats,file=stats.out,ncolumns=statsnum,sep=",",append=FALSE)
    }

    model <-NULL
    if (useLasso)
    {
	library ("glmnet")
	suspect.out  <- rejectsFile
	suspect.list<- vector(length=2)
	suspect.list[1]="Sample"
	suspect.list[2]="Pred"
	write (suspect.list,file=suspect.out,ncolumns=2,sep=",",append=FALSE)

	if (modelFile != "")
	{
	    load (file=modelFile)
	    print(paste("Loaded model file ",modelFile,sep=""))
	}
	else
	{
	    data (model)
	    pkgDir<-system.file(package="DivMelt")
	    load (file=paste(pkgDir,"/data/model.rda",sep=""))
	    print("Loaded default model file")
	}
    }


    # --------------------------------------------------------------------------- #
    # open our data file and format our data objects
    # --------------------------------------------------------------------------- #

    if (singleFile == "")
	files<-list.files(data.path,pattern="*.FLO")
    else
	files=list(gsub("ABT","FLO",singleFile))

    firstFile=1

    for (fileind in 1:length(files))
    {
	file.flor <- files[fileind]
	file.desc <- gsub("FLO","ABT",files[fileind])

	print(paste("Processing file ",files[fileind],sep=""))

	if (format == "pdf" && selectList=="" && sampleMask=="")
	{
	    plot.name  <- paste(data.dir,"_",gsub(".FLO","",file.flor),sep="")
	    plot.name  <- paste(plot.name,"pdf",sep=".")
	    pdf(plot.name,width=8.5)
	    if (plot.asp > 1)
	    {
		#par (mfrow = c(2, 2),pin=c(3.2,(3.2/plot.asp)),oma=c(0,5,0,0))
		## divide the device into three rows and two columns
		## allocate figure 1 all of row 1
		## allocate figure 2 column 1 of row 2
		## allocate figure 3 column 2 of row 2
		## allocate figure 4 column 1 of row 3
		## allocate figure 5 column 2 of row 3
		layout(matrix(c(1,1,2,3,4,5),3,2,byrow = TRUE),
		       widths=c(4,2), heights=c(1,2,2),TRUE)
		par (oma=c(1,2,1,1))

	    } else {
		par (mfrow = c(2, 2),pin=c((3.2/plot.asp),3.2),oma=c(0,5,0,0))
	    }
	}

	# load our data
	data.flor <- .load.flo (file.path (data.path, file.flor))
	data.desc <- .load.abt (file.path (data.path, file.desc))

	if (firstFile==1)
	{
	    sample.tally <- data.frame (data.desc[,1], as.integer(0))
	    names (sample.tally) <- c ("Name", "count")
	    firstFile=0
	}
	else
	{
	    for (sname in 1:length(data.desc[,1]))
	    {
		sampFound=0
		for (i in 1:nrow(sample.tally))
		{
		    if (sample.tally[i,1]==data.desc[sname,1])
		    {
			sampFound=1
		    }
		}
		if (sampFound == 0)
		{
		    tallyTemp <- data.frame(data.desc[sname,1], as.integer(0))
		    names (tallyTemp) <- c ("Name", "count")
		    sample.tally=rbind(sample.tally,tallyTemp)
		}
	    }
	}

	# let's summon our sample name
	for (samp in 1:length(data.desc[,2]))
	{
	    sample.well <- data.desc[samp,2]
	    sample.row = chr(asc('A') + rowFromWell(sample.well))
	    sample.col = as.character(colFromWell(sample.well))

	    sample.name <- data.desc[samp, 1]
	    print(sample.name)

	    if (sample.name==""||
		substring(sample.name,1,6)=="Sample" ||
		substring(sample.name,1,3)=="NTC" ||
		substring(sample.name,1,3)=="ntc" ||
		substring(sample.name,1,1)==" ")
	    {
		print("skipped..")
		if (! is.null(progressFunc)) {
		    fileProg = (as.double(samp)/length(data.desc[,2]))
		    progressFunc((as.double(fileind-1) / length(files))+(fileProg / length(files)))
		}
		next
	    }
	    # skip Pos/Neg controls unless generating training data 
	    # i.e. not using a modle and generating (stats file)
	    if ((statsFile == "" || useLasso) &&
	        (sample.name=="Pos" || sample.name=="Neg"|| 
		substring(sample.name,1,3)=="NEG" ||
		substring(sample.name,1,3)=="POS" ||
		substring(sample.name,1,3)=="neg" ||
		substring(sample.name,1,3)=="pos"))
	    {
		print("skipped..")
		if (! is.null(progressFunc)) {
		    fileProg = (as.double(samp)/length(data.desc[,2]))
		    progressFunc((as.double(fileind-1) / length(files))+(fileProg / length(files)))
		}
		next
	    }

	    # increment sample counts (typically there are 2 per name)
	    sample.copy=0
	    for (i in 1:nrow(sample.tally))
	    {
		if (sample.tally[i,1]==sample.name)
		{
		    sample.tally[i,2]=sample.tally[i,2]+1
		    sample.copy <- sample.tally[i,2]
		    print(sprintf("sample: %s, copy %d",sample.name,sample.copy))
		}
	    }
	    # note that selected samples are printed one to a file
	    if (format == "png" || selectList != "" || sampleMask != "")
	    {
		sample.selected<-0
		if (selectList != "")
		{
		    for (i in 1:nrow(selected))
		    {
			if ( selected[i,1] == paste(sample.copy,sample.name,sep="_"))
			{
			    sample.selected=1
			}
		    }
		}
		if (sampleMask != "" && 
		    regexpr(sampleMask,c(paste(sample.copy,sample.name,sep="_")))[1] >= 0)
		{
		    sample.selected=1
		}

		if ((selectList != "" || sampleMask != "") && sample.selected==0)
		{
		    print("skipped..")
		    if (! is.null(progressFunc)) {
			fileProg = (as.double(samp)/length(data.desc[,2]))
			progressFunc((as.double(fileind-1) / length(files))+(fileProg / length(files)))
		    }
		    next
		}
		
		plot.name=paste(sample.copy,gsub("/","_",sample.name),sep="_") 
		if (format == "pdf")
		{
		    pdf(paste(plot.name,"pdf",sep="."),width=8.5)
		} else if (format == "png")
		{
		    png(paste(plot.name,"png",sep="."),width = 8.5,height=8.5,units="in",res=300)
		}
		if (format != "")
		{
		    if (plot.asp > 1)
		    {
			layout(matrix(c(1,1,2,3,4,5),3,2,byrow = TRUE),
			   widths=c(4,2), heights=c(1,2,2),TRUE)

			par (oma=c(1,2,1,1))
			#par (mfrow = c(2, 2),pin=c(3.2,(3.2/plot.asp)),oma=c(0,5,0,0))
		    } else {
			par (mfrow = c(2, 2),pin=c((3.2/plot.asp),3.2),oma=c(0,5,0,0))
		    }
		}
	    }

	    # find nth copy in DivMelt scores (if available)
	    DivMelt.index = 0
	    if (file.DivMelt != "")
	    {
		for (sname in 1:length(newdat[,1]))
		{
			if (newdat[sname,1] == paste(sample.copy,sample.name,sep="_"))
			{
			    if ( newdat[sname,6]== sample.row &&
				newdat[sname,7] == sample.col )
			    {
				DivMelt.index = sname
				break
			    }
			    else
				print(paste("Failed to match row/col in DivMelt score file for sample ",
				    newdat[DivMelt.index,1]," (Row=",newdat[sname,6]," vs. ",
				    sample.row,",Col=",newdat[sname,7]," vs. ",sample.col,")",sep=""))
			}
		}
		if (DivMelt.index != 0)
		    print(sprintf("Matched against %s copy number %d (index %d)",
			sample.name, sample.copy, DivMelt.index))
		else
		    print(paste("Failed to find",sample.name,", copy",sample.copy,"in DivMelt score file"))
	    }

	    print("Reading data...")
	    # create our data of interest
	    data.sample <- data.frame (data.flor[data.flor[,1] == data.desc[samp,2], 6], 
				       data.flor[data.flor[,1] == data.desc[samp,2], 7])
	    names (data.sample) <- c ("Temperature", "Flo1")

	    # --------------------------------------------------------------------------- #
	    # set up and derive our working data
	    # --------------------------------------------------------------------------- #
	    # set up our working data
	    data.orig <- data.frame (data.sample[,1], data.sample[,2])

	    # plot
	    if (format != "")
	    {
		print("Plotting Fluorescence...")

		par(mar=c(0,0,0,0))  
		plot.new()
		fileString=paste("File Name:",files[fileind],sep=" ")
		if (DivMelt.index != 0)
		    fileString = paste(fileString,
			"\n\nDivMelt Score file:",scoreFile,sep=" ")
		fileString= paste(fileString,"\n\nSample Well#:",sample.well,
		    "- Row:",sample.row,",",
		    "Col:",sample.col,sep=" ")

		text(0.5,0.5,fileString,cex=1.5)

		par(mar=c(4,4,4,1))  
		if (xrangeOverride)
		    plot (data.orig, type='l', xlab = "Temperature (C)", ylab = "Fluorescence",
			xlim=range(c(xrangemin,xrangemax)),
			main=paste(sample.copy,"_",sample.name,sep=""))
		else
		    plot (data.orig, type='l', xlab = "Temperature (C)", ylab = "Fluorescence",
			main=paste(sample.copy,"_",sample.name,sep=""))

		if (showFlorSmooth)
		    lines(lowess(data.orig[,1],data.orig[,2]),col="red")   #fit with lowess function
		par(mar=c(0,0,0,0))  
		plot.new()
		if (format == "pdf")
		{
		    settingStr= paste("Theta 1: ",as.character(theta1),
			"\nT1 hold interval: ",as.character(t1Cutoff), 
			"\nT1 slope window: ",as.character(t1SlopeWindow),"\n", sep="")
		    settingStr= paste(settingStr,"Theta 2: ",as.character(theta2),
			"\nT2 hold interval: ",as.character(t2Cutoff), 
			"\nT2 slope window: ",as.character(t2SlopeWindow),"\n", sep="")
		    if (excludeEarly)
		    {
			settingStr=paste(settingStr,
			    sprintf("Exclude early pk. CO: %.1f\n",earlyCutoff),sep="")
		    }
		    else
		    {
			settingStr=paste(settingStr, "Exclude early peaks: Off\n",sep="")
		    }
		    if (excludeLate)
		    {
			settingStr=paste(settingStr,
			    sprintf("Exclude late pk. CO: %.1f\n",lateCutoff),sep="")
		    }
		    else
		    {
			settingStr=paste(settingStr, "Exclude late peaks: Off\n",sep="")
		    }
		    if (includeShoulders)
		    {
			settingStr=paste(settingStr,
			    sprintf("Shoulder ht. threshold: %.1f %%\n",
				shoulderCutoff * 100.),sep="")
		    }
		    else
		    {
			settingStr=paste(settingStr, "Include shoulders: Off\n",sep="")
		    }
		    text(0.2,0.9,"Settings",adj=0,cex=1.5)
		    text(0.2,0.5,settingStr,adj=0,cex=1.2)
		}
	    }
	    print("Gathering stats on data...")

	    model1 = lm(data.orig[,2]~data.orig[,1])  #simple regression 
	    sigma=summary(model1)$sigma   # basic statistical summary

	    minx=min(data.orig[,1])
	    maxx=max(data.orig[,1])
	    miny=min(data.orig[,2])
	    maxy=max(data.orig[,2])
	    meany=mean(data.orig[,2])
	    stddevy=sd(data.orig[,2])

	    # generate the negative derivative values
	    y.scaled <-  .scale_linear (data.orig[,2])
	    dy.orig <- .derive_flo (y.scaled,data.orig[,1])
	    data.derived <- data.frame (data.orig[,1], dy.orig)

	    if (format != "")
	    {
		print("Plotting dFluorescence...")
		par(mar=c(4,4,4,1))  
		if (xrangeOverride)
		    plot (x=data.orig[,1], y=dy.orig, type='l', 
			xlim=range(c(xrangemin,xrangemax)),
			xlab="Temperature (C)", ylab="-dFluorescence/dTemp")
		else
		    plot (x=data.orig[,1], y=dy.orig, type='l', xlab="Temperature (C)", ylab="-dFluorescence/dTemp")
	    }

	    # smooth the derivative curve
	    x.work=data.orig[,1]
	    dy.loess <- loess (y ~ x, span=loess.span, data.frame (x=data.orig[,1], y=dy.orig))
	    dy.predict <- stats::predict (dy.loess, data.frame (x=x.work))

	    if (format != "" && showDFlorSmooth)
	    {
		lines (x.work, dy.predict, col="red")
	    }

	    dminy=min(dy.predict)
	    dmaxy=max(dy.predict)
	    dmeany=mean(dy.predict)
	    dstddevy=sd(dy.predict)

	    scaleF=(maxx-minx)/(dmaxy-dminy)
	    cutoff1=(tan(theta.low.rad) * plot.asp)/ scaleF
	    altT1s<- vector(length=4)
	    for (i in 1:4)
		altT1s[i]=0
	    altT1Cnt=0

	    cutoff2=(tan(theta.high.rad) * plot.asp)/ scaleF
	    altT2s<- vector(length=4)
	    for (i in 1:4)
		altT2s[i]=0
	    altT2Cnt=0


	    # find "peak" intensity (uses 20 sample window)
	    scan.upi=0
	    scan.upstop=0
	    scan.downi=0
	    scan.maxspread=0
	    for (i in 1:length(data.orig[,1]))
	    {
		if (dmaxy==dy.predict[i])
			scan.peaki=i
	    }
	    for (i in 11:(length(data.orig[,1])-11))
	    {
		if (dy.predict[i] > dy.predict[i-10] && 
		    dy.predict[i] > dy.predict[i+10] &&
		    (dy.predict[i] - dy.predict[i-10]) + 
		    (dy.predict[i] - dy.predict[i+10]) > scan.maxspread)
		{
		    scan.upstop=i
		    scan.maxspread=(dy.predict[i] - dy.predict[i-10]) +
				(dy.predict[i] - dy.predict[i+10])
		}
	    }
	    if (scan.upstop==0)
		scan.upstop=scan.peaki

	    # big peak error means sloping graoph (peak not at max dFlores)
	    else if (abs(data.orig[scan.upstop,1] - data.orig[scan.peaki,1]) > 1.)
		scan.peaki=scan.upstop

	    # fill histograms on each side of peak
	    scan.hist1=vector(length=plot.bins)
	    scan.hist2=vector(length=plot.bins)
	    for (j in 1:plot.bins)
	    {
		scan.hist1[j]=0
		scan.hist2[j]=0
	    }
	    scan.start=1
	    for (i in 1:(length(data.orig[,1])))
	    {
		if (data.orig[i,1] < plot.minTemp)
		    scan.start = i
	    }

	    dmaxy1=max(dy.predict[1:scan.upstop])
	    dmaxy2=max(dy.predict[scan.upstop:length(dy.predict)])
	    dminy1=min(dy.predict[scan.start:scan.upstop])
	    dminy2=min(dy.predict[scan.upstop:length(dy.predict)])
	    scan.range1=dmaxy-dminy1
	    scan.range2=dmaxy-dminy2

	    dminy1i=1;

	    scan.high1=dy.predict[1]

	    # this will be set to index of 1st value in lowest downslope bin
	    scan.defaultDn = length(data.orig[,1])

	    # this will be set to index of 1st value in 2nd upslope bin
	    scan.defaultUp = scan.start

	    for (i in scan.start:(length(data.orig[,1])))
	    {
		if (i <= scan.upstop)
		{
		    if (dy.predict[i] == dminy1)
			dminy1i=i;

		    if (scan.high1<dy.predict[i])
			scan.high1=dy.predict[i]

		    #if (scan.high1-dy.predict[i] > 0.1*scan.range1 && 
		    #	dy.predict[i] > dminy1+(scan.range1/plot.bins) &&
		    #	scan.upstop==scan.peaki)
		    #	scan.upstop=i

		    for (j in 1:plot.bins)
			if (dy.predict[i] >= dminy1+(j-1)*(scan.range1/plot.bins) && dy.predict[i] <= dminy1+j*(scan.range1/plot.bins))
			{
			    if (j==2 && scan.hist1[j] == 0)
				scan.defaultUp = i
			    scan.hist1[j] = scan.hist1[j] + 1
			}
			
		}
		else
		{
		    for (j in 1:plot.bins)
			if ((dy.predict[i] >= dminy2+(j-1)*(scan.range2/plot.bins)) && (dy.predict[i] <= dminy2+j*(scan.range2/plot.bins)))
			{
			    if (j==1 && scan.hist2[j] == 0)
				scan.defaultDn = i
			    scan.hist2[j] = scan.hist2[j]+1
			}
		}
	    }

	    # identify biggest histogram bins on either side of peak
	    scan.maxbin1=0
	    scan.maxbin2=0
	    scan.maxcnt1=0
	    scan.maxcnt2=0
	    for (j in 1:plot.bins)
	    {
		if (scan.hist1[j] > scan.maxcnt1)
		{
		    scan.maxbin1=j
		    scan.maxcnt1=scan.hist1[j]
		}
		if (scan.hist2[j] > scan.maxcnt2)
		{
		    scan.maxbin2=j
		    scan.maxcnt2=scan.hist2[j]
		}
	    }
	    scan.secondbin=1
	    if (format == "pdf")
	    {
		segments(data.orig[scan.peaki,1], dy.predict[scan.peaki]-2*cutoff1,
		    data.orig[scan.peaki,1], dy.predict[scan.peaki]+2*cutoff1,
		   col = "purple")
		segments(data.orig[scan.peaki,1]-2, dy.predict[scan.peaki],
		    data.orig[scan.peaki,1]+2, dy.predict[scan.peaki],
		   col = "purple")
	    }

	    scan.upHeight = 0.  # record "Shoulder" height in case it is to be ignored
	    scan.settingHeight = FALSE
	    scan.slopeThis = 0.
	    scan.intrough = FALSE
	    scan.dnlock = FALSE
	    scan.lastBelow = 0
	    scan.currentT1 = 0
	    scan.preT2=0

	    print(t1SlopeWS/2)
	    print(t2SlopeWS/2)
	    side1=max(1,t1SlopeWS/2)
	    side2=max(1,t2SlopeWS/2)
	    for (i in max(scan.start,10):(length(data.orig[,1])-side2))
	    {

		# this section is for "rising" analysis
		if (i < scan.upstop)
		{
		    slope=(dy.predict[i+side1]-dy.predict[i-side1])/(data.orig[i+side1,1]-data.orig[i-side1,1])
		    # record first sample in 2nd histogram bucket (used
		    # if for fallback T1 if no theta1 angle is found
		    if (scan.secondbin==1 && 
			dy.predict[i] > dminy1+(scan.range1/plot.bins)
			&& i > dminy1i)
		    {
			scan.secondbin=i
		    }

		    if (slope < cutoff1) {
			scan.lastBelow=i
		    }

		    if ((scan.lastBelow != 0) && 
			(data.orig[i,1] - data.orig[scan.lastBelow,1] > t1Cutoff) &&
			(scan.currentT1 !=scan.lastBelow+1) &&
			(slope > cutoff1))
		    {
			scan.currentT1=scan.lastBelow+1
			if (altT1Cnt < 4)
			{
			    altT1Cnt=altT1Cnt+1
			    altT1s[altT1Cnt]=scan.currentT1
			}
			if (format == "pdf")
			    segments(data.orig[scan.currentT1,1], dy.predict[scan.currentT1]-2*cutoff1,
				data.orig[scan.currentT1,1], dy.predict[scan.currentT1]+2*cutoff1,
			       col = "green")

			if (!scan.settingHeight && includeShoulders)
			    print(paste("current Shoulder height ",scan.upHeight,
			    "vs. ",shoulderCutoff * scan.range1))

			# if already watching shoulder height, then hold off on reset
			if (!scan.settingHeight &&
			    (scan.upi == 0 || 
			    !includeShoulders || 
			    scan.upHeight < (shoulderCutoff * scan.range1)))
			{
				scan.upi=scan.currentT1
				scan.upHeight = 0
				scan.settingHeight = TRUE
			   print(paste("Shoulder start at t=",data.orig[scan.currentT1,1],"start val=",
				dy.predict[scan.currentT1]))
			}
			scan.uptrough = 0
		    }

		    # watch slope to know when to stop checking height
		    if (scan.upi != 0 && scan.upi != i && scan.settingHeight)
		    {
			scan.slopeThis=(dy.predict[i+t1SlopeWS/2]-dy.predict[i-t1SlopeWS/2])/(data.orig[i+t1SlopeWS/2,1]-data.orig[i-t1SlopeWS/2,1])

			#record height of peak
			if ( scan.slopeThis > 0. )
			{
			    scan.upHeight = max(scan.upHeight, dy.predict[i] - dminy1)
			}
			else 
			{
			   scan.settingHeight = FALSE
			   print(paste("Shoulder end at t=",data.orig[i,1]," height=",
				scan.upHeight,"(",
				(scan.upHeight/scan.range1)*100.,"%)",sep=""))
			}
		    }
		    # if enabled,
		    # watch for trough after initial peak 
		    # (resets after specified degree trough)
		    if (excludeEarly && scan.upi != 0 )
		    {
			if (dy.predict[i] <= dminy1+(0.1)*scan.range1)
			{
			    # grow the trough
			    scan.uptrough = scan.uptrough + data.orig[i,1]-data.orig[i-1,1]
			}
			else
			    scan.uptrough = 0

			# erase previous T1 if we get a "large" trough
			if (scan.uptrough > scan.uptrough_limit)
			{
			    scan.upi = 0
			    scan.upHeight = 0
			}
		    }
		}
		# This section is for "falling" analysis
		else
		{
		    slope=(dy.predict[i+side2]-dy.predict[i-side2])/(data.orig[i+side2,1]-data.orig[i-side2,1])
		    # if enabled,
		    # watch for trough after initial peak 
		    # (resets after specified degree trough)
		    if (excludeLate && scan.downi != 0 && ! scan.dnlock )
		    {
			if (dy.predict[i] <= dminy2+(0.1)*scan.range2)
			{
			    # grow the trough
			    if (scan.intrough)
				scan.dntrough = scan.dntrough + data.orig[i,1]-data.orig[i-1,1]
			    else
			    {
				scan.dntrough = 0
			        scan.intrough = TRUE
			    }
			}

			# erase previous T2 if we get a "large" trough
			# followed by another shoulder
			if (scan.dntrough > scan.dntrough_limit)
			{
			    scan.dnlock = TRUE
			}

		    }
		    if (slope < cutoff2) {
			scan.preT2=i
		    }

		    if ((scan.preT2 != 0) && 
			(data.orig[i,1] - data.orig[scan.preT2,1] > t2Cutoff) &&
			(slope > cutoff2))
		    {
			if (altT2Cnt < 4 &&
			    (altT2Cnt==0 || altT2s[altT2Cnt] != scan.preT2+1))
			{
			    altT2Cnt=altT2Cnt+1
			    altT2s[altT2Cnt]=scan.preT2+1
			}
			if (format == "pdf")
			    segments(data.orig[scan.preT2+1,1], dy.predict[scan.preT2+1]-2*cutoff2,
				data.orig[scan.preT2+1,1], dy.predict[scan.preT2+1]+2*cutoff2,
			       col = "green")
			if ((dy.predict[scan.preT2+1] >= dminy2+(scan.maxbin2-1)*(scan.range2/plot.bins)) && 
				(dy.predict[scan.preT2+1] <= dminy2+(scan.maxbin2+1)*(scan.range2/plot.bins)))
			{
				if (!scan.dnlock)
				    scan.downi=scan.preT2+1
			}
		    }
		}
	    }

	    scan.error=""
	    errorCode=0
	    if (scan.upi == 0) {
		if (scan.secondbin == 1)
		    scan.secondbin=dminy1i

		scan.upi = scan.defaultUp
		scan.error = paste("Warning: Could not find T1 theta = ",as.character(theta1)," for sample ",sample.name)
		print(scan.error)
		errorCode=1
	    }
	    if (format != "")
	    {
		abline(v=data.orig[scan.upi,1])
		segments(data.orig[scan.upi,1]-2, dy.predict[scan.upi]-2*cutoff1,
			data.orig[scan.upi,1]+2, dy.predict[scan.upi]+2*cutoff1,
		       col = "blue")
	    }

	    if (scan.downi == 0) {
		scan.downi = scan.defaultDn
		scan.error=paste("Warning: could not find T2 theta = ",as.character(-theta2)," for sample ",sample.name)
		print(scan.error)
		errorCode=2
	    }
	    if (format != "")
	    {
		abline(v=data.orig[scan.downi,1])
		segments(data.orig[scan.downi,1]-2, dy.predict[scan.downi]-2*cutoff2,
			data.orig[scan.downi,1]+2, dy.predict[scan.downi]+2*cutoff2,
		       col = "blue")

		if (DivMelt.index != 0)
		{
		    abline(v=as.double(newdat[DivMelt.index,2]),col="gray")
		    abline(v=as.double(newdat[DivMelt.index,3]),col="gray")
		}
	    }

	    print(sprintf("Score = %.2f",data.orig[scan.downi,1]-data.orig[scan.upi,1]))

	    statsnum=1
	    stats[1]=paste(sample.copy,"_",sample.name,sep="")
	    statsnum=statsnum+1
	    stats[statsnum]=sample.row
	    statsnum=statsnum+1
	    stats[statsnum]=sample.col
	    statsnum=statsnum+1
	    stats[statsnum]=as.character(miny)
	    statsnum=statsnum+1
	    stats[statsnum]=as.character(maxy)
	    statsnum=statsnum+1
	    stats[statsnum]=as.character(meany)
	    statsnum=statsnum+1
	    stats[statsnum]=as.character(stddevy)
	    statsnum=statsnum+1
	    stats[statsnum]=as.character(sigma)
	    statsnum=statsnum+1
	    stats[statsnum]=as.character(dminy)
	    statsnum=statsnum+1
	    stats[statsnum]=as.character(dmaxy)
	    statsnum=statsnum+1
	    stats[statsnum]=as.character(dmeany)
	    statsnum=statsnum+1
	    stats[statsnum]=as.character(dstddevy)
	    statsnum=statsnum+1

	    if (statsFile != "")
	    {
		write (stats,file=stats.out,ncolumns=statsnum,sep=",",append=TRUE)
	    }

	    if (useLasso)
	    {
		pred.data<-matrix(nrow=1,ncol=9)

		# set up our working data
		for(i in 4:12)
		{
		    pred.data[1,i-3]=as.numeric(stats[i])
		}

		print("running prediction...")
		print(pred.data)

		preds=predict(model,pred.data,s="lambda.min")

		print("ran prediction...")

		print(sprintf("Pred = %.2f",preds[1]))

		if (preds[1] < 0) {
		    suspect.list[1]=as.character(stats[1])
		    suspect.list[2]=as.character(preds[1])
		    write (suspect.list,file=suspect.out,ncolumns=2,sep=",",append=TRUE)
		}
	    }

	    if ((useLasso) && (preds[1] < 0))
	    {
		title.suffix="- Warning: curve characteristics rejected by Lasso regression"
		print(title.suffix)
		errorCode=3
	    }
	    else if (scan.error != "")
	    {
		title.suffix=scan.error
	    }
	    else if (DivMelt.index != 0)
	    {
		title.suffix= sprintf("(vs. %.2f)",as.double(newdat[DivMelt.index,5]))
	    }
	    else
		title.suffix=""

	    if (format != "")
	    {
		title(sprintf("Score = %.2f%s",
		    data.orig[scan.downi,1]-data.orig[scan.upi,1],
		    title.suffix))

		par(mar=c(0,0,0,0))  
		plot.new()
		if (format=="pdf" && scan.error == "")
		{
		    resStatStr1 = sprintf("T1:    %.2f",data.orig[scan.upi,1])
		    if (DivMelt.index != 0)
			resStatStr1 = paste(resStatStr1,
			    sprintf("(vs. %.2f)",as.double(newdat[DivMelt.index,2])),sep=" ")
		    resStatStr1 = paste(resStatStr1,
			    sprintf("\nT2:    %.2f",data.orig[scan.downi,1]))
		    if (DivMelt.index != 0)
			resStatStr1 = paste(resStatStr1,
			    sprintf("(vs. %.2f)",as.double(newdat[DivMelt.index,3])),sep=" ")

		    resStatStr2 = sprintf("Peak:  %.2f",data.orig[scan.peaki,1])
		    if (DivMelt.index != 0)
			resStatStr2 = paste(resStatStr2,
			    sprintf("(vs. %.2f)",as.double(newdat[DivMelt.index,4])),sep=" ")

		    trueT1Alts=0
		    resStatStr3=""
		    if (altT1Cnt >= 1)
		    {
			printed=0
			for(i in 1: altT1Cnt)
			{
			    if (altT1s[i] != scan.upi && printed < 3) #true for all but 1
			    {
				resStatStr3=paste(resStatStr3,
				    sprintf("\nAltT1%s:  %.2f",chr(asc('a') +trueT1Alts),
					data.orig[altT1s[i],1]))
				trueT1Alts=trueT1Alts+1
				printed=printed+1
			    }
			}
		    }
		    trueT2Alts=0
		    if (altT2Cnt >= 1)
		    {
			printed=0
			for(i in 1: altT2Cnt)
			{
			    if (altT2s[i] != scan.downi && printed < 3) #true for all but 1
			    {
				resStatStr3=paste(resStatStr3,
				    sprintf("\nAltT2%s:  %.2f",chr(asc('a') +trueT2Alts),
					data.orig[altT2s[i],1]))
				trueT2Alts=trueT2Alts+1
				printed=printed+1
			    }
			}
		    }

		    if (DivMelt.index != 0)
			text(0.2,0.9,"Results (vs. score file)", adj=0,cex=1.5,col="black")
		    else
			text(0.2,0.9,"Results", adj=0,cex=1.5,col="black")

		    text(0.2,0.75,resStatStr1, adj=0,cex=1.2,col="blue")
		    text(0.2,0.6,resStatStr2, adj=0,cex=1.2,col="purple")
		    if (resStatStr3 != "")
			text(0.2,0.6-(0.05 * (trueT1Alts+trueT2Alts)),resStatStr3, adj=0,cex=1.2,col="green")

		}

		if (format == "png" || selectList != "" || sampleMask != "")
		    dev.off()
	    }

	    outmeasure=1
	    out[outmeasure]=gsub(".FLO","",files[fileind])
	    outmeasure=outmeasure+1
	    out[outmeasure]=paste(sample.copy,"_",sample.name,sep="")
	    outmeasure=outmeasure+1
	    out[outmeasure]=sample.row
	    outmeasure=outmeasure+1
	    out[outmeasure]=sample.col
	    outmeasure=outmeasure+1
	    out[outmeasure]=as.character(data.orig[scan.upi,1])
	    outmeasure=outmeasure+1
	    out[outmeasure]=as.character(data.orig[scan.downi,1])
	    outmeasure=outmeasure+1
	    out[outmeasure]=as.character(data.orig[scan.peaki,1])
	    outmeasure=outmeasure+1
	    out[outmeasure]=as.character(data.orig[scan.downi,1]-data.orig[scan.upi,1])
	    outmeasure=outmeasure+1
	    out[outmeasure]=theta1
	    outmeasure=outmeasure+1
	    out[outmeasure]=theta2
	    outmeasure=outmeasure+1
	    out[outmeasure]=includeShoulders
	    outmeasure=outmeasure+1
	    out[outmeasure]=shoulderCutoff
	    outmeasure=outmeasure+1
	    out[outmeasure]=excludeEarly
	    outmeasure=outmeasure+1
	    out[outmeasure]=earlyCutoff
	    outmeasure=outmeasure+1
	    out[outmeasure]=excludeLate
	    outmeasure=outmeasure+1
	    out[outmeasure]=lateCutoff
	    outmeasure=outmeasure+1
	    out[outmeasure]=t1Cutoff
	    outmeasure=outmeasure+1
	    out[outmeasure]=t2Cutoff
	    outmeasure=outmeasure+1
	    out[outmeasure]=t1SlopeWindow
	    outmeasure=outmeasure+1
	    out[outmeasure]=t2SlopeWindow
	    outmeasure=outmeasure+1
	    if (scan.upi==scan.defaultUp)
	    {
		for(i in 1: 3)
		{
		    out[outmeasure]=0
		    outmeasure=outmeasure+1
		}
	    }
	    else
	    {
		printed=0
		for(i in 1: 4)
		{
		    if (altT1s[i] != scan.upi && printed < 3) #true for all but 1
		    {
			if (altT1s[i] != 0)
			    out[outmeasure]=as.character(data.orig[altT1s[i],1])
			else
			    out[outmeasure]=0
			outmeasure=outmeasure+1
			printed=printed+1
		    }
		}
	    }
	    if (scan.downi==scan.defaultDn)
	    {
		for(i in 1: 3)
		{
		    out[outmeasure]=0
		    outmeasure=outmeasure+1
		}
	    }
	    else
	    {
		printed=0
		for(i in 1: 4)
		{
		    if (altT2s[i] != scan.downi && printed < 3) #true for all but 1
		    {
			if (altT2s[i] != 0)
			    out[outmeasure]=as.character(data.orig[altT2s[i],1])
			else
			    out[outmeasure]=0
			outmeasure=outmeasure+1
			printed=printed+1
		    }
		}
	    }

	    out[outmeasure]=errorCode
	    outmeasure=outmeasure+1

	    if (file.DivMelt != "")
	    {
		if (DivMelt.index != 0)
		{
		    out[outmeasure]=newdat[DivMelt.index,1]
		    outmeasure=outmeasure+1
		    out[outmeasure]=newdat[DivMelt.index,2]
		    outmeasure=outmeasure+1
		    out[outmeasure]=newdat[DivMelt.index,3]
		    outmeasure=outmeasure+1
		    out[outmeasure]=newdat[DivMelt.index,4]
		    outmeasure=outmeasure+1
		    out[outmeasure]=newdat[DivMelt.index,5]
		    outmeasure=outmeasure+1
		    out[outmeasure]=newdat[DivMelt.index,6]
		    outmeasure=outmeasure+1
		    out[outmeasure]=newdat[DivMelt.index,7]
		    outmeasure=outmeasure+1
		    out[outmeasure]=as.character(data.orig[scan.upi,1]-as.double(newdat[DivMelt.index,2]))
		    outmeasure=outmeasure+1
		    out[outmeasure]=as.character(data.orig[scan.downi,1]-as.double(newdat[DivMelt.index,3]))
		    outmeasure=outmeasure+1
		    out[outmeasure]=as.character(data.orig[scan.peaki,1]-as.double(newdat[DivMelt.index,4]))
		    outmeasure=outmeasure+1
		    out[outmeasure]=as.character(abs((data.orig[scan.downi,1]-data.orig[scan.upi,1])-as.double(newdat[DivMelt.index,5])))
		    outmeasure=outmeasure+1
		}
		else
		{
		    for (m in 1:extraplaces)
		    {
			out[outmeasure]=""
			outmeasure=outmeasure+1
		    }
		}
	    }
	    if (useLasso)
	    {
		out[outmeasure] = preds[1] >= 0
		outmeasure=outmeasure+1
	    }
	    write (out,file=data.out,ncolumns=outmeasure,sep=",",append=TRUE)

	    if (! is.null(progressFunc)) {
		fileProg = (as.double(samp)/length(data.desc[,2]))
		progressFunc((as.double(fileind-1) / length(files))+(fileProg / length(files)))
	    }
	}
	if (format == "pdf" && selectList=="" && sampleMask=="")
	{
	    dev.off()
	}
    }
}
