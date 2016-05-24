# Error warnings from pam
cor.options <- list(debug=TRUE, #whether to turn on debugging or not
                    err.file=ifelse(.Platform$OS.type=="windows",
                      paste(tempdir(), "cortrace.txt", sep=.Platform$file.sep), "cortrace.txt"),
                    image.file=ifelse(.Platform$OS.type=="windows",
                      paste(tempdir(), "corimage.Rdata", sep=.Platform$file.sep),"corimage.Rdata"),                     
                    reserved.class.label="Unspecified")

##
## Our error handler
##
cor.xl.error.trace <- function() {
  err.message <- geterrmessage()
  if (!is.null(cor.options$image.file)) {
    save.image(cor.options$image.file)
  }
  if (!is.null(cor.options$err.file)) {
    sink(cor.options$err.file)
    print(err.message)
    traceback()
    sink()
  }
  winDialog(type="ok", message=err.message)
}

##
## Upon loading, if we are in a windows environment, we use the windows
## dialog mechanism to display errors. Useful for debugging COM apps
##
.onLoad <- function(lib, pkg) {

# Rob changed this next line on  apr 10, 2005, requested by Uwe Ligges

#  if ( .Platform$OS.type == "windows"  ) {
  if ( .Platform$OS.type == "windows" && interactive() ) {
    options(error=cor.xl.error.trace)
  }

}

##
## Upon unload, we set things back the way they were...
##
.onUnload <- function(libpath){
  if ( .Platform$OS.type == "windows") {
    options(error=NULL)
  }
}






#Functions that get called:
#
#cor.xl.LoadData
#cor.xl.RunAnalysis
#cor.xl.Plot1, cor.xl.Plot2, cor.xl.Plot3
#


# So I dont have to worry about scope or user choices for values
cor.xl.glob.assign = function(x, val){
	assign(paste("cor.", x, sep=""), val, pos=.GlobalEnv)
}

cor.xl.glob.get = function(x){
	return(get(paste("cor.", x, sep="")))
}	


# Takes the row labels and the data and preps it into blocks
cor.xl.GetBlocks = function(data, labels, which){
	blockcol = labels[,3]
	labels = labels[,-3]
	
	#save data and labels as blocks
	blocks = unique(blockcol)
	for(block in blocks){
		cor.xl.glob.assign(paste("block",which,block,sep="."), data[blockcol == block,])
		cor.xl.glob.assign(paste("labels",which,block,sep="."), labels[blockcol == block,1:2])
	}
	
	#save list of blocks
	cor.xl.glob.assign(paste("blocklist",which,sep="."), sapply(blocks, function(x){paste("block",which,x,sep=".")}, USE.NAMES=F))
	
	#save blockcol ordered by blocks
	cor.xl.glob.assign(paste("blockcol",which,sep="."), Reduce(function(x, y){c(x, rep(y, length(blockcol[blockcol==y])))}, blocks, init=NULL))
	
}

# Converts block names form blocklist to label names
cor.xl.Blocks2Labels = function(blocks){
	labels = vector(mode="character", length=length(blocks))
	
	which = substr(blocks[1], 7, 7)
	for(i in 1:length(blocks)){
		
		# 9 is for the letter after "block.#."
		name = substr(blocks[i], 9, 1000000L)
		labels[i] = paste("labels",which,name,sep=".")
	}
	
	return(labels)
}

# Used to combine blocks (of data or labels)
cor.xl.CombineBlocks = function(blocks){
	full = cor.xl.glob.get(blocks[1])
	for(block in blocks[-1]){
		full = rbind(full, cor.xl.glob.get(block))
	}
	return(full)
}

cor.xl.ImputeData = function(which){
	cor.xl.glob.assign(paste("datawasimputed",which,sep="."), F)
	blocks = cor.xl.glob.get(paste("blocklist",which,sep="."))
	for(block in blocks){
		bl.data = cor.xl.glob.get(block)
		if(sum(is.na(bl.data)) > 0){
			
			#check if enough rows to use knn 
			if(length(bl.data[,1]) > 49) {
				bl.data = impute.knn(bl.data)$data
			}
			
			#else use row mean
			else {
				tmp = rowMeans(bl.data, na.rm=T) %*% 
					  matrix(1,nrow=1,ncol=ncol(bl.data))
				bl.data[is.na(bl.data)] = tmp[is.na(bl.data)]
			}
			
			cor.xl.glob.assign(block, bl.data)
			cor.xl.glob.assign(paste("datawasimputed.",which,sep=""), T)
		}
	}	
}

cor.xl.LoadData = function(bl.col, smp.row, which){

	data = cor.xl.glob.get(paste("rawdata",which,sep="."))
	labels = cor.xl.glob.get(paste("rawlabels",which,sep="."))
	
	#If it is the first dataset, load in the sample labels.
	# Else, arrange the cols to match the first dataset
	smplabels = cor.xl.glob.get(paste("header",which,sep="."))[smp.row, -(1:(if(bl.col){3}else{2}))]
	if(which == 1){
		cor.xl.glob.assign("samplelabels", smplabels)
	}
	else{
		matched = match(cor.xl.glob.get("samplelabels"), smplabels)
		if(length(matched) != length(smplabels)){
			return("Sample labels don't match up.  Sample labels must be the same for both datasets. You might have misspecified a sample label row
or neglected to indicate the presence of a  Block Column.")	
		}
		data = data[,matched]
	}
	
	
	# If there is a blocking col split data/labels up into blocks
	if(bl.col) {
		cor.xl.GetBlocks(data, labels, which)
	}
	
	# Otherwise just assign it all to one big block
	else {
		cor.xl.glob.assign(paste("blocklist.",which,sep=""), paste("block",which,"full dataset",sep="."))
		cor.xl.glob.assign(paste("block",which,"full dataset",sep="."), data)
		cor.xl.glob.assign(paste("labels",which,"full dataset",sep="."), labels)
		cor.xl.glob.assign(paste("blockcol",which,sep="."), F)
	}
	
	
	#Impute the data if need be.  On error, return a message
	a = try(cor.xl.ImputeData(which), TRUE)
	if(class(a) == "try-error") {
		return("An error occured while trying to impute data.  This may be caused by too much missing data.")
	}	
		
		
	#If data was imputed, then prep something so that we can paste into excel
	if(cor.xl.glob.get(paste("datawasimputed.",which,sep=""))){
		imputeddata = cor.xl.CombineBlocks(cor.xl.glob.get(paste("blocklist.",which,sep="")))
		labels = cor.xl.CombineBlocks(cor.xl.Blocks2Labels(cor.xl.glob.get(paste("blocklist.",which,sep=""))))
		if(bl.col){
			labels = cbind(labels, cor.xl.glob.get(paste("blockcol",which,sep=".")))
		}
		cor.xl.glob.assign(paste("imputeddata",which,"1",sep="."), labels)
		cor.xl.glob.assign(paste("imputeddata",which,"2",sep="."), imputeddata)
	}
	return("")
}

cor.xl.PrepforRun = function(boollst, which) {
	datalst = cor.xl.glob.get(paste("blocklist",which,sep="."))[boollst]
	x = cor.xl.CombineBlocks(datalst)
	labels = cor.xl.CombineBlocks(cor.xl.Blocks2Labels(datalst))
	# This is probably the easiest place to get back the blocking info... lets try
	bl.col = cor.xl.glob.get(paste("blockcol",which,sep="."))
	if(bl.col != F){
		allblocks = cor.xl.glob.get(paste("blocklist",which,sep="."))
		for(i in 1:length(boollst)){
			if(!boollst[i]){
				
				# 9 is for "block.#."  The other arg is just to get the rest of 
				# the string
				bl.col = bl.col[bl.col!=substr(allblocks[i], 9, 1000000L)]
			}
		}
		cor.xl.glob.assign(paste("final.blockcol",which,sep="."), bl.col)
		labels = cbind(labels,bl.col)	
	}

	#PMA needs transposed matrix
	x = t(x)
	return(list("x"=x, "labels"=labels))
}
	
	
# Might be worth having cor.xl.PrepforRun also get the row labels.  
#Also, why not have it take a boollst instead of blocklst?
cor.xl.RunAnalysis = function(boollst1, boollst2, K)	{
	prep1 = cor.xl.PrepforRun(boollst1, 1)
	prep2 = cor.xl.PrepforRun(boollst2, 2)
	mat1 = prep1$x
	mat2 = prep2$x
#        stop(ncol(mat2)) # this is fine....
	
	# Save nrows (remember that the data is transposed)
	cor.xl.glob.assign("nrows.1", ncol(mat1))
	cor.xl.glob.assign("nrows.2", ncol(mat2))
	
	
	# Set parameters for CCA/CCA.permute
	typex = cor.xl.glob.get("type.1")
	typez = cor.xl.glob.get("type.2")

        if(ncol(mat1)>nrow(mat1) && cor.xl.glob.get("unpenalized.1")) return("Do not use type Unpenalized if number of rows exceeds the number of columns.") # remember it's transposed.
        
        if(ncol(mat2)>nrow(mat2) && cor.xl.glob.get("unpenalized.2")) return("Do not use type Unpenalized if number of rows exceeds the number of columns.") # remember it's transposed.
        
	penaltyxs = if(cor.xl.glob.get("unpenalized.1")){1}#sqrt(ncol(mat1))}
				else{cor.xl.glob.get("penalty.1")}
	penaltyzs = if(cor.xl.glob.get("unpenalized.2")){1}#sqrt(ncol(mat2))}
				else{cor.xl.glob.get("penalty.2")}
				
	# Check penalties
	if(!is.null(penaltyxs)){
		if(typex=="standard" && (penaltyxs < 0 || penaltyxs >1)){
			return(paste("penalty on dataset1 must be between 0 and 1")) }
		else if(typex=="ordered" && penaltyxs<0){
			return("penalty on dataset1 must be nonnegative")}
	}
	if(!is.null(penaltyzs)){
		if(typez=="standard" && (penaltyzs < 0 || penaltyzs > 1)){
			return(paste("penalty on dataset2 must be between 0 and 1")) }
		else if(typez=="ordered" && penaltyzs<0){
			return("penalty on dataset2 must be nonnegative")}
	}
				
	const = cor.xl.glob.get("constraint.1")
	if(const == "any sign") 
		{uneg = F; upos = F}
	else if (const == "positive")
		{uneg = F; upos = T}
	else 
		{uneg = T; upos = F}
		
	const = cor.xl.glob.get("constraint.2")
	if(const == "any sign") 
		{vneg = F; vpos = F}
	else if (const == "positive")
		{vneg = F; vpos = T}
	else 
		{vneg = T; vpos = F}
		
		

	if(typex == "ordered" & cor.xl.glob.get("blockcol.1") != F)
		{chromx = prep1$labels[,3]}
	else
		{chromx = NULL}
		
	if(typez == "ordered" & cor.xl.glob.get("blockcol.2") != F)
		{chromz = prep2$labels[,3]}
	else
		{chromz = NULL}	
		
	
	set.seed(cor.xl.glob.get("randomseed"))
	
	nperms = cor.xl.glob.get("nperms")
	
	if(cor.xl.glob.get("permute")){
#          stop(paste("My name is ", sep="", paste(ncol(mat1), ncol(mat2))))
	step1 = CCA.permute(mat1, mat2, typex=typex, typez=typez, penaltyxs=penaltyxs, penaltyzs=penaltyzs, upos=upos, uneg=uneg, vpos=vpos, vneg=vneg, chromx=chromx, chromz=chromz, nperms=nperms)
	
	
	step2 = CCA(mat1, mat2, typex=typex, typez=typez, penaltyx=step1$bestpenaltyx, penaltyz = step1$bestpenaltyz, v=step1$v.init, K=K, upos=upos, uneg=uneg, vpos=vpos, vneg=vneg, chromx=chromx, chromz=chromz)
	
	
	cor.xl.glob.assign("permute.output", step1)
	cor.xl.glob.assign("penalty.1", step1$bestpenaltyx)
	cor.xl.glob.assign("penalty.2", step1$bestpenaltyz)
	}
	
	else{
	step2 = CCA(mat1, mat2, typex=typex, typez=typez, penaltyx=penaltyxs, penaltyz = penaltyzs, K=K, upos=upos, uneg=uneg, vpos=vpos, vneg=vneg, chromx=chromx, chromz=chromz)	
	cor.xl.glob.assign("penalty.1", penaltyxs)
	cor.xl.glob.assign("penalty.2", penaltyzs)
	}
	
	
	# Save values to output to excel 
	cor.xl.glob.assign("output", step2)
	cor.xl.glob.assign("numsamples", nrow(mat1))
	
	
	cor.xl.OrderRows(step2$u, prep1$labels, 1)
	cor.xl.OrderRows(step2$v, prep2$labels, 2)
	
	
	nonzero1 = vector(length=length(step2$u[1,]))
	for (i in 1:length(nonzero1)){
		nonzero1[i] = sum(step2$u[,i] != 0)	
	}
	nonzero2 = vector(length=length(step2$v[1,]))
	for (i in 1:length(nonzero2)){
		nonzero2[i] = sum(step2$v[,i] != 0)	
	}
	
	cor.xl.glob.assign("cors", step2$cors)
	cor.xl.glob.assign("nonzero.1", nonzero1)
	cor.xl.glob.assign("nonzero.2", nonzero2)
	
	return("")
}


cor.xl.Plot1 = function() {
  fileName = tempfile(fileext=".png")
  cor.xl.glob.assign("plot1File", fileName)
  png(fileName)
  plot(cor.xl.glob.get("permute.output"))
  dev.off()
		
}



cor.xl.Plot2 = function() {
	
	if(cor.xl.glob.get("blockcol.1") != F)
		{chrom = cor.xl.glob.get("final.blockcol.1")}
			#tmp = cor.xl.glob.get("final.blockcol.1")
		#chrom = match(tmp, sort(unique(cor.xl.glob.get("blockcol.1"))))
	else
		{chrom = NULL}
	#!!!!!
        fileName = tempfile(fileext=".png")
        cor.xl.glob.assign("plot2File", fileName)

          png(fileName)
	PlotCGH(cor.xl.glob.get("output")$u[,1], chrom=chrom)
        title("Component 1, Data Set 1")
	dev.off()
}

cor.xl.Plot2a = function() {

        if(cor.xl.glob.get("blockcol.1") != F)
                {chrom = cor.xl.glob.get("final.blockcol.1")}
                        #tmp = cor.xl.glob.get("final.blockcol.1")
                #chrom = match(tmp, sort(unique(cor.xl.glob.get("blockcol.1"))))
        else
                {chrom = NULL}
        #!!!!!
        fileName = tempfile(fileext=".png")
        cor.xl.glob.assign("plot2aFile", fileName)
        png(fileName)
        PlotCGH(cor.xl.glob.get("output")$u[,2], chrom=chrom)
        title("Component 2, Data Set 1")
        dev.off()
}

cor.xl.Plot2b = function() {

        if(cor.xl.glob.get("blockcol.1") != F)
                {chrom = cor.xl.glob.get("final.blockcol.1")}
                        #tmp = cor.xl.glob.get("final.blockcol.1")
                #chrom = match(tmp, sort(unique(cor.xl.glob.get("blockcol.1"))))
        else
                {chrom = NULL}
        #!!!!!

        fileName = tempfile(fileext=".png")
        cor.xl.glob.assign("plot2bFile", fileName)
          png(fileName)
        PlotCGH(cor.xl.glob.get("output")$u[,3], chrom=chrom)
        title("Component 3, Data Set 1")
        dev.off()
}


cor.xl.Plot3 = function() {
	
	if(cor.xl.glob.get("blockcol.2") != F)
		{chrom = cor.xl.glob.get("final.blockcol.2")}
			#tmp = cor.xl.glob.get("final.blockcol.2")
		#chrom = match(tmp, sort(unique(cor.xl.glob.get("blockcol.2"))))
	else
		{chrom = NULL}
	#!!!!!
        fileName = tempfile(fileext=".png")
        cor.xl.glob.assign("plot3File", fileName)
        png(fileName)
	PlotCGH(cor.xl.glob.get("output")$v[,1], chrom=chrom)
        title("Component 1, Data Set 2")
	dev.off()
		
}

cor.xl.Plot3a = function() {

        if(cor.xl.glob.get("blockcol.2") != F)
                {chrom = cor.xl.glob.get("final.blockcol.2")}
                        #tmp = cor.xl.glob.get("final.blockcol.2")
                #chrom = match(tmp, sort(unique(cor.xl.glob.get("blockcol.2"))))
        else
                {chrom = NULL}
        #!!!!!
        fileName = tempfile(fileext=".png")
        cor.xl.glob.assign("plot3aFile", fileName)
        png(fileName)
        PlotCGH(cor.xl.glob.get("output")$v[,2], chrom=chrom)
        title("Component 2, Data Set 2")
        dev.off()
}

cor.xl.Plot3b = function() {

        if(cor.xl.glob.get("blockcol.2") != F)
                {chrom = cor.xl.glob.get("final.blockcol.2")}
                        #tmp = cor.xl.glob.get("final.blockcol.2")
                #chrom = match(tmp, sort(unique(cor.xl.glob.get("blockcol.2"))))
        else
                {chrom = NULL}
        #!!!!!
        fileName = tempfile(fileext=".png")
        cor.xl.glob.assign("plot3bFile", fileName)
        png(fileName)
        PlotCGH(cor.xl.glob.get("output")$v[,3], chrom=chrom)
        title("Component 3, Data Set 2")
        dev.off()
}

# Returns a vector of bools which are T if that row of mat is all 0
cor.xl.Get0Rows = function(mat){
	all0 = vector(length=nrow(mat))
	for(i in 1:nrow(mat)){
		reduced = unique(mat[i,])
		if(length(reduced)==1 && reduced==0){all0[i] = T}
	}
	return(all0)
}

# Sorts mat and labels and removes all rows for which mat is all 0
cor.xl.OrderRows = function(mat, labels, which) {
	non01 = !cor.xl.Get0Rows(mat)
	non0u = matrix(mat[non01,], ncol=ncol(mat))
	non0l1 = matrix(labels[non01,], ncol=ncol(labels))
	
	pos = non0u[,1] > 0
	neg = non0u[,1] <= 0
	
	
	# calls to as.matrix are to avoid issues when ncol(mat) == 1 or nrow(mat) == 1
	cor.xl.glob.assign(paste("weights",which,sep="."), rbind(matrix(matrix(non0u[pos,], ncol=ncol(non0u))[order(matrix(non0u[pos,], ncol=ncol(non0u))[,1], decreasing=T),], ncol=ncol(non0u)), matrix(matrix(non0u[neg,], ncol=ncol(non0u))[order(matrix(non0u[neg,], ncol=ncol(non0u))[,1]),], ncol=ncol(non0u))))
	cor.xl.glob.assign(paste("labels",which,sep="."), rbind(matrix(non0l1[pos,], ncol=ncol(non0l1))[order(matrix(non0u[pos,], ncol=ncol(non0u))[,1], decreasing=T),],
								  matrix(non0l1[neg,], ncol=ncol(non0l1))[order(matrix(non0u[neg,],ncol=ncol(non0u))[,1]),]))

}
