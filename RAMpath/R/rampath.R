## Original RAMpath input
## By Johnny Zhang, Jack McArdle, & Aki Hamagami
## version 0.1, May 20, 2012
ramIndex<-function(input){
	sptInput<-unlist(strsplit(input, "="))
	value<-as.numeric(sptInput[2])
	left<-sptInput[1]
	#if (substr(left,1,5)!="arrow" || substr(left,1,5)!="sling") warning('The keyword arrow or sling may not be input correctly!')
	leftLength<-nchar(left)
	leftInd<-substr(left,7,leftLength-1)
	leftSplit<-unlist(strsplit(leftInd, ","))
	ind1<-as.numeric(leftSplit[1])
	ind2<-as.numeric(leftSplit[2])
	c(ind1,ind2,value)
}

ramMatrix<-function(model){
	# check for empty syntax
    if(length(model) == 0) {
        stop("ERROR: empty model syntax")
    }

    # break up in lines 
    model <- unlist( strsplit(model, "\n") )

    # remove comments starting with '#'
    model <- gsub("#.*","", model)

    # replace semicolons by newlines and split in lines again
    model <- gsub(";","\n", model)
    model <- unlist( strsplit(model, "\n") )

    # strip all white space
    model <- gsub("[[:space:]]+", "", model)

    # keep non-empty lines only
    idx <- which(nzchar(model))
    model <- model[idx]
    
    ## find the number of manifest variables
    manifestInd<-grep('manifest', model)
    if (length(manifestInd)>1) warning('Multiple manifest keywords found and the first was used.')
    manifestLine<-model[manifestInd]
    manifestSplit<-unlist( strsplit(manifestLine, "=") )
    manifest<-manifestSplit[2]
    ## find the number of latent variables 
    latentInd<-grep('latent', model)
    if (length(latentInd)>0){
    	if (length(latentInd)>1) warning('Multiple latent keywords found and the first was used.')
    	latentLine<-model[latentInd]
    	latentSplit<-unlist( strsplit(latentLine, "=") )
    	latent<-latentSplit[2]
    }else{
    	latent<-0
    }
    
    ## total number of variables and the dimenstion of A and S matrix
    manifest<-as.numeric(manifest)
    latent<-as.numeric(latent)
    nrow<-manifest+latent
    
    A<-S<-matrix(0, nrow, nrow)
    
    ## Construct the matrix for the arrows
    arrowInd<-grep('arrow', model)
    arrowLine<-model[arrowInd]
    for (i in 1:length(arrowLine)){
    	temp<-ramIndex(arrowLine[i])
    	A[temp[1],temp[2]]<-temp[3]
    }
    
    ## Construct the matrix for the slings
    slingInd<-grep('sling', model)
    slingLine<-model[slingInd]
    for (i in 1:length(slingLine)){
    	temp<-ramIndex(slingLine[i])
    	S[temp[1],temp[2]]<-temp[3]
    	S[temp[2],temp[1]]<-temp[3]
    }

	## The filter matrix
	F<-diag(manifest)
	if (latent>0) FLatent<-matrix(0,manifest,latent)
	F<-cbind(F,FLatent)
	
	## variable names
	labelInd<-grep('label', model)
	varname<-character(0)
	if (length(labelInd)>0){
		varnameLine<-model[labelInd]
		varname<-substr(varnameLine,7,nchar(varnameLine))
		varname<-unlist(strsplit(varname,",")) 
	}
	if (length(varname)!=nrow){
		if (length(varname)>0){
			warning('The provided number of varible labels is not equal to the total number of variables. Default variables names will be used')
		}
		varname<-c(paste('V',1:manifest,sep=''),paste('F',1:latent,sep=''))		
	} 
	rownames(S)<-colnames(S)<-varname
	rownames(A)<-colnames(A)<-varname
	rownames(F)<-varname[1:manifest]
	colnames(F)<-varname
	
	cat("Matrix A\n")
	print(A)
	cat("Matrix S\n")
	print(S)
	cat("Matrix F\n")
	print(F)
  lname<-NULL
  if (nrow>manifest)  lname=varname[(manifest+1):nrow]
	invisible(return(list(F=F, A=A, S=S, nvar=nrow, manifest=manifest,latent=latent,lname=lname,varname=varname)))
}

ramFlip<-function(input){
	input<-unlist(strsplit(input,">"))
	n<-length(input)
	input<-input[n:1]
	paste(input[1:(n-1)], collapse="<")
}
ramRmOne<-function(input){
	input<-unlist(strsplit(input,">"))
	n<-length(input)
	paste(input[2:n], collapse=">")
}
## Modified from Boker's Rampath2000 article and S scripts
## makePathList, makeSpanList, makeBridgeList

makePathList <- function(AMatrix, Ase, indirect=TRUE) {
	k <- 0
	tIndex <- 1:10000
	tFrom <- tTo <- tStartID <- tFromID <- tLength <- tValue <- rep(0, 10000)
	tPath <- rep(NULL,10000) ## Save paths
	tPathName<-tPathSig<-rep(NULL,10000)
		
	tNames<-rownames(AMatrix)
	if (is.null(tNames)) tNames<-1:nrow(AMatrix)	
	for (i in 1:nrow(AMatrix)) {
		for (j in 1:ncol(AMatrix)) {
			if (AMatrix[i,j] != 0) {
				k <- k + 1
				tFrom[k] <- j
				tTo[k] <- i
				tStartID[k] <- tIndex[k]
				tFromID[k] <- tIndex[k]
				tLength[k] <- 1
				tValue[k] <- AMatrix[i,j]
				tPath[k]<-paste(j,' > ',i,sep='')
				tPathName[k]<-paste(tNames[j],' > ',tNames[i],sep='')
				if (!missing(Ase)){
					if (Ase[i,j] != 0){
						tPathSig[k]<-abs(AMatrix[i,j]/Ase[i,j])
					}
				}				
			}
		}
	}
	
	m<-k
	if (indirect){
	t1 <- 0
	maxLength <- 0
	while (t1 != k) {
		t1 <- k
		maxLength <- maxLength + 1
		tSelect1 <- tLength == 1
		tSelect2 <- tLength == maxLength 
		if (length(tIndex[tSelect1]) == 0 | length(tIndex[tSelect2]) == 0) 
			break
		for (i in tIndex[tSelect2]) {
			for (j in tIndex[tSelect1]) {
				if (tTo[i] == tFrom[j]) {
					k <- k + 1
					tFrom[k] <- tFrom[i]
					tTo[k] <- tTo[j]
					tStartID[k] <- tIndex[i]
					tFromID[k] <- tIndex[j]
					tLength[k] <- 1 + maxLength
					tValue[k] <- tValue[i] * tValue[j]
					tPath[k]<-paste(tPath[i],' > ',tTo[j],sep='')
					tPathName[k]<-paste(tPathName[i],' > ',tNames[tTo[j]],sep='')
				}
			}
		}		
	}
	}
	return(list(index=tIndex[1:k], fromVar=tFrom[1:k], 
	            toVar=tTo[1:k], startID=tStartID[1:k], 
	            fromID=tFromID[1:k], length=tLength[1:k], 
				value=tValue[1:k], tPath=tPath[1:k], tPathName=tPathName[1:k], nPath=m, tPathSig=tPathSig[1:m]))
}


makeSpanList <- function(SMatrix, Sse) {
	k <- 0
	tIndex <- c(1:10000)
	tVarA <- tVarB <- tValue <- rep(0, 10000)
	tPath <- rep(NULL,10000) ## Save paths
	tPathName <- tPathSig <- tPathSE <- rep(NULL,10000)
	
	tNames<-rownames(SMatrix)
	if (is.null(tNames)) tNames<-1:nrow(SMatrix)
	
	for (i in 1:dim(SMatrix)[1]) {
		for (j in 1:dim(SMatrix)[2]) {
			if (SMatrix[i,j] != 0) {
				k <- k + 1
				tVarA[k] <- j
				tVarB[k] <- i
				tValue[k] <- SMatrix[i,j]
				tPath[k]<-paste(j,' <> ',i,sep='')
				tPathName[k]<-paste(tNames[j],' <> ',tNames[i],sep='')
				
				if (!missing(Sse)){
					tPathSE[k] <- Sse[i,j]
					if (tPathSE[k] != 0){
						tPathSig[k]<-abs(tValue[k]/tPathSE[k])												
						#cat(tValue[k], ' ', tPathSE[k], ' ', tPathSig[k],  " \n")
					}
				}
			}
		}
	}
	return(list(index=tIndex[1:k], varA=tVarA[1:k], varB=tVarB[1:k], value=tValue[1:k], tPath=tPath[1:k], tPathName=tPathName[1:k], tPathSig=tPathSig[1:k]))
}

makeBridgeList <- function(pathList, spanList) {
	k <- 0
	tIndex <- c(1:10000)
	tVarA <- tVarB <- tSpanID <- tPath1ID <- tPath2ID <- tValue <- rep(0, 10000)
	tPath <- rep(NULL,10000) ## Save paths
	tPathName <- rep(NULL,10000)
	
	for (i in 1:length(spanList$index)) {
		k <- k + 1
		tVarA[k] <- spanList$varA[i]
		tVarB[k] <- spanList$varB[i]
		tSpanID[k] <- spanList$index[i]
		tPath1ID[k] <- 0
		tPath2ID[k] <- 0
		tValue[k] <- spanList$value[i]
		tPath[k]<-spanList$tPath[k]
		tPathName[k]<-spanList$tPathName[k]
	}
	
	for (i in 1:length(spanList$index)) {
		for (j in 1:length(pathList$index)) {
			
			if (spanList$varA[i] == pathList$fromVar[j]) {
				k <- k + 1
				tVarA[k] <- pathList$toVar[j]
				tVarB[k] <- spanList$varB[i]
				tSpanID[k] <- spanList$index[i]
				tPath1ID[k] <- pathList$index[j]
				tPath2ID[k] <- 0
				tValue[k] <- spanList$value[i] * pathList$value[j]
				tPath[k]<-paste(ramFlip(pathList$tPath[j]), ' < ', tPath[i], sep='')
				tPathName[k]<-paste(ramFlip(pathList$tPathName[j]), ' < ', tPathName[i], sep='')
			}
			
			if (spanList$varB[i] == pathList$fromVar[j]) {
				k <- k + 1
				tVarA[k] <- spanList$varB[i]
				tVarB[k] <- pathList$toVar[j]
				tSpanID[k] <- spanList$index[i]
				tPath1ID[k] <- 0
				tPath2ID[k] <- pathList$index[j]
				tValue[k] <- spanList$value[i] * pathList$value[j]
				tPath[k]<-paste(ramFlip(pathList$tPath[j]), ' < ', tPath[i], sep='')
				tPathName[k]<-paste(ramFlip(pathList$tPathName[j]), ' < ', tPathName[i], sep='')
			}
			
		}
	}
	
	for (i in 1:length(spanList$index)) {
		for (j in 1:length(pathList$index)) {
			if (spanList$varA[i] == pathList$fromVar[j]) {
				for (h in 1:length(pathList$index)) {
					if (spanList$varB[i] == pathList$fromVar[h]) {
						k <- k + 1
						tVarA[k] <- pathList$toVar[j]
						tVarB[k] <- pathList$toVar[h]
						tSpanID[k] <- spanList$index[i]
						tPath1ID[k] <- pathList$index[j]
						tPath2ID[k] <- pathList$index[h]
						tValue[k] <- spanList$value[i] * pathList$value[j] * pathList$value[h]
						tPath[k]<-paste(ramFlip(pathList$tPath[h]), " < ", tPath[i],  substring(pathList$tPath[j],2), sep='')
						tPathName[k]<-paste(ramFlip(pathList$tPathName[h]), " < ", tPathName[i], ">", ramRmOne(pathList$tPathName[j]), sep='')
					}
				}
			}
		}
	}
	
	return(list(index=tIndex[1:k], varA=tVarA[1:k], 
	            varB=tVarB[1:k], spanID=tSpanID[1:k], 
	            path1ID=tPath1ID[1:k], path2ID=tPath2ID[1:k], 
				value=tValue[1:k], path=tPath[1:k],pathName=tPathName[1:k], tPathSig=pathList$tPathSig, tSpanSig=spanList$tPathSig))
}

ramPathBridge<-function(rammatrix, allbridge=FALSE, indirect=TRUE){
	Amatrix<-rammatrix$A
	Smatrix<-rammatrix$S
	if (is.null(rammatrix$Ase)) {
		tPathlist <- makePathList(Amatrix,indirect=indirect)
	}else{
		tPathlist <- makePathList(Amatrix,rammatrix$Ase,indirect=indirect)
	}
    if (is.null(rammatrix$Sse)) {
    	tSpanlist <- makeSpanList(Smatrix)
    }else{
    	tSpanlist <- makeSpanList(Smatrix,rammatrix$Sse)
    }
    tBridgelist<-NULL
    if (allbridge) tBridgelist <- makeBridgeList(tPathlist, tSpanlist)
    ramobject<-list(path=tPathlist, bridge=tBridgelist, span=tSpanlist, ram=rammatrix)
	class(ramobject)<-'RAMpath'
	ramobject 
}

plot.RAMpath <- function (x, file, from, to, type=c("path","bridge"), size = c(8, 8), node.font = c("Helvetica", 14), edge.font = c("Helvetica", 10), rank.direction = c("LR","TB"), digits = 2, output.type=c("graphics", "dot"), graphics.fmt="pdf", dot.options=NULL, ...)
{
    pathbridge<-x
	if (length(type)>1) type<-type[1]
	tPathlist<-pathbridge$path
	tBridgelist<-pathbridge$bridge
	tSpanlist<-pathbridge$span
	rammatrix<-pathbridge$ram
	varname<-rammatrix$varname
	nvar<-length(varname)
	latent <- rammatrix$lname
	npath<-tPathlist$nPath
	nspan<-length(tSpanlist$index)
	pathValue<-round(tPathlist$value[1:npath],digits=digits)
	spanValue<-round(tSpanlist$value)
	tPathlist$value<-round(tPathlist$value,digits)
	tSpanlist$value<-round(tSpanlist$value,digits)
	if (!is.null(tBridgelist)) tBridgelist$value<-round(tBridgelist$value,digits)
	tPathSig<-pathbridge$path$tPathSig
	tSpanSig<-pathbridge$span$tPathSig
	tBridgelist$pathName <- gsub("[[:space:]]+", "", tBridgelist$pathName)
	tPathlist$tPathName <- gsub("[[:space:]]+", "", tPathlist$tPathName)
	
	output.type <- match.arg(output.type)
	
	if (missing(from)|missing(to)){
		## plot the solid path diagram if missing from or to
		if (!missing(file)) {
			dot.file <- paste(file, ".dot", sep="")
			handle <- file(dot.file, "w")
			on.exit(close(handle))
			if (output.type == "graphics") graph.file <- paste(file, ".", graphics.fmt, sep="")
		}else handle <- stdout()
	
		rank.direction <- match.arg(rank.direction)
		cat(file = handle, paste("digraph \"pathdiagram\" {\n", sep = ""))
		cat(file = handle, paste("  rankdir=", rank.direction, ";\n",sep = ""))
		cat(file = handle, paste("  size=\"", size[1], ",", size[2],"\";\n", sep = ""))
		cat(file = handle, paste("  node [fontname=\"", node.font[1], "\" fontsize=", node.font[2], " shape=box];\n", sep = ""))
		cat(file = handle, paste("  edge [fontname=\"", edge.font[1],"\" fontsize=", edge.font[2], "];\n", sep = ""))
		cat(file = handle, "  center=1;\n")
		
    	if (!is.null(latent)){
			for (lat in latent) {
				cat(file = handle, paste("  \"", lat, "\" [shape=ellipse]\n",sep = ""))
			}
    	}
    	
    	## plot the mean if necessary
    	if (max(abs(rammatrix$M))>0){
    		cat(file = handle, paste("  \"1\" [shape=triangle]\n",sep = ""))
    		cat(file = handle, paste("  \"1\" -> \"1\" [label=\"1\"   dir=both]\n",sep = ""))
    	}
    	
		## single headed arrows	
		for (i in 1:npath){
			sig<-''
			if (!is.null(tPathSig)){
				if (!is.na(tPathSig[i])){
					if (tPathSig[i]>1.96) sig<-'*'
				}
			}
			
			cat(file = handle, paste("  \"", varname[tPathlist$fromVar[i]], "\" -> \"", varname[tPathlist$toVar[i]], "\" [label=\"",tPathlist$value[i], sig, "\"];\n", sep = ""))
		}
		
		## for intercept / means
		if (max(abs(rammatrix$M))>0){
			Mpath<-which(abs(rammatrix$M)>0)
    		for (i in Mpath){
    			sig<-''
    			if (rammatrix$Mse[i]!=0){
    				temp<-rammatrix$M[i]/rammatrix$Mse[i]
    				if (abs(temp)>1.96) sig<-'*'
    			}
    			
    			cat(file = handle, paste("  \"1\" -> \"", rownames(rammatrix$M)[i], "\" [label=\"",round(rammatrix$M[i],digits=digits), sig, "\"];\n", sep = ""))
    		}
    	}
    	
		## double headed arrows	
		for (i in 1:nspan){
			sig<-''
			if (!is.null(tSpanSig)){
				if (!is.na(tSpanSig[i])){
					if (tSpanSig[i]>1.96) sig<-'*'
				}
			}
			if (tSpanlist$varA[i] >= tSpanlist$varB[i]) cat(file = handle, paste("  \"", varname[tSpanlist$varA[i]], "\" -> \"", varname[tSpanlist$varB[i]], "\" [label=\"", tSpanlist$value[i], sig, "\"  dir=both];\n", sep = ""))
		}
		
	
		cat(file = handle, "}\n")
		if (output.type == "graphics" && !missing(file)){
			cmd <- paste("dot -T", graphics.fmt, " -o ", graph.file, " ", dot.options, " ", dot.file, sep="")
			cat("Running ", cmd, "\n")
			result <- try(system(cmd))
		}
	}else{
		## plot certain path or bridge
		if (!is.numeric(from)){ varA<-which(varname==from) }else varA<-from
		if (!is.numeric(to)){ varB<-which(varname==to) } else varB<-to 
		
		if (length(varA)==0|varA>nvar|varA<0) stop("The from variable does not exist")
		if (length(varB)==0|varB>nvar|varB<0) stop("The to variable does not exist")
		
		if (type=="path"){
		## plot the paths
		index<-which(tPathlist$fromVar==varA & tPathlist$toVar==varB)
		if (length(index)==0) stop('No such path exists')
		for (j in 1:length(index)){
			## plot the solid path diagram if missing from or to
			if (!missing(file)) {
				dot.file <- paste(file, j, ".dot", sep="")
				handle <- file(dot.file, "w")
				on.exit(close(handle))
				if (output.type == "graphics") graph.file <- paste(file, j, ".", graphics.fmt, sep="")
			}else handle <- stdout()
	
			rank.direction <- match.arg(rank.direction)
			cat(file = handle, paste("strict digraph \"pathdiagram\" {\n", sep = ""))
			cat(file = handle, paste("  rankdir=", rank.direction, ";\n",sep = ""))
			cat(file = handle, paste("  size=\"", size[1], ",", size[2],"\";\n", sep = ""))
			cat(file = handle, paste("  node [fontname=\"", node.font[1], "\" fontsize=", node.font[2], " shape=box  style=dashed];\n", sep = ""))
			cat(file = handle, paste("  edge [fontname=\"", edge.font[1],"\" fontsize=", edge.font[2], " style=dashed];\n", sep = ""))
			cat(file = handle, "  center=1;\n")			
			
			## plot specific solid paths and shapes
			sPath<-tPathlist$tPathName[index[j]]
			sPathInd<-unlist(strsplit(sPath,">"))
			nPathToPlot<-length(sPathInd)
			for (k in 1:nPathToPlot){
				if (sPathInd[k] %in% latent){
					cat(file = handle, paste("  \"", sPathInd[k], "\" [shape=ellipse style=solid];\n",sep = ""))
				}else{
					cat(file = handle, paste("  \"", sPathInd[k], "\" [style=solid];\n",sep = ""))
				}
			}
			
			for (k in 1:(nPathToPlot-1)){
				value<- pathValue[(tPathlist$fromVar[1:npath]==which(sPathInd[k]==varname) ) & (tPathlist$toVar[1:npath]==which(sPathInd[k+1]==varname))]
				cat(file = handle, paste("  \"", sPathInd[k], "\" -> \"", sPathInd[k+1], "\" [label=\"", value, "\" style=solid];\n", sep = ""))
			}
			
			## the rest of the path diagram
			if (!is.null(latent)){
			for (lat in latent) {
				cat(file = handle, paste("  \"", lat, "\" [shape=ellipse];\n",sep = ""))
			}
			}
	
			## single headed arrows
			for (i in 1:npath){
				cat(file = handle, paste("  \"", varname[tPathlist$fromVar[i]], "\" -> \"", varname[tPathlist$toVar[i]], "\" [label=\"",tPathlist$value[i], "\"];\n", sep = ""))
			}
			## double headed arrows
			for (i in 1:nspan){
				if (tSpanlist$varA[i] >= tSpanlist$varB[i]) cat(file = handle, paste("  \"", varname[tSpanlist$varA[i]], "\" -> \"", varname[tSpanlist$varB[i]], "\" [label=\"", tSpanlist$value[i], "\"  dir=both];\n", sep = ""))
			}
			cat(file = handle, "label=\"Effect ",tPathlist$tPathName[index[j]],"=",tPathlist$value[index[j]] ,"\";\n labelloc = \"top\";\n")
	
			cat(file = handle, "}\n")
			if (output.type == "graphics" && !missing(file)){
				cmd <- paste("dot -T", graphics.fmt, " -o ", graph.file, " ", dot.options, " ", dot.file, sep="")
				cat("Running ", cmd, "\n")
				result <- try(system(cmd))
			}
		}
		}
		
		if (type=="bridge"){
		
		## plot the bridge
		index<-which(tBridgelist$varA==varA & tBridgelist$varB==varB)
		for (j in 1:length(index)){
			## plot the solid path diagram if missing from or to
			if (!missing(file)) {
				dot.file <- paste(file, j, ".dot", sep="")
				handle <- file(dot.file, "w")
				on.exit(close(handle))
				if (output.type == "graphics") graph.file <- paste(file, j, ".", graphics.fmt, sep="")
			}else handle <- stdout()
	
			rank.direction <- match.arg(rank.direction)
			cat(file = handle, paste("strict digraph \"pathdiagram\" {\n", sep = ""))
			cat(file = handle, paste("  rankdir=", rank.direction, ";\n",sep = ""))
			cat(file = handle, paste("  size=\"", size[1], ",", size[2],"\";\n", sep = ""))
			cat(file = handle, paste("  node [fontname=\"", node.font[1], "\" fontsize=", node.font[2], " shape=box  style=dashed];\n", sep = ""))
			cat(file = handle, paste("  edge [fontname=\"", edge.font[1],"\" fontsize=", edge.font[2], " style=dashed];\n", sep = ""))
			cat(file = handle, "  center=1;\n")			
			
			## plot specific solid paths and shapes
			sPath<-tBridgelist$pathName[index[j]]
			sBridge<-unlist(strsplit(sPath,"<>"))
			sLeft<-unlist(strsplit(sBridge[1],"<"))
			sRight<-unlist(strsplit(sBridge[2],">"))
			sPathName<-unique(c(sLeft,sRight))
			nPathToPlot<-length(sPathName)
			
			## Plot the shapes
			for (k in 1:nPathToPlot){
				if (sPathName[k] %in% latent){
					cat(file = handle, paste("  \"", sPathName[k], "\" [shape=ellipse style=solid];\n",sep = ""))
				}else{
					cat(file = handle, paste("  \"", sPathName[k], "\" [style=solid];\n",sep = ""))
				}
			}
			
			## plot the paths
			## single headed arrows
			
			## left side
			nLeft<-length(sLeft)
			if (nLeft>1){
				for (k in 1:(nLeft-1)){
					value<- pathValue[(tPathlist$fromVar[1:npath]==which(sLeft[k+1]==varname) ) & (tPathlist$toVar[1:npath]==which(sLeft[k]==varname))]
					cat(file = handle, paste("  \"", sLeft[k+1], "\" -> \"", sLeft[k], "\" [label=\"", value, "\" style=solid];\n", sep = ""))
				}
			}
			
			## right side
			nRight<-length(sRight)
			if (nRight>1){
				for (k in 1:(nRight-1)){
					value<- pathValue[(tPathlist$fromVar[1:npath]==which(sRight[k]==varname) ) & (tPathlist$toVar[1:npath]==which(sRight[k+1]==varname))]
					cat(file = handle, paste("  \"", sRight[k], "\" -> \"", sRight[k+1], "\" [label=\"", value, "\" style=solid];\n", sep = ""))
				}
			}
			
			## plot the double headed arrows
			value<- spanValue[(tSpanlist$varA==which(sLeft[nLeft]==varname) ) & (tSpanlist$varB==which(sRight[1]==varname))]
			cat(file = handle, paste("  \"", sLeft[nLeft], "\" -> \"", sRight[1], "\" [label=\"", value, "\" style=solid   dir=both];\n", sep = ""))
			
			
			## the full path diagram
			if (!is.null(latent)){
			for (lat in latent) {
				cat(file = handle, paste("  \"", lat, "\" [shape=ellipse];\n",sep = ""))
			}
			}
	
			## single headed arrows
			for (i in 1:npath){
				cat(file = handle, paste("  \"", varname[tPathlist$fromVar[i]], "\" -> \"", varname[tPathlist$toVar[i]], "\" [label=\"",tPathlist$value[i], "\"];\n", sep = ""))
			}
			## double headed arrows
			for (i in 1:nspan){
				if (!(tSpanlist$varA[i]==which(sLeft[nLeft]==varname) & tSpanlist$varB[i]==which(sRight[1]==varname)) & !(tSpanlist$varB[i]==which(sLeft[nLeft]==varname) & tSpanlist$varA[i]==which(sRight[1]==varname))){
					if (tSpanlist$varA[i] >= tSpanlist$varB[i]) cat(file = handle, paste("  \"", varname[tSpanlist$varA[i]], "\" -> \"", varname[tSpanlist$varB[i]], "\" [label=\"", tSpanlist$value[i], "\"  dir=both];\n", sep = ""))
				}				
			}
			
			cat(file = handle, "label=\"Bridge ",tBridgelist$pathName[index[j]],"=",tBridgelist$value[index[j]] ,"\";\n labelloc = \"top\";\n")
	
			cat(file = handle, "}\n")
			if (output.type == "graphics" && !missing(file)){
				cmd <- paste("dot -T", graphics.fmt, " -o ", graph.file, " ", dot.options, " ", dot.file, sep="")
				cat("Running ", cmd, "\n")
				result <- try(system(cmd))
			}
		}
		}
	}
	invisible(NULL)
}

ramUniquePath<-function(tPathlist){
	tUniquePath<-cbind(tPathlist$fromVar, tPathlist$toVar)
	name<-tPathlist$tPathName[1]
	k<-1
	for (i in 2:length(tPathlist$index)){
		check<-cbind(tUniquePath[1:(i-1), 1]==tPathlist$fromVar[i], tUniquePath[1:(i-1), 2]==tPathlist$toVar[i])
		maxT<-max(apply(check,1,sum))
		if (maxT<2){
			k<-k+1
			tUniquePath[k,]<-c(tPathlist$fromVar[i], tPathlist$toVar[i])
			name<-c(name, tPathlist$tPathName[i])
		}
	}
	return(list(tUniquePath=tUniquePath[1:k, ], tUniqueName=name))
}

summary.RAMpath<-function(object, from, to, type=c("path","bridge"), se=FALSE, ...){
  pathbridge<-object
	if (length(type)>1) type<-type[1]
	tPathlist<-pathbridge$path
	tBridgelist<-pathbridge$bridge
	tSpanlist<-pathbridge$span
	rammatrix<-pathbridge$ram
	varname<-rammatrix$varname	
	nvar<-length(varname)
	
	if (missing(from)|missing(to)){	
	## print the paths
	if (type=="path"){
	## col1: path col2: value col3: percent
	tUniquePath<-ramUniquePath(tPathlist)
	npath<-length(tUniquePath$tUniqueName)
	cat('Path and its decomposions:\n\n')
	nString<-max(nchar(tPathlist$tPathName))

	txt<-sprintf(paste("%-",nString+8,"s %4s %12s %9s\n",sep=""), "Path Name", "", "Value", "Pecent")
	if (se) txt<-sprintf(paste("%-",nString+8,"s %4s %12s %12s %9s\n",sep=""), "Path Name", "", "Value", "se", "Pecent")
	cat(txt)
	for (i in 1:npath){
		index<-which(tPathlist$fromVar==tUniquePath$tUniquePath[i,1] & tPathlist$toVar==tUniquePath$tUniquePath[i,2])
		if (length(index)>0){
			name<-paste(varname[tUniquePath$tUniquePath[i,1]]," > ",varname[tUniquePath$tUniquePath[i,2]],sep="")
			path<-length(index)
			value<-sum(tPathlist$value[index])
			percent<-100
			txt<-sprintf(paste("%-",nString+8,"s %4.0f %12.3f %9.2f\n",sep=""), name, path, value, percent)
			if (se) txt<-sprintf(paste("%-",nString+8,"s %4.0f %12.3f %12s %9.2f\n",sep=""), name, path, value, " ", percent)
			cat(txt)
			for (j in 1:length(index)){
				name<-tPathlist$tPathName[index[j]]
				path<-j
				value.j<-tPathlist$value[index[j]]
				percent.j<-value.j/value*100
				if (se){ 
					se.j<-ramEffectSE(object, name)
					txt<-sprintf(paste("  %-",nString+6,"s   %4.0f %12.3f %12.3f %9.2f\n",sep=""), name, path, value.j, se.j, percent.j)
				}else{
					txt<-sprintf(paste("  %-",nString+6,"s   %4.0f %12.3f %9.2f\n",sep=""), name, path, value.j, percent.j)
				}
					
				cat(txt)
			}
		}
	}
	}
	if (type=="bridge"){
		## col1: path col2: value col3: percent
		nString<-max(nchar(tBridgelist$pathName))
		cat('Covariance and its bridges:\n\n')
		txt<-sprintf(paste("%-",nString+8,"s %4s %12s %9s\n",sep=""), "Path Name", "", "Value", "Percent")
		cat(txt)
		
	
		for (i in 1:nvar){
			for (k in 1:i){
				index<-which(tBridgelist$varA==i & tBridgelist$varB==k)
				if (length(index)>0){
					name<-paste(varname[i]," <> ",varname[k],sep='')
					path<-length(index)
					value<-sum(tBridgelist$value[index])
					percent<-100
					txt<-sprintf(paste("%-",nString+8,"s %4.0f %12.3f %9.2f\n",sep=""), name, path, value, percent)
					cat(txt)
					for (j in 1:length(index)){
						name<-tBridgelist$pathName[index[j]]
						path<-j
						value.j<-tBridgelist$value[index[j]]
						percent.j<-value.j/value*100
						txt<-sprintf(paste("  %-",nString+6,"s   %4.0f %12.3f %9.2f\n",sep=""), name, path, value.j, percent.j)
						cat(txt)
					}## end of j loop
		
				}else{
					name<-paste(varname[i],"<>",varname[k],sep='')
					path<-0
					value<-NA
					percent<-NA
					txt<-sprintf(paste("%-",nString+8,"s   %4.0f %12.3f %9.2f\n",sep=""), name, path, value, percent)
					cat(txt)
				}
			}
		}
	}
	}else{
		if (!is.numeric(from)){ varA<-which(varname==from) }else varA<-from
		if (!is.numeric(to)){ varB<-which(varname==to) } else varB<-to 
		
		if (length(varA)==0|varA>nvar|varA<0) stop("The from variable does not exist")
		if (length(varB)==0|varB>nvar|varB<0) stop("The to variable does not exist")
		if (type=="path"){						
			cat('Path and its decomposions:\n\n')
			nString<-max(nchar(tPathlist$tPathName))
			txt<-sprintf(paste("%-",nString+8,"s %3s %12s %9s\n",sep=""), "Path Name", "", "Value", "Percent")
			cat(txt)

			index<-which(tPathlist$fromVar==varA & tPathlist$toVar==varB)
			if (length(index)==0) stop(paste("No such path from ",varname[varA]," to ",varname[varB], sep=""))
			if (length(index)>0){
				name<-paste(varname[varA]," > ", varname[varB],sep="")
				path<-length(index)
				value<-sum(tPathlist$value[index])
				percent<-100
				txt<-sprintf(paste("%-",nString+8,"s %3.0f %12.3f %9.2f\n",sep=""), name, path, value, percent)
				cat(txt)
				for (j in 1:length(index)){
					name<-tPathlist$tPathName[index[j]]
					path<-j
					value.j<-tPathlist$value[index[j]]
					percent.j<-value.j/value*100
					txt<-sprintf(paste("  %-",nString+6,"s %3.0f %12.3f %9.2f\n",sep=""), name, path, value.j, percent.j)
					cat(txt)
				}
			}
		}
		if (type=="bridge"){
			## col1: path col2: value col3: percent
			nString<-max(nchar(tBridgelist$pathName))
			cat("Covariance and its bridges:\n\n")
			txt<-sprintf(paste("%-",nString+8,"s %3s %12s %9s\n",sep=""), "Path Name", "", "Value", "Percent")
			cat(txt)
					
			index<-which(tBridgelist$varA==varA & tBridgelist$varB==varB)
			if (length(index)==0) stop(paste("No covariance between ", varname[varA], " and ", varname[varB], sep=""))
			if (length(index)>0){
				name<-paste(varname[varA]," <> ",varname[varB],sep='')
				path<-length(index)
				value<-sum(tBridgelist$value[index])
				percent<-100
				txt<-sprintf(paste("%-",nString+8,"s %3.0f %12.3f %9.2f\n",sep=""), name, path, value, percent)
				cat(txt)
				for (j in 1:length(index)){
					name<-tBridgelist$pathName[index[j]]
					path<-j
					value.j<-tBridgelist$value[index[j]]
					percent.j<-value.j/value*100
					txt<-sprintf(paste("  %-",nString+6,"s %3.0f %12.3f %9.2f\n",sep=""), name, path, value.j, percent.j)						
					cat(txt)
				}## end of j loop
		
			}		
		}
	}
}

isNumeric <- function(constant){
	save <- options(warn = -1)
	on.exit(save)
	!(is.na(as.numeric(constant)))
}

ramParseLavaan<-function(input, manifest, type=0){
	sptInput<-unlist(strsplit(input, "="))
	sLeft<-sptInput[1]
	sRight<-sptInput[2]
	
	## the left side
	sLeft1<-unlist(strsplit(sLeft, "(", fixed=TRUE))
	sLeft2<-unlist(strsplit(sLeft1[2], ")", fixed=TRUE))
	leftSplit<-unlist(strsplit(sLeft2, ","))
	toVar<-as.numeric(leftSplit[1])
	fromVar<-as.numeric(leftSplit[2])
	
	## the right side
	## find ? or @
	findAt<-grep("@", sRight)
	if (length(findAt)>0){
		fixed<-0
		if (nchar(sRight)==1){
			fixedAt<-1
			label<-NA
		}else{
			splitRight<-unlist(strsplit(sRight, "@", fixed=TRUE))
			
			if (length(splitRight)==1){
				label<-splitRight[1]
				fixedAt<-1
			}else{
				if (splitRight[1]==""){
					label<-NA
				     fixedAt<-as.numeric(splitRight[2])
				}else{
					label<-splitRight[1]
				     fixedAt<-as.numeric(splitRight[2])
				}
			}
		} 
	}else{
		findQM<-grep("?", sRight, fixed=TRUE)
		if (length(findQM>0)){
			## freely estimated parameter
			fixed<-1
			if (nchar(sRight)==1){
				fixedAt<-NA
				label<-NA
			}else{
				splitRight<-unlist(strsplit(sRight, "?", fixed=TRUE))
				
				if (length(splitRight)==1){
					label<-splitRight[1]
					fixedAt<-NA
				}else{
					if (splitRight[1]==""){
						label<-NA
					     fixedAt<-as.numeric(splitRight[2])
					}else{
						label<-splitRight[1]
				   	  fixedAt<-as.numeric(splitRight[2])
					}
				}
			} 
		}else{
			if (isNumeric(sRight)){
				fixed<-1
				label<-NA
				fixedAt<-as.numeric(sRight)
			}else{
				fixed<-1
				label<-sRight
				fixedAt<-NA
			}
		}
	}
	
	if (type==0){
		if (toVar<=manifest & fromVar>manifest){arrowType<-2}else{arrowType<-1} 
	}else{
		arrowType<-3
	}
	
	list(ramValue=c(fromVar, toVar, arrowType, fixed, fixedAt), ramLabel=label)
}

ram2lavaan<-function(model){
	# check for empty syntax
    if(length(model) == 0) {
        stop("ERROR: empty model syntax")
    }

    # break up in lines 
    model <- unlist( strsplit(model, "\n") )

    # remove comments starting with '#'
    model <- gsub("#.*","", model)

    # replace semicolons by newlines and split in lines again
    model <- gsub(";","\n", model)
    model <- unlist( strsplit(model, "\n") )

    # strip all white space
    model <- gsub("[[:space:]]+", "", model)

    # keep non-empty lines only
    idx <- which(nzchar(model))
    model <- model[idx]
    
    ## find the number of manifest variables
    manifestInd<-grep('manifest', model)
    if (length(manifestInd)>1) warning('Multiple manifest keywords found and the first was used.')
    manifestLine<-model[manifestInd]
    manifestSplit<-unlist( strsplit(manifestLine, "=") )
    manifest<-manifestSplit[2]
    
    ## find the number of latent variables 
    latentInd<-grep('latent', model)
    if (length(latentInd)>0){
    	if (length(latentInd)>1) warning('Multiple latent keywords found and the first was used.')
    	latentLine<-model[latentInd]
    	latentSplit<-unlist( strsplit(latentLine, "=") )
    	latent<-latentSplit[2]
    }else{
    	latent<-0
    }
    
    ## total number of variables and the dimenstion of A and S matrix
    manifest<-as.numeric(manifest)
    latent<-as.numeric(latent)
    nrow<-manifest+latent
    
    ## construct a matrix with from variable, to variable, arrows, fixed/free parameters, starting values, labels
    ## fixed: 0 Free: 1
    ## obs to obs/obs 1; latent to obs 2; double 3
    
    ramMatrix<-NULL
    ramLabel<-NULL
    
    ## Construct the matrix for the arrows
    arrowInd<-grep('arrow', model)
    arrowLine<-model[arrowInd]
    for (i in 1:length(arrowLine)){
    		temp<-ramParseLavaan(arrowLine[i], manifest, 0)
    		ramMatrix<-rbind(ramMatrix, temp$ramValue)
    		ramLabel<-c(ramLabel, temp$ramLabel)
    }
    
    ## Construct the matrix for the slings
    slingInd<-grep('sling', model)
    
    slingLine<-model[slingInd]
    for (i in 1:length(slingLine)){
    		temp<-ramParseLavaan(slingLine[i], manifest, 1)
    		ramMatrix<-rbind(ramMatrix, temp$ramValue)
    		ramLabel<-c(ramLabel, temp$ramLabel)
    }
    
    ## variable names
	labelInd<-grep('label', model)
	varname<-character(0)
	if (length(labelInd)>0){
		varnameLine<-model[labelInd]
		varname<-substr(varnameLine,7,nchar(varnameLine))
		varname<-unlist(strsplit(varname,",")) 
	}
	if (length(varname)!=nrow){
		if (length(varname)>0){
			warning('The provided number of varible labels is not equal to the total number of variables. Default variables names will be used')
		}
		varname<-c(paste('V',1:manifest,sep=''),paste('F',1:latent,sep=''))		
	} 

	## Generate a lavaan model
	nEq<-nrow(ramMatrix)
	lavaanModel<-NULL
	for (i in 1:nEq){
		arrowType<-"~"
		if (ramMatrix[i,3]==2) arrowType="=~"
		if (ramMatrix[i,3]==3) arrowType="~~"
		
		if (ramMatrix[i,4]==1){
			if (is.na(ramMatrix[i,5]) & is.na(ramLabel[i])) preLabel<-NULL
			if (is.na(ramMatrix[i,5]) & !is.na(ramLabel[i])) preLabel<-paste(ramLabel[i], "*", sep="")
			if (!is.na(ramMatrix[i,5]) & is.na(ramLabel[i])) preLabel<-paste("start(",ramMatrix[i,5],")", "*", sep="")
			if (!is.na(ramMatrix[i,5]) & !is.na(ramLabel[i])) {
				multMod<-varname[ramMatrix[i,1]]
				if (ramMatrix[i,3]==2) {multMod<-varname[ramMatrix[i,2]]}
				preLabel<-paste(ramLabel[i], "*",multMod,"+start(",ramMatrix[i,5],")", "*", sep="")
			}
		}else{
			preLabel<-paste(ramMatrix[i,5],"*",sep="")
		}
		
		temp<-paste(varname[ramMatrix[i,2]], arrowType, preLabel, varname[ramMatrix[i,1]], "\n", sep="")
		if (ramMatrix[i,3]==2) temp<-paste(varname[ramMatrix[i,1]], arrowType, preLabel, varname[ramMatrix[i,2]], "\n", sep="")

		lavaanModel<-paste(lavaanModel, temp)
	}
	return(lavaanModel)
}


## Fit a RAM model using Lavaan and return to the RAM matrix
ramFit<-function(ramModel, data, type=c('ram','lavaan'), digits=3, zero.print="0", ...){
  if (missing(type)) type = 'ram'
  if (type=='ram') {
    lavaanModel<-ram2lavaan(ramModel)
  }else{    
    lavaanModel<-ramModel
  }
	fitModel<-lavaan(model= lavaanModel, data=data, fixed.x=FALSE, warn=FALSE, ...)
	parTable<-fitModel@ParTable
	parEst<-fitModel@Fit@est
	parSE<-fitModel@Fit@se
	fitInd<-fitMeasures(fitModel)
	
	varALL<-unique(c(parTable$rhs, parTable$lhs))
	varData<-names(data)
	
	obsVar<-latVar<-NULL
	for (i in 1:length(varALL)){
		if (varALL[i] %in% varData){ 
			obsVar<-c(obsVar, varALL[i])
		}else{
			latVar<-c(latVar, varALL[i])
		}
	}
	varName<-c(obsVar, latVar)
	manifest<-length(obsVar)
	latent<-length(latVar)
	
	nrow<-length(varName)
	A<-S<-Ase<-Sse<-matrix(0,nrow,nrow,dimnames=list(varName, varName))
	M<-Mse<-matrix(0,nrow,1,dimnames=list(varName, 'M'))
	for (j in parTable$id){
		if (parTable$op[j]=="~"){
			A[parTable$lhs[j], parTable$rhs[j]]<-parEst[j]
			Ase[parTable$lhs[j], parTable$rhs[j]]<-parSE[j]
		}
		if (parTable$op[j]=="=~"){
			A[parTable$rhs[j], parTable$lhs[j]]<-parEst[j]
			Ase[parTable$rhs[j], parTable$lhs[j]]<-parSE[j]
		}
		if (parTable$op[j]=="~~"){
			S[parTable$lhs[j], parTable$rhs[j]]<-parEst[j]
			Sse[parTable$lhs[j], parTable$rhs[j]]<-parSE[j]
			S[parTable$rhs[j], parTable$lhs[j]]<-parEst[j]
			Sse[parTable$rhs[j], parTable$lhs[j]]<-parSE[j]
		}
	}
	
	## Print some results
	## Print the Model fit
	#cat("--------------------------------\n")
	#cat("Model fit statistics and indices\n")
	#cat("--------------------------------\n")
	#print(fitInd)
	## Print the parameter estimates
	A.na<-A
	A.na[A==0]<-NA
	S.na<-S
	S.na[S==0]<-NA
	Ase.na<-Ase
	Ase.na[Ase==0]<-NA
	Sse.na<-Sse
	Sse.na[Sse==0]<-NA
	
	cat("\n--------------------\n")
	cat("Parameter estimates:\n")
	cat("--------------------\n")
	cat("\nMatrix A\n\n")
	print(A.na, digits=digits,na.print = zero.print)
	cat("\nMatrix S\n\n")
	print(S.na,digits=digits,na.print = zero.print)
  
	cat("\n----------------------------------------\n")
	cat("Standard errors for parameter estimates:\n")
	cat("----------------------------------------\n")
	cat("\nMatrix A\n\n")
	print(Ase.na,digits=digits,na.print = zero.print)
	cat("\nMatrix S\n\n")
	print(Sse.na,digits=digits,na.print = zero.print)
  lname<-NULL
	if (nrow>manifest) lname=varName[(manifest+1):nrow]
	invisible((list(A=A, S=S, Ase=Ase, Sse=Sse, M=M, Mse=Mse, fit=fitInd, lavaan=fitModel, nvar=nrow, manifest=manifest,latent=latent,lname=lname,varname=varName, model=lavaanModel, data=data)))
}


