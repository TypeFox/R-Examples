############
############  All these functions are from the no-longer-maintainted phyext package,
############  and reproduced here to support the SigTree package,
############  with permission of phyext package author Conrad Stack.
############

# data information:
noneData <- function() { return("none") }
discData <- function() { return("discrete") }
contData <- function() { return("cont") }


guess.datatype <- function(datvals)
{
	if(is.null(dim(datvals)))
		return(character(0))
	
	ndatcols = ncol(datvals)
	datatypes = rep("",ndatcols)
	
	if(ndatcols > 0)
	{
		datatypes[sapply(seq(ndatcols),function(i) is.factor(datvals[,i]))] = discData()
		datatypes[sapply(seq(ndatcols),function(i) is.numeric(datvals[,i]))] = contData()
	}
	
	return(datatypes)
}


write.characters2 <- function(xdf,blockname="CHARACTERS",dtype=c(contData(),discData()),missing.char="?")
{	
	
	dtype = match.arg(dtype)
	# set up state labels:
	use.state.labels=F
	state.labels=character(0)
	
	# For disc data, figure out what the symbols should be 
	discover.symbols<-function(dat,zero.based=T)
	{			
		syms = integer(0)
		if(ncol(dat)!=0){
			for(ii in seq(ncol(dat)))
			{
				syms = c(syms,unique((as.integer(dat[,ii]))))
			}
		}
		syms = sort(unique(syms))
		
		# assume that syms are in order
		if(zero.based)
		{
			if(syms[1] != 0){
				offs = -syms[1]
				syms = syms + offs
			}
		}
		
		return(syms)
	}
	
	# convert factors to zero-based ints
	# CAUTION: this function might not work in next NCL version
	convert.to.int <- function(dat,zero.based=T)
	{
		offset = 0
		if(zero.based)
			offset = -1
		
		for(ii in seq(ncol(dat)))
		{
			dat[,ii] = as.integer(dat[,ii])+offset
		}
		return(dat)
	}
	
	# generate state labels:
	make.state.labels <- function(dat)
	{
		cnames = colnames(dat)
		outlabs = character(ncol(dat))
		for(ii in seq(ncol(dat)))
		{
			# order by
			snames = unique(as.character(dat[,ii])[order(as.integer(dat[,ii]))])
			outlabs[ii] = sprintf("%d %s / %s",ii,cnames[ii],paste(snames,collapse=" "))
		}
		
		return(outlabs)
	}
	
	
	if(!is.data.frame(xdf))
		stop("Internal function .write.characters.block needs a data.frame as the first argument")
	
	header = paste("BEGIN ",blockname,";",sep="")
	header.title = paste("TITLE ",blockname,"_matrix;",sep="")
	header.dims = sprintf("DIMENSIONS NTAX=%d NCHAR=%d;",nrow(xdf),ncol(xdf))
	header.format = sprintf("FORMAT DATATYPE=%s MISSING=%s",ifelse(dtype==contData(),"CONTINUOUS","STANDARD"),missing.char)   # TODO: add GAP, SYMBOLS 
	if(dtype == discData()){
		# This check is any of the levels are NOT integers
		# if they are integers, then it is assumed that they do not 
		# need to be writen as state labels
		#
		use.state.labels = any(is.na(as.integer(levels(xdf[,1]))))
		if(use.state.labels){
			state.labels = make.state.labels(xdf)
			xdf = convert.to.int(xdf) # this does the zero-based conversion
		}
		header.format = sprintf("%s SYMBOLS=\"%s\";",header.format,paste(discover.symbols(xdf),collapse=" "))
	} else {
		header.format = paste(header.format,";",sep="")
	}
	
	if(use.state.labels){
		header.labels = sprintf("CHARSTATELABELS\n\t%s;", paste(state.labels,collapse=","))
	}else{
		header.labels = sprintf("CHARSTATELABELS\n\t%s;", paste(paste(seq(ncol(xdf)),colnames(xdf)),collapse=","))
	}

	mmatrix = "MATRIX"
	#mmatrix.data = unname(cbind(rownames(xdf),apply(xdf,2,as.character)))
	if(dtype==contData()){
		mmatrix.data = apply(xdf,1:2,function(i) sprintf("%0.15f",i))
	} else {
		mmatrix.data = apply(xdf,2,as.character)
	}
	if(any(is.na(mmatrix.data)))
		mmatrix.data[which(is.na(mmatrix.data),T)] <- missing.char
	
	mmatrix.data = apply(mmatrix.data,1,paste,collapse=" ")
	mmatrix.data = unname(cbind(rownames(xdf),mmatrix.data))
	mmatrix.data = unname(apply(mmatrix.data,1,paste,collapse="\t"))
	mmatrix.end = ";\n\nEND;"
	

	return(c(header,
			header.title,
			header.dims,
			header.format,
			header.labels,
			mmatrix,
			mmatrix.data,
			mmatrix.end))

}




# Check for evidence that this file contains simmap-formatted trees
# @param 'vers' is the simmap version that should be checked for
# @param 'quick' if TRUE, then the text string containing a tree is 
#		 only checked to see if the formatting looks right.  If FALSE,
#		 the method attempts to create a tree from the string.
# @param 'vers' can take the values (1.0,1.1,1.5); a meaningless result
#		 will be returned otherewise
#
is.simmap <- function(finput="",text=NULL,vers=c(1.1),quick=TRUE)
{
	if(!is.null(text)){
		# TODO: check text for newlines and split on them if they exist.
		rawtext=text
	} else {
		if(!file.exists(finput))
			stop("Assumed 'finput' was a filename and could not be fond")
		
		rawtext = scan(finput,what=character(0),strip.white=T,sep="\n")		
	}

	## TODO: split individual strings on ';' character
	
	# if there is only one string, then just check that,
	# otherwise assume the file is text from a nexus file
	if(length(rawtext)==1)
	{
		treesblock = rawtext
		treelines = 1
	} else {
		treesblock = read.nexus.block(txt=rawtext,block="trees",rm.comments=T)
		treelines = which(tolower(substr(treesblock,1,4))=="tree")
	}
	
	if(length(treesblock)==0)
	{
		warning("No trees in this file...")
		return(FALSE)
	}
	
	potentialsimmap = logical(length(treelines))
	count = 1
	for(linenumb in treelines)
	{
		# check for simmap styles
		if(vers == 1.1 || vers == 1.0)
		{
			# remove all comments from the string
			junk = gsub("\\[(.*?)\\]","",treesblock[linenumb])
			potentialsimmap[count] = grepl(":\\{.*?\\}",junk)
			if(potentialsimmap[count] && !quick)
			{
				tr = try(read.simmap(text=junk,vers=vers),silent=T)
				if(is(tr,"try-error"))
					potentialsimmap[count] = FALSE
			}
			
		} else {
			
			potentialsimmap[count] = grepl("\\[(.*?)\\]",treesblock[linenumb])
			if(potentialsimmap[count] && !quick)
			{
				junk = gsub("\\[(.*?)\\]","",treesblock[linenumb])
				
				# TODO: change this to read.simmap.new when that function is finished:
				tr = try(read.tree(text=junk),silent=T)
				if(is(tr,"try-error"))
					potentialsimmap[count] = FALSE
			}
		}
		
		count = count + 1
	}
	
	return(any(potentialsimmap))
}






#-----------------------
# Phylo4d class extended
# @ author Conrad Stack
#
# This extension will support "sub-nodes" so that branches can be fragmented to 
# accomadate the results of ancestral state reconstructions (state changes along 
# branches).
#-----------------------

# @param subnode.id  	ID for this node 
# @param subnode.data  	data.frame the should look like superclass data.frame
# @param subnode.branch The branch that this subnode is located on. (Duplicate information - the reason this uses node id's instead of edge indices is due to the fact that edge indices seem unstable)
# @param subnode.pos  	The distance to the subnode from the ancestor (root) of the branch. Can be a non-zero range
setClass("phylo4d_ext", 
			representation(	subnode.id="integer",
							subnode.data="data.frame",
							subnode.branch="matrix",  
							subnode.pos="matrix",
							weight="numeric" ),
			,contains="phylo4d",
			prototype=prototype(
				subnode.id=integer(0),
				subnode.data=data.frame(NULL),
				subnode.branch=matrix(0,nrow=0,ncol=2),
				subnode.pos=matrix(0,nrow=0,ncol=2),
				weight=numeric(0)
				)
			)


#------------------
# check validity
#------------------
validPhylo4d_ext <-function(object)
{
	
	retval = TRUE

	# number of ids and branches are the same
	retval = retval && (length(object@subnode.id) == nrow(object@subnode.branch))

	# number of branches and number of entries in the data matrix are the same
	if(length(object@subnode.data)!=0)
		retval = retval && (nrow(object@subnode.branch) == nrow(object@subnode.data))
	
	# number of ids and number of rows in position matrix are similar
	retval = retval && (length(object@subnode.id) == nrow(object@subnode.pos))
	
	# Check that branch ids exist in superclass 
	if(nrow(object@subnode.branch)!=0)
		for(br.ind in seq(nrow(object@subnode.branch)))
			retval = retval && any(apply(edges(object),1,function(i) all(object@subnode.branch[br.ind,] == i)))
		
	# check that data.frames in 'data' slot and 'subnode data' slot are structurally similar
	# retval = retval && all( names(object@data) == names(object@subnode.data) )
	
	return (retval)
}
setValidity("phylo4d_ext",validPhylo4d_ext)


#----------------------------------------------
#  Constructors: extended phylo4d  
#  (overloaded)
#----------------------------------------------

# make it a generic function with return class of phylo4d_ext (or list of these)
setGeneric("phyext", function(x, ...) { standardGeneric("phyext")}, valueClass=c("phylo4d_ext","list") )

# return itself (otherwise it will use the function below):
setMethod("phyext", "phylo4d_ext",
    function(x, ...) {
	    return(x)
})
	
# if first arg is a phylo4d object:
setMethod("phyext", "phylo4d",
    function(x, snode.data=NULL, snode.branch=NULL, snode.pos=NULL, ...) {
	
	# default values:
	if(missing(snode.data) || is.null(snode.data))
		snode.data = data.frame(x@data[0,])
	
	if(missing(snode.branch) || is.null(snode.branch))
		snode.branch = matrix(0,nrow=0,ncol=2)
	
	if(missing(snode.pos) || is.null(snode.pos))
		snode.pos = matrix(0,nrow=0,ncol=2)
	
	# convert to matrix is it isn't already
	if(is.numeric(snode.pos))
		snode.pos = matrix(rep(snode.pos,2),ncol=2)
	
	if(is.numeric(snode.branch))
		snode.branch = matrix(snode.branch,ncol=2)
		
	# Convert any singletons into subnodes:
	#
	if(hasSingle(x))
	{
		
		###########
		# Changed (7/11) by conrad to process more than one subnode:
		#
		# process subnodes:		
		
		
		# Get singleton node ids:
		tab=table(edges(x)[,1])
		snodeid = as.integer(names(tab)[which(tab==1)])
		snodeid = snodeid[snodeid!=0]  # 0 is a dummy node
		
		
		# Get singleton indices: (edges where singletons are the descendant)
		#sings = sapply(snodeid,function(i) which(edges(x)[,2] == i))
		sings = snodeid
		
		# Get singleton data: (why are edge indices used here(sings) instead of node ids?)
		#decs.data = data.frame(tdata(x,"all")[sings,])
		#dnames <- colnames(tdata(x))
		#names(decs.data) <- dnames
		decs.data = tdata(x,"all")[sings,,drop=F]
		
		
		########################################################################################
		# Setup other temporary data containers for transferring singleton data to subnodes
		# 
		# Idea:
		# 1.)	einds gets a first singletons:
		# 		new.ancestor <<----------->> first.singleton
		# 2.)	For loop the processes all subsequent singletons, looping until non-singleton is found
		# 3.) 	Put data back together
		#
		
		# 1.)
		einds = which(apply(edges(x),1,function(i) !(i[1]%in%snodeid) && i[2]%in%snodeid  )) # get all first singletons-on-a-branch
		alle = edges(x)[einds,]
		if(!is.matrix(alle))  alle = matrix(alle,nrow=1)
		elens = edgeLength(x)[einds]
		#edata = data.frame()
		#colnames(edata) <- dnames
		edata = tdata(x,"all")[alle[,2],,drop=F]
		ee = cbind(alle,elens)
		newnodes = matrix(0,ncol=3,nrow=0)
		
		# 2.)
		for(ii in seq(nrow(ee)))
		{
			anc = ee[ii,1]
			dec = ee[ii,2]
			
			# if decendent is still a singleton:
			if(dec %in% snodeid)
			{
				newlens = ee[ii,3]
				newdata = data.frame(NULL)
	
				# travel back until a non-singleton is found:
				while(dec %in% snodeid)
				{
					#cat("dec=",dec," - ")
					
					# length to new subnode:
					tmplen = unname(edgeLength(x)[which(edges(x)[,1]==dec)])
					newlens = append(newlens, (tail(newlens,1)+tmplen) )
					#cat(levels(tmpdata[1,])[as.integer(tmpdata[1,])],"\n")
					
					# new subnode data:
					#tmpdata = data.frame(decs.data[which(sings==dec),])
					#colnames(tmpdata) <- dnames
					tmpdata = decs.data[which(sings==dec),,drop=FALSE]
					newdata = rbind( newdata, tmpdata )						
					
					dec = edges(x)[which(edges(x)[,1]==dec),2]
				}
				
				# use any
				newlens = head(tail(newlens,-1),-1)
				newdata = tail(newdata,-1)
				if(length(newlens)!=0)
				{
					for(kk in seq(length(newlens)))
					{
						newnodes = rbind(newnodes,c(anc,dec,newlens[kk]))
					}
					edata = rbind(edata,newdata)
				}
				ee[ii,2] <- dec
			}
		}
		
		# 3.)
		eetmp = rbind(ee,newnodes)
		ancs = unname(eetmp[,1])
		decs = unname(eetmp[,2])
		snode.length = unname(eetmp[,3])
		snode.data = edata
		#names(snode.data) <- dnames
		
		
		###############################################################################################################
		
		# new labels:
		ancs.label = labels(x)[ancs]
		decs.label = labels(x)[decs]
		
		# Remove singletons from original tree
		x = collapse.singletons(x)
		
		# condition the branches:
		for(ii in seq(length(ancs)))
		{
			# NOTE: This function assumes that labels are unique (might not always be true...)
			snode.branch = rbind(snode.branch, edges(x)[.edge.index(x, ancs.label[ii], decs.label[ii]),])
			# Conrad: make subnode position a fraction of the parent branch length:
			#snode.pos = rbind(snode.pos, rep(snode.length[ii],2))
			snode.pos = rbind(snode.pos, rep(snode.length[ii],2)/edgeLength(x)[.edge.index(x,snode.branch[ii,1],snode.branch[ii,2])]) 
		}
	}
	
	# create dummy ids for subnodes 
	# NOTE: subnode.ids not currently being used
	#
	idvals = integer(0)
	if(nrow(snode.branch)!=0)
	{
		idstart = (nTips(x) + nNodes(x) + 1)
		idvals = as.integer(seq(from = idstart, length.out=nrow(snode.branch)))
	} 
	
	# Annotate collapsed tree with subnode data:
	retobj = new("phylo4d_ext",
		x,
		subnode.id=idvals, 
		subnode.data = snode.data,
		subnode.branch = snode.branch, 
		subnode.pos = snode.pos)
	
	return(retobj)
})



setMethod("phyext", "phylo4",
    function(x, snode.data=NULL, snode.branch=NULL, snode.pos=NULL, ...) {

	# create new phylo4 object (allow arguments to be passed via ellipsis)
	phyd = phylo4d(x,...)
	
	phyext(phyd,snode.data,snode.branch,snode.pos)
})



setMethod("phyext", "phylo",
    function(x, snode.data=NULL, snode.branch=NULL, snode.pos=NULL, ...) {
	phyd = phylo4d(x,...)
	phyext(phyd,snode.data,snode.branch,snode.pos)
	
})


setMethod("phyext","list",
	function(x,...){
	
	x = sapply(x,phyext,...)
	return(x)
})


# assume character points to a file or is a piece of text
# TODO: check if trees are simmap formatted!
setMethod("phyext", "character",
    function(x, snode.data=NULL, snode.branch=NULL, snode.pos=NULL, ...) {
	 
	# Figure out what x's format is:
	tr = NULL
	if(file.exists(x))
	{
		# x points to a file  
		if(is.simmap(finput=x)){
			tr = try(read.nexus.simmap(x)[[1]],silent=T)
		} else {
			tr = try(read.tree(x),silent=T)
		}
		
		if(is(tr,"try-error"))
		{
			tr = try(read.nexus(x),silent=T)
			if(is(tr,"try-error"))
				stop("If first argument is a file then it must contains trees in newick or nexus format")
		}
	} else {
		# x is a character string:
		if(is.simmap(text=x)){
			tr = try(read.simmap(text=x),silent=T)
		} else {
			tr = try(read.tree(text=x),silent=T)
		}
		if(is(tr,"try-error"))
			stop("If first argument is a string, then it must be newick-formated.")
	}
	
	phyd = phylo4d(tr,...)
	phyext(phyd,snode.data,snode.branch,snode.pos)
})








# get coalescent intervals from tip dated trees


tipdate.ci<-function(tr,show.plot=F)
{
	
	if ( !(class(tr) %in% c("multiPhylo","phylo")) )
	  	stop("object \"tr\" is not of class \"phylo\"")
	if(!is.binary.tree(tr))
		stop("Need a binary tree")
	if(!is.rooted(tr))
		stop("Need a rooted tree")

	nleaves = length(tr$tip.label)
	ileaves = seq(1,nleaves)
	#rnode = nleaves + 1 # root node
	rnode =  which(tabulate(tr$edge[,2])==0)[1]
	cipos <- numeric(0)  #!# cipos<<-numeric(0)
	citype <- numeric(0) #!# citype<<-numeric(0)
	cinode <- numeric(0) #!# cinode<<-numeric(0)
	
	isleaf<-function(n)
	{
		return( (length(which(tr$edge[,1]==n)) == 0) )
	}
	
#!#	addlen<-function(len)
#!#	{
#!#		cipos<<-c(cipos,len)
#!#	}
	
	
#!#	addtype<-function(type)
#!#	{
#!#		citype<<-c(citype,type)
#!#	}	
	
#!#	addnode<-function(nnode)
#!#	{
#!#		cinode<<-c(cinode,nnode)
#!#	}
	
	# USE SAPPLY for recursion from now on!
	# Start at root node and trace back from there....
	# node = node number,
	# tr = tree (should stay constant)
	# len = interval lengths
	#
	chknode<-function(node,len)
	{
		#Sys.sleep(0.25)
		#cat("node: ",node," length: ",len,"\n")
		# check for leaf node
		if(isleaf(node))
		{
			cipos <- c(cipos,len) #!# addlen(len)
			citype <- c(citype,0) #!# addtype(0)		
			cinode <- c(cinode,node) #!# addnode(node)
		} else {
			cipos <- c(cipos,len) #!# addlen(len)
			citype <- c(citype,1) #!# addtype(1)	
			cinode <- c(cinode, node) #!# addnode(node)
			# tr$edge 	 col1 = FROM,
			#		 col2 = TO
			inds = which(tr$edge[,1]==node)
			if(show.plot) edgelabels(as.character(tr$edge.length[inds]),inds,adj=-0.5)			
			Recall(tr$edge[inds[1],2],(len + tr$edge.length[inds[1]]))
			Recall(tr$edge[inds[2],2],(len + tr$edge.length[inds[2]]))
		}
	}
	
	if(show.plot) plot.phylo(tr,show.tip.label=F)

	chknode(rnode,0)
	sorder = order(cipos, decreasing=TRUE)
	tmppos = cipos[sorder]
	tmptype = citype[sorder]
	#rm(list=c("cipos","citype","cinode"))
	
	lin = numeric(0) #lineages
	il = abs(diff(tmppos)) #interval.length
	ic = length(cipos)-1 # interval.count
	td = diff(range(cipos)) #total.depth
	
	count=0;
	for (i in tmptype)
	{
		if(i==0) count = count + 1 	# sampling event
		else count= count - 1		# coalescent event
		lin = append(lin,count)
	}
	
	return(list(lineages=lin[1:(length(lin)-1)],interval.length=il,interval.count=ic,total.depth=td,I=tmptype[-1]))
}


setGeneric("treeHeight", function(x) { standardGeneric("treeHeight") })
setMethod("treeHeight", signature(x="phylo4d_ext"),
  function(x) {	
	options(warn=-1)  # suppress warnings for this sections
	retval = tipdate.ci(as(x,"phylo"))$total.depth
	options(warn=0)  
	return(retval)
})

setMethod("treeHeight", signature(x="phylo4d"),
  function(x) {
	options(warn=-1)  # suppress warnings for this sections
	retval = tipdate.ci(as(x,"phylo"))$total.depth
	options(warn=0)  
	return(retval)
})

setMethod("treeHeight", signature(x="phylo4"),
  function(x) {
	options(warn=-1)  # suppress warnings for this sections
	retval = tipdate.ci(as(x,"phylo"))$total.depth
	options(warn=0)  
	return(retval)
})










#-------------------------------------------------------
#  Methods 
#  -Are concerned with subnode stuff only.  
#   Higher-level tree modifications should be done
#  	with phylobase functions.
#-------------------------------------------------------

# Table of Contents:

## Generic Methods: (and overloads)
#

## SIMMAP Methods:
# get.nodenames -- internal; get names of nodes
# is.simmap	 --  check whether file or text string contains simmap formatted trees
# read.simmap  --  read simmap text (v1.0-1.1) to phylo4d object (with singletons)
# read.simmap.new  --  read simmap text (v1.5) to phylo4d_ext
# read.nexus.simmap  -- reads nexus file which contains simmap (all versions) trees
# expand.singles  --  experimental; expand singleton nodes into a bifurcation with one 0-length branch terminating at a dummy node
# collapse.to.singles  --  experimental;  collapse bifurcations which contain a 0-length branch connected to a leaf node into singletons
# collapse.singletons  --  remove singleton nodes from a tree
# collapse.subnodes  --  experimental;  collapse subnodes to singletons
# write.simmap (w/ two subroutines)  -- write a tree (and parts of it's data) to simmap format (vers=c(1.0,1.1,1.5))
# write.simmap.old  --  wrapper; calls write.simmap with v1.0 or v1.1 parameters
# write.simmap.new  --  wrapper; calls write.simmap with v1.5 parameters
# write.nexus.simmap  --  uses write.simmap and other internal functions to write a nexus file using write.simmap to output each tree.

## Subnode Methods:
# .edge.index  --  internal;
# nSubNodes  --  get number of subnodes
# hasSubNodes  --  does the tree have subnodes?
# hasData  --  does the tree contain any data? (is the data slot empty?)
# getSubNodeData  --  returns subnode data.frame
# getSubNodePosish  -- returns subnode positions
# getSubNodeEdgeInds  --  returns subnode edge indices
# getEmptyDataFrame  --  returns empty data frame in the same format as the data frame from the data slot
# addSubNode  --  add a subnode to a tree
# showSubNodes  --  print out a representation of the subnodes to stdout


## Generics
setGeneric("sndata", function(x, ...) { standardGeneric("sndata") })
setGeneric("sndata<-", function(x,datnames=NULL, value) { standardGeneric("sndata<-") })
setGeneric("snid", function(x, ...) { standardGeneric("snid") })
setGeneric("snposition", function(x, ...) { standardGeneric("snposition") })
setGeneric("snbranch", function(x, ...) { standardGeneric("snbranch") })
setGeneric("rmdata", function(x,index,subindex) { standardGeneric("rmdata")} )
setGeneric("weight", function(x) { standardGeneric("weight")} )
setGeneric("weight<-", function(x,value) { standardGeneric("weight<-")} )
setGeneric("hasWeight",function(x,strict=TRUE) { standardGeneric("hasWeight")} )

# subnode
setGeneric("hasSubNodes", function(x) { standardGeneric("hasSubNodes") })
setGeneric("hasData", function(x) { standardGeneric("hasData") })
setGeneric("hasDataColumn",function(x,index) { standardGeneric("hasDataColumn") })

# overloading phylobase methods:

#--------------------------
# SIMMAP Processing Methods		
#--------------------------

# internal:
# get node names from newick string
get.nodenames<-function(newick.txt)
{
	nnames = character(0)
	ttmp = newick.txt
	
	if(grepl("^.*;$",ttmp)) {
		# remove trailing semi-colon, otherwise it will mess up the regular expression matching
		ttmp = substring(ttmp,first=1,last=(nchar(ttmp)-1))
	}
	nname.pat="(\\(|,|\\))([a-zA-Z0-9'_\\.]{1,})(:|;)"
	junk = regexpr(nname.pat,ttmp)
	count = 0
	while(junk[1] != -1)
	{
		tmpname = substr(ttmp,(junk[1]+1),(junk[1]+attr(junk,"match.length")-2))
		nnames = append(nnames,tmpname)
		ttmp = sub(nname.pat,"",ttmp)	
		junk = regexpr(nname.pat,ttmp)
		count = count + 1
		stopifnot(count < 100000)
	}
	return(nnames)
}




# read modified newick file
# citation: 
# Bollback J. P. (2006) SIMMAP: Stochastic character mapping of 
# discrete traits on phylogenies. BMC Bioinformatics. 7:88
# 
# Assume that branch lengths are present (otherwise, why use SIMMAP?)
# Assume that the root node can only have one simmap state (or, only use the first)
#
# @param as.num - convert the data to numeric (by default, a factor data type is used)
# @param text - use text string instead of a file name
# @param vers - which version of simmap to use
# @param as.num - no idea
# @param add.root - should the tree be treated as unrooted?
#
read.simmap <- function(file="",text=NULL, vers=1.1, as.num=FALSE, add.root=TRUE, ...)
{
	
	if(is.null(text))
	{
		#stop("Need to have text for now")
		if(!file.exists(file))
			stop("Assumed 'file' was a filename and could not be fond.")
		
		text = scan(file,what=character(0),strip.white=T,sep="\n")
	}

	# clear whitespace	
	text = gsub("\\s","",text)
	
	# add semicolon to end
	if(substring(text,nchar(text))!=";")
		text = paste(text,";",sep="")
	
	if(TRUE)
	{
		# add root node and internal node names
		count = 1
		while(regexpr("):",text)[1] != -1)
		{
			text = sub( "):", paste(")Internal",sprintf("%0.7d",count),":",sep=""), text)
			count = count + 1
		}
			
		# Root
		if(add.root && regexpr(");",text)[1] != -1)
			text = sub( ");", ")Root:0;" , text)
		
	}
	
	## Poor replacement for regular expressions
	#
	edge.ind = integer(0)
	edge.state = character(0)
	edge.len = numeric(0)
	junk = strsplit(text,"")[[1]]
	sub.branches = cbind( which(junk=="{"), which(junk=="}") )
	br.lens = numeric(nrow(sub.branches))
	
	for(ii in seq(nrow(sub.branches)))
	{
		br.total = 0
		# get the internal part
		within = paste(junk[seq( sub.branches[ii,1]+1 , sub.branches[ii,2]-1 )],collapse="")
		splitchar = ifelse(vers==1.0,";",":")
		within.sub = strsplit(within,splitchar)[[1]]
		for(jj in within.sub)
		{
			slptmp = strsplit(jj,",")[[1]]
			edge.ind = append(edge.ind,ii)
			edge.state = append(edge.state,slptmp[1])
			edge.len = append(edge.len,slptmp[2])
			br.total = br.total + as.numeric(slptmp[2])
		}
		br.lens[ii] = br.total
	}

	
	# horrible way to put it back together:
	# switch to regular expressions
	#
	newick.str = paste(junk,collapse="")
	for(ii in seq(nrow(sub.branches)))
	{
		replen = diff(sub.branches[ii,])+1
		substr(newick.str, sub.branches[ii,1], sub.branches[ii,2]) <- paste(rep(" ",replen),collapse="")
		substr(newick.str, sub.branches[ii,1], sub.branches[ii,2]) <- as.character(br.lens[ii])
	}
	newick.str = gsub("\\s","",newick.str) # strip whitespace
	
	# convert to ape format:
	tr = read.tree(text=newick.str)
	
	# make singleton nodes
	#
	edge.len = as.numeric(edge.len)
	all.labels = c(tr$tip.label,tr$node.label)
	
	# all names except the root node (which ought to be named, otherwise this will break)
	#
	#nnames = head(get.nodenames(newick.str),-1)  # This should match up with edge.ind
	nnames = get.nodenames(newick.str)
	scount = 1
	nTips = length(tr$tip.label)
	nInts = length(tr$node.label)
	
	# node information
	dataVal = character(0)
	dataNode = integer(0)
		
	# loop through split branches
	for(kk in unique(edge.ind))
	{
		# this information shouldn't change:
		is.root = F
		edgeind = which(edge.ind == kk)
		trind = which(nnames[kk] == all.labels)
		esplice = which(tr$edge[,2] == trind) # which edge to chop
		anc = tr$edge[esplice,1]
		dec = tr$edge[esplice,2]
	
		# assume that this is the root (no way to check for it otherwise..)
		if( length(esplice) == 0 )
		{
			dec = trind
			is.root = T
		}
		
		# store data information:
		dataNode = append(dataNode,dec)
		dataVal = append(dataVal, edge.state[edgeind[1]])

		newsingles = length(edgeind)-1
		
		if(newsingles==0 || is.root)
			next
		
		# add new singleton nodes
		for(mm in seq(newsingles))
		{
			#cat("- adding a new single\n")
			### add new internal nodes to 'phylo' object:
			#
			tr$Nnode = tr$Nnode + 1; nInts = nInts + 1
			tr$node.label = append(tr$node.label,sprintf("Singleton%0.7d",scount))
			nodeid = nTips + nInts # should be the last one available
			tr$edge[esplice,1] = nodeid
			tr$edge = rbind(tr$edge,c(anc,nodeid))
			tr$edge.length[esplice] = edge.len[edgeind[mm]]
			tr$edge.length = append(tr$edge.length,edge.len[edgeind[(mm+1)]])
			
			# store data information:
			dataNode = append(dataNode,nodeid)
			dataVal = append(dataVal, edge.state[edgeind[(mm+1)]])
			
			
			# update info:
			esplice = length(tr$edge.length) # should be the last one
			dec = nodeid
			scount = scount + 1
		}		
	}
	

	# create phylo4d object
	new.labels = c(tr$tip.label,tr$node.label)
	if(as.num)
		dataVal = as.numeric(dataVal)
	tmpdf = data.frame("simmap_state"=dataVal, row.names=new.labels[dataNode])
	rettree = phylo4d(tr,all.data=tmpdf,missing.data="OK",rownamesAsLabels=TRUE)
	
	#write.nexus(tr,file="written.tree")
	return(rettree)
}


# helper functions for read.simmap.new (below)
#
match.order <- function(basevec,cmpvec)
{
	neworder = order(basevec)
	neworder[order(basevec)] <- order(cmpvec)
	stopifnot(all(cmpvec[neworder] == basevec))
	return(neworder)		
}


strip<-function(str,left=TRUE,right=TRUE)
{
	if(left)
		str = sub("^[ \t\n]+","",str)
	if(right)
		str = sub("[ \t\n]+$","",str)
	return(str)
}


# Read simmap v1.5 files
# ----------------------
# Notes:
# 1. all mutational (mapping) information is contained with a comment block (e.g., [&map={0,0.144863,1,0.453803,0}] ),
# 2. within this block the mutational information is contained within curly braces, { }, preceded by a map command, &map= { },
# 3. the mutational information is from the descendant node to the ancestor (i.e., tips to root assignment),
# 4. the information is STATE, LENGTH, STATE,... excluding the last length, e.g., two changes (0=>1=>0) might look like this, [&map={0,0.144863,1,0.453803,0}], with state 0 having a length of 0.144863, followed by state 1 with a length of 0.453803, with the final state 0 having a length of BRLEN - (0.144863 + 0.453803); this approach helps reduce the size of the file by not including redundant information,
# 5. if a branch has no changes then it will look like this, [&map={STATE}]LENGTH
#
# @param specialpatt indicates which variables should be treated in the SIMMAP way
#		 any variables found which are not in specialpatt are treated as containing 
#		 vector data.  By default, all variable are considered to be in the special
#		 format.
#
read.simmap.new <- function(file="",text=NULL, specialpatt=character(0), add.root=TRUE)
{
	
	str.has <- function(patt,token,lowercase=T){
		len = length(grep(sprintf("%s",patt),tolower(token)))
		if(len == 1)
			return(TRUE)
		return(FALSE)
	}
	
	if(is.null(text))
	{
		if(file.exists(file)){
			text = paste(readLines(file),collapse="")
		}else{
			stop("Could not find file ",file)
		}
	}

	# clear whitespace
	txt = gsub("\\s","",text)
	
	# add semicolon to end
	if(substring(txt,nchar(txt))!=";")
		txt = paste(txt,";",sep="")
	
	# check to see if it is actually a SIMMAP v1.5 formatted
	if(!is.simmap(text=txt,vers=1.5))
		stop("'text' does not seem to contain any simmap-formatted trees")
	
	
	# Extract comments from the tree string (which should contain the mapped data)
	comments.pos = gregexpr("\\[&.*?\\]",txt)[[1]]
	comments.len = attr(comments.pos,"match.length")
	comments = character(length(comments.pos))
	for(ii in seq(length(comments)))
	{
		comments[ii] = substr(txt,(comments.pos[ii]+2), ((comments.pos[ii] + comments.len[ii] - 2)))
	}
	
	length(comments)
	
	
	####
	# Process Trees
	# 1.) strip SIMMAP comments from the tree string
	#
	cleaned = gsub("\\[&.*?\\]","",txt)
	textstr = cleaned
	
	# 2.) add root node and internal node names
	count = 1
	while(regexpr("):",textstr)[1] != -1)
	{
		textstr = sub( "):", paste(")Internal",sprintf("%0.7d",count),":",sep=""), textstr)
		count = count + 1
	}
	
	# check if root has a value, if it does not then 
	# the program needs to know so that it can be set
	# when processing the comments (check if anc == root).
	#
	norootval=FALSE
	if(add.root && regexpr(");",textstr)[1] != -1){
		textstr = sub( ");", ")Root:0;" , textstr)
		norootval=TRUE
	}
	
	# 3.) create phylo4d object from tree string where comments have been removed (data is added next)		
	nodenames = get.nodenames(textstr)
	phy = as(read.tree(text=textstr),'phylo4d')
	#
	# End Process Trees
	####
	
	# Sanity check: (leave out root, if it did not come with an explicit state)
	#if(length(comments)!= (length(nodenames) - as.integer(norootval)) )
	#	stop("Node names does not match up with comments (and they should).")
	
	
	####
	# Process Comments:
	#
	# NOTE: Only one branch mapping is allowed at this point (TODO: fix this)
	#
	subnode.names = c("map","simmap_state",specialpatt)
	subnode.patt = sprintf("^(%s)$",paste(subnode.names,collapse="|"))
	subnode.pos = matrix(NA,nrow=0,ncol=2)
	subnode.branch = matrix(NA,nrow=0,ncol=2)
	subnode.data = character(0)
	subnode.commentid = integer(0)
	subnode.subind = integer(0)	

	phy.data = NULL
	for (kk in seq(length(comments)))
	{
		## Tokenize command comment:
		comment.tokens = strsplit(comments[kk],",")[[1]]
		
		# check if there is an assignment (shouldn't there always be?)
		ex.inds = (grepl("=",comment.tokens))
	
		# reconstitute tokens if they were broken up in comment.tokens
		curr = 1
		using = rep(TRUE,length(ex.inds))
		count = 1
		while(count <= length(ex.inds)){
			if(ex.inds[count]){
				curr=count
			} else {
				using[count] = FALSE
				comment.tokens[curr] = paste(comment.tokens[curr],comment.tokens[count],sep=",")
			}
			count = count + 1
		}
		comment.tokens = comment.tokens[using]
		## end Tokenize command comment
		
		# strip beginning / ending whitespace from each token
		# comment.tokens = strip(comment.tokens)
	
		## Process name/value pairs from comment tokens
		## NOTE: branch mapping is treated differently that other attributes
		nnames = character(0)
		vals = character(0)
		
		for(attr.str in comment.tokens)
		{
			# strip whitespace from tokens
			#
			tmp=strsplit(attr.str,"=")[[1]]
			metric.basename = strip(tmp[1])
			metric.value = strip(tmp[2])
			if(str.has("\\{",metric.value) && str.has("\\}",metric.value))
			{
				metric.value = sub("^\\{(.*?)\\}$","\\1",metric.value)
				metric.value = strip(strsplit(metric.value,",")[[1]])
			}
			
			if(!grepl(subnode.patt,metric.basename))
			{
				# process regular values
				# Expand basenames:
				if(length(metric.value) != 1)
					metric.basename = sprintf("%s_%04d",metric.basename,seq(length(metric.value)))
				
				nnames = append(nnames,metric.basename)
				vals = append(vals,metric.value)
				
			} else {
				
				# process map nodes:
				# values format: STATE,LENGTH,STATE,....
				# FROM desc to ancs (so, 1st STATE and LENGTH pair are the values for the current node)
				#
				if(length(metric.value) == 1){
					nnames = append(nnames,metric.basename)
					vals = append(vals,metric.value[1])
				} else {
					
					desc = unname(which(labels(phy)==nodenames[kk]))
					# check if desc is the root node (if it is then the ancestor is node 0 per phylobase's standard, I think).
					if(isRooted(phy)){
						anc = ifelse(rootNode(phy) == desc ,0 ,unname(ancestor(phy,desc)))
					} else {
						anc = unname(ancestor(phy,desc))
					}
					eind = .edge.index(phy,anc,desc)
					elen = unname(edgeLength(phy)[eind])
	
					# if one more length is needed:
					if( length(metric.value) %% 2 != 0 )
					{
						tonow = sum(as.numeric(metric.value[seq(2,length(metric.value),by=2)]))
						metric.value = c(metric.value,(elen-tonow))
					}
	
					
					# Separate states from lengths and make lengths fractions of BRLEN 
					# NOTE: fractions start at ANC(!) not DESC in the phyext format.  This may cause some confusion
					#
					states = rev(metric.value[seq(1,length(metric.value),by=2)])
					lens = rev(as.numeric(metric.value[seq(2,length(metric.value),by=2)]))
					lens = cumsum(lens) / elen
					
					# Assign data from DESC node
					nnames = append(nnames,metric.basename)
					vals = append(vals,tail(states,1)) # last one the DESC node state
					
					# Cache data from subnodes on this branch:
					for (subind in seq(length(states)-1))
					{
						subnode.branch = rbind(subnode.branch,c(anc,desc)) # which branch
						subnode.pos = rbind(subnode.pos,c(lens[subind],lens[subind]))  # which position on the branch
						tmpstate=states[subind]
						names(tmpstate) <- metric.basename
						subnode.data = append(subnode.data,tmpstate)  # TODO: make sure this works well
						subnode.commentid = append(subnode.commentid,kk)
						subnode.subind = append(subnode.subind, subind)
					}
				}
			}
		}
		
		stopifnot(length(nnames) == length(vals))
			
		if(is.null(phy.data))
		{
			phy.data = matrix(NA,nrow=length(nodenames),ncol=length(nnames))
			rownames(phy.data) <- nodenames
			colnames(phy.data) <- nnames
		} else {
	
			newcols = setdiff(nnames,colnames(phy.data))
			if(length(newcols) != 0){
				
				for(cname in newcols)
					phy.data = cbind(phy.data,NA)
				
				colnames(phy.data) <- c(head(colnames(phy.data),-length(newcols)),newcols)
			}
			
			unusedcols = setdiff(colnames(phy.data),nnames)
			if(length(unusedcols) != 0)
			{
				nnames <- append(nnames,unusedcols)
				vals <- append(vals,rep(NA,length(unusedcols)))
			}
		}
		
		phy.data[kk,] = vals[match.order(colnames(phy.data),nnames)]
		
	} # end kk comments loop
	
	# Guess as phy.data types:
	options(warn=-1)  # suppress warnings for this sections
	phy.data = data.frame(apply(phy.data,2,I),stringsAsFactors=FALSE)
	for(jj in seq(ncol(phy.data)))
	{
		exclude = is.na(phy.data[,jj])
		if( all(!is.na(as.numeric(phy.data[(!exclude),jj]))) )
			phy.data[,jj] = as.numeric(phy.data[,jj])
	}
	options(warn=0)
	#
	# End Process Comments
	#####
	

	###
	# Combine processed tree and processed comments:
	if(isRooted(phy)){
		phyd = addData(phy,all.data=phy.data,match.data=TRUE,rownamesAsLabels=TRUE)
	} else {
		phyd = addData(phy,all.data=phy.data,match.data=TRUE,rownamesAsLabels=TRUE,missing.data="OK")
	}
	
	if(length(subnode.data)!=0){
		# add subnode stuff:
		phyd = phyext(phyd)
		usubbranches = unique(subnode.commentid)
		
#		for(ii in seq(length(subnode.data)))
#		{
		for(ii in seq(length(usubbranches)))
		{
			cominds = which(subnode.commentid == usubbranches[ii])
			usubnodes = unique(subnode.subind[cominds])
			for(jj in seq(length(usubnodes)))
			{
				thisinds = cominds[which(subnode.subind[cominds] == usubnodes[jj])]
				tmpdat = subnode.data[thisinds]
				# sanity check:
				check = apply(subnode.pos[thisinds,,drop=FALSE],1,mean)
				stopifnot(all(check == check[1]))
				phyd = addSubNode(phyd,subnode.branch[thisinds[1],1],subnode.branch[thisinds[1],2],subnode.pos[thisinds[1]],tmpdat,pos.is.fraction=TRUE)
			}
		}
	}
	
	return(phyd)
}




# Read trees from a nexus file.  This function is only really necessary for 
# nexus files where trees have simmap formatting.  If they don't, then
# readNexus should be used instead.
#
# note: if type=="all", data is read in from the "characters" block
#		and it's rownames are assumed to match the taxa names.
#
read.nexus.simmap <- function(finput="", text=NULL, vers=NULL, type=c("all", "tree"),...)
{
	
	type <- match.arg(type)

	if (type == "all") {
		returnData <- TRUE
	} else {
		returnData <- FALSE
	}


	outtrees = NULL
	if(!is.null(text)){
		# TODO: check text for newlines and split on them if they exist.
		rawtext=text
	} else {
		if(!file.exists(finput))
			stop("Assuming finput is a file and could not find it")
		
		rawtext = scan(finput,what=character(0),strip.white=T,sep="\n")		
	}

	treesblock = read.nexus.block(txt=rawtext,block="trees",rm.comments=T)
	
	if(length(treesblock)==0)
	{
		warning("No trees in this file...")
		return(FALSE)
	}

	## TODO: split individual strings on ';' character into newlines
	
	treelines = which(tolower(substr(treesblock,1,4))=="tree")
	if(length(treelines)==0)
		return(NULL)

	
	# make sure the the translated table is used if found
	transstart = which(sapply(treesblock,function(i) tolower(i)=="translate",USE.NAMES=F))
	translated = length(transstart)!=0
	transend=integer(0)
	taxatrans=""
	index=integer(0)
	taxname=""
	
	# assume that trees are not translated otherwise:
	# these should be in order numerically.
	if(translated)
	{
		transend = seq(transstart,min(treelines)-1)[min(grep(";",treesblock[seq(transstart,min(treelines)-1)]))]
		taxatrans = treesblock[(transstart+1):transend]
		
		# should ALWAYS be the last token:
		if(any(taxatrans==";"))
			taxatrans = taxatrans[-which(taxatrans==";")]
		
		index = as.integer(sapply(strsplit(taxatrans,"\\s"),function(i) i[1]))  # I'm assuming that they are integers
		taxname = sapply(strsplit(taxatrans,"\\s"),function(i) sub("(,|;)","",i[2]))
	}

	count = 0
	for(linenumb in treelines)
	{
		
		cat("tree number",linenumb,"\n")
		# clearly, this will not work for non-simmap, non-newick trees (or simmap v1.5 strings)
		#tmpstr = tail(strsplit(treesblock[linenumb],"=")[[1]],1)
		# so, trying this instead:
		#tmpstr = sub("^tree.+?=(.+)$","\\1",treesblock[linenumb])
		
		# This might not be necessary:
		if(has.weights(text=treesblock[linenumb])){
			tmpstr = sub("\\[(.*?)\\]","",treesblock[linenumb])
		} else {
			tmpstr = treesblock[linenumb]
		}
		#tmpstr = sub("^tree.+?=(.+)$","\\1",tmpstr)
		#tmpstr = sub("^tree(.+?)=\\s{0,}\\[(.*?)\\](.+)$","\\3",tmpstr)
		#tmpstr = sub("^tree(.+?)=(.*?)(\\(.+)$","\\3",tmpstr)
		roottype = sub("^(TREE|tree)(.+?)=(.*?)(\\(.+)$","\\3",tmpstr)
		roottype = gsub("\\s","",roottype)
		tmpstr = sub("^(TREE|tree)(.+?)=(.*?)(\\(.+)$","\\4",tmpstr)
		
		addroot = !grepl("&U",roottype)
		# check for simmap style
		# remove first comment (should others be removed?
		if(is.simmap(text=tmpstr)){
			trtmp = phyext(unname(read.simmap(text=tmpstr,add.root=addroot)))
		} else if (is.simmap(text=tmpstr,vers=1.5)) {
			trtmp = phyext(unname(read.simmap.new(text=tmpstr,add.root=addroot,...)))
		} else if (is.simmap(text=tmpstr,vers=1.0)){
			trtmp = phyext(unname(read.simmap(text=tmpstr,vers=1.0,add.root=addroot)))
		} else {
			trtmp = phyext(read.tree(text=tmpstr))
		}
		
		
		if(!is(trtmp,'phylo4'))
			trtmp = as(trtmp,'phylo4')
		
		if(!inherits(trtmp,"phylo4")){
			cat("Trouble parsing line for trees:\n",treesblock[linenumb],"\n")
			stop("Wha?")
		}
		
		outtrees = append(outtrees,trtmp)
		count = count + 1
	}
	
	# if the trees have been translated, then
	# add the names back in.
	#
	if(translated)
	{
		for(treeind in seq(length(outtrees)))
		{
			tnames = as.integer(tipLabels(outtrees[[treeind]]))
			if(!any(is.na(tnames)))
			{
				tipLabels(outtrees[[treeind]]) <- taxname[tnames]
			} else {
				warning("Tree translations found, but could not be mapped")
			}
		}
	}
	
	# Find and add data:
	if(returnData && has.block(txt=rawtext,blockname="characters"))  #!# block changed to blockname 7/1/15
	{
		data.part = read.characters2(txt=rawtext,blockname="characters")
		if(length(data.part) > 0)
		{
			for(treeind in seq(length(outtrees)))
			{
				outtrees[[treeind]] <- addData(outtrees[[treeind]],tip.data=data.part,rownamesAsLabels=TRUE)
			}
		}
	}
	
	
	return(outtrees)
}




# expand singleton nodes into bifurcating nodes with one junk node
# this function is mainly for plotting and saving files
# TODO: convert all explicit calls to @'slot' to their accessor counterpart (e.g. tree@edge.length goes to edgeLength(tree)
#
expand.singles <- function(tree,keep.data=FALSE)
{
	# note: tips should be indexed 1...N, where N is the number of tips
	if(!is(tree,"phylo4"))
		stop("tree needs to be of class phylo4")
	
	if(hasSingle(tree))
	{
		has.data = hasData(tree) && keep.data		

		# Need to add labels to trees if any are NA:
		# TODO: make sure labels are unique
		if(any( is.na(labels(tree)) ))
		{
			count = round(runif(1) * 1e6)
			for(ii in which(is.na(labels(tree))))
			{
				labels(tree)[ii] <- sprintf("Internal%07d",count)  # warning, these labels might not be unique
				count = round(runif(1) * 1e6)
			}
		}
		
		# rename 'Singleton' labels
		#labels(tree) = sub("^Singleton(.*)$","Internal\\1",labels(tree))

		tmptable=table(tree@edge[,1])
		snodes = as.integer(names(tmptable)[which(tmptable==1)])
		snodes = snodes[snodes!=0] # 0 is the root node (check on this...)
		
		nold = nrow(tree@edge)
		nnew = length(snodes)
		count = 1
		for(ii in snodes)
		{
			tree@label = append(sprintf("JUNK%0.7d",count),tree@label)
			tree@edge = rbind(tree@edge, c(ii, count))
			tree@edge.length = append(tree@edge.length,0)
			count = count + 1
		}
		
		# rearrange
		#tree@order = "unknown"
		rootind = which(tree@edge[,1] == 0)
		tree@edge[seq(nold),] = tree@edge[seq(nold),] + nnew
		tree@edge[seq(from=nold+1,to=nold+nnew),1] = tree@edge[seq(from=nold+1,to=nold+nnew),1] + nnew
		tree@edge[rootind,1] = 0
		
		# rename
		nnodes = nTips(tree) + nNodes(tree)
		names(tree@label) <- as.character(seq(1,nnodes))
		names(tree@edge.length) <- apply(tree@edge,1,paste,collapse="-")
		if(has.data)
		{
			rownames(tree@data) <- as.integer(rownames(tree@data)) + nnew
		}
	}
	
	if(!keep.data)
		tree = as(tree,"phylo4")
	
	return(tree)
}


# 
# sister function of expand.singles -> converts internal nodes with
# zero-length branches into singletons (in phylo4 format).  Will only work if
# the zero-len branch is connected to a tip
#
collapse.to.singles <- function(tree,by.name=NULL)
{
	if(is(tree,"phylo"))
		tree = as(tree,"phylo4")
	
	if(!is(tree,"phylo4"))
		stop("tree argument needs to be of class phylo4")

	# don't use the root node
	if(any(edgeLength(tree, seq(1,(nNodes(tree)+nTips(tree)))[nodeType(tree)!='root']) == 0))
	{
		# Remove zero-length branches and their tips
		zinds = which(edgeLength(tree)==0)
		torem = edges(tree)[zinds,]
		torem = torem[edges(tree)[zinds,1]!=0,][,2] # don't include root
		zinds = zinds[edges(tree)[zinds,1]!=0] 		# don't include root
		torem = torem[nodeType(tree)[torem]=="tip"]
		
		if(length(torem)!=0)
		{
			# remove excess tips
			rm.count = length(torem)
			tip.count = length(tipLabels(tree))
			total.count = length(labels(tree))
			
			tree@edge <- edges(tree)[-zinds,]  # TODO: add 'edges<-' to phylobase
			tree@edge.length = tree@edge.length[-zinds]
			tree@label = tree@label[-torem]
			
			## begin reindexing nodes
			# tips
			t.oldseq = seq(1, (tip.count))[-torem]
			t.newseq = seq(1,(tip.count-rm.count))
			t.replace.inds = which(tree@edge[,2] %in% t.oldseq)
			stopifnot(length(t.replace.inds) == length(t.oldseq))  # it should be a 1-1 relationship for tips
			tree@edge[,2] = replace(tree@edge[,2],t.replace.inds ,t.newseq)
			
			# internal
			n.oldseq = seq(tip.count+1,total.count)
			n.newseq = seq(tip.count+1-rm.count, total.count-rm.count)
			for(kk in seq(length(n.oldseq)))
			{
				replaceinds = which(tree@edge[,1] == n.oldseq[kk])
				if(length(replaceinds)!=0)
					tree@edge[replaceinds,1] = rep(n.newseq[kk],length(replaceinds))
					
				replaceinds = which(tree@edge[,2] == n.oldseq[kk])
				if(length(replaceinds)!=0)
					tree@edge[replaceinds,2] = rep(n.newseq[kk],length(replaceinds))
			}
			
			
			## rename
			names(tree@label) <- as.character(seq(1,length(labels(tree))))
			names(tree@edge.length) <- apply(tree@edge,1,paste,collapse="-")
		}
	} else {
		warning("No zero-length branches found to be removed")
	}
	return(tree)
}


# collapse singleton nodes using ape<->phylo4 conversion functions:
# TODO: make this less reliant on ape functions
#
collapse.singletons <- function(phy)
{
	rettree = phy
	
	# if singletons exist
	if(hasSingle(rettree))
	{
		if(is(rettree,"phylo4d"))
		{
			tab=table(edges(rettree)[,1])
			snodeid = as.integer(names(tab)[which(tab==1)])
			snodeid = snodeid[snodeid!=0]  # 0 is a dummy node
			
			# not sure why this is used:
			#ancs =  sapply(snodeid,function(i) which(edges(rettree)[,1] == i)) # where the subnode is the ancestor
			#decs =  sapply(snodeid,function(i) which(edges(rettree)[,2] == i)) # where the subnode is the descendant
			#newdata = data.frame(tdata(rettree,"all")[-decs,],row.names=labels(rettree)[-decs])
			newdata = tdata(rettree)[-snodeid,,drop=F]
			#colnames(newdata) <- colnames(tdata(rettree))
			
			# hack it for now:
			rettree = as(rettree,"phylo") # convert to ape
			rettree <- collapse.singles(rettree) # collapse singles
			if(is.rooted(rettree)) {
				rettree = phylo4d(rettree,all.data=newdata,rownamesAsLabels=TRUE) # create new phylo4d object
			} else {
				rettree = phylo4d(rettree,all.data=newdata,rownamesAsLabels=TRUE,missing.data="OK") # 
			}
		} else {
			rettree = as(rettree,"phylo") # convert to ape
			rettree <- collapse.singles(rettree) # collapse singles
			rettree = phylo4(rettree)
		}
	}
	
	return(rettree)
}


# turn subnodes into singletons
# x - phylobase tree
# rm.ex.data - "remove 'extra' data":
#  should data frames unique to tdat be removed (TRUE)
#  or should extra columns be added to sdat (FALSE)
#  WARNING: rm.ex.data = TRUE could lead to a loss of data (obviously)
#
collapse.subnodes <- function(x,rm.ex.data = TRUE)
{
	
	if(!inherits(x,"phylo4d_ext") || !hasSubNodes(x))
	{
		warning("No subnodes to collapse")
		return(x)
	}
	
	# Step 1: cache data from x and setup variables to hold info for the new tree
	
	## a.) Handle data:
	tdat = tdata(x)
	sdat = sndata(x)
	col.diffs = setdiff(colnames(tdat),colnames(sdat))
	if(length(col.diffs) > 0)
	{
		if(rm.ex.data){
			tdat = tdat[,-which(colnames(tdat) %in% col.diffs)]
		} else {
			tmpnames = colnames(sdat)
			sdat = cbind(sdat,matrix(NA,nrow=nrow(sdat),ncol=length(col.diffs)))
			colnames(sdat) <- c(tmpnames,col.diffs)
		}
	}
	newdat = tdat
	
	# TODO: might need to reorder columns of tdat and sdat to match up....
	
	## b.) Subnode positions
	posmeans = snposition(x)
	if(!is.matrix(posmeans))
		posmeans = matrix(posmeans,nrow=1,ncol=2)
	posmeans = apply(posmeans,1,mean)
	
	## c.) Handle edges
	el = edges(x)
	newel = el
	elens = edgeLength(x)
	newelens = elens
	snedges.inds = getSubNodeEdgeInds(x)
	u.sninds = unique(snedges.inds)
	
	## d.) Handle node ids and labels
	nids = nodeId(x) # node ids
	lastnodeid = max(nids) # base from which to create new node ids
	nodelabs = nodeLabels(x)
	
	# Step 2.) Systematically convert subnodes into new sets of edges (with singletons)
	for(ii in seq(length(u.sninds)))
	{
		# get edge and subnode indices:
		eid = u.sninds[ii]
		sninds = which(snedges.inds == eid)
		curr.anc = el[eid,1]
		curr.desc = el[eid,2]	
		elen = elens[eid]
		newnodeids = lastnodeid + seq(length(sninds))
		lastnodeid = max(newnodeids)
		tmp = sprintf("Singleton%07d",newnodeids)
		names(tmp) <- newnodeids
	
		# Create new nodes and edges	
		# positions FROM anc TO desc:
		newlens = sort(posmeans[sninds] * elen) # absolute positions along the branch
		newlens = unname(c(newlens[1],diff(newlens),elen - max(newlens))) # convert to spacing btw subnodes
		newedges = matrix(NA,ncol=2,nrow=length(newlens))
		newedges[,1] <- unname(c(curr.anc,newnodeids))
		newedges[,2] <- unname(c(newnodeids,curr.desc))
	
		tmpdat = sdat[sninds,,drop=F]
		rownames(tmpdat) <- c(tmp)
		newdat = rbind(newdat,tmpdat)
		
		# Append data:	
		newel = rbind(newel,newedges)
		newelens = c(newelens,newlens)
	
		nodelabs = c(nodelabs,tmp)
	}
	
	# Step 3.) Create new tree with new edgelist,labels, and data
	
	# remove old indices:
	newel = newel[-u.sninds,]
	newelens = newelens[-u.sninds]
	
	# create new tree:
	newphy = phylo4(newel,newelens,tipLabels(x),nodelabs)
	# add data:
	newphy = addData(newphy,all.data=newdat)
	
	return(newphy)
}



# subroutine of write.simmap:
newlabels.v1x <- function(x,usestate,splitchar,write.nas=TRUE) 
{
	
	# Note: this version of simmap only handle one state
	if(is.null(usestate))
		usestate = colnames(tdata(x))
			
	usestate = usestate[1]
	tdat = tdata(x)[,usestate,drop=F]
	snodes.present = hasSubNodes(x)
	sdat = NULL; snedges.inds = NULL
	
	if(snodes.present){
		sdat = sndata(x)[,usestate,drop=F]
		snedges.inds = apply(snbranch(x),1,function(i) .edge.index(x,i[1],i[2]))
	}
		
	es = edges(x)[,2]
	if(!isRooted(x))
		es = union(es,edges(x)[,1])  # added 3/8 to accomodate unrooted trees:
	
	elens = edgeLength(x); elens[is.na(elens)] <- 0.0
	newlenlab=character(length(es))
	
	for(ii in seq(length(es)))
	{
		nodeid=es[ii]
		if(!(ii %in% snedges.inds))
		{
			# added 3/8 to accomodate unrooted trees:
			if(is.na(elens[ii])){
				newlenlab[ii] = ""
			} else if(!is.na(tdat[nodeid,1]) || write.nas){
				newlenlab[ii] = paste("{",tdat[nodeid,1] ,",", elens[ii],"}",sep="")
			} else {
				newlenlab[ii] = as.character(elens[ii])
			}
		} else {
			
			snind = which(snedges.inds == ii)
			snlens = snposition(x)[snind,] * elens[ii]
			snstates = sdat[snind,1]
			if(!is.matrix(snlens))
				snlens = matrix(snlens,nrow=1)
			
			# reorder
			snpos = apply(snlens,1,mean)
			neword = order(snpos,decreasing=T)
			snpos = snpos[neword]
			snstates = snstates[neword]
			
			newlenlab[ii] = paste("{",tdat[nodeid,1],",",(elens[ii]-max(snpos)),sep="")
			
			if(length(snpos)>1)
				for(jj in seq(length(snpos)-1))
					newlenlab[ii] = paste(newlenlab[ii],splitchar,snstates[jj],",",(snpos[jj]-snpos[(jj+1)]),sep="")
			
			newlenlab[ii] = paste(newlenlab[ii],splitchar,tail(snstates,1),",",tail(snpos,1),"}",sep="")
			
		}
	}
	
	oldlabs = labels(x)[es]
	names(newlenlab) <- oldlabs
	oldlabs[which(is.na(oldlabs))] <- ""
	newlab = paste(oldlabs,":", newlenlab ,sep="")
	newlab[which(newlenlab=="")] <- ""  # remove any <NA> data

	return(newlab)
}


# subroutine of write.simmap:
newlabels.v15 <- function(x,usestate,splitchar)
{
	
	if(is.null(usestate))
		usestate = colnames(tdata(x))
	else {
		if(is.numeric(usestate))
		{
			usestate = colnames(tdata(x))[usestate]
		}
	}
	
	# get tip and subnode data:
	tdat = tdata(x)[,usestate,drop=F]
	snodes.present = hasSubNodes(x)
	sdat = NULL; snedges.inds = NULL
	
	if(snodes.present){
		sdat = sndata(x)[,usestate,drop=F]
		snedges.inds = apply(snbranch(x),1,function(i) .edge.index(x,i[1],i[2]))
	}
	
	es = edges(x)[,2]
	if(!isRooted(x))
		es = union(es,edges(x)[,1])  # added 3/8 to accomodate unrooted trees:

	elens = edgeLength(x); elens[is.na(elens)] <- 0.0
	newlenlab=character(length(es)) # new "length" labels for each node
	
	# Setup SIMMAP comments to hold data values
	for(ii in seq(length(es)))
	{
		nodeid=es[ii]
		
		if(!(ii %in% snedges.inds))
		{ # if there is no subnode on this branch:
			if(is.na(elens[ii])) {
				newlenlab[ii] = ""  # added to deal with unrooted trees
			} else {		
				topper = "[&"
				backer = paste("]",elens[ii],sep="")
				middle = ""
				for(colind in seq(length(usestate)))
					if(!is.na(tdat[nodeid,usestate[colind]]))
						middle = paste(middle,sprintf("%s={%s}%s",usestate[colind],tdat[nodeid,usestate[colind]],ifelse(colind==length(usestate),"",",")),sep="")
			
				if(middle==""){
					newlenlab[ii] = as.character(elens[ii])
				} else {
					newlenlab[ii] = paste(topper,middle,backer,sep="")
				}
			}
		} else {
			
			topper = "[&"
			backer = paste("]",elens[ii],sep="")
			middle = ""

			# retrieve subnode data:
			snind = which(snedges.inds == ii)
			snlens = snposition(x)[snind,] * elens[ii]
			
			for(colind in seq(length(usestate)))
			{
				tmp = ""
				snstates = sdat[snind,usestate[colind]]
				if(!is.matrix(snlens))
					snlens = matrix(snlens,nrow=1)
				
				snpos = apply(snlens,1,mean)
				neword = order(snpos,decreasing=T)
				
				snpos = snpos[neword]
				snstates = snstates[neword]
				
				# first, write the DESC state:
				tmp = paste(sprintf("%s={%s", usestate[colind],tdat[nodeid,usestate[colind]]),(elens[ii]-max(snpos)),sep=",")
				#newlenlab[ii] = paste(sprintf("[&%s={%s",usestate[colind],tdat[nodeid,usestate[colind]]),(elens[ii]-max(snpos)),sep=",")
				
				# intermediates, from DESC -> ANC
				if(length(snpos)>1)
					for(jj in seq(length(snpos)-1))
						tmp = paste(tmp,splitchar,snstates[jj],",",(snpos[jj]-snpos[(jj+1)]),sep="")
						#newlenlab[ii] = paste(newlenlab[ii],splitchar,snstates[jj],",",(snpos[jj]-snpos[(jj+1)]),sep="")
				
				# last subnode
				#newlenlab[ii] = paste(newlenlab[ii],splitchar,tail(snstates,1),ifelse(colind==length(usestate),"",","),sep="")
				tmp = paste(tmp,splitchar,tail(snstates,1),"}",sep="")
				middle = paste(middle,tmp,ifelse(colind==length(usestate),"",","),sep="")
			}
			
			if(middle==""){
				newlenlab[ii] = as.character(elens[ii])
			} else {
				newlenlab[ii] = paste(topper,middle,backer,sep="")
			}
			
		}
	}
	# End setup SIMMAP comments
	
	
	# Substitute SIMMAP comments into tree
	oldlabs = labels(x)[es]
	names(newlenlab) <- oldlabs
	oldlabs[which(is.na(oldlabs))] <- ""
	#newlab = paste(oldlabs,":{", newlenlab ,"}",sep="")
	newlab = paste(oldlabs,newlenlab,sep=":")
	newlab[which(newlenlab=="")] <- ""  # remove any <NA> data
	
	# DEBUG: at this point, oldlabs and newlabs look like they have been lined up correctly
	
	return(newlab)		
}


# write modified newick file:
# -This is kind of a hack: it basically writes subnodes and lengths to a new character label and then
#  creates a new 'phylo' class using only those labels (without edge lengths).  Then it uses
#  APEs algorithm to write those labels to a newick string.  
#
# -If any of the SIMMAP datatypes are is.na, then they are left out!  TODO:  Look into changing this in the future
#
write.simmap <- function(x,usestate=NULL,file="",vers=1.1,...)
{
	splitchar = ifelse(vers==1.0,";",":")
	
	if(is.null(usestate))
		usestate = colnames(tdata(x))
	
	#if(hasSubNodes(x))
	if(hasData(x) && hasDataColumn(x,usestate))
	{
		newlab=NULL
		if(vers == 1.5){
			newlab <- newlabels.v15(x,usestate,",")
		} else {
			newlab <- newlabels.v1x(x,usestate,splitchar)
		}
		
		# reorder:
		es = edges(x)[,2]
		if(!isRooted(x))
			es = union(es,edges(x)[,1])  # added 3/8 to accomodate unrooted trees:
	
		altorder = order(es)
		newlab = newlab[altorder]
		ntype = nodeType(x)[es[altorder]]
		phy = as(x,'phylo')
		newedges = edges(x)
		newphy = list(edge=newedges,
						tip.label=newlab[which(ntype=="tip")],
						node.label=newlab[which(ntype!="tip")],
						Nnode=length(which(ntype!="tip"))
					)
		
		class(newphy) <- "phylo"
		
		################################################################
		# borrowed code from APE (write.tree.R):
		#
		output.tree.names=FALSE
		append = FALSE
		digits = 10
		brl <- !is.null(newphy$edge.length)
		nodelab <- !is.null(newphy$node.label)
		f.d <- paste("%.", digits, "g", sep = "")
		
		## Helper functions
		cp <- function(s) STRING <<- paste(STRING, s, sep = "")

		add.internal <- function(i) {
		    cp("(")
		    br <- which(newphy$edge[, 1] == i)
		    for (j in br) {
		        desc <- newphy$edge[j, 2]
		        if (desc > n) add.internal(desc)
		        else add.terminal(j)
		        if (j != br[length(br)])  cp(",")
		    }
		    cp(")")
		    if (nodelab) cp(newphy$node.label[i - n])
		    if (brl) {
		        cp(":")
		        cp(sprintf(f.d, newphy$edge.length[which(newphy$edge[, 2] == i)]))
		    }
		}
		
		add.terminal <- function(i) {
		    cp(newphy$tip.label[newphy$edge[i, 2]])
		    if (brl) {
		        cp(":")
		        cp(sprintf(f.d, newphy$edge.length[i]))
		    }
		}
		## End helper functions

		n <- length(newphy$tip.label)
		STRING <- if (output.tree.names) paste("ERROR!", "(", sep = "") else "("
		br <- which(newphy$edge[, 1] == n + 1)
		for (j in br) {
		    desc <- newphy$edge[j, 2]
		    if (desc > n) add.internal(desc)
		    else add.terminal(j)
		    if (j != br[length(br)]) cp(",")
		}
		if (is.null(newphy$root.edge)) {
		    cp(")")
		    if (nodelab) cp(newphy$node.label[1])
		    cp(";")
		} else {
		    cp(")")
		    if (nodelab) cp(newphy$node.label[1])
		    cp(":")
		    cp(sprintf(f.d, newphy$root.edge))
		    cp(";")
		}
		if (file == "") return(STRING)
	    cat(STRING, file = file, append = append, sep = "\n")
	    
		################################################################		
	} else {
		phy = as(x,'phylo')
		write.tree(phy,file,...)
	}
}

	
# write simmap version 1.1 strings
#
write.simmap.old <- function(x,usestate=NULL,file="",...)
{
	write.simmap(x,usestate,file=file,vers=1.1,...)
}

# write simmap version 1.5 strings
#
write.simmap.new <- function(x,usestate=NULL,file="",...)
{
	write.simmap(x,usestate,file=file,vers=1.5,...)
}



# This is the main write function for phylo4d_ext
# Mainly ripped from APE write.nexus function
#
# @param obj a list or single phylo4d_ext objects
# @param file the file to write to
# @param translate should the tree names be moved to a special section
# @param vers which SIMMAP version should be used (can be vers=c(1.1, 1.0, 1.5))
# @param usestates if NULL, then all data columns are used
# @param dtype for compatibility with older simmap versions
#  	   anything in 'usestate' will be written into simmap comments
#	   anything left over that is compatible with dtype will be written as 
#	   a character.
#				
#
write.nexus.simmap <- function(obj, file = "", translate = TRUE, dtype=c(noneData(),contData(),discData()), ...)
{
	
	if(!is.list(obj))
	{
		if(!is(obj,"phylo4d_ext"))
			stop("This function is only made to work with phylo4d_ext objects. Use write.nexus() instead.")
		
		obj <- list(obj)
	} else {
		if(!all(sapply(obj,is,'phylo4d_ext')))
			stop("This function is only made to work with phylo4d_ext objects or lists of such. Use write.nexus() instead.")
	}
	ntree <- length(obj)

	dtype = match.arg(dtype)
	dat = NULL
	if(hasData(obj[[1]]))
		dat = tdata(obj[[1]],"tip")
	
    cat("#NEXUS\n", file = file)
    cat(paste("[R-package APE, ", date(), "]\n\n", sep = ""),file = file, append = TRUE)
	
    #N <- length(obj[[1]]$tip.label)
	N <- length(tipLabels(obj[[1]]))

        cat("BEGIN TAXA;\n", file = file, append = TRUE)
        cat(paste("\tDIMENSIONS NTAX = ", N, ";\n", sep = ""),file = file, append = TRUE)
        cat("\tTAXLABELS\n", file = file, append = TRUE)
        cat(paste("\t\t", tipLabels(obj[[1]]), sep = ""),sep = "\n", file = file, append = TRUE)
        cat("\t;\n", file = file, append = TRUE)
        cat("END;\n", file = file, append = TRUE)
   	
    cat("BEGIN TREES;\n", file = file, append = TRUE)
    if (translate) {
        ## We take arbitrarily the labels of the first tree, and
        ## translate them as "1", "2", "3", ...
        cat("\tTRANSLATE\n", file = file, append = TRUE)
        tmp <- checkLabel(tipLabels(obj[[1]]))
        X <- paste("\t\t", 1:N, "\t", tmp, ",", sep = "")
        ## We remove the last comma:
        X[length(X)] <- gsub(",", "", X[length(X)])
        cat(X, file = file, append = TRUE, sep = "\n")
        cat("\t;\n", file = file, append = TRUE)
        token <- as.character(1:N)
        names(token) <- tipLabels(obj[[1]])
        tipLabels(obj[[1]]) <- token
        if (ntree > 1) {
            for (i in 2:ntree)
                tipLabels(obj[[i]]) <- token[tipLabels(obj[[i]])]
            class(obj) <- NULL
        }
    } else {
        for (i in 1:ntree)
          tipLabels(obj[[i]]) <- checkLabel(tipLabels(obj[[i]]))
    }
    for (i in 1:ntree) {
	    tprefix = "\tTREE * UNTITLED"
	    #weights
	    if(hasWeight(obj[[i]])){
		    tprefix = sprintf("%s [&W %f]",tprefix,weight(obj[[i]]))
	    }
	    #is rooted
        if (isRooted(obj[[i]])){
        	#cat("\tTREE * UNTITLED = [&R] ", file = file, append = TRUE)
        	tprefix = paste(tprefix," = [&R] ",sep="")
    	}else{
	    	#cat("\tTREE * UNTITLED = [&U] ", file = file, append = TRUE)
	    	tprefix = paste(tprefix," = [&U] ",sep="")
		}
		cat(tprefix, file = file, append = TRUE)
		
		# NOTE (2/13) - write.simmap does this step now, and better:
		#if(is.null(usestates))
		#	usestates = colnames(tdata(obj[[i]]))
		
        cat(write.simmap(obj[[i]], file="",...),"\n", sep = "", file = file, append = TRUE)
    }
    cat("END;\n", file = file, append = TRUE)

	
	# begin characters block write:
	if( dtype != noneData() && !is.null(dat) )
	{
		datnames = colnames(dat)

		treeargs <- list(...)
		if( "usestate" %in% names(treeargs) ){
			
			exclude = integer(0)
			if(!is.character(treeargs$usestate)) {
				exclude = as.integer(treeargs$usestate)
			} else {
				exclude = which(datnames == treeargs$usestate)
			}
			
			if(length(exclude) > 0)
				dat = dat[,-exclude,drop=FALSE]
		}
		
		dat = dat[,which(guess.datatype(dat) == dtype),drop=FALSE]
		
		if(ncol(dat) != 0)
		{
			cat( write.characters2(dat,dtype=dtype), sep = "\n", file= file, append=TRUE)
		} else {
			warning("Could not find any suitable data to write.")
		}
	}
 	# end characters block write

	
}








#------------------------------------------
# Modifying phylo4d-extension objects
# 
#------------------------------------------

# Internal:
# NOTE: labels are case-sensitive
# NOTE: this function needs to deal with NAs
.edge.index <- function(tree,anc,dec)
{
	if(!is(tree,'phylo4'))
		stop("Tree needs to inherit from class phylo4")
	
	# convert from characters to indices
	if(is.character(c(anc,dec))){
		
		if(is.na(anc)){
			anc = which(is.na(labels(tree)))
		} else {
			anc = which(labels(tree) == anc)
		}
		
		if(is.na(dec)){
			dec = which(is.na(labels(tree)))
		} else {
			dec = which(labels(tree) == dec)
		}
		
	}
			
	if(is.na(anc) || is.na(dec) || is.null(anc) || is.null(dec))
		stop("Must specify valid ancestor and decendent: ",anc,"-",dec)
	
	#eind = which(edges(tree)[,1] == anc & edges(tree)[,2] == dec)
	eind = which(edges(tree)[,1] %in% anc & edges(tree)[,2] %in% dec)

	if(length(eind) != 1){
		warning("No connection between ",anc," and ",dec,"\n")
		return(NA)
	}
	
	return(eind)
}


nSubNodes <- function(x)
{
	return (length(x@subnode.id))
}

setMethod("hasSubNodes", signature(x="list"),
  function(x) {
	return(sapply(x,hasSubNodes))
})


setMethod("hasSubNodes", signature(x="phylo4d_ext"),
  function(x) {
	return(nSubNodes(x)!=0)
})


setMethod("hasSubNodes", signature(x="phylo4"),
  function(x) {
	return(FALSE)
})

setMethod("hasSubNodes", signature(x="phylo"),
  function(x) {
	return(FALSE)
})

setMethod("hasData", signature(x="phylo"),
  function(x) {
	return(FALSE)
})

setMethod("hasData", signature(x="phylo4"),
  function(x) {
	return(FALSE)
})

setMethod("hasData", signature(x="phylo4d"),
  function(x) {
	return(ncol(tdata(x))!=0)
})

setMethod("hasData", signature(x="list"),
	function(x) {
		return(sapply(x,hasData))
})

setMethod("hasDataColumn", signature("phylo",index="ANY"),
	function(x,index){
	return(FALSE)
})

setMethod("hasDataColumn", signature("phylo4",index="ANY"),
	function(x,index){
	return(FALSE)
})

setMethod("hasDataColumn", signature("phylo4d",index="numeric"),
	function(x,index){
	retval = FALSE
	if( hasData(x) && length(index) != 0 && all(index <= ncol(tdata(x))) )
	{
		retval=TRUE
	}
	return(retval)
})

setMethod("hasDataColumn", signature("phylo4d",index="character"),
	function(x,index){

	uinds = unique(index)
	ind = which(colnames(tdata(x)) %in% uinds)
	
	if(length(ind) != length(uinds))
	{
		return(FALSE)
	} else {
		return(hasDataColumn(x,ind))
	}
})



getSubNodeData <- function(x,colname)
{
	if(missing(colname))
		return(x@subnode.data)
	
	return(x@subnode.data[colname])
}

getSubNodePosish <- function(x)
{
	return(x@subnode.pos)
}


getSubNodeEdgeInds <- function(x)
{
	edgeinds = integer(nSubNodes(x))
	if(length(edgeinds)!=0)
	{
		for(ii in seq(nrow(x@subnode.branch)))
			edgeinds[ii] = .edge.index(x,x@subnode.branch[ii,1],x@subnode.branch[ii,2])
	}
	return (edgeinds)
}

# return empty data.frame styled like
# @data slot
getEmptyDataFrame <- function(x,...)
{
	tmpdf = data.frame(x@data[0,],...)
	colnames(tmpdf) <- colnames(tdata(x))
	return(tmpdf)
}


# TODO: overload tdata here
# 
addSubNode <- function(x,anc,dec,position,dataf,pos.is.fraction=FALSE)
{
	
	if(!is(x,'phylo4d_ext')){
		warning("Converting x to an extended phylo4d object")
		x = phyext(x)
	}
	
	if(missing(dataf))
		stop("dataf needs to contain something")

	eind = .edge.index(x,anc,dec)
	if(is.na(eind))
		stop("Failure to find edge from ",anc, " to ",dec,"\n")
	
	elen = edgeLength(x)[eind]
	if(elen == 0)
		stop("Cannot place subnode on a zero length branch")
	
	if(position > ifelse(pos.is.fraction,1,elen))
		stop("Position: ",position,", is greater than allowed value of ",ifelse(pos.is.fraction,1,elen),"\n")
	
	# TODO: check for overlapping subnodes
	# TODO: check to see if every col in data.frame has 
	#		an associated col in tdata() AND that col is  
	#		defined for internal and tip nodes

	# construct data frame:
	newdf = getEmptyDataFrame(x)
	if(is.data.frame(dataf))
	{
		if( all(names(dataf) == names(newdf)) ){
			newdf = rbind(newdf, dataf)
		} else {
			# only use cols which exist in tdata(x)
			colinds = which(names(dataf) %in% names(newdf)) 
			if( length(colinds) == 0 )
				stop("Could not find any common columns between dataf and tdata(x)")
			
			newdf = merge(newdf,dataf[,colinds,drop=F],all=T)
		}
		
	} else {
		


		newdf = getEmptyDataFrame(x,stringsAsFactors=FALSE)  # ??

		# special case (if lengths are the same, then assume that they are in order):
		if( length(dataf) == ncol(newdf) && is.null(names(dataf)) )
			names(dataf) <- names(newdf)

		newdf[1,] <- rep(NA,ncol(newdf))
		# try to match up names
		ndf = names(dataf)
		count = 1
		
		for(nm in ndf){
			nmind = which(names(newdf)==nm)
			if(!is.null(nmind) && length(nmind) != 0)
			{
				if(is.factor(newdf[,nmind])  && !(dataf[count] %in% levels(newdf[,nmind])) )
				{
					# need to add more levels if this column is a factor
					newlevs = c(levels(newdf[,nmind]),unname(dataf[count]))
					levels(newdf[,nmind]) <- newlevs
					levels(tdata(x)[,nmind]) <- newlevs	# TODO: TEST THIS!
				}
				newdf[1,nmind] = dataf[count]
			}
			count = count + 1
		}

		
	}
	
	if(!pos.is.fraction) 
		position = position / elen
	
	x@subnode.id = append(x@subnode.id, as.integer(nTips(x) + nNodes(x) + nSubNodes(x) + 1))
	x@subnode.data = rbind(x@subnode.data, newdf)
	x@subnode.branch = rbind(x@subnode.branch, edges(x)[eind,])
	x@subnode.pos = rbind(x@subnode.pos, rep(position,2))
	
	return (x)
}


#
showSubNodes <- function(x)
{
	if(hasSubNodes(x))
	{
		charlen = 80
		brchar = "-"
		terminalchar = "-"
		regchar="0"
		overlapchar = "*"
		# use them all by default:
		
		einds = apply(x@subnode.branch,1,function(i) .edge.index(x,i[1],i[2]))
		
		for(eind in unique(einds))
		{
			snid = which(einds == eind)
			anc = x@subnode.branch[snid,1]
			dec = x@subnode.branch[snid,2]
			
			# setup output string
			tmpstr = rep(brchar,charlen)
			tmpstr[1] = terminalchar
			tmpstr[charlen] = terminalchar
			
			# get relative positions of nodes:
			elen = edgeLength(x)[eind]
			breaksize = elen / charlen
			snpos = x@subnode.pos[snid,] * elen 
			if(!is.matrix(snpos)) 
				snpos = matrix(snpos,nrow=1)
			
			snmeans = apply(snpos,1,mean)
			
			for(xx in seq(length(snid)))
			{
				nbreaks = max(1, floor(diff(snpos[xx,]) / breaksize) )
				from = floor(snpos[xx,] / breaksize)
				tmpstr[seq(from,length.out=nbreaks)][tmpstr[seq(from,length.out=nbreaks)]==regchar] = overlapchar
				tmpstr[seq(from,length.out=nbreaks)][tmpstr[seq(from,length.out=nbreaks)]==brchar] = regchar
				#tmpstr[seq(from,length.out=nbreaks)] = rep(regchar,nbreaks)
			}
			
			# collapse and print:
			tmpstr = paste(tmpstr,collapse="")
			if(length(snid)==1){
				names(tmpstr) <- sprintf("Subnode at ~%0.2f on branch: '%s' to '%s' (brlen=%0.2f)",snpos[1]+diff(snpos[1,])/2,labels(x)[anc],labels(x)[dec],elen)
			} else {
				names(tmpstr) <- sprintf("%d Subnodes on branch: %d to %d (bren=%0.2f)",length(snid),anc[1],dec[1],elen)
			}
			print(tmpstr)
			cat("\n")
		}
		print("0 indicates a subnode; * indicates 2+ subnodes overlapping. Positions are relative.")
	}
}



#-----------------------------
# Set Generics
#-----------------------------

# show
#!#setMethod("show","phylo4d_ext", function(object){ printphylo4(object); showSubNodes(object)}) # printphylo4 changed to print 07.27.15 to adapt to new phylobase package
setMethod("show","phylo4d_ext", function(object){ print(object); showSubNodes(object)})

setMethod("snid", signature(x="phylo4d_ext"),
  function(x) {
	return(x@subnode.id)
})

setMethod("snid", signature(x="list"),
	function(x) {
		return(x[[1]]@subnode.id)
})

setMethod("sndata", signature(x="phylo4d_ext"),
  function(x) {
	return(x@subnode.data)
})

setMethod("sndata", signature(x="list"),
  function(x) {
	return(x[[1]]@subnode.data)
})


# this adds data
# TODO: Change 'adds' behavior to 'replaces'
setReplaceMethod("sndata", signature(x="phylo4d_ext"),
  function(x,datnames=NULL,value) {
	
	if(!is.data.frame(value) && (length(value) != nSubNodes(x)))
	{
	  warning("Can only add vectors of data if they are the same length as the number of subnodes")
	  return(x)
	}
	if(!is.data.frame(value)){
		value = data.frame(value)
	}
	if(!is.null(datnames))
	{
		if(is.character(datnames) && length(datnames)==ncol(value))
		{
			names(value) <- datnames
		} else {
			warning("Not using specified names")
		}
	}
	x@subnode.data = cbind(x@subnode.data,value)
	return(x)
})


setReplaceMethod("sndata",signature(x="list"),
	function(x,datnames=NULL,value) {
		for(ii in seq(length(x)))
			sndata(x[[ii]],datnames=datnames) <- value
			
	return(x)
})


# begin rmdata:
# 'rmdata' method removes data columns from the 'data' and 'subnode.data' slots
#

setMethod("rmdata", signature(x="phylo4d_ext",index="numeric",subindex="numeric"),
  function(x,index,subindex) {
	
	if(length(index)>0 && abs(index) <= ncol(tdata(x)) && index != 0)
	{
		x@data = x@data[,-index,drop=F]
	} else {
		warning("The data column to be removed could not be found in tdata(x)")
	}

	if( hasSubNodes(x) && length(subindex) >0 && abs(subindex) <= ncol(sndata(x)) && subindex != 0)
	{
		x@subnode.data = x@subnode.data[,-index,drop=F]
	}
	
	return(x)
})
	

setMethod("rmdata", signature(x="phylo4d_ext",index="numeric",subindex="missing"),
  function(x,index) {
	return( rmdata(x,index=index,subindex=index) )
})

setMethod("rmdata", signature(x="phylo4d_ext",index="character",subindex="missing"),
  function(x,index) {
	ind = which(colnames(tdata(x)) %in% index)
	if(hasSubNodes(x)){
		subind = which(colnames(sndata(x)) %in% index)
	} else {
		subind = numeric(0)
	}
	return( rmdata(x,index=ind,subindex=subind) )
})


setMethod("rmdata", signature(x="phylo4d",index="numeric",subindex="missing"),
  function(x,index) {
	  if(length(index)>0 && abs(index) <= ncol(tdata(x)) && index != 0)
		x@data = x@data[,-index,drop=F]
	
	return(x)
})

setMethod("rmdata", signature(x="phylo4d",index="character",subindex="missing"),
  function(x,index) {
	inds = which(names(tdata(x)) %in% index)
	return(rmdata(x,inds))
})


setMethod("rmdata", signature(x="list",index="ANY",subindex="missing"),
	function(x,index) {
		x = sapply(x,rmdata,index)
	return(x)
})
# end rmdata


setMethod("snposition", signature(x="phylo4d_ext"),
  function(x) {
	return(x@subnode.pos)
})

setMethod("snposition", signature(x="list"),
  function(x) {
	return(x[[1]]@subnode.pos)
})

setMethod("snbranch", signature(x="phylo4d_ext"),
  function(x) {
	return(x@subnode.branch)
})

setMethod("snbranch", signature(x="list"),
  function(x) {
	return(x[[1]]@subnode.branch)
})


#---------------
## WEIGHT:
#
#---------------
setMethod("weight", signature(x="phylo4d_ext"),
  function(x) {
	return(x@weight)
})

setMethod("weight", signature(x="list"),
	function(x) {
		if(hasWeight(x))
			return( sapply(x,weight) )
		
		return(numeric(0))
})

setReplaceMethod("weight", signature(x="phylo4d_ext"),
  function(x,value) {
	x@weight = value
	return(x)
})

setReplaceMethod("weight",signature(x="list"),
	function(x,value) {
	
	if(length(x) != length(value))
		stop("Replacement values need to be same length as the list")
	
	for(ii in seq(length(x)))
		weight(x[[ii]]) <- value[ii]
	
	return(x)
})

setMethod("hasWeight",signature(x="phylo4d_ext"),
	function(x,strict=TRUE){
		return( length(x@weight)!=0 )
})


setMethod("hasWeight",signature(x="list"),
	function(x,strict=TRUE){
		retbool=ifelse(strict,TRUE,FALSE)
		for(ii in seq(length(x))){
			tmpbool = hasWeight(x[[ii]])
			if(strict && !tmpbool)
				return(FALSE)
			
			if(!strict && tmpbool)
				return(TRUE)
			
		}
		return(retbool)
})

# 
# setMethod("addData", signature(x="phylo4d_ext"),
# 	function(x,...,snode.data=NULL) {
# 		
# 		oldcols = names(tdata(x))
# 		# Add data the normal way:
# 		x = getMethod("addData","phylo4d")(x,...)
# 		
# 		newcols = setdiff(names(tdata(x)),names(sndata(x)))
# 		supercols = setdiff(oldcols,names(tdata(x)))
# 		
# 		cat("newcols: ", newcols,"\n")
# 		
# 		# Add data to subedges:
# 		if(is.null(snode.data))
# 		{
# 			if(length(newcols)!=0)
# 			{
# 				newdat = data.frame(matrix(NA,nrow=nSubNodes(x),ncol=length(newcols)))
# 				names(newdat) <- newcols
# 				x@subnode.data = cbind(x@subnode.data,newdat)
# 			}
# 			
# 			# Remove superfluous columns
# 			if(length(supercols)!=0)
# 			{
# 				for(scol in supercols)
# 					x = rmdata(x,scol)
# 			}
# 		}else{
# 			# add this data:
# 			if(length(newcols)!=0)
# 			{
# 				print("adding this way")
# 				sndata(x,datnames=newcols) <- snode.data
# 			}
# 		}
# 		
# 		return(x)
# })
# 
# 








#----------------------------------------
#  Phylo4d extension plots
#----------------------------------------

# 
# @param usestate specifies which column index for the @data slot should be accessed and plotted.  
#		 It can be either an integer (between 1...ncol) or a character string (which is a column 
#		 name.)
#
phyextPlot <- function(x,states,states.col,
						states.na="none", 
						usestate=1,
						plot.subnodes=T,
						plot.points=T,
						line.widths,line.types, ... )
{

	# TODO: sanity check to make sure usestate is a proper data column index
	
	# plot base phylogeny using phylobase functions:
	junk <- x
		
	gtree <- extractTree(junk)  # TODO: overload this so that it works on phyext class (keeps subnode info)
	posi = phyloXXYY(gtree)
	grid.newpage()
	pushViewport(viewport(layout=grid.layout(nrow=1, ncol=1), name="base"))
	pushViewport(viewport(layout.pos.col=1, name="plot1"))
		
	treePlot(gtree, newpage=FALSE,...)

	if(hasTipData(junk) || hasNodeData(junk))
	{
		if(missing(states))
		{
			tmp = as.character(unique(tdata(junk,'all')[,usestate,drop=F])[,1])
			tmp = tmp[!is.na(tmp)]
			if(hasSubNodes(junk) && plot.subnodes)
			{
				tmp = c(tmp,as.character(unique(sndata(junk)[,usestate,drop=F])[,1]))
				tmp = tmp[!is.na(tmp)]
			}
			tmp = unique(tmp)
			states = c(states.na,tmp)
			states.col = c(1,seq(from=2,length.out=length(states)-1))
		}
		
		if(missing(line.widths)){
			line.widths = rep(1,length(states))
		} else {
			if(length(states) != length(line.widths))
				stop("line.widths need to be the same length as states")
		}
		
		if(missing(line.types)){
			line.types = rep(1,length(states))
		} else {
			if(length(states) != length(line.types))
				stop("line.types need to be the same length as states")
		}
				
		seekViewport("tree")
		
		eord = edges(gtree)[posi$eorder,] # this is the order used
		treedata = tdata(junk,"all")[eord[,2],usestate,drop=T]  
		datamap = sapply(treedata,function(i) which(states == i),simplify=T)
		if(is.list(datamap))
			datamap = unlist(lapply(datamap, function(i) ifelse(length(i)==0,1,i[1])))
		
		# replot edges:
		grid.segments(posi$segs$h0x,posi$segs$h0y, posi$segs$h1x,posi$segs$h1y,gp=gpar(col=states.col[datamap],lwd=line.widths[datamap],lty=line.types[datamap]))
		grid.segments(posi$segs$v0x,posi$segs$v0y, posi$segs$v1x,posi$segs$v1y,gp=gpar(col=states.col[datamap],lwd=line.widths[datamap],lty=line.types[datamap]))
		if(plot.points) grid.points(posi$xx,posi$yy,pch=20,gp=gpar(col=states.col[datamap],cex=0.5))
		
		
		# plot sub nodes:
		if(hasSubNodes(junk) && plot.subnodes)
		{
			esub = matrix(edges(junk)[getSubNodeEdgeInds(junk),],ncol=2)
			posi.inds = apply(esub,1,function(i) which(i[1] == eord[,1] & i[2] == eord[,2]))
			subdata = getSubNodeData(junk,usestate)
			submapping = sapply(subdata[,1],function(i) which(states == i))
			subposi = getSubNodePosish(junk)
			
			get.x.offset <- function(xxyy,inds)
			{
				apply(cbind(xxyy$segs$h0x[inds],xxyy$segs$h1x[inds]),1,diff)
			}
			
			# reorder (so that lines don't completely cover each other):
			neword = order(rowMeans(subposi),decreasing=T)
			subposi = matrix(subposi[neword,],ncol=2)
			submapping = submapping[neword]
			posi.inds = posi.inds[neword]
			
			# subbranch positions:
			subposi.x0 = posi$segs$h0x[posi.inds]
			subposi.y0 = posi$segs$h0y[posi.inds]
			subposi.y1 = posi$segs$h1y[posi.inds]
			subposi.x1 = posi$segs$h0x[posi.inds] + (get.x.offset(posi,posi.inds) * rowMeans(subposi))
			subposi.vx0 = posi$segs$v0x[posi.inds]
			subposi.vy0 = posi$segs$v0y[posi.inds]
			subposi.vy1 = posi$segs$v1y[posi.inds]
			subposi.vx1 = posi$segs$v1x[posi.inds]
			
			# plot subnodes:
			grid.segments(subposi.x0,subposi.y0,subposi.x1,subposi.y1,gp=gpar(col=states.col[submapping],lwd=line.widths[submapping],lty=line.types[submapping]))
			grid.segments(subposi.vx0,subposi.vy0,subposi.vx1,subposi.vy1,gp=gpar(col=states.col[submapping],lwd=line.widths[submapping],lty=line.types[submapping]))
			if(plot.points) grid.points(subposi.x1,subposi.y1,pch=20,gp=gpar(col=states.col[submapping],cex=0.5))
			
		}
		
		upViewport(2)
	}
}

# The line below is unnecessary since phylobase sets 'plot' as a generic function
# (through graphics I think)
#setGeneric('plot')
setMethod('plot', signature(x='phylo4d_ext', y='missing'), function(x, y, ...) {
    phyextPlot(x, ...)
})







#---------------------------------------------
# Process Nexus files  
# -	Extra methods for extracting different information
# 	from nexus-formatted files.
#
#---------------------------------------------


# Strip comments from a tree string
#
.strip.tree.comments <- function(text=NULL)
{
	
}


# Method to read the first comment in a line in the format '[&...]'
# This is for reading tree weights chiefly
#
# Example:
# get.nexus.comments("example.txt")->lala
#
get.nexus.comments<-function(finput,text=NULL)
{
	
	if(!is.null(text)){
		# TODO: check text for newlines and split on them if they exist.
		rawtext=text
	} else {
		if(!file.exists(finput))
			stop("Assuming finput is a file and could not find it")
		
		rawtext = scan(finput,what=character(0),strip.white=T,sep="\n")		
	}
	
	# TODO: return named pair {treename, comment}
	comments = character(0)
	for(ii in seq(length(rawtext)))
	{
		junk =  gsub("^.*?\\[(.*?)\\].*$","\\1",rawtext[ii])
		if(length( grep("^&(.*)$",junk) ) != 0)
			comments = append(comments,junk)
	}
	
	return(comments)	
}


has.weights <- function(finput,text=NULL)
{
	if(!is.null(text)){
		# TODO: check text for newlines and split on them if they exist.
		rawtext=text
	} else {
		if(!file.exists(finput))
			stop("Assuming finput is a file and could not find it")
		
		rawtext = scan(finput,what=character(0),strip.white=T,sep="\n")		
		rawtext = read.nexus.block(txt=rawtext,block="trees",silent=T)
	}
	
	comments = grep("(&lnP|&W)( |=)",rawtext)
	return( length(comments)>0 )
}


# get tree weights from file or tree string
# Assuming this format: 
# [.... &W -122235 ........]
#
get.tree.weights <- function(finput,text=NULL,starters=c("&W","&lnP"),splitchar="( |=)")
{
	if(!is.null(text)){
		# TODO: check text for newlines and split on them if they exist.
		rawtext=text
	} else {
		if(!file.exists(finput))
			stop("Assuming finput is a file and could not find it")
		
		rawtext = scan(finput,what=character(0),strip.white=T,sep="\n")		
		rawtext = read.nexus.block(txt=rawtext,block="trees",silent=T)
	}
	
	comments = get.nexus.comments(text=rawtext)
	tokens = strsplit(comments,splitchar)
	find.weight <- function(ii,ww)
	{
		return(ii[(which(toupper(ii) %in% toupper(ww)))+1])
	}
	weighs = unlist(lapply(tokens,find.weight,starters))
	
	return(as.numeric(weighs))
}



# Internal function - split tokens by a certain character
.split.tokens <- function(txt,char)
{
	nlines = length(txt)
	newtokens = character(nlines)
	charlines = unname(sapply(txt,function(i) length(grep(char,i))))
	curline=1
	for(ll in seq(nlines))
	{
		if(charlines[ll] == 0){
			newtokens[curline] = txt[ll]
			curline = curline + 1
		}else{
			tmp = strsplit(txt[ll],char)[[1]]
			for(tmpline in tmp){
				#if(tmpline!=""){
					#newtokens[curline] = paste(tmpline,char,sep="")
					newtokens[curline] = tmpline
					curline = curline + 1
				#}
			}
		}
	}
	return (newtokens)
}


# Internal function - get all content within a nexus block
.get.nexus.block.inds <- function(filename,blockname,text=NULL)
{
	# choose character vector
	if(!is.null(text))
	{
		filetext = text
	} else {
		filetext = scan(filename,what="character",sep="\n",strip.white=T)
	}
	filetext = tolower(filetext)
	
	search.for=paste(paste("begin",tolower(blockname),sep="\\s"),"[\\s]{0,};",sep="")
	start.ind = grep( search.for ,filetext,ignore.case=T)
	if(length(start.ind) == 0)
		return (integer(0))
	
	end.ind = grep("end;",filetext,ignore.case=T)
	end.ind = end.ind[head(which(end.ind > start.ind),1)]

	return( c((start.ind),(end.ind)) )
}



# Internal:
# Get text of a read nexus block
read.nexus.block <- function(finput,txt=NULL,block,rm.comments=F,silent=F)
{
	if(!is.null(txt)){
		# Using the text argument is not recommended
		rawtext=txt
	} else {
		
		if(!file.exists(finput))
			stop("Assuming finput is a file and could not find it")
		
		rawtext = scan(finput,what=character(0),strip.white=T,sep="\n",quiet=silent)		
	}
	
	inds = .get.nexus.block.inds(blockname=block,text=rawtext)
	if(length(inds)==0)
	{
		if(!silent) warning(paste("This file has",block, "no block"))
		return (character(0))
	}
	
	rawtext = rawtext[(inds[1]+1):(inds[2]-1)]
	# TODO: split up newlines if they exist
	
	if(rm.comments)
	{
		# assume comments are start a line with [ 
		# and end a line with ]
		#
		comment.starts = grep("^\\[",rawtext)
		comment.ending = grep("\\]$",rawtext)
		if(length(comment.starts) == length(comment.ending)){
			comment.pairs = cbind(comment.starts,comment.ending)
		} else {
			comment.pairs = cbind(character(0),character(0))
		}
		
		if(nrow(comment.pairs)>0)
		{
			for(pair in seq(nrow(comment.pairs)))
			{
				print(pair)
				rawtext = rawtext[-(comment.pairs[pair,1]:comment.pairs[pair,2])]
			}
		}
	}
	
	return (rawtext)
}


has.block <- function(finput,txt=NULL,blockname="characters2")
{
	retbool = (length(read.nexus.block(finput,txt,block=blockname,silent=T))!=0)
	return(retbool)	
}

has.characters2 <- function(finput,txt=NULL)
{
	return(has.block(finput,txt,"characters2"))
}



# Alternative way to read in characters
# Also works as a work around to readNexus' annoying habit
# or crashing when using type="data" but where trees block
# contains simmap formatted trees
read.characters2 <- function(finput,txt=NULL,blockname="characters2")
{
	rawtext = NULL
	if(!is.null(txt))
	{
		rawtext = txt
	} else {
			if(!file.exists(finput)){
				stop("Assuming finput is a file and could not find it")
			} else {
				rawtext = scan(finput,what=character(0),strip.white=T,sep="\n")						
			}
	} 
	
	tmpfile = tempfile()
	tmphead = "#NEXUS\n"
	tmptaxa = c("BEGIN TAXA;",read.nexus.block(txt=rawtext,block="taxa"),"END;")
	tmptext = read.nexus.block(txt=rawtext,block=blockname)
	tmptext = c("BEGIN CHARACTERS;",tmptext,"END;")
	writeLines(c(tmphead,tmptaxa,tmptext),con=tmpfile)
	data2.part = readNexus(tmpfile,type="data",levels.uniform=F)
	rownames(data2.part) <- checkLabel(rownames(data2.part))
	return(data2.part)
}










#------------------------
# Coersing functions	|
#------------------------

# pulled from implicit(?) 'coerce<-' function
deep.phy.copy <- function (from, to, value,usedata=TRUE) 
{
	useslots = c("data", "metadata", "edge", "edge.length","label", "edge.label", "annote")	 # leaving 'order' out...
	if(!usedata) 
		useslots = useslots[-c(1,2)]
    for (what in useslots) 
    	slot(from, what) <- slot(value, what)
    from
}

# up convert:
# NOTE: not sure if it is necessary to define all of these...
# 		there might be a more 'S4' way to do this.
#
setAs('phylo','phylo4d_ext', function(from, to) {
	return(phyext(from))
})

setAs('phylo4','phylo4d_ext', function(from, to) {
	return(phyext(from))
})

setAs('phylo4d','phylo4d_ext', function(from, to) {
	return(phyext(from))
})

# 
# # down convert:
# # TODO: figure out if these functions are necessary:
# #
# setAs('phylo4d_ext','phylo4d',function(from,to) {
# 
# 	newphy = deep.phy.copy(new(to,order="unknown"),to,from,TRUE)
# 	newphy
# })
# 
# 
# setAs('phylo4d_ext','phylo4',function(from,to) {
# 	
# 	newphy = deep.phy.copy(new(to,order="unknown"),to,from,FALSE)
# 	newphy
# })
# 
# 
# 
# setAs('phylo4d_ext','phylo',function(from,to) {
# 	
# 	# TODO: Need to convert subnodes back into singlestons, 
# 	#		and then collapse the singletons.
# 	return(as(as(from,"phylo4d"),"phylo"))
# })
# 
# 






#!# This section commented out on inclusion in SigTree package, 7/1/15
#!## use this with a namespace
#!#.onLoad <- function(lib, pkg) {
#!#    require(phylobase)
#!#    #!#require(methods) # not sure if this is needed
#!#    phylobase.options(singleton="ok")
#!#}

