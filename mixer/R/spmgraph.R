setClass("spmgraph",representation(name="character",edgesrcs="numeric",edgedests="numeric",nodenames="character",nodeindexes="numeric"),
                    prototype(name='Unnamed',edgesrcs=vector(mode="numeric",length=0),edgedests=vector(mode="numeric",length=0),
                              nodenames=vector(mode="character",length=0),nodeindexes=vector(mode="numeric",length=0)))

addedge.spmgraph<-function(dest,srcindex,initenv) {
	g=get("instance",envir=initenv)
	destindex=which(g@nodenames==dest)
	if (length(destindex)==0) {
		destindex=length(g@nodenames)+1
		g@nodenames[destindex]<-dest
		g@nodeindexes[dest]<-destindex
	}
	g@edgesrcs<-c(g@edgesrcs,srcindex)
	g@edgedests<-c(g@edgedests,destindex)
	assign("instance",g,envir=initenv)
}

parseline.spmgraph<-function(line,initenv) {
	fields=unlist(strsplit(line,"[[:space:]]+",perl=TRUE))
	if (length(fields)>1) {
		src=fields[1]
		dest=fields[2:length(fields)]
		g=get("instance",envir=initenv)
		nextindex=length(g@nodeindexes)+1	
		srcindex=which(g@nodenames==src)
		if (length(srcindex)==0) {
			srcindex=length(g@nodenames)+1
			g@nodenames[srcindex]<-src
			g@nodeindexes[src]<-srcindex
		}
		assign("instance",g,envir=initenv)
		lapply(dest,'addedge.spmgraph',srcindex,initenv)
	}
}

readfromfile.spmgraph<-function(filename,initenv) {
	lines=readLines(filename)
	lapply(lines,'parseline.spmgraph',initenv)
	return(TRUE)
}


addadjmatedge.spmgraph<-function(index,initenv,nnodes)
{
	index=index-1
	g=get("instance",envir=initenv)
	col=as.integer(index/nnodes)
	row=index-col*nnodes+1
	col=col+1
	srcindex=which(g@nodenames==row)
	if (length(srcindex)==0) {
		srcindex=length(g@nodenames)+1
		g@nodenames[srcindex]<-row
		g@nodeindexes[row]<-srcindex
	}
	destindex=which(g@nodenames==col)
	if (length(destindex)==0) {
		destindex=length(g@nodenames)+1
		g@nodenames[destindex]<-col
		g@nodeindexes[col]<-destindex
	}
	g@edgesrcs<-c(g@edgesrcs,srcindex)
	g@edgedests<-c(g@edgedests,destindex)
	assign("instance",g,envir=initenv)
}


buildfromadjacencymatrix.spmgraph<-function(adjmat,initenv) {
	nonnulls=which(adjmat!=0)
	nnodes=dim(adjmat)[1]
	lapply(nonnulls,'addadjmatedge.spmgraph',initenv,nnodes)
	return(TRUE)
}
	

setMethod("initialize",
	"spmgraph",
	function(.Object,initialdata=NULL) {
                initialized=FALSE
		if (is.character(initialdata) == TRUE) {
			.Object@name=initialdata
			initenv=new.env()
			assign("instance",.Object,envir=initenv)
			initialized=readfromfile.spmgraph(initialdata,initenv)
			.Object=get("instance",envir=initenv)
                }
		if (is.matrix(initialdata) == TRUE) {
			initenv=new.env()
			assign("instance",.Object,envir=initenv)
			initialized=buildfromadjacencymatrix.spmgraph(initialdata,initenv)
			.Object=get("instance",envir=initenv)
		}
		if (initialized == FALSE) {
			warning("Could not correctly initialize spmgraph object.")
		}
		.Object
        })

setGeneric("getEdges",function(g,...) { standardGeneric("getEdges")})

indextolabel.spmgraph<-function(index,g) {
	g@nodenames[index]
}

setMethod("getEdges",
	"spmgraph",
	function(g,withlabels=FALSE) {
		edgelist=c(g@edgesrcs,g@edgedests)
		if (withlabels==TRUE)
			edgelist=lapply(edgelist,'indextolabel.spmgraph',g)
		matrix(edgelist,nrow=2,byrow=T)
	})

setGeneric("getAdjacencyMatrix",function(g,...) { standardGeneric("getAdjacencyMatrix")})

addedgetoadjacencymatrix.spmgraph<-function(index,g,adjmatenv)
{
	adjmat=get("adjmat",envir=adjmatenv)
	src=g@edgesrcs[index]
	dest=g@edgedests[index]
	adjmat[src,dest]=1
	adjmat[dest,src]=1
	assign("adjmat",adjmat,envir=adjmatenv)
}

setMethod("getAdjacencyMatrix",
	"spmgraph",
	function(g,withlabels=FALSE) {
		nnodes=length(g@nodenames)
		adjmat=matrix(rep(0,nnodes*nnodes),nrow=nnodes)
		adjmatenv=new.env()
		assign("adjmat",adjmat,envir=adjmatenv)
		lapply(seq(length(g@edgesrcs)),'addedgetoadjacencymatrix.spmgraph',g,adjmatenv)
		adjmat=get("adjmat",envir=adjmatenv)
		if (withlabels==TRUE) {
			nodenames=lapply(1:length(g@nodenames),'indextolabel.spmgraph',g)
			colnames(adjmat)<-nodenames
			rownames(adjmat)<-nodenames
		}
		return(adjmat)
	})

