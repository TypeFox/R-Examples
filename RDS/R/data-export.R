


#' Export an rds.data.frame to file
#' @param x The rds.data.frame to export
#' @param file The name of the file to create.
#' @export
write.rdsobj <- function(x,file){
	if(!inherits(x,"rds.data.frame"))
		stop("data is not an rds.data.frame object")
	saveRDS(x,file)
}


#' writes an rds.data.frame recruitment tree as a GraphViz file
#' @param x An rds.data.frame.
#' @param file A character vector representing the file
#' @export
write.graphviz <- function(x,file){
	x <- as.rds.data.frame(x)
	id <- as.char(get.id(x))
	rid <- as.char(get.rid(x))
	sid <- as.char(get.seed.rid(x))
	w <- get.wave(x)
	el <- cbind(rid,id)
	el <- el[order(w),]
	el <- el[rid!=sid,]
	
	cat("digraph rds {\nsize=\"8.5,8.5\";\nlayout=neato;\nmode=\"hier\";",file=file)
	cat(paste0("\"",el[,1],"\" -> \"",el[,2],"\"",collapse=";\n"),
			file=file,append=TRUE)
	cat("}",file=file,append=TRUE)
}




#' Writes out the RDS tree in NetDraw format
#' @param x An rds.data.frame.
#' @param file a character vector representing a file.
#' @param by.seed If true, seperate files will be created for each seed.
#' @details If by.seed is false, two files are created using 'file' as a base name.
#' \code{paste0(file,".DL")} contains the edge information, and \code{paste0(file,".vna")}
#' contains the nodal attributes
#' @export
write.netdraw <- function(x,
	file=NULL,
	by.seed=FALSE){
rds.data <- as.rds.data.frame(x)

id <- get.id(rds.data)
recruiter.id <- get.rid(rds.data)
seed.id <- get.seed.id(rds.data)
seed.rid <- get.seed.rid(rds.data)


##################################################################################################
# function to write the nodes and edges in the DL language.
write.parent.child.connections <- function(idx,participants){
	children <- get.children(idx)
	children <- children[!is.na(children)]
	
	for(child.idx in children[!is.na(children)]){            
		if(length(unlist(children)) > 0){
			cat(sprintf("%s %s\n",
							match(idx,participants),
							match(child.idx,participants)),
					file=DL.file)
		}
	}
	
}


##################################################################################################
# Functions to navigate the tree structure. 

get.children <- function(idx){
	id[recruiter.id == idx]
}

get.recruiter.row <- function(idx,rid){
	rid[id == idx]
}  

get.child.idx <- function(child,rid=get.rid(rds.data)){
	parent <- get.recruiter.row(child,rid)
	if(parent != seed.rid)
		siblings <- get.children(parent)
	return((1:length(siblings))[siblings == child])
}





############################################################################
#   #  Now we create the temporary DL file.
if(is.null(file)){
	file.base<-paste(getwd(),substitute(rds.data),sep="/")
}else{
	file.base <- file
}


#Write the DL header
header.to.DL <- paste(
		'DL n = ',nrow(rds.data),', format = edgelist1\n','labels:\n',
		sep="")

# Now we write the nodes and the edges.    
if(by.seed){
	seeds <- sort(unique(seed.id))
	for(seed in seq(along=seeds)){
		n <- sum(seed.id == seeds[seed])
		header.to.DL <- paste(
				'DL n = ',n,', format = edgelist1\n','labels:\n',
				sep="")
		if(n<2)
			next
		#Write edges
		DL.file <- file(sprintf("%s_Seed%s.DL",file.base,seeds[seed]),"wt")
		cat(header.to.DL,file=DL.file,append=TRUE,fill=FALSE)
		participants <- id[seed.id == seeds[seed]]

		cat(paste("'",paste(participants,collapse="','"),"'",sep=""),file=DL.file,
				append=TRUE,fill=FALSE)
		cat(paste('\n','data:\n',sep=""),file=DL.file,append=TRUE,fill=FALSE)
		sapply(participants,write.parent.child.connections,participants)
		close(DL.file)
		#Write attribute file
		ADL.file <- sprintf("%s_Seed%s.vna",file.base,seeds[seed])
		rd1 <- cbind(ID=id,rds.data)[seed.id == seeds[seed],]
		cat("*node data\n",file=ADL.file,append=FALSE,fill=FALSE)
		suppressWarnings(
				utils::write.table(rd1,file=ADL.file,append=TRUE,row.names=FALSE,
						eol="\r\n",sep="\t")
		)
	}
}else{
	#Write edges
	DL.file <- file(sprintf("%s.DL",file.base),"wt")
	cat(header.to.DL,file=DL.file,append=FALSE,fill=FALSE)
	participants <- id

	cat(paste("'",paste(participants,collapse="','"),"'",sep=""),
			file=DL.file,append=TRUE,fill=FALSE)
	cat(paste('\n','data:\n',sep=""),file=DL.file,append=TRUE,fill=FALSE)
	sapply(participants,write.parent.child.connections,participants)
	cat("\n",file=DL.file,append=TRUE,fill=FALSE)
	close(DL.file)
	#Write attribute file
	ADL.file <- sprintf("%s.vna",file.base)
	rd1 <- cbind(ID=id,rds.data)

	cat("*node data\n",file=ADL.file,append=FALSE,fill=FALSE)
	suppressWarnings(
			utils::write.table(rd1,file=ADL.file,append=TRUE,row.names=FALSE,
					eol="\r\n",sep="\t")
	)
}


invisible()   
}     


#' Writes out the RDS tree in RDSAT format
#' @param x An rds.data.frame.
#' @param file a character vector representing a file.
#' @export
write.rdsat <- function(x,
	file=NULL){
rds.data <- as.rds.data.frame(x)

id <- get.id(rds.data)
recruiter.id <- get.rid(rds.data)
        max.coupons <- length(tabulate(table(recruiter.id)))

############################################################################
#   #  Now we create the temporary RDSAT file.
if(is.null(file)){
	file.base<-paste(getwd(),substitute(rds.data),sep="/")
}else{
	file.base <- file
}

get.children <- function(idx){
	id[recruiter.id == idx]
}

        coupons <- matrix("",ncol=max.coupons,nrow=nrow(rds.data))
        for(idx in id){
 children <- get.children(idx)
 children <- children[!is.na(children)]
         if(length(children)>0){
          coupons[match(idx,id),seq_along(children)] <- children
         }
        }
        coupons[coupons==""] <- paste(nrow(rds.data)+(1:sum(coupons=="")))
        network.size <- attr(rds.data,"network.size")
        full.rds <- cbind(id,rds.data[,network.size],id,coupons,rds.data)
        colnames(full.rds)[2] <- "network.size"

#Write the RDSAT header
header.to.RDSAT <- paste(
		'RDS\n',nrow(rds.data)," ",max.coupons," 0 ",
                        paste(colnames(rds.data),collapse=" "),'\n',
		sep="")

# Now we write the nodes and the edges.    
	#Write edges
RDSAT.file <- file(sprintf("%s.rdsat",file.base),"wt")
cat(header.to.RDSAT,file=RDSAT.file,append=FALSE,fill=FALSE)

utils::write.table(full.rds,
	file=RDSAT.file,quote=FALSE,append=TRUE,row.names=FALSE,col.names=FALSE)
close(RDSAT.file)

invisible()   
}     
