"export.merlin" <- 
		function(data,pedfile="merlin.ped",datafile="merlin.dat",
				mapfile="merlin.map",format="merlin",fixstrand="no",
				extendedmap=TRUE,traits=1, order = TRUE, stepids = 100, dpieceFun="new") {
	if (!is(data,"gwaa.data")) stop("Data argument should be of gwaa.data-class")
	dpieceFunOptions <- c("old","new")
	dpieceFunInt <- match(dpieceFun,dpieceFunOptions,nomatch=0)
	if (dpieceFunInt<=0) {
		out <- paste("dpieceFun argument should be one of",dpieceFunOptions,"\n")
		stop(out)
	} else {
		if (dpieceFunInt==1) dump.piece=dump.piece
		else if (dpieceFunInt==2) dump.piece=dump.piece.New
		else stop("weird stuff!")
	}
	formats <- c("merlin","plink")
	if (!(match(format,formats,nomatch=0)>0)) {
		out <- paste("format argument should be one of",formats,"\n")
		stop(out)
	}
	fixes <- c("no","+","-")
	if (!(match(fixstrand,fixes,nomatch=0)>0)) {
		out <- paste("fixstrand argument should be one of",fixes,"\n")
		stop(out)
	}
	if (order) {
		ord <- sortmap.internal(chromosome(data),map(data))
		if (any(ord$ix != c(1:length(ord$ix)))) {data <- data[,ord$ix]; gc()}
	}
	if (fixstrand != "no") {
		###### as.raw(1) == "+"
		###### as.raw(2) == "-"
		if (fixstrand == "-") {tos <- as.raw(1);fos<-as.raw(2);} 
		else if (fixstrand == "+") {tos <- as.raw(2);fos <- as.raw(1)}
		else {stop("strange strand...")}
		stf <- which(as.raw(data@gtdata@strand) == tos) 
		if (length(stf)>0) {
			revs <- alleleID.revstrand()
			savnam <- names(data@gtdata@strand)
#
#		data@gtdata@coding[stf] <- new("snp.coding",as.raw(revs[as.character(data@gtdata@coding[stf])]))
#
#		ugly fix -- 
#
			data@gtdata@coding[stf] <- new("snp.coding",as.raw(revs[as.character(as.raw(data@gtdata@coding[stf]))]))
			names(data@gtdata@coding) <- savnam
			data@gtdata@strand[stf] <- new("snp.strand",rep(fos,length(stf)))
			names(data@gtdata@strand) <- savnam
		}
	}
	bstp <- stepids #100
	if (data@gtdata@nids>(bstp*1.5)) {
		steps <- seq(from=0,to=data@gtdata@nids,by=bstp)
		if (data@gtdata@nids != steps[length(steps)]) steps[length(steps)+1] <- data@gtdata@nids
		dump.piece(data=data,from=1,to=steps[2],traits=traits,
				pedfile=pedfile,append=FALSE,format=format)
		cat("dumped ids",1,"to",steps[2],"\n")
		for (jjj in c(2:(length(steps)-1))) {
			dump.piece(data=data,from=(steps[jjj]+1),to=(steps[jjj+1]),traits=traits,
					pedfile=pedfile,append=TRUE,format=format)
			cat("dumped ids",steps[jjj]+1,"to",steps[jjj+1],"\n")
		}
	} else {
		dump.piece(data=data,from=1,to=data@gtdata@nids,traits=traits,
				pedfile=pedfile,append=F,format=format)
	}
	snps <- data@gtdata@snpnames
	inf <- data.frame(a=rep("M",length(snps)),b=snps)
	if (traits>0) {
		tran <- paste("faket",seq(1:traits),sep="")
		inf0 <- data.frame(a=rep("T",traits),b=tran)
		inf <- rbind(inf0,inf)
	}
	if (!is.null(datafile)) {
		write.table(inf,file=datafile,col.names=FALSE,row.names=FALSE,quote=FALSE)
	}
	if (format=="merlin") {
		map <- data.frame(chromosome=as.character(data@gtdata@chromosome),markername=data@gtdata@snpnames,position=data@gtdata@map)
		write.table(map,file=mapfile,col.names=TRUE,row.names=FALSE,quote=FALSE)
	} else if (format=="plink") {
		map <- data.frame(chromosome=as.character(chromosome(data)),
				markername=snpnames(data),gpos=0,position=map(data))
		write.table(map,file=mapfile,col.names=FALSE,row.names=FALSE,quote=FALSE)
	} else {
		stop("non-formalized format")
	}
	if (extendedmap) {
		map <- data.frame(chromosome=as.character(data@gtdata@chromosome),markername=data@gtdata@snpnames,position=data@gtdata@map,strand=as.character(data@gtdata@strand),coding=as.character(data@gtdata@coding))
		write.table(map,file=paste(mapfile,".ext",sep=""),col.names=TRUE,row.names=FALSE,quote=FALSE)
	}
}

dump.piece <- function(data,fromid,toid,traits,pedfile,append,format="merlin") {
	if (toid < fromid) stop("toid<fromid")
	x <- as.character(data@gtdata[c(fromid:toid),])
	x <- replace(x,(x==""),"0/0")
	x <- replace(x,is.na(x),"0/0")
	if (format=="plink") {
		x <- sub("/"," ",x)
	}
	ids <- rownames(x)
	nids <- length(ids)
	sx <- data@phdata$sex[c(fromid:toid)]
	sx <- replace(sx,(sx==0),2)
	if (traits > 0) {
		tr <- matrix(rep(0,nids*traits),ncol=traits)
		x <- data.frame(seq(fromid,toid),ids,rep(0,nids),rep(0,nids),sx,tr,x)
	} else {
		x <- data.frame(seq(fromid,toid),ids,rep(0,nids),rep(0,nids),sx,x)
	}
	write.table(x,file=pedfile,col.names=FALSE,row.names=FALSE,quote=FALSE,append=append,sep=" ")
}

dump.piece.New <- function(data,fromid,toid,traits,pedfile,append,format="merlin") {
	if (toid < fromid) stop("toid<fromid")
	# ugly line -- should be in the C++ code
	if (!append) unlink(pedfile);
	if (format=="plink") {plink <- TRUE;} else {plink <- FALSE;}
	#data <- data[c(fromid:toid),]
	result <- .Call("export_plink",as.character(idnames(data)[c(fromid:toid)]),
			as.raw(gtdata(data)@gtps), as.integer(nsnps(data)), as.integer(nids(data)),
			as.character(coding(data)), as.integer(fromid), as.integer(toid),
			as.integer(male(data)), as.integer(traits), as.character(pedfile), 
			as.logical(plink), as.logical(append))
	#
	return(1)
}
