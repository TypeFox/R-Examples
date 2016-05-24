buildSNPannotation<-function(pkg,rs=TRUE,allele=TRUE,gene=TRUE,chromosome=FALSE,
		position=FALSE,strand=FALSE,cytoband=FALSE,max.genes=0,lib.loc=NULL,
		others=NULL,subset=NULL,pattern=NULL,na.rm=TRUE){
	require(oligoClasses)
	require(pkg,character.only=TRUE,lib.loc=lib.loc) || stop(paste("Package",pkg,
		"not available."))
	conn<-db(get(pkg))
	what<-c("man_fsetid","dbsnp_rs_id","chrom","physical_pos","strand","cytoband",
		"allele_a","allele_b","gene_assoc")
	cn<-c("Probe-Set-ID","RefSNP","Chromosome","Position","Strand","Cytoband",
		"Allele_A","Allele_B","Gene")
	interest<-c(TRUE,rs,chromosome,position,strand,cytoband,allele,allele,gene)
	what<-what[interest]
	cn<-cn[interest]
	if(!is.null(others)){
		what<-unique(c(what,others))
		cn<-unique(c(cn,others))
	}
	if(any(!what%in%dbListFields(conn,"featureSet")))
		stop("Some of the specified annotations seem to be not available.")
	sql<-paste("SELECT", paste(what,collapse=", "), "FROM featureSet")
	if(!is.null(pattern))
		sql<-paste(sql," WHERE man_fsetid LIKE '",pattern,"'",sep="")
	out<-dbGetQuery(conn,sql)
	rn<-out[,1]
	rownames(out)<-rn
	colnames(out)<-cn
	str<-out$Strand
	if(any(str%in%(0:1))){
		str[str==0]<-"+"
		str[str==1]<-"-"
		out$Strand<-str
	}
	if(max.genes>0)
		out<-shortenGeneDescription(out,max.length=max.genes)
	if(!is.null(subset)){
		ids<-match(subset,rn)
		if(any(is.na(ids))){
			warning(sum(is.na(ids))," of the ",length(ids)," Probe-Set-IDs specified",
				" by 'subset'\n","are not available in the annotation package.")
			if(na.rm){
				subset<-subset[!is.na(ids)]
		 		ids<-ids[!is.na(ids)]
			}
		}
		out<-out[ids,]
		rownames(out)<-subset
	}
	out
}



shortenGeneDescription<-function(dat,colname="Gene",max.length=2,sep="///",add.ldots=TRUE){
	if(max.length<1)
		stop("max.length must be at least 1.")
	if(!is.data.frame(dat))
		stop("dat must be a data frame.")
	ids<-which(colnames(dat)==colname)
	if(length(ids)!=1)
		stop("Exactly one column of 'dat' must be called ",colname,".")
	genes<-strsplit(dat[,ids],sep)
	shorten<-function(x,limit,sep,add.ldots){
		le<-length(x)
		tmp<-paste(x[1:min(le,limit)],collapse=sep)
		if(add.ldots && le>limit)
			tmp<-paste(tmp,"...",sep=sep)
		tmp
	}
	dat[,ids]<-sapply(genes,shorten,max.length,sep,add.ldots)
	dat
}

	