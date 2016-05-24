#.packageName<- "DSviaDRM"


  #####################################################################################################
  ## a gene filtering function: Genes which have a Between-Experiment Mean  Expression Signal (BEMES) lower than the median of BEMES's 
  ## of all genes will be filtered out. 
  ## 'exprs' is the expression data for two conditions.
  # output: A data frame or matrix with a reduced number of rows.
  #####################################################################################################
"expressionBasedfilter" <- function(exprs) {
 	all.mean <- apply(exprs,1,mean)
  	median <- median(all.mean)
  	exprs.filter <- exprs[all.mean>median,]
  	return(exprs.filter)
}

  #################################################################################################### 
  ## a gene filtering function: Those genes not significantly more variable than the median gene are filtered out.
  ## 'exprs' is the expression data for two conditions.
  ## 'p' is the probability  cut-off of Chi Squared distribution.
  # output: A data frame or matrix with a reduced number of rows.
  ####################################################################################################
"varianceBasedfilter" <- function(exprs,p) {
 	n <- ncol(exprs)
 	Var.i <- apply(exprs,1,var)
 	Var.median <- median(Var.i)
 	quantity <- (n-1)*Var.i/Var.median
 	degree <- ncol(exprs)-1
 	prob <- pchisq(quantity, degree, ncp=0, lower.tail = F, log.p = FALSE)
 	exprs.filter <- exprs[prob<=p,]
 	return(exprs.filter)
}


###############################################################################
## DCgene&link: for unveil differential coexpression genes and links based on DCp and DCe in DCGL;
## input: disease and control expression profile;
## output: differential coexpression genes and links;
## 
###############################################################################
"DCEA" <- function(exprs.1,exprs.2,r.method=c('pearson','spearman')[1],link.method=c('qth','rth','percent')[1],cutoff=0.25,N=0,N.type=c('pooled','gene_by_gene')[1],q.method=c("BH","holm", "hochberg", "hommel", "bonferroni", "BY","fdr")[1],nbins=20,p=0.1) {
	cor.filtered.1<-cor.filtered.2<-NULL
	if (nrow(exprs.1)!=nrow(exprs.2)) stop('the two expression matrices must have the same number of rows (genes).')
	if (length(rownames(exprs.1))==0 | length(rownames(exprs.2))==0) stop('the expression matrices must have row names specifying the gene names.')
	if ( min(ncol(exprs.1),ncol(exprs.2))<3 ){
		stop('each expression matrix must have at least three or more columns.')
	} else if (min(ncol(exprs.1),ncol(exprs.2))<5 ) {
		warning('the minimum number of columns is less than five and the result may not be reliable.')
	}
	m <- nrow(exprs.1) # exprs.1, exprs.2 is the expression data for different conditions.
  if (m>5000) warning('the number of genes exceeds 5000 and the program may takes long time to run.')
	genes = rownames(exprs.1)
  	
	cor.filtered = switch(link.method,
		rth = rLinkfilter(exprs.1,exprs.2,r.method=r.method,cutoff=cutoff),
	 	qth =  	qLinkfilter(exprs.1,exprs.2,r.method=r.method,cutoff=cutoff),
		percent = percentLinkfilter(exprs.1,exprs.2,r.method=r.method,cutoff=cutoff)
	)
	rth.1 = cor.filtered$rth.1
	rth.2 = cor.filtered$rth.2
	cor.filtered.1 = cor.filtered$cor.filtered.1
	cor.filtered.2 = cor.filtered$cor.filtered.2
 
  #####################################################################################################################################################
  ## For one gene, there are n-1 correlation value pairs also q value pairs(From the two conditions). For a correlation pair for this gene,
  ## if one of the q values is less than the threshold, the pair will be retained. Then there are 'number.i.uniq' unique pairs retained that is there 
  ## are two vectors of correlation values. 
  ## Then a length normalized Euclidean Distance for the two vectors will be calculated (LNED). 
  ##################################################################################################################################################### 

	dC.length = calc.dC(cor.filtered.1,cor.filtered.2,genes)	
	dC = dC.length$dC
	number_uniq = dC.length$length
	
 ########################################################################################################################
 ## Disturb the sample labels for the two conditions and re-assign the samples to two datasets,then calculate the 'dC0' for 
 ## N times and then pool all the dC0 together to construct a 'NULL' distribution.
 #########################################################################################################################
	if(N>0){
		dC0 <- matrix(nrow=length(genes),ncol=N)
		rownames(dC0) <- genes
		exprs <- cbind(exprs.1,exprs.2)
		expSamples <- colnames(exprs)
		n.1 = ncol(exprs.1)
		n.2 = ncol(exprs.2)
		cat.j = 0
		for(j in 1:N) {
			if ( (j*100/N)%/%10>cat.j) {
				cat.j = cat.j+1
				cat(cat.j*10,'%','\n')
			}
			seq <- sample(n.1+n.2)
			exprs.1 <- exprs[,seq[1:n.1]];
			exprs.2 <- exprs[,seq[(n.1+1):(n.1+n.2)]];
			rownames(exprs.1) <- rownames(exprs.2) <- genes

			cor.filtered = switch(link.method,
                		rth = rLinkfilter(exprs.1,exprs.2,r.method=r.method,cutoff=cutoff),
                		qth =   qLinkfilter(exprs.1,exprs.2,r.method=r.method,cutoff=cutoff),
                		percent = percentLinkfilter(exprs.1,exprs.2,r.method=r.method,cutoff=cutoff)
        		)
			dC0.j <- calc.dC(cor.filtered$cor.filtered.1,cor.filtered$cor.filtered.2,cor.filtered$genes)
  		dC0[,j] = dC0.j$dC
		}

		p.value.DCp = switch(N.type,
			gene_by_gene = apply(cbind(dC0,dC),1,function(x) sum(x[1:(length(x)-1)]>x[length(x)],na.rm=T)/length(!is.na(x[1:(length(x)-1)]))),
			pooled = sapply(dC,function(x) sum(as.vector(dC0)>x,na.rm=T)/length(!is.na(as.vector(dC0))))
		)
  		q.value.DCp <- p.adjust(p.value.DCp,method=q.method)
  		Result.DCp <- data.frame(dC=dC,links=number_uniq,p.value=p.value.DCp,q.value=q.value.DCp);
  		row.names(Result.DCp) <- genes;
	} else { 
  		Result.DCp<- data.frame(dC=dC,links=number_uniq,p.value=rep(NA,length(dC)),q.value=rep(NA,length(dC)));
  		row.names(Result.DCp) <- genes;
  	}
  	
#############################################################
## decide three sets of correlation pairs and organize them into two-columned matrices.
#############################################################  	
  	
  	idx.same = (cor.filtered.1*cor.filtered.2)>0; idx.same[is.na(idx.same)] <- TRUE  ##fixing special cases where cor = NA (caused by at least one constant gene expression vector)
  	idx.diff = (cor.filtered.1*cor.filtered.2)<0; idx.diff[is.na(idx.diff)] <- FALSE
  	idx.switched = (cor.filtered.1*cor.filtered.2<0)& ( abs(cor.filtered.1)>=rth.1 & abs(cor.filtered.2)>=rth.2 ); idx.switched[is.na(idx.switched)] <- FALSE
  	
  	cor.same = cbind(cor.filtered.1[idx.same],cor.filtered.2[idx.same])
  	rownames(cor.same) <- names(cor.filtered.1)[idx.same]
  	cor.switched = cbind(cor.filtered.1[idx.switched],cor.filtered.2[idx.switched])
  	rownames(cor.switched) <- names(cor.filtered.1)[idx.switched]  	
  	cor.diff = cbind(cor.filtered.1[idx.diff & (!idx.switched)],cor.filtered.2[idx.diff & (!idx.switched)])
  	rownames(cor.diff) <- names(cor.filtered.1)[idx.diff & (!idx.switched)]

	n.switchedDCL = nrow(cor.switched)
	if ( is.null( rownames(cor.same) ) ) {name.same = NULL}
	if ( !is.null( rownames(cor.same) ) ) {
		name.same = strsplit(rownames(cor.same),',')
		name.same = matrix(unlist(name.same),length(name.same),2,byrow=T)
	}
	if ( is.null( rownames(cor.switched) ) ) {name.switched = NULL}
	if ( !is.null( rownames(cor.switched) ) ) {
		name.switched = strsplit(rownames(cor.switched),',')
		name.switched = matrix(unlist(name.switched),length(name.switched),2,byrow=T)
	}
	if ( is.null( rownames(cor.diff) ) ) {name.diff = NULL}
	if ( !is.null( rownames(cor.diff) ) ) {
		name.diff = strsplit(rownames(cor.diff),',')
		name.diff = matrix(unlist(name.diff),length(name.diff),2,byrow=T)
	}
	name.all = rbind(name.same,name.switched,name.diff)

  #############################################################
  ## Determine DCLs from same sign correlation pairs
  #############################################################
  	if(nrow(cor.same)>1){
		de.s = LFC(cor.same,nbins,p,sign="same")
		DCL.same = cor.same[de.s,]
		name.same = name.same[de.s,]
		n.sameDCL = nrow(DCL.same)
		DCL.same <- data.frame(name.same,DCL.same);
		colnames(DCL.same) <- c("Gene.1","Gene.2","cor.1","cor.2")
 	} else stop("only one or no same-signed pair in all!")

  #############################################################
  ## Determine DCLs from different sign correlation pairs
  #############################################################	
  	if(nrow(cor.diff)>1){
		de.d = LFC(cor.diff,nbins,p,sign="diff")
		DCL.diff = cor.diff[de.d,]
		name.diff = name.diff[de.d,]
  	n.diffDCL = nrow(DCL.diff)
  	DCL.diff <- data.frame(name.diff,DCL.diff);
		colnames(DCL.diff) <- c("Gene.1","Gene.2","cor.1","cor.2")
	} else stop("only one or no differently-signed pair in all!")

	
################################################################################################
## Determine Switched DCLs if they exist 
################################################################################################
	
	pairs=rbind(name.same,name.diff,name.switched);
	if(n.switchedDCL>0) {
		DCL.switched <- data.frame(name.switched,cor.switched);
		colnames(DCL.switched) <- c("Gene.1","Gene.2","cor.1","cor.2")
		cor.max <- apply(abs(cor.switched),1,max)
		middle <-sort(cor.max,method = "quick", index.return=TRUE,decreasing=TRUE)$ix ########
		DCL.switched<- DCL.switched[middle,]
	}

####################################
## All links
####################################
	g.all <- graph.data.frame(name.all);
	gene.all <- as.matrix(V(g.all)$name);
	de.all <- degree(g.all);  
#####################################
## DCLs
#####################################
	g <- graph.data.frame(pairs);
	gene.1 <- as.matrix(V(g)$name);
	de <- degree(g);  
######################################
##DCLs of same sign
######################################
	g.same <- graph.data.frame(name.same);
	g.same.name <- as.matrix(V(g.same)$name);
	degree.same <- as.matrix(degree(g.same));

########################################
## DCLs of different sign
########################################

	g.diff <- graph.data.frame(name.diff);
	g.diff.name <- as.matrix(V(g.diff)$name);
	degree.diff <- as.matrix(degree(g.diff));

#######################################
## DCLs of switched correlation
#######################################
	if(n.switchedDCL>0) {
		g.switch <- graph.data.frame(name.switched);
		g.switch.name <- as.matrix(V(g.switch)$name);
		degree.switch <- as.matrix(degree(g.switch));
	} else { degree.switch = matrix(0,1,1)
	DCL.switched = matrix("NULL",1,1)
	}

#######################################
## Numbers for DCLs of different type. 
#######################################

	degree.bind <- matrix(0,m,5)
	row.names(degree.bind) <- genes
	colnames(degree.bind) <- c("All.links","DC.links","DCL.same","DCL.diff","DCL.switched")

	degree.bind[gene.all,1]=de.all
	degree.bind[gene.1,2]=de
	degree.bind[g.same.name,3]=degree.same
	degree.bind[g.diff.name,4]=degree.diff
	if(n.switchedDCL>0) {
		degree.bind[g.switch.name,5]=degree.switch
	}


########################################################
## DCGs Identification
########################################################
#
# 	prob <- nrow(pairs)/nrow(name.all)
#	p.value.DCe <- pbinom(degree.bind[,'DC.links']-1, degree.bind[,'All.links'], prob, lower.tail = F, log.p = FALSE);
# 	q.value.DCe <- p.adjust(p.value,method=q.method);
# 
# 	degree.bind <- cbind(degree.bind,p.value.DCe,q.value.DCe)
# 	colnames(degree.bind) <- c("All.links","DC.links","DCL_same","DCL_diff","DCL_switch","p","q")
#
# 	middle <-sort(as.numeric(degree.bind[,'q']),method = "quick", decreasing=FALSE,index.return=TRUE)$ix 
# 	DCGs <- degree.bind[middle,]
# 
##########################################################
 
	DCLs <- rbind(data.frame(DCL.same,type='same signed'),data.frame(DCL.diff,type='diff signed'))
	if (n.switchedDCL>0) DCLs <- rbind(DCLs,data.frame(DCL.switched,type='switched opposites'))
	DCLs <- data.frame(DCLs,cor.diff=FALSE)
	DCLs[,'cor.diff'] <- abs(DCLs[,'cor.1']-DCLs[,'cor.2'])
 	Result.DCe <- DCLs

	Result<-list(genes=Result.DCp,links=Result.DCe)
	return(Result)
}

"DCpathway" <- function(DCEA.res=DCEA.res, DisName="COPD",pathways){
	##pathways: pathID geneID
	gene.dC<-merge(pathways, DCEA.res$genes, by.x="geneID", by.y="row.names", all.x=T)
#	gene.dC<-gene.dC[!is.na(gene.dC$dC),]
	pathway.factor<-factor(gene.dC$pathID, levels=unique(gene.dC$pathID),ordered=T)
	pathway.dC<-as.numeric(unlist(by(gene.dC$dC, pathway.factor, mean, na.rm=T)))
#	pathway.dC<-as.numeric(unlist(by(gene.dC$dC, pathway.factor, function(x) mean(x,na.rm=T))))
	Result<-data.frame(pathID=levels(pathway.factor),dC=pathway.dC)
	colnames(Result)[2]<-paste(DisName,"dC",sep=".")
#	rownames(Result)<-Result[,"pathID"]
#	Result
#	Result<-as.matrix(Result)
#	Result
#	Result<-Result[,-1]
#	Result
#	Result<-as.matrix(Result)
#	Result
#	colnames(Result)<-paste(DisName,"dC",sep=".")
#	Result
	return(Result)
}

"DS"<- function(DCpathway.disn,Ndis=3,DCEA.disn,
								DisNames=c("AA","IgA","T2D"), pathways, cutoff=0.05, Nper=0, 
								FigName="DisNetwork.pdf",vsize=5,lcex=0.5,ewidth=1.5 ) {
	DCpathway.disn<-DCpathway.disn[,!duplicated(colnames(DCpathway.disn))]
	colnames(DCpathway.disn)[2:ncol(DCpathway.disn)]<-DisNames
	if (ncol(DCpathway.disn)-1!=Ndis | length(names(DCEA.disn))!=Ndis | length(DisNames)!=Ndis) {
		stop (' "DCpathway.disn", "DCEA.disn", "DisNames" and "Ndis" must have the same number.')
	}
#	DCpathway.disn<-DCpathway.disn[,!duplicated(colnames(DCpathway.disn))]
	if (Nper>0){
		parcor0<-PartialCor(DCpathway.disn=DCpathway.disn)
		parcor.N<-matrix(,nrow(parcor0),0)
		for (i in 1:Nper){
#			pathway.i<-pathway2gene.random(pathways)
			pathway.i<-data.frame(pathID=pathways[,"pathID"],geneID=sample(pathways[,"geneID"]))
			DCpathway.disn.j<-matrix(,nrow(DCpathway.disn),0)
			for (j in 1:Ndis){
				Result.j<-DCpathway(DCEA.res=DCEA.disn[[j]],DisName=DisNames[j],pathway.i)
				DCpathway.disn.j<-cbind(DCpathway.disn.j,Result.j)
			}
			DCpathway.disn.j<-DCpathway.disn.j[,!duplicated(colnames(DCpathway.disn.j))]
			parcor.i<-PartialCor(DCpathway.disn=DCpathway.disn.j)
			parcor.N<-cbind(parcor.N,parcor.i)
		}
		p.value=apply(abs(cbind(parcor0,parcor.N)),1,function(x) sum(x[1:(length(x)-1)]>x[length(x)],na.rm=T)/length( x[1:(length(x)-1)] ))
		q.value=p.adjust(p.value,method="BH")
		Result<-data.frame(parcor0,p.value=p.value,q.value=q.value)
	} else {
		parcor0<-PartialCor(DCpathway.disn=DCpathway.disn)
		Result<-data.frame(parcor0,p.value=rep(NA,length(nrow(parcor0))),q.vlaue=rep(NA,length(nrow(parcor0))))
	}
#	return(Result)
###########################################################
## vis the disease associations.
## red line for positive relationship, green line for negative relationship, gray line for un-significant relationship
###########################################################	
	a<-t(matrix(unlist(strsplit(rownames(Result),",")),2,nrow(Result)))
	relation<-cbind(a,Result)
	colnames(relation)[1:2]<-c("dis1","dis2")
	if (nrow(relation)>1000) {
		warning ("the  number of significant disease pairs(>1000) is too large to display clearly and the program maybe need long time to run.\n")
	}
	node<-unique(c(as.character(relation$dis1),as.character(relation$dis2)))
	f <- graph.data.frame(relation)
#	if (length(V(f)) != nrow (node.classes)) {
#		stop('oops - why rows of node classes not equal to nodes of graph?\n')
#		} else {
#		node.classes = node.classes[V(f)$name,]
		V(f)$shape <- 'circle'
		V(f)$color <- 'skyblue'
		V(f)$label.color <- "black";
		V(f)$label.cex <- lcex;
		V(f)$size <- vsize
		V(f)$label <- V(f)$name
		E(f)$lty <- 1; 
		E(f)$arrow.mode <- '-'
		E(f)$color <- "gray"
		E(f)$color[relation$p.value<cutoff & relation$parcor>0] <- "red"
		E(f)$color[relation$p.value<cutoff & relation$parcor<0] <- "green"
		E(f)$width <- ewidth
		pdf(FigName)
		plot(f,layout=layout.fruchterman.reingold)
		dev.off()
		cat("The graph of", FigName, "has been completed and saved in your working directory.\n")
#	}
#	RIF_reg_rank<-RIF_reg_rank[order(-abs(as.numeric(RIF_reg_rank[,'RIFscore']))),]
	Result<-Result[order(as.numeric(Result[,'p.value'])),]
	return(Result)
}

"comDCGL"<-function(Ndis=3,DCEA.disn,DisNames=c("AA","CKD","T2D"),cutoff=0.05, tf2target ) {
	if ( length(names(DCEA.disn))!=Ndis | length(DisNames)!=Ndis) {
		stop (' "DCEA.disn", "DisNames" and "Ndis" must have the same number.')
	}
	
	for (i in 1:Ndis){
		DCEA.disn[[i]]$genes<-data.frame(DCGs=rownames(DCEA.disn[[i]]$genes),DCEA.disn[[i]]$genes)
		colnames(DCEA.disn[[i]]$genes)[2:ncol(DCEA.disn[[i]]$genes)]<-paste(colnames(DCEA.disn[[i]]$genes)[2:ncol(DCEA.disn[[i]]$genes)],DisNames[i],sep=".")
		colnames(DCEA.disn[[i]]$links)[3:ncol(DCEA.disn[[i]]$links)]<-paste(colnames(DCEA.disn[[i]]$links)[3:ncol(DCEA.disn[[i]]$links)],DisNames[i],sep=".")
	}
	comDCGs<-data.frame(DCGs=DCEA.disn[[1]]$genes[,1])
	comDCLs<-DCEA.disn[[1]]$links[,c("Gene.1","Gene.2")]
	for (i in 1:Ndis){
		genes.i<-DCEA.disn[[i]]$genes
#		colnames(genes.i)<-paste(colnames(genes.i),DisNames[i],sep=".")
		links.i<-DCEA.disn[[i]]$links
#		colnames(genes.i)<-paste(colnames(genes.i),DisNames[i],sep=".")
		if (sum(is.na(genes.i$p.value))==nrow(genes.i)){
			dCname<-paste("dC",DisNames[i],sep=".")
			genes.i <- genes.i[order(-genes.i[,dCname]),]
			DCGs.i <- genes.i[1:ceiling(nrow(genes.i)*cutoff),]
		} else {
			DCGs.i<-genes.i[genes.i$p.value<cutoff,]
		}
		comDCLs<-merge(links.i,comDCLs,by.x=c("Gene.1","Gene.2"),by.y=c("Gene.1","Gene.2"))
		comDCGs<-merge(DCGs.i,comDCGs,by.x="DCGs",by.y="DCGs")
#		colnames(comDCGs)[1]<-"row.names"
	}
	if (nrow(comDCGs)<1) warning('there is no intersection DCGs between', Ndis, 'diseases.')
	if (nrow(comDCLs)<1) warning('there is no intersection DCLs between', Ndis, 'diseases.')
	
	
	DCLs<-comDCLs
	DCG<-comDCGs$DCGs
	tf.target <- paste(tf2target[,'TF'],tf2target[,'gene'],sep='; ')
	DCL.pair1 <- paste(DCLs$Gene.1,DCLs$Gene.2, sep='; ')
	TF1.idx <- DCL.pair1 %in% tf.target
	DCL.pair2 <- paste(DCLs$Gene.2,DCLs$Gene.1, sep='; ')
	TF2.idx <- DCL.pair2 %in% tf.target
	DCLs <- data.frame(TF=NA,DCLs)
	DCLs[TF1.idx,'TF'] <- as.character(DCLs[TF1.idx,'Gene.1'])
	DCLs[TF2.idx,'TF'] <- as.character(DCLs[TF2.idx,'Gene.2'])
	DCLs[TF1.idx & TF2.idx, 'TF'] <- paste(as.character(DCLs[TF1.idx & TF2.idx,'Gene.1']),as.character(DCLs[TF1.idx & TF2.idx,'Gene.2']),sep='; ')
	colnames(DCLs)[1] <- 'internal.TF'
	DCG1.idx<-DCLs$Gene.1 %in% DCG
	DCG2.idx<-DCLs$Gene.2 %in% DCG
	DCLs<-data.frame(DCG=NA,DCLs)
	DCLs[DCG1.idx,'DCG']<-as.character(DCLs[DCG1.idx,'Gene.1'])
	DCLs[DCG2.idx,'DCG']<-as.character(DCLs[DCG2.idx,'Gene.2'])
	DCLs[DCG1.idx & DCG2.idx,'DCG']<-paste(as.character(DCLs[DCG1.idx & DCG2.idx,'Gene.1']),as.character(DCLs[DCG1.idx & DCG2.idx,'Gene.2']),sep=';')
	DCLs<-DCLs[,c(3,4,2,1,5:ncol(DCLs))]
	TF<-unique(sort(tf2target$TF))
	comDCGs<-data.frame(DCGisTF=FALSE,comDCGs)
	comDCGs[,'DCGisTF']<-comDCGs[,'DCGs'] %in% TF
	comDCGs<-comDCGs[,c(2,1,3:ncol(comDCGs))]
	Result<-list(comDCGs=comDCGs,comDCLs=DCLs)
	return(Result)
}

"comDCGLplot"<-function(comDCGL.res,FigName="comDCGL.pdf",tf2target,
												vsize=5,asize=0.3,lcex=0.3,ewidth=1.5){
#Gene.1  Gene.2 internal.TF  DCG    cor.1.T2D  cor.2.T2D    type.T2D cor.diff.T2D   cor.1.IgA  cor.2.IgA    type.IgA cor.diff.IgA    cor.1.AA    cor.2.AA
#1 AKR1A1  IGFBP7        <NA> <NA> -0.937724906  0.0780705 diff signed    1.0157954 -0.83962129  0.1085241 diff signed    0.9481454 -0.05800985 -0.78363613
#2   CCNC SH3BGRL        <NA> <NA>  0.427680424 -0.8203428 diff signed    1.2480232  0.29612477 -0.7735561 diff signed    1.0696809  0.91564968  0.13898127
#3  JOSD1   AP2M1        <NA> <NA> -0.033814698 -0.7297648 same signed    0.6959501 -0.07225289 -0.6868257 same signed    0.6145728  0.91416922 -0.17348650

	DCGs<-comDCGL.res$comDCGs
	DCLs<-comDCGL.res$comDCLs
	nodes<-unique(sort(c(as.character(DCLs$Gene.1),as.character(DCLs$Gene.2))))
	DCG<-DCGs$DCGs
	TF<-unique(sort(tf2target$TF))
	nodes.classes<-data.frame(TF=nodes %in% TF, DCG= nodes %in% DCG)
	rownames(nodes.classes) <- as.character(nodes)
	
	relation<-DCLs[,1:4]
	if (nrow(relation)>1000) {
		warning ("the edge number of common DCL(>1000) is too large to display clearly and the program maybe need long time to run.\n")
	}
	g <- graph.data.frame(relation)
	if (length(V(g))!=nrow(nodes.classes)) {
		stop('oops - why rows of node classes not equal to nodes of graph?\n')
	} else {
		nodes.classes = nodes.classes[V(g)$name,]
		V(g)$shape <- 'circle'; V(g)$shape[nodes.classes$TF] <- 'square'
		V(g)$color <- 'skyblue'; V(g)$color[nodes.classes$DCG] <- 'pink'
		E(g)$color <- "black"
		E(g)$width <- ewidth
		E(g)$arrow.mode <- '-'; E(g)$arrow.mode[!is.na(relation$internal.TF)] <- '>'
		E(g)$arrow.size <- asize
		V(g)$size <- vsize
		V(g)$label.color <- "black";
		V(g)$label.cex <- lcex; 
		V(g)$label <- V(g)$name
		pdf(FigName)
		plot(g,layout=layout.fruchterman.reingold)
		dev.off()
		cat("The graph of", FigName, "has been completed and saved in your working directory.\n")
	}
}

"qLinkfilter" <-function(exprs.1,exprs.2,cutoff=0.25,r.method=c('pearson','spearman')[1],q.method=c("BH","holm", "hochberg", "hommel", "bonferroni", "BY","fdr")[1]) {
        # m <- nrow(exprs.1) # exprs.1, exprs.2 is the expression data for different conditions.
        degree.1 <- ncol(exprs.1)-2
        degree.2 <- ncol(exprs.2)-2

        genes <- rownames(exprs.1)
        exprs.1 <- as.matrix(exprs.1)
        exprs.2 <- as.matrix(exprs.2)
        cor.1 <- cor(t(exprs.1),method=r.method,use="pairwise.complete.obs")
        cor.2 <- cor(t(exprs.2),method=r.method,use="pairwise.complete.obs")
        cor.1 <- cor.1[lower.tri(cor.1,diag=F)]
        cor.2 <- cor.2[lower.tri(cor.2,diag=F)]

        rm(exprs.1); rm(exprs.2)

	
        t.1 <- cor.1*sqrt(degree.1)/sqrt(1-cor.1*cor.1)
        t.2 <- cor.2*sqrt(degree.2)/sqrt(1-cor.2*cor.2)

        p0.1 <- 2*pt(-abs(t.1), degree.1, lower.tail = TRUE, log.p = FALSE)
        p0.2 <- 2*pt(-abs(t.2), degree.2, lower.tail = TRUE, log.p = FALSE)
        #diag(p0.1) <- NA
        #diag(p0.2) <- NA
        rm(t.1); rm(t.2)

        q.1<- p.adjust(p0.1, method = q.method)
        q.2<- p.adjust(p0.2, method = q.method)
        #dim(q.1)<- dim(p0.1); 
        #dim(q.2)<- dim(p0.2)
        #diag(q.1) <- 1
        #diag(q.2) <- 1
        rm(p0.1); rm(p0.2)

        rth.1 <- abs(cor.1[which.min(abs(q.1-cutoff))])
        rth.2 <- abs(cor.2[which.min(abs(q.2-cutoff))])
        cor.1[q.1>=cutoff & q.2>=cutoff] <- cor.2[q.1>=cutoff & q.2>=cutoff] <- 0
        
        #cor.1 <- cor.2 <- diag(rep(0,length(genes)))
        #cor.1[lower.tri(cor.1,diag=F)] = cor.1; cor.1 = cor.1+t(cor.1)
        #cor.2[lower.tri(cor.2,diag=F)] = cor.2; cor.2 = cor.2+t(cor.2)
        #rownames(cor.1) <- rownames(cor.2) <- colnames(cor.1) <- colnames(cor.2) <- genes


	name.row <- matrix(rep(genes,length(genes)),length(genes),length(genes))
	name.col <- matrix(rep(genes,length(genes)),length(genes),length(genes),byrow=T)
	name.pairs <- matrix(paste(name.row,name.col,sep=','),length(genes),length(genes))
	rm(list=c('name.row','name.col'))
	name.pairs <- name.pairs[lower.tri(name.pairs,diag=F)]
	names(cor.1) <- names(cor.2) <- name.pairs


        cor.filtered <- list(rth.1 = rth.1, rth.2 = rth.2, cor.filtered.1 = cor.1, cor.filtered.2 = cor.2, genes=genes)
        return(cor.filtered)
}


"rLinkfilter" <- function(exprs.1,exprs.2,cutoff=0.8,r.method=c('pearson','spearman')[1]) {

 	genes <- rownames(exprs.1)
        exprs.1 <- as.matrix(exprs.1)
        exprs.2 <- as.matrix(exprs.2)
        cor.filtered.1 <- cor(t(exprs.1),method=r.method,use="pairwise.complete.obs")
        cor.filtered.2 <- cor(t(exprs.2),method=r.method,use="pairwise.complete.obs")
        cor.filtered.1 <- cor.filtered.1[lower.tri(cor.filtered.1,diag=F)]
        cor.filtered.2 <- cor.filtered.2[lower.tri(cor.filtered.2,diag=F)]
        rm(exprs.1); rm(exprs.2)

        cor.filtered.1[abs(cor.filtered.1)<=cutoff & abs(cor.filtered.2)<=cutoff] <- 0
        cor.filtered.2[abs(cor.filtered.1)<=cutoff & abs(cor.filtered.2)<=cutoff] <- 0

	name.row <- matrix(rep(genes,length(genes)),length(genes),length(genes))
	name.col <- matrix(rep(genes,length(genes)),length(genes),length(genes),byrow=T)
	name.pairs <- matrix(paste(name.row,name.col,sep=','),length(genes),length(genes))
	rm(list=c('name.row','name.col'))
	name.pairs <- name.pairs[lower.tri(name.pairs,diag=F)]
	names(cor.filtered.1) <- names(cor.filtered.2) <- name.pairs

        list(rth.1=cutoff, rth.2=cutoff,cor.filtered.1 = cor.filtered.1, cor.filtered.2 = cor.filtered.2, genes=genes)

}
  ###############################################################################################################
  ## Select part of the correlation pairs, the max(abs(cor.1),abs(cor.2))
  ## exprs.1 a data frame or matrix for condition A, with rows as variables (genes) and columns as samples.  
  ## exprs.2 a data frame or matrix for condition B, with rows as variables (genes) and columns as samples.  
  ## percent percent of links to be reserved. 
  # output: A list with two components of lists: one lists the rth (thresholds of correlation values) for both conditions, the other lists the two matrices of filtered pairwise correlation values.
  ###############################################################################################################
"percentLinkfilter" <-function(exprs.1,exprs.2,cutoff=0.25,r.method=c('pearson','spearman')[1]) {
  	# m <- nrow(exprs.1) # exprs.1, exprs.2 is the expression data for different conditions.
  	n.1 <- ncol(exprs.1)
  	n.2 <- ncol(exprs.2)
  
  	#degree.1 <- n.1-2
  	#degree.2 <- n.2-2
  
  	genes <- rownames(exprs.1)
  	exprs.1 <- as.matrix(exprs.1)
  	exprs.2 <- as.matrix(exprs.2)
  	cor.filtered.1 <- cor(t(exprs.1),method=r.method,use="pairwise.complete.obs")
  	cor.filtered.2 <- cor(t(exprs.2),method=r.method,use="pairwise.complete.obs")
  	cor.filtered.1 <- cor.filtered.1[lower.tri(cor.filtered.1,diag=F)]
	cor.filtered.2 <- cor.filtered.2[lower.tri(cor.filtered.2,diag=F)]
	rm(exprs.1); rm(exprs.2)
	
	#diag(cor.filtered.1) <- 0
  	#diag(cor.filtered.2) <- 0
	cor.1.2.max <- pmax(abs(cor.filtered.1),abs(cor.filtered.2))
	#cor.1.2.max <- cor.1.2.max[lower.tri(cor.1.2.max,diag=F)]
  	#cor.1.2.vector <- cbind(cor.pairwise.1[lower.tri(cor.pairwise.1, diag = FALSE)],cor.pairwise.2[lower.tri(cor.pairwise.2, diag = FALSE)])
  	#cor.1.2.max <- apply(abs(cor.1.2.vector),1,max)  #max of the two cor
  	#cat('passed lower tri extraction \n')
	cor.1.2.max <- sort(cor.1.2.max,decreasing = TRUE);
	#cat('passed max corr operation \n')
  	# cor.1.2.max.sort <- as.matrix(cor.1.2.max.sort)
  	Rth <- cor.1.2.max[as.integer(length(cor.1.2.max)*cutoff)];
  	rm(cor.1.2.max)
  	#cor.filtered.1 <- cor.pairwise.1
  	#cor.filtered.2 <- cor.pairwise.2
	#rm(cor.pairwise.1); rm(cor.pairwise.2)
  	cor.filtered.1[abs(cor.filtered.1)<=Rth & abs(cor.filtered.2)<=Rth] <- 0
  	cor.filtered.2[abs(cor.filtered.1)<=Rth & abs(cor.filtered.2)<=Rth] <- 0
	#cat('passed 0 substitution \n')
  	
	#cor.2 <- cor.1 <- diag(rep(0,length(genes)))
	#cor.1[lower.tri(cor.1,diag=F)] = cor.filtered.1; cor.1 = cor.1+t(cor.1)
	#cor.2[lower.tri(cor.2,diag=F)] = cor.filtered.2; cor.2 = cor.2+t(cor.2)
	#rownames(cor.1) <- rownames(cor.2) <- colnames(cor.1) <- colnames(cor.2) <- genes
	#rm(cor.filtered.1); rm(cor.filtered.2)
	#cat('reform to two matrices\n')

	name.row <- matrix(rep(genes,length(genes)),length(genes),length(genes))
	name.col <- matrix(rep(genes,length(genes)),length(genes),length(genes),byrow=T)
	name.pairs <- matrix(paste(name.row,name.col,sep=','),length(genes),length(genes))
	name.pairs <- name.pairs[lower.tri(name.pairs,diag=F)]
	rm(list=c('name.row','name.col'))
	names(cor.filtered.1) <- names(cor.filtered.2) <- name.pairs

	cor.filtered <- list(rth.1=Rth,rth.2=Rth,cor.filtered.1 = cor.filtered.1, cor.filtered.2 = cor.filtered.2, genes=genes)
	return(cor.filtered)
}
	calc.dC <- function(cor.filtered.1,cor.filtered.2,genes) {
		nzero.vec <- (cor.filtered.1 != 0)|(cor.filtered.2 != 0 )
		nzero.sm <- diag(rep(0,length(genes)))
		nzero.sm[lower.tri(nzero.sm,diag=F)] <- nzero.vec; nzero.sm = nzero.sm+t(nzero.sm)
		number_uniq <- apply(nzero.sm,1,sum)
 
 		squares = (cor.filtered.1-cor.filtered.2)^2
  #	number_uniq = apply(cor.filtered.1!=0 | cor.filtered.2!=0,1,sum)
#  	ss = apply(squares,1,sum)
		squares.sm <- diag(rep(0,length(genes)))
		squares.sm[lower.tri(squares.sm,diag=F)] <- squares; squares.sm = squares.sm+t(squares.sm)
		ss = apply(squares.sm,1,sum)  
		LNED.result = as.vector(matrix(NA,length(genes),1))
  		LNED.result[number_uniq!=0] = sqrt(ss[number_uniq!=0])/sqrt(number_uniq[number_uniq!=0])
  		names(LNED.result) <- genes
		list(dC=LNED.result,length=number_uniq)
	}
"LFC" <- function(exprs,nbins=20,p=0.1,sign) {
 if(sign=='same'){
	exprs.min <- apply(abs(exprs),1,min)
  	exprs.max <- apply(abs(exprs),1,max)
  	exprs.diff <- exprs.max/exprs.min
 }else {
	 		exprs.min <- apply(abs(exprs),1,min)
  			exprs.max <- apply(abs(exprs),1,max)
  			exprs.diff <- exprs.max/exprs.min
  			a <- exprs.max
	 		exprs.max <- exprs.diff
			exprs.diff <- a	
	 	}
  n.tol = length(exprs.max)
	num.bin = ceiling(n.tol/nbins)
	exprs.max.sorted = sort(exprs.max)
	steps = min(exprs.max)
	for (i in 1:nbins) {
		if (i == nbins) {
			steps = c(steps, exprs.max.sorted[n.tol]+1)
		} else {
			steps = c(steps, exprs.max.sorted[i*num.bin])
		}	
	}
	exprs.bin<-rep(1,n.tol);
	for (i in 1:nbins) {
		exprs.bin[exprs.max>=steps[i] & exprs.max<steps[i+1]]<-i
	}
	
	steps.x<-(steps[1:(length(steps)-1)]+steps[2:length(steps)])/2
	# For bin 1 to 2 to 3, each time get the decisive point s x and y coordinates - the gradually elongating vectors exprs.x.sub and exprs.y.sub.
	exprs.x.sub<-exprs.y.sub<-NULL
	for (i in 1:nbins) {
		if (sum(exprs.bin==i)>0) {
			exprs.diff.sub<-exprs.diff[exprs.bin==i]
			diff.rank<-sort(exprs.diff.sub,decreasing=T,index.return=T)$ix
			exprs.x.sub<-c(exprs.x.sub,exprs.max[exprs.bin==i][diff.rank[ceiling(length(exprs.diff.sub)*p)]])
			exprs.y.sub<-c(exprs.y.sub,exprs.diff[exprs.bin==i][diff.rank[ceiling(length(exprs.diff.sub)*p)]])
			
		}
	}
	# important changement: steps.x substitutes for exprs.x.sub.
	mm<-glm(y~x,data=data.frame(y=exprs.y.sub,x=1/steps.x))
	exprs.diff.threshold<-predict(mm,newdata=data.frame(x=1/exprs.max))
	delink<- (exprs.diff>exprs.diff.threshold)
	delink
  }
  
"PartialCor"<-function(DCpathway.disn=DCpathway.disn){
#	DCpathway.disn<-as.matrix(DCpathway.disn)
	rownames(DCpathway.disn)<-DCpathway.disn[,1]
	DCpathway.disn<-DCpathway.disn[,-1]
	DCpathway.disn<-na.omit(DCpathway.disn)
	parcor<-pcor(DCpathway.disn)
	parcor<-parcor$estimate
	gse<-colnames(DCpathway.disn)
	parcor<-parcor[lower.tri(parcor,diag=F)]
	name.row <- matrix(rep(gse,length(gse)),length(gse),length(gse))
	name.col <- matrix(rep(gse,length(gse)),length(gse),length(gse),byrow=T)
	name.pairs <- matrix(paste(name.row,name.col,sep=','),length(gse),length(gse))
	name.pairs <- name.pairs[lower.tri(name.pairs,diag=F)]
	names(parcor)<-name.pairs
	parcor<-as.data.frame(parcor)
	return(parcor)
}	
 
 
"pathway2gene.random"<-function(pathways) {
	pathways <- data.frame(pathways)
	colnames(pathways) <- c('pathID','geneID')
	genes <- unique(as.character(pathways$geneID))
	nPath <- sort(table(pathways$pathID))
	factor.nPath <- factor(nPath,levels=unique(nPath),ordered=T)
	Paths <- by(names(nPath),factor.nPath,paste,sep='')
	Genes <- sapply(levels(factor.nPath),function(x) {sample(genes,x)})
	comb <- mapply(expand.grid,Paths,genes,SIMPLIFY=FALSE)
	Pathways.r <- matrix(unlist(lapply(comb,function(x) unlist(t(x)))),ncol=2,byrow=T)
	Pathways.r		
}
