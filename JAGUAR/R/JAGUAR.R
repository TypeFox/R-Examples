################################################################################################################
# 
# jaguar_slice - R function to slice pre-processed gene expression data
#
# Author: Chaitanya Acharya
# Updated on: Sep 15, 2015
#
# Creates subdirectories to store partitioned gene expression data
#
# Arguments: 
#			geneexp		->	Matrix of preprocessed gene expression data
#							Individuals with missing tissue gene expression data will have NAs in their place
#
#			size		->	Indicates the size of the partition. Default value is 100 genes i.e. each
#							gene expression slice will consist of 100 genes
#			
#			path		->	Indicates the "full" path to the subdirectories. Default value is the current working directory
#
#
# Returns: Partitioned gene expression data
################################################################################################################

"jaguar_slice" <- function(geneexp,size=100,path=getwd()){
	
	ngenes = size;
	samples = colnames(geneexp);
	tot.size = dim(geneexp)[1];
	cat("\nTotal number of genes: ",tot.size," Partition size: ",size,"\n\n")
	part = round(tot.size/ngenes);
	if(part == 1) stop("The size of the partition must be smaller than the data size")
	cat( "Slicing the gene expression data in to ", part ," partitions", "\n\n" );
	output = matrix(1:((part-1)*ngenes),nrow=part-1,byrow=TRUE);
	dataSplit = split(output, 1:nrow(output));
	dataSplit[[length(dataSplit)+1]]= (length(dataSplit)*ngenes+1):tot.size;
	for(i in 1:length(dataSplit)){
		dir.create(paste(path,"/dir",i,sep=""));
		setwd(paste(path,"/dir",i,sep=""));
		print(getwd())
		write.table(geneexp[dataSplit[[i]],],"GeneExp_matrix.txt",sep="\t",col.names=NA,quote=F);
		setwd("../");
	}
	cat( "\nPlease check for the subdirectories under the listed path: ",path, "\n" );
}

################################################################################################################
# 
# jaguar_gwa - JAGUAR R function for genome-wide analysis
#
# Author: Chaitanya Acharya
# Updated on: Sep 15, 2015
#
# 
# Returns a matrix of p-values indicating the significance of association between a gene-SNP pair
# Our joint score test statistic is computed as --
#
# U_\psi = Y^t . V^{-1} . (0.5 a_\gamma XX^t + a_\beta GG^t) . V^{-1} . Y
# 
# Arguments: 
#			geneexp		->	Matrix of preprocessed gene expression data
#							Individuals with missing tissue gene expression data will have NAs in their place
#
#			geno		->	Matrix of genotype data in allele dosage format (0/1/2)
#
#
# Returns: Matrix of p-values (computed from the Satterthwaite method) with genes on rows and SNPs on columns
################################################################################################################

"jaguar_gwa" <- function(geneexp,genomat,write=FALSE){
	
	cat("\nRunning Genome-Wide eQTL analysis on ",nrow(geneexp)," genes and ",nrow(genomat)," SNPs","\n\n")
	geneID = rownames(geneexp);
	snpID = rownames(genomat);
	nobs = ncol(genomat);
	k = ncol(geneexp) / ncol(genomat)
	ind = as.factor(rep(1:nobs,k))
	tissue = as.factor(rep(1:k,each=nobs))
	
	temp = llply(suppressWarnings( data.frame(t(geneexp),check.names=F) ),as.numeric);
	k_new = k - (apply(matrix(temp[[1]],nobs,k),1,function(i) sum(is.na(i))))
	R = 1 - t(apply(is.na(matrix(temp[[1]],nobs,k)),1,as.numeric));
	
	GeneOBJ = llply(temp,function(x){
					fit = lmer(x ~ 0+tissue+(1|ind),REML=F)
					x = matrix(x,nobs,k)
					est.eps = attr(VarCorr(fit),"sc")^2;
					est.tau = VarCorr(fit)[[1]][1];
					Ynew = apply(x,1,function(xx) xx-colMeans(x,na.rm=T));
					Ynew[is.na(Ynew)]<-0
					return(list("Eps"=est.eps,"Tau"=est.tau,"Y"=Ynew));
					})
	Eps = laply(GeneOBJ,function(x){x$Eps})
	Tau = laply(GeneOBJ,function(x){x$Tau})
	Y = llply(GeneOBJ,function(x){x$Y})
	U = GENEapply(as.matrix(genomat),Y,Eps,Tau,k_new,R)
	rownames(U)=geneID;
	colnames(U)=snpID;	
	if(write) write.table(U,"jaguar_out.txt",sep="\t",col.names=NA,quote=F)
	return(U)
}

################################################################################################################
# 
# jaguar_cis - JAGUAR R function for cis analysis
#
# Author: Chaitanya Acharya
# Updated on: Sep 15, 2015
#
# 
# If permutations are used to determine the null distribution of the score test statistic for a single gene,
# 
#
# Our joint score test statistic is computed as --
#
# U_\psi = Y^t . V^{-1} . (0.5 a_\gamma XX^t + a_\beta GG^t) . V^{-1} . Y
# 
# Arguments: 
#			geneexp		->	Matrix of preprocessed gene expression data
#							Individuals with missing tissue gene expression data will have NAs in their place
#
#			geno		->	Matrix of genotype data in allele dosage format (0/1/2)
#		
#			snp.bed		->  BED file format of SNP description
#
#			gene.bed	->  BED file format of gene description
#
#			cisDist		->  cis distance is defined as the maximum absolute distance between the gene and a SNP. Default value is 100Kb
#
#			nperm		->  Number of permutations. Default value is 10,000. Note that if it is 0, no permutations are performed
#
#			seed		->  Seed value for permutations. Default value is 123.
#
# Returns:	In case of permutations, returns a vector of p-values for every gene
#			In case of no permutations, returns a list of p-values for every gene
################################################################################################################

"jaguar_cis" <- function(geneexp,genomat,snp.bed,gene.bed,cisDist=100000,nperm=10000,seed=100){
	
	if(nperm<0) nperm = 0
	geneID = rownames(geneexp);
	snpID = rownames(genomat);
	nobs = ncol(genomat);
	k = ncol(geneexp) / ncol(genomat);
	ind = as.factor(rep(1:nobs,k));
	tissue = as.factor(rep(1:k,each=nobs));
	
	temp = llply(suppressWarnings( data.frame(t(geneexp),check.names=F) ),as.numeric);
	k_new = k - (apply(matrix(temp[[1]],nobs,k),1,function(i) sum(is.na(i))));
	R = 1 - t(apply(is.na(matrix(temp[[1]],nobs,k)),1,as.numeric));
	
	GeneOBJ = llply(temp,function(x){
					fit = lmer(x ~ 0+tissue+(1|ind),REML=F);
					x = matrix(x,nobs,k);
					est.eps = attr(VarCorr(fit),"sc")^2;
					est.tau = VarCorr(fit)[[1]][1];
					Ynew = apply(x,1,function(xx) xx-colMeans(x,na.rm=T));
					Ynew[is.na(Ynew)]<-0;
					return(list("Eps"=est.eps,"Tau"=est.tau,"Y"=Ynew));
					})
	Eps = laply(GeneOBJ,function(x){x$Eps});
	Tau = laply(GeneOBJ,function(x){x$Tau});
	Y = llply(GeneOBJ,function(x){x$Y});
	snpC.split = split(snp.bed,snp.bed[,1]);
	geneC.split = split(gene.bed,gene.bed[,1]);
	newGENO = apply(gene.bed,1,function(x){
					gene_chr = x[1]; gene_tss = as.numeric(x[2]);
					snp_chr = match(as.character(gene_chr),names(snpC.split));
					snp_sites = snpC.split[[snp_chr]][,2];
					maxDist = abs(gene_tss-snp_sites);
					keep = which(maxDist<=cisDist);
					cis_snps = snpC.split[[snp_chr]][keep,4];
					geno_mat = as.matrix(genomat[rownames(genomat) %in% cis_snps,]);
					})
	geno_dim = ldply(newGENO,dim)[,1];
	out = which(geno_dim==0);
	if(length(out)>0){
		Y=Y[-out]; Eps=Eps[-out]; Tau=Tau[-out]; newGENO=newGENO[-out]; geno_dim=geno_dim[-out];
	}
	if(nperm>0){
	  cat("\nRunning permutation-based cis-eQTL analysis on",length(Y), "genes with",nperm,"permutations", "\n\n")
		nperm=nperm+1;
		permMAT = matrix(0,nperm,nobs);
		permMAT[1,] = 0:(nobs-1);
		set.seed(seed);
		for(i in 2:nrow(permMAT)){
			permMAT[i,] = sample(0:(nobs-1),replace=F);
		}
		U_perm = cis_eqtl(newGENO,Y,permMAT,Eps,Tau,k_new,R,TRUE);
		snps = lapply(newGENO,rownames); names(snps) = names(Y)
		U_out = list( "results"=data.frame("Genes"=names(Y),"PermPval"=U_perm,"NbCisSNPS"=geno_dim,"NbPerm"=nperm-1),"cisSNPs"=snps )
					 
	}else{
	  cat("\nRunning cis-eQTL analysis on",length(Y),"genes with no permutation-resampling", "\n\n")
		permMAT = matrix(0,1,1);
		U_perm = cis_eqtl(newGENO,Y,permMAT,Eps,Tau,k_new,R,FALSE);
		names(U_perm) = names(Y);
		U_out = melt(U_perm)[,-1]; U_out[,1] = unlist(lapply(newGENO,rownames)); U_out = U_out[,c(1,3,2)]
		colnames(U_out)=c("SNP","Gene","pvalue")
	}
	return(U_out);
}

################################################################################################################
# 
# jaguar_process - R function to process results after running JAGUAR based on a predetermined threshold
#
# Author: Chaitanya Acharya
# Updated on: Sep 15, 2015
#
# NOTE: Output the files after running JAGUAR for this to work. 
#
# Arguments: 
#			jaguar.out		->	Matrix containing p-values outputted by gwa() function.
#
#			theshold		->	Predetermined threshold value to call for significant gene-SNP pairs
#
#			plot			->	Boolean value for a QQplot. Default value if FALSE
#								NOTE that the QQplot is generated from a randomly sampled set of p-values (50,000 values)
#								due to memory issues
#
# Returns: Matrix with significant gene-SNP pairs and their associated p-values along with an optional QQplot
################################################################################################################

"jaguar_process" <- function(jaguar.out,threshold,plot=FALSE){
  
  cat("\nProcessing JAGUAR results\n\n")
  if(class(jaguar.out)!="matrix") stop("Wrong data type passed!");
  if(plot){
    U.vec = as.vector(jaguar.out);
    if(length(U.vec)>50000){
      U.vec.temp = U.vec[sample(1:length(U.vec),50000,replace=F)];
    }else{
      U.vec.temp = U.vec[sample(1:length(U.vec),length(U.vec),replace=F)];
    }
    theoretical.quantiles = (1:length(U.vec.temp))/(1+length(U.vec.temp)); 
    data.quantiles = sort(U.vec.temp);
    plot(-log10(theoretical.quantiles),-log10(data.quantiles),main="QQ-plot of the score test pvals",xlab="Theoritical Quantiles",ylab="Data Quantiles");
    abline(0,1,col="red")
  }
  w = which(jaguar.out<threshold,arr.ind=T);
  val = jaguar.out[jaguar.out<threshold];
  genes = rownames(jaguar.out)[w[,1]];
  snps = colnames(jaguar.out)[w[,2]];
  eQTL = data.frame("Gene"=genes,"SNP"=snps,"pvalue"=val);
  return(eQTL);
}

################################################################################################################
# 
# jaguar_sim - R function to run simulations to test the association between one gene-SNP pair in
#					multiples tissues
#
# Author: Chaitanya Acharya
# Updated on: Sep 15, 2015
#
# We test the "Global" null "H_0: \beta=0; \gamma=0" using our joint score test, and
# test the "Local" null "H_0: \gamma=0" using a variance component score test.
#
# Date is generated under the following model --
# 
# y_{i,t} = \mu_{t} + \beta_{i}*g_{i} + \tau_{i} + b_{t}*g_{i} + \epsilon_{i,t}
# where i = individual and t = tissue
#
# Arguments: 
#			nobs		->	Number of observations. Default value is 500.
#			k			->  Number of tissues. Default value is 5.
#			tau			->  Subject-specific random intercept. Default value is 1.
#			eps			->	Random error. Default value is 1.
#			PVEg		->  Proportion of variance explained by \gamma (tissue-specific random effect). Default value is 0.
#			bta			->  Main genotypic effect indicating the overall shift in mean. Default value is 0.
#			maf			->  Minor allele frequency. Default value is 0.10
#
# NOTE: Running the simulations using the default values will give the type-I error of our model
#
# Returns: P-values from our joint score test and the variance component score test 
################################################################################################################

"jaguar_sim" = function(nobs = 500, k = 5, tau = 1, eps = 1,PVEg = 0,bta = 0,maf = 0.10){
	
	if(maf<0.05){
		warning("Low minor allele frequencies (MAF) lead to loss of power. MAF adjusted to 0.10");
		maf = 0.10;
	}
	gamma = ((eps+tau)*PVEg) / (100-PVEg);
	mu = rep(0.5,k);
	repeat{
		snp = rbinom(nobs,2,maf)
		if(sum(snp)>1){ break }
	}
	bkg=rnorm(k,0,gamma)
	Y = do.call("rbind",lapply(1:nobs,function(i) mu + bta*snp[i] + rnorm(1,0,tau) + bkg*snp[i] + rnorm(k,0,eps)))
	data = data.frame("IND"=as.factor(rep(1:nobs,each=k)),"Gene"=as.vector(t(Y)),"Geno" = rep(snp,each=k),"Tissue"=as.factor(rep(1:k,nobs)))
  
	# Fit the model at global null
	fit0 = lmer(Gene~0+Tissue+(1|IND),data,REML=F)
	est0.eps = attr(VarCorr(fit0),"sc")^2; est0.tau = VarCorr(fit0)[[1]][1];
	Y0new = apply(Y,1,function(x) x-colMeans(Y))
	joint_ST = jagSIM(est0.eps,est0.tau,k,Y0new,snp)
  
	# Fit the model at local null for U_gamma
	fit = lmer(Gene~0+Geno+Tissue+(1|IND),data,REML=F)
	est.mu = as.numeric(fixef(fit)); est.eps = attr(VarCorr(fit),"sc")^2;
	est.tau = VarCorr(fit)[[1]][1];	J = model.matrix(fit);
	Ynew = matrix(data$Gene-(J%*%est.mu),k,nobs)
	gamma_ST = vcSIM(est.eps,est.tau,k,Ynew,snp)
  
	return(c("VCScoreTest"=gamma_ST,"JointScoreTest"=joint_ST))
}

################################################################################################################
# 
# jaguar_plotqtl - R function to plot eQTL
#
# Author: Chaitanya Acharya
# Updated on: Sep 15, 2015
#
# NOTE that this is a slightly modified plotting function originally written by Wei Sun as a part of eMap R-Package
# Source: http://www.bios.unc.edu/~weisun/software/eMap.pdf
#
# Arguments: 
#			geneID		-> A vector indicating the row IDs of genes to be mapped
#			snpID		-> A vector indicating the row IDs of SNPs to be mapped
#			gene.chr	-> A vector indicating the chromosomal location of the genes to be mapped
#			gene.pos	-> A vector indicating the start site of all the genes on the Gene Chip
#			snp.chr		-> A vector indicating the chromosomal location of the SNPs to be mapped
#			snp.pos		-> A vector indicating the chromosomal location of all the SNPs on the SNP Chip
#			scores		-> A vector of p-values of each gene-SNP pair
#			chroms		-> A vector indicating the number of chromosomes to me mapped. Default value is 1 to 22
#
################################################################################################################


"jaguar_plotqtl" = function(geneID,snpID,gene.chr,gene.pos,snp.chr,snp.pos,scores,chroms=1:22){
  
	if(length(geneID) == 0){
		stop("length(geneID)=0\n")
	}
	if(length(geneID) != length(snpID)){
		stop("length(geneID) != length(snpID)\n")
	}
	if(length(geneID) != length(scores)){
		stop("length(geneID) != length(scores)\n")
	}
	valid.chroms = c(1:90, "X", "Y")
	wrongChr = chroms[!(chroms %in% valid.chroms)]
	if(length(wrongChr)>0){
		stop(wrongChr, " are not valid chromosome labels\n")
	}
  
	todrop = which(!(gene.chr %in% chroms))
	gene.chr[todrop] = NA
	gene.pos[todrop] = NA
	gene.chr[gene.chr=="X"] = 99
	gene.chr[gene.chr=="Y"] = 100
	gene.chr = as.integer(gene.chr)
	gene.pos = as.numeric(gene.pos)
  
	todrop = which(!(snp.chr %in% chroms))
	snp.chr[todrop] = NA
	snp.pos[todrop] = NA
	snp.chr[snp.chr=="X"] = 99
	snp.chr[snp.chr=="Y"] = 100
	snp.chr = as.integer(snp.chr)
	snp.pos = as.numeric(snp.pos)
  
	chrs = union(unique(gene.chr), unique(snp.chr))
	chrs = sort(as.numeric(chrs))
	num.chrs = chrs[chrs <= 90]
	chr.chrs = c(num.chrs, "X", "Y")
	max.num = max(num.chrs)
	if(max.num > length(num.chrs)){
		str = "there is no SNP/gene location information for some chromosomes\n"
		stop(str)
	}
	gene.chr[gene.chr==99]  = max.num+1
	gene.chr[gene.chr==100] = max.num+2
	snp.chr[snp.chr==99]  = max.num+1
	snp.chr[snp.chr==100] = max.num+2
	chreMax = tapply(gene.pos, gene.chr, max, na.rm=TRUE, simplify = FALSE)
	chrmMax = tapply(gene.pos, gene.chr, max, na.rm=TRUE, simplify = FALSE)
	chrMax = numeric(length(chrs))
	for(i in 1:length(chrs)){
		ch = as.character(i)
		chrMax[i] = max(chreMax[[ch]], chrmMax[[ch]])
	}
	temp = data.frame(Gene_ID=geneID, Marker_ID=snpID)
	nChr = length(chrMax)
	chrLen = c(0, cumsum(chrMax))
	ep = gene.pos + chrLen[gene.chr]
	mp = snp.pos + chrLen[snp.chr]
	ymax = chrLen[nChr+1]
	bdr1 = -0.016*ymax
	bdr2 = -0.006*ymax
  
	par(mar=c(3,4,0,0))
	plot(c(bdr1,ymax*1.05), c(bdr1,ymax*1.05), type="n", xlab="",ylab="", main="", xaxt="n", yaxt="n", bty="n")
	mtext("SNP Location", side=1, line=1)
	mtext("Transcript Location", side=2, line=1)
	gpos = ep[temp$Gene_ID]
	mpos = mp[temp$Marker_ID]
	points(mpos, gpos, col="black",pch=20, cex=1)
  
	if(nChr>1){
		nchr.plot = floor(nChr/2)
		rect(rep(bdr1,nchr.plot), chrLen[seq(1,nChr,by=2)],rep(bdr2,nchr.plot), chrLen[seq(2,nChr+1,by=2)],border=NA, col="red")
		rect(chrLen[seq(1,nChr,by=2)], rep(bdr1,nchr.plot),chrLen[seq(2,nChr+1,by=2)], rep(bdr2,nchr.plot),border=NA, col="red")
	}
	kp = seq(1,nChr,by=2)
	ats = 0.5*(chrLen[-1] + chrLen[-length(chrLen)])
	mtext(chr.chrs[kp], at=ats[kp], side=1, line=-0.5, cex=0.8)
	mtext(chr.chrs[kp], at=ats[kp], side=2, line=-0.5, cex=0.8)
	lines(c(0,ymax*1.01), rep(ymax*1.01,2), lty=2)
	lines(rep(ymax*1.01,2), c(0,ymax*1.01), lty=2)
}
