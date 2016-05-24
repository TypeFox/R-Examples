## Functions for haplotype initialization for the 'g' option


# Use an EM approach from the package haplo.stats to estimate
# haplotype frequencies. Then (optionally) use the estimated
# frequencies to set up the initial haplotype configuration
# and a list of "known" haplotypes for the dictionary rephaser
estimateHap=function(args, HaploFreqFile, InitialHaplos=TRUE,
		InitialHaploFile="initialhaps", HaploList=TRUE, HaploListFile="initialhaplist",
		tol=0.00001){

	#require(haplo.stats)
	
	if (args$DataType=="h"){
		stop("Datatype is haplotype; haplotype frequencies do not need to be estimated\n")
	} else if (args$DataType!="g"){
		stop("Datatype must be 'g' to estimate haplotype frequencies\n")
	} else {
		
		args=changeArgs(args,HaploFreqFile=HaploFreqFile)
		
		genos=read.table(args$DataFile)
		
		# Convert genotype format
		genos2col=geno1to2(genos)
		
		# Estimate haplotype frequencies
		haps=haplo.em(genos2col)
		
		haplomat=as.matrix(haps$haplotype)
		haplomat=matrix(as.numeric(haplomat),byrow=F,ncol=ncol(haplomat))-1
		freqmat=cbind(haplomat,haps$hap.prob)
		
		
		# Set up the two locus haplotype frequency matrix
		numloci=ncol(haplomat)
		nsam=nrow(genos)
		hapfreqs=matrix(ncol=4, nrow=(numloci-1))
		
		for (i in 1:(numloci-1)){
			colcount=1
			for (j in 0:1){
				for (k in 0:1){
					rows=(freqmat[,i]==j)&(freqmat[,(i+1)]==k)
					hapfreqs[i,colcount]=(sum(haps$hap.prob[rows])*nsam)+1
					colcount=colcount+1
				}
				
			}
		}
		
		hapfreqs=round(hapfreqs,5)
		hapfreqs=hapfreqs/rowSums(hapfreqs)
		write.table(hapfreqs,file=HaploFreqFile,quote=F,row.names=F,col.names=F)
		
		
		# Sample initial configuration for each individual using the
		# results from haplo.em
		if (InitialHaplos==TRUE){
			
			args=changeArgs(args,InitialHaploFile=InitialHaploFile)
			
			res=cbind(haps$indx.subj, haps$hap1code,haps$hap2code,haps$post)
			outdat=NULL 
			nid=nrow(genos)
			
			for (j in 1:nid){
				sub=matrix(res[res[,1]==j,],ncol=ncol(res))
				config=sample(1:nrow(sub),size=1,prob=sub[,4])
				hap1=as.numeric(haps$haplotype[sub[config,2],])-1
				hap2=as.numeric(haps$haplotype[sub[config,3],])-1
				if ( sum( (hap1+hap2)!=genos[j,] )>0 ){
					stop("Estimated haplotypes don't match genotypes\n")
				}
				outdat=c(outdat,paste(hap1,collapse=""),paste(hap2,collapse=""))
			}
			write(outdat,file=InitialHaploFile,ncolumns=1)
		}
		
		
		# Write out a list of haplotypes for the dictionary rephaser
		if (HaploList==TRUE){		
			args=changeArgs(args,HaploListFile=HaploListFile)
			keep=haplomat[haps$hap.prob>tol,]
			haplolist=apply(keep,1,paste,collapse="")
			write(haplolist,file=HaploListFile,ncolumns=1,append=F)
		}
		
		
	}
	return(args)
	
}