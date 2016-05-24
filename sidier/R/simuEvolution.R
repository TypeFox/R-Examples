simuEvolution <-
function(input,seqL,iLength,nReplicates){
for(reps in 1:nReplicates)
{
nodes<-max(input[,c(1,2)])
OutSeqs<-matrix(ncol=seqL,nrow=nodes)

nt<-c("A","C","T","G")
ancestor<-sample(nt,seqL,replace=T)
OutSeqs[1,]<-ancestor

	for(k in 1:nrow(input))
	{
	f0<-OutSeqs[input[k,1],]
	f1<-f0
	#substitutions
	S1<-round(input[k,3]*length(which(f0!="-")))
	mutations<-sample(which(f0!="-"),S1)
	for(l in mutations)
	f1[l]<-nt[sample(which(f0[l]!=nt),1)]
	OutSeqs[input[k,2],]<-f1
	#indels
	I1<-input[k,4]*length(which(f0!="-"))
	i1<-round(I1*input[k,5])
	d1<-round(I1*(1-input[k,5]))
	#insertions
	if(i1>0)
		{
		mutations2<-sample(which(f0!="-"),(i1))
		for(l in mutations2)
		#f1[l]<-nt[sample(which(f0[l]!=nt),1)]
			{
			OutSeqs2<-matrix("-",ncol=(iLength + ncol(OutSeqs)),nrow=nrow(OutSeqs))
			OutSeqs2[,c(1:l)]<-OutSeqs[,c(1:l)]
			OutSeqs2[,c((l+iLength):ncol(OutSeqs2))]<-OutSeqs[,c(l:ncol(OutSeqs))]
			OutSeqs2[k,c((l+1):((l+1)+(iLength-1)))]<-sample(nt,iLength,replace=T)
	
			}
		}
	#delections
	if(d1>0)
		{
		mutationSites<-length(f0)-iLength
		mutations2<-sample(which(f0[1:mutationSites]!="-"),(d1))
		for(l in mutations2)
		#f1[l]<-nt[sample(which(f0[l]!=nt),1)]
			{
		OutSeqs2<-OutSeqs
		OutSeqs2[k,c((l+1):((l+1)+(iLength-1)))]<-rep("-",iLength)
			}
		}
	OutSeqs<-OutSeqs2
	}
	

gapMAL<-c()
for (l3 in 1:ncol(OutSeqs))
if(paste(OutSeqs[,l3],collapse="")==paste(rep("-",nrow(OutSeqs)),collapse=""))
gapMAL<-c(gapMAL,l3)
if (is.null(gapMAL)==F) 
OutSeqsn<-OutSeqs[,-c(gapMAL)] #removes positions full of gaps


for(l2 in 1:nrow(OutSeqs))
write.table(file=paste("simulated_sequences_rep",reps,".fas",sep=""),paste(">",l2,"\n",paste(OutSeqs[l2,],collapse=""),sep=""),col.names=F, row.names=F,quote=F,append=T)
	
#kk<-read.dna(file="simulated_sequences_rep",reps,".fas",format="fasta")
#dist.dna(kk,model="N")
#plot(nj(dist.dna(kk,model="N")))
	
tips<-input[which(is.na(match(input[,2],input[,1]))),2]
for(l2 in tips)
write.table(file=paste("simulated_sequences_tips_rep",reps,".fas",sep=""),paste(">",l2,"\n",paste(OutSeqs[l2,],collapse=""),sep=""),col.names=F, row.names=F,quote=F,append=T)

#kk<-read.dna(file="simulated_sequences_tips_rep",reps,".fas",format="fasta")
#dist.dna(kk,model="N")
#plot(nj(dist.dna(kk,model="N")))
}
}
