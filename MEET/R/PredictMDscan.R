PredictMDscan <-function(iicc){
 require("seqinr")
 require("MEET")
           
 write.fasta <- get("write.fasta",pos="package:seqinr")
 read.fasta <- get("read.fasta",pos="package:seqinr")
 call.MDscan<-iicc$MDscan
 
 pvalor<-iicc$pvalor
 k<-length(iicc$DNA[[1]])
 len_motif<-iicc$lenmotif
 num_motif<-iicc$nummotif
 factor<-as.matrix(iicc$Transcriptionfactor)
 listfactor<-lapply(c(1:nrow(factor)),function(x){factor[x,]})
 output<-lapply(seq(1,length(iicc$DNA), 1), function(x){})
 
 for(i in c(1:length(iicc$DNA))){   	
    Sequence <- iicc$DNA[[i]]
    write.fasta(listfactor, names=c(1:(nrow(factor))), nbchar = ncol(factor), file.out="background.fa",open="w")
    resultats_MDscan<-run.read.MDscan(Sequence,k, len_motif, num_motif, call.MDscan)
    Score<-scoreMDscan(resultats_MDscan,k,factor,direction=NULL)
    DetectedFactors<-cbind(which(Score>=pvalor),Score[Score>=pvalor])
    if(nrow(DetectedFactors)==0){
      output<-"No Binding Sites Found"
    }else{
      output[[i]]<-lapply(c(1:nrow(DetectedFactors)),function(x){cbind(Sequence=paste(Sequence[DetectedFactors[x,1]:(DetectedFactors[x]+len_motif)],sep="",collapse=""),pvalue=DetectedFactors[x,2], position=DetectedFactors[x,1])})
      }
    }
    return(output) 

}

