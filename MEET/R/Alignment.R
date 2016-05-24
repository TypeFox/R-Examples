Alignment<-function(TF,iicc){
require("MEET")
alignedSequences<-switch(iicc$alignment, "CLUSTALW"=align.clustalw(filein=TF, fileout="SetTF.fa", call=iicc$clustalw), "MUSCLE"=align.muscle(filein=TF, fileout="SetTF.fa", gapopen=iicc$gapopen, maxiters=iicc$maxiters, gapextend=iicc$gapextend, call=iicc$muscle),"MEME"=align.MEME(filein=TF,fileout="SetTF.fa",iicc),"NONE"=Read.aligned(TF), stop("Alignment method not included"))
return(alignedSequences)
}