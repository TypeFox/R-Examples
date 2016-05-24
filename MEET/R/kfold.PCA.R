kfold.PCA <-
function(iicc, TF){

require("MEET")
require("ROCR")
require("pcaMethods")


  write.fasta <- get("write.fasta",pos="package:seqinr")
    read.fasta <- get("read.fasta",pos="package:seqinr")
    
    x<-read.fasta(file=TF,forceDNAtolower=FALSE)
    y<-x[-iicc$outsequence]
       
  
  position<-iicc$position
  Prob<-iicc$background
  Factortrans<-matriu<-factor
  missing<-iicc$missing
  vectorP<-iicc$vector

  DNA<-iicc$DNAreal
 
  NumericalMatrix<-numericalDNA(Prob)
 
  matriuROC<-vector(mode="list", length=length(vectorP))
  proba<-vector(mode="list", length=length(vectorP))
  
  for(i in 1:length(vectorP)){
  iicc$parameters<-vectorP[i]
 
   TotalResidus <- vector(mode="list", length=length(y))
 
    for(m in 1:length(y)){

      
	    z<-x[-m]
	    write.fasta(sequences=z,names=c(1:length(z)),file.out="background.fa",open="w") 
	    Factortrans<-switch(iicc$alignment, "CLUSTALW"=align.clustalw(filein="background.fa", fileout="Sq.fa", call=iicc$clustalw), "MUSCLE"=align.muscle(filein="background.fa", fileout="Sq.fa", gapopen=iicc$gapopen, maxiters=iicc$maxiters,gapextend=iicc$gapextend, call=iicc$muscle),"MEME"=align.MEME(filein="background.fa",fileout="Sq.fa",iicc),"NONE"=Read.aligned("background.fa"), stop("Alignment method not included"))
      Factortrans=as.matrix(Factortrans)
     
      window<-floor(ncol(Factortrans)) 

      suma<-apply(Factortrans,2,function(r){sum(r=="-")})
      threshold<-floor(nrow(Factortrans)*missing/100)
      Factortrans<-Factortrans[, suma<=threshold]
         
      training.set<-apply(Factortrans,1, function(x){as.vector(t(NumericalMatrix[x,]))})
      training.set<-t(training.set)
       Background<-sapply(X=c(1:(length(DNA)-ncol(Factortrans)+1)), FUN=function(X, Factortrans){DNA[X:(X+ncol(Factortrans)-1)]}, Factortrans=Factortrans)
      Background<-t(Background)
      Backgroundnum<-apply(Background,1,function(x){as.vector(t(NumericalMatrix[x,]))})
      Backgroundnum<-t(Backgroundnum)
 

      
       residus<- (PCanalysis(TFBS=training.set, nPCs=iicc$parameters, Sequences=Backgroundnum))
       #residus<-residus/ncol(training.set)
       parametersQ<-JacksonParameters(nPCs=iicc$parameters, TFBS=training.set)
      residus<-sapply(residus,QtoJackson, h0=parametersQ$h0, x1=parametersQ$x1,x2=parametersQ$x2,x3=parametersQ$x3 )
	TotalResidus[[m]]<-residus
    }
     
  matriuROC[[i]]<-unlist(TotalResidus)

 
  }

  matriuROC
}

