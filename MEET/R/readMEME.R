readMEME<-function(resultat, num_motif){
  j<-0
  mi<-0
  iii<-0

  Sequencies<-vector(mode="list", length=num_motif)
  MEME.aligned<-vector(mode="list", length=(num_motif+1))
  Sites<-vector(mode="logical", length=num_motif)
  width<-vector(mode="logical", length=num_motif)
  consensus<-vector(mode="logical", length=num_motif)

  for(i in 1:length(resultat)){
    if(substr(resultat[i], start=1, stop=5)=="MOTIF"){
    #print(i)
    mi<-mi+1
    width[mi]<-as.numeric(substr(resultat[i], start=20, stop=21))
    Sites[mi]<-as.numeric(substr(resultat[i], start=34, stop=35))
    #print("Sites=")
    #print(substr(resultat[i], start=34, stop=35))
   
    #print("mi=")
    #print(mi)
    }
    else if(substr(resultat[i], start=1, stop=13)=="Sequence name" && substr(resultat[i], start=27, stop=31)== "Start"){
	j<-j+1
	#print("j=")
	#print(j)
	#print("Sites=")
	#print(Sites[j])
	m<-i+2
	Sequencies[[j]]<-resultat[m:(m+Sites[j]-1)]
	print(Sequencies[[j]])
	
      } else if(substr(resultat[i], start=1, stop=10)=="Multilevel"){
	  iii<-iii+1
	    #print("iii=")
	    #print(iii)
	    consensus[iii]=substr(resultat[i], start=22, stop=(22+width[iii]))
	    #print("consensus=")
	    #print(consensus[iii])
	    
      }
    
  }
  
  for(i in 1:length(Sequencies)){
      #print("alineacio")
      #print(i)
      alineacio<-Sequencies[[i]]
      if(is.null(alineacio)){
      print("null")
      }else{
      alineacio<-strsplit(alineacio, NULL)
      alineacio<-sapply(alineacio,function(X){X<-X[42:length(X)]} )
      alineacio<-t(alineacio)
      suma<-apply(alineacio,2,function(X){sum(X==" ")})
      alineacio<-alineacio[,suma<dim(alineacio)[1]]
      suma<-apply(alineacio,2,function(X){sum(X==".")})
      alineacio<-alineacio[,suma<dim(alineacio)[1]]
      alineacio[alineacio==" "]<-"-"
      alineacio[alineacio=="."]<-"-"
      }
      MEME.aligned[[i]]<-alineacio

      
  }
  
  print(consensus)
  MEME.aligned[[(num_motif+1)]]<-consensus
  
return(MEME.aligned)
}
