VocIndex<-function(base, num.text, Fmin=5, sep.punctuation=TRUE, separ=NULL){

if (!is.null(num.text)) {
 if (is.character(num.text)) 
 num.text<- which(colnames(base) %in% num.text)
  if (is.numeric(num.text)) 
  num.text<- num.text
 if(length(num.text)==1)
  num.text<-num.text 
 if(length(num.text)>1){
   for(i in 1:length(num.text)){
    if(i==1)
    text1<-base[,num.text[1]]
    else text1<-paste(text1,base[,num.text[i]],sep=".") 
   }
    base[,(ncol(base)+1)]<-text1
    num.text<-ncol(base)
  }
}
	
    r<-length(rownames(base))
     base[,num.text]<-str_replace_all( base[,num.text],"[']", " ")
     p1 <- unlist(strsplit(as.character(base[1:r,num.text]),split=" "))  
    if(sep.punctuation){
                 filt="(['?\n<U+202F><U+2009>]|[[:punct:]]|[[:space:]]|[[:cntrl:]])+"
                 if (!is.null(separ)) filt=paste(separ,filt,sep="|") 
		      p2<-str_replace_all(p1,filt, "")
                      p3<-str_replace_all(p2,"[,;.:0123456789*]", "")
                    p1<-p3
	}   
	sel <- which(p1=="")      
	if (length(sel)==0){
          p1<-tolower(p1) 
	}
	 if (length(sel)!=0){	 		
    	p1 <- p1[-sel]           
	pOK <- tolower(p1)       
	}
		Lpal <- length(pOK)       
		uLpal <- unique(pOK)      

	indices <- vector()        
	freq<-vector() 
	Posit<-vector("list",length(uLpal)) 
	Distancias<-vector("list",length(uLpal))
      Exponencial<-vector()
   
	   for(i in 1:length(uLpal)){  
 		pos <- which(pOK==uLpal[i]) 
		F <- length(pos)          
 			freq[i]<-F
			T <- Lpal/F               
		posIni <- pos        
		posEnd<-c(pos[2:length(pos)],Lpal+pos[1])    
		D <- posEnd-posIni    
			k <- which(D>T)       
		Distancias[[i]]<-D

		if(length(k)==0){     
 				Nprim <- Lpal
 			}

 			if(length(k)!=0){   
 					Di <- D[k]     
 				 Nprim <- Lpal - sum(Di-T) 
 			}

 	R <- (Nprim-F)/(Lpal-F) 
 		indices[i] <- round(R,2) 
		Posit[[i]]<-pos       

	}
IndVocab<-cbind.data.frame(uLpal, indices, freq)
IndVocab<-subset(IndVocab, IndVocab$freq>=Fmin)
colnames(IndVocab)<-c("Word","Index","Freq")
IndComp<-cbind.data.frame(uLpal, indices, freq) 
IndVocab1<-subset(IndComp, IndComp$freq>=Fmin)
colnames(IndVocab1)<-c("Word","Index","Freq")
IndVocab<-with(IndVocab,IndVocab[order(Word),])
Voc.Regular<-with(IndVocab1, subset(IndVocab1, Index >=0.60))
Voc.Regular<-with(Voc.Regular,Voc.Regular[order(Index, decreasing=TRUE),])
Voc.Local<-with(IndVocab1, subset(IndVocab1, Index <= 0.40))
Voc.Local<-with(Voc.Local,Voc.Local[order(Index,Freq),])
 
res<-vector(mode='list')
		 res$VocIndex<-IndVocab
  		 res$RegVoc<-Voc.Regular
		 res$LocalVoc<-Voc.Local
         class(res)<-c("VocIndex","list")
return(res)

}
