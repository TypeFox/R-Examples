uSentences <-
function(base, num.text, SentLength){		
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
	p1 <- unlist(strsplit(as.character(base[1:r,num.text]),split=" ")) 
		sel <- which(p1=="")      
		if (length(sel)==0){
         		 p1<-tolower(p1) 
		}
		 if (length(sel)!=0){	 		
    			p1 <- p1[-sel]           
		p1 <- tolower(p1)       
		}
       Tfrase<-SentLength		
	listagrupos<-list() 		
	ngrupos<-length(p1)%/%Tfrase 	
        for(i in 1:ngrupos){ 
		if (i!=ngrupos){ 
			listagrupos[[i]]<-p1[((i-1)*Tfrase+1):(i*Tfrase)]
		      	}else{
			listagrupos[[i]]<-p1[((i-1)*Tfrase+1):length(p1)]
			listagrupos[[i]]
			
		}
	}

	listita<-listagrupos
	  for (i in 1:length(listita)){
 		 names(listita[[i]])<-c(1:length(listita[[i]]))
	  }
		numeracio<-paste(c(1:length(p1[((i-1)*Tfrase+1):length(p1)])))
		listita2 <- lapply(listita, FUN=function(x) x[numeracio])
		Frase<- do.call("cbind", listita2)
		row.names(Frase) <-numeracio
		Fr<-t(Frase)
		Frse<-vector("list",nrow(Fr))
		for(i in 1:nrow(Fr)){
		for(j in 1:length(numeracio)){
			 if (!is.na(Fr[i,j]))Frse[[i]]<-paste(Frse[[i]],Fr[i,j])
			}
		}
		
	for (i in 1:length(Frse)){
		names(Frse[[i]])<-c(1:length(Frse[[i]]))
	}
		numeracio<-paste(c(1:1))
		listita<- lapply(Frse, FUN=function(x) x[numeracio])
		Frse<- do.call("cbind", listita)
		row.names(Frse) <-"Sentence"
		Frse<-as.data.frame(t(Frse))
return (Sentences=Frse)
}
