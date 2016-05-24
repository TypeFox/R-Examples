FindElement<-function(What,Object,Element=list()){
	#str(Object)
	if(class(Object)=="data.frame"){
		#search in columns
		if(What %in% colnames(Object)){			
			Element[[length(Element)+1]]<-Object[,What]
			names(Element)[length(Element)]=paste(What,"_",length(Element),sep="")
		}
		else if(What %in% rownames(Object)){
			Element[[length(Element)+1]]<-Object[What,]
			names(Element)[length(Element)]=paste(What,"_",length(Element),sep="")
		}	
	}
	if(class(Object)=="list"){
		#Element=list()
		
		for(i in 0:length(Object)){
			if(i==0){
				Names=names(Object)
				if(What%in%Names){
					for(j in which(What==Names)){
						Element[length(Element)+1]=Object[j]
						names(Element)[length(Element)]=paste(What,"_",length(Element),sep="")
						return(Element)
					}
				}
			}
			else if(class(Object[[i]])[1]=="list"){
				#Names=names(Object[[i]])
				#if(What%in%Names){
				#	for(j in which(What==Names)){
				#			Element[length(Element)+1]=Object[[i]][j]
				#		names(Element)[length(Element)]=paste(What,"_",length(Element),sep="")
				#
				#	}
				#}
				Element=FindElement(What,Object[[i]],Element=Element)
				#for(j in 1:length(temp)){
				#	Element[length(Element)+1]=temp[j]
				#	names(Element)[length(Element)]=paste(What,"_",length(Element),sep="")
				#}
				
			}
			else if(class(Object[[i]])[1]=="data.frame"){
				Element=FindElement(What,Object[[i]],Element=Element)
				
			}
		}
	}	
	return(Element)
}
