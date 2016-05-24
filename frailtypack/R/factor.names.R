#
#  names without 'factor' or 'as.factor'
#

factor.names <- function(x){
 	x <- matrix(x,nrow=1)
	Names <- apply(x,MARGIN=2,FUN=function(x){
 		if(length(grep("factor",x))!= 0 && length(grep(":",x))==0){
			pos1 <- grep("\\(",unlist(strsplit(x,split="")))+1
			pos2 <- grep("\\)",unlist(strsplit(x,split="")))-1
			compris.factor <- substr(x,start=pos1,stop=pos2)
			after.factor <- substr(x,start=(pos2+2),stop=length(unlist(strsplit(x,split=""))))
			paste(compris.factor,after.factor,sep="")
		}else if(length(grep("factor",x))!= 0 && length(grep(":",x))>0){
       if(grep('\\(',unlist(strsplit(x,split="")))[1]<grep(":",unlist(strsplit(x,split="")))[1] && length(grep('\\(',unlist(strsplit(x,split=""))))==1){
        pos1 <- grep("r",unlist(strsplit(x,split="")))[1]+2
		    pos2 <- grep("\\)",unlist(strsplit(x,split="")))[1]-1
		    compris.factor <- substr(x,start=pos1,stop=pos2)
		    after.factor <- substr(x,start=(pos2+2),stop= grep(":",unlist(strsplit(x,split="")))[1]-1)
		    pos3 <- grep(":",unlist(strsplit(x,split="")))[1]
		    pos4 <- length(unlist(strsplit(x,split="")))
		    compris.reste <- substr(x,start=pos3,stop=pos4)
		   
		    paste(compris.factor,after.factor,compris.reste,sep="")
		  }else if(grep("\\(",unlist(strsplit(x,split="")))[1]>grep(":",unlist(strsplit(x,split="")))[1] && length(grep('\\(',unlist(strsplit(x,split=""))))==1){
		    pos2 <- grep(":",unlist(strsplit(x,split="")))[1]
		    compris.reste <- substr(x,start=1,stop=pos2)
		    pos3 <- grep("\\(",unlist(strsplit(x,split="")))[1]+1
		    pos4 <- grep("\\)",unlist(strsplit(x,split="")))[1]-1
		    
		    compris.factor <- substr(x,start=pos3,stop=pos4)
		    after.factor <- substr(x,start=(pos4+2),stop= length(unlist(strsplit(x,split=""))))
		  
		    paste(compris.reste,compris.factor,after.factor,sep="")
		  }else{
		   
		    pos1 <- grep("\\(",unlist(strsplit(x,split="")))[1]+1
		    pos2 <- grep("\\)",unlist(strsplit(x,split="")))[1]-1
		    compris.factor1 <- substr(x,start=pos1,stop=pos2)
		    after.factor1 <- substr(x,start=(pos2+2),stop= grep(":",unlist(strsplit(x,split="")))[1]-1)
		    pos3 <- grep("\\(",unlist(strsplit(x,split="")))[2]+1
		    pos4 <- grep("\\)",unlist(strsplit(x,split="")))[2]-1
		    compris.factor2<- substr(x,start=pos3,stop=pos4)
		    after.factor2 <- substr(x,start=(pos4+2),stop= length(unlist(strsplit(x,split=""))))
		    		    paste(compris.factor1,after.factor1,":",compris.factor2,after.factor2,sep="")
		    
		  }
      }else{
			x
		}
	}
	)
	return(Names)
}

