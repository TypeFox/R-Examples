at.level<-function(x,level){
          if(is.factor(x)==FALSE){stop("x in at.level is not a factor")}
	  if(is.numeric(level)){
	  	 M<-outer(x, levels(x)[level], "==")
	  	}else{
	     M<-outer(x, level, "==")
		}
	    mode(M)<-"numeric"
	    M
	}
	
at.set<-function(x,level){
         if(is.factor(x)==FALSE){stop("x in at.set is not a factor")}
	  if(is.numeric(level)){
	    M<-x%in%(levels(x)[level])
	  	}else{
	    M<-x%in%level
		}
	    mode(M)<-"numeric"
	   as.matrix(M)
	}
		
#leg<-function(x,degree, normalized=TRUE){
#             if(requireNamespace("orthopolynom", quietly = TRUE)==FALSE){
#               stop("orthopolynom not loaded")
#             } 			
#	     lp<-orthopolynom::legendre.polynomials(n=abs(degree), normalized=normalized)
#             if(degree<0){
#               lp<-lp[-1]
#             }			
#             M<-sapply(lp,function(lp){as.function(lp)(x)})
# 	     colnames(M)<-paste(as.character(as.list(substitute(list(x)))[[2]]),  (0+1*(degree<0)):abs(degree), sep=".")
#             if(degree>=0){
#               M[,1][which(is.na(M[,1]))]<-as.function(lp[[1]])(0)
#             }
#	     M
# }
