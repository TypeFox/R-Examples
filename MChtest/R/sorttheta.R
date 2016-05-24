"sorttheta" <-
function(theta){
	## produce a ``standard ordering'' of theta
	w<-as.vector(theta[,1])
	m<-as.vector(theta[,2])
	v<-m-w
	s<-length(w)
	
	sorted.theta<-matrix(NA,s,2)
	for (i in 1:s){
		if (i==1){
			sorted.theta[1,]<-c(w[v==0],m[v==0])
			current.w<-sorted.theta[1,1]
			current.m<-sorted.theta[1,2]
			top<-1
			inn<-v==0
			}
		else {
			### else number 1
			index<-w==current.w & m==(current.m+1) & !inn
			if (any(index)==TRUE){
				sorted.theta[i,]<-c(w[index],m[index])
				current.w<-sorted.theta[i,1]
			    current.m<-sorted.theta[i,2]
			    inn<-inn | index
			    top<-i
			}
			else {
			### else number 2
				index<-w==(current.w+1) & m==(current.m+2) & !inn
				if (any(index)==TRUE){
				sorted.theta[i,]<-c(w[index],m[index])
				current.w<-sorted.theta[i,1]
			    current.m<-sorted.theta[i,2]
			    inn<-inn | index
			    top<-i
			    }
				else {
					### else number 3
					index<-w==(current.w-1) & m==(current.m) & !inn
					if (any(index)==TRUE){
					sorted.theta[i,]<-c(w[index],m[index])
					current.w<-sorted.theta[i,1]
			    	current.m<-sorted.theta[i,2]
			    	inn<-inn | index
			    	}
					else {
						### else number 4
						index<-w==(current.w-1) & m==(current.m-1) & !inn
						if (any(index)==TRUE){
						sorted.theta[i,]<-c(w[index],m[index])
						current.w<-sorted.theta[i,1]
			    		current.m<-sorted.theta[i,2]
			    		inn<-inn | index
			    		}
						else {
						### else number 5
						index<-w==(current.w-1) & (!inn) 
						if (any(index)==TRUE){
						if (length(w[index])>1){ stop("error in sort function")}	
						sorted.theta[i,]<-c(w[index],m[index])
						current.w<-sorted.theta[i,1]
			    		current.m<-sorted.theta[i,2]
			    		inn<-inn | index
			    		}
						else{
								### else number 6
						 		stop("error 2 in sorting function")	
						} ### end else 6
                     } ### end else 5
                 } ### end else 4
				}	### end else 3
			}   ### end else 2
		}  ### end else 1
	} ### end for loop
out<-list(theta=sorted.theta,top=top)
out
}

