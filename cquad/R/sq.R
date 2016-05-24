sq <-
function(J,s=NULL){

# provides the matrix of all possible binary sequences of length J with sum equal s
# if s is omitted all possibile binary sequences of length J

	M = NULL
	if(J>0){
		if(!is.null(s)){                  # s is not provided
			if(s == J) M = matrix(1,1,J)  # alle elements equal 1
            if(s>1 & s<J) for(i in 1:(J-s+1)){
				S = sq(J-i,s-1)
				if(is.null(S)) r = 0
				else r = dim(S)[1]
	        	M1 = cbind(matrix(0,r,i-1),matrix(1,r,1),S)
        		M = rbind(M1,M)
        	}
    		if(s==1){                    # with only one element equal to 1       
				M = matrix(0,J,J)
      			for(j in 1:J) M[j,J-j+1] = 1
      		}
    		if(s==0) M = matrix(0,1,J)   # without element equal to 1
    	}
		else{                            # if s is given
	    	if(J==1){
	    		M = matrix(c(0,1),2,1)
	    	}else{
				M1 = sq(J-1)
	      		M = rbind(cbind(0,M1),cbind(1,M1))
	    	}
      	}
	}
	return(M)
}
