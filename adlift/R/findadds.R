`findadds` <-
function (rem, neighbrs, po, index = 1:(length(rem) + length(po))) 
{

l <- length(rem)
n <- l + length(po)			#data length
adds <- NULL
wh <- init <- init2<-NULL
init <- match(index, rem)		#finds where index[i] is lifted
init[is.na(init)] <- l + 1		#must be in po, so don't already a scaling coeff.
init<-l-(init-1)			#finds out how many steps are needed for this part.

if (l == 0) {
	init2 <- rep(1,times=length(index))
}
else{
	for (i in 1:length(index)) {
	        wh <- NULL
       		for (j in 1:l) {
               		wh[j] <- any(neighbrs[[j]] == index[i])
       		}
       		if (sum(wh) == 0) {		#index[i] isn't ever a neighbour, so no minimum steps needed.
               		init2[i] <- l + 1
      		 }
       		else {
              		init2[i] <- min(which(wh==1))
       		}
       	}
}
init2<-l-(init2-1)              #finds out how many adding steps to perform.

tmp<-rbind(init,init2)
adds<-apply(tmp,2,max)

adds
}

