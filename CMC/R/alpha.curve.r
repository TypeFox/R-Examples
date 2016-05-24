alpha.curve = function(x){
	data = x	
	n = nrow(data) #n. of obs
	k = ncol(data) #n. of item
	max.vec 			= c()
	which.max.vec = c()
	
	#Compute alpha for the complete dataset
	complete.alpha = alpha.cronbach(data)#,n,k)

	#Compute alpha max removing j elements from the complete dataset
	j = 1 #(1 element to remove)
	elements.to.remove.1	= seq(1:k)
	alpha.1				 		= c()
	for(i in 1:length(elements.to.remove.1)){
		data.reduced 	= data[,-elements.to.remove.1[i]]
		#n.reduced 		= nrow(data.reduced)
		#k.reduced		= ncol(data.reduced)
		alpha.1[i]		= alpha.cronbach(data.reduced)#,n.reduced,k.reduced)
	}
	max.vec[j]	 			= max(alpha.1)
	which.max.vec[j]   	= which.max(alpha.1)
	
	for (j in 2:(k-2)){
		elements.to.remove	= matrix(0, nrow=(k-j+1), ncol=(j-1))
		
		for(r in 1:length(which.max.vec)){
			elements.to.remove[,r] = matrix(rep(which.max.vec[r], k-j+1), ncol=1)
			}
		elements.to.remove = cbind(elements.to.remove,seq(1:k)[-which.max.vec[1:(j-1)]])

		alpha				 	= c()
		for(i in 1:nrow(elements.to.remove)){
			data.reduced 	= data[,-elements.to.remove[i,]]
			#n.reduced 		= nrow(data.reduced)
			#k.reduced		= ncol(data.reduced)
			alpha[i]			= alpha.cronbach(data.reduced)#,n.reduced,k.reduced)
		}
		
		max.vec[j]			= max(alpha)
		which.max.vec[j]	= elements.to.remove[which.max(alpha),j]
	}
	
	output = data.frame("N.Item" = seq(2,k),
								 "Alpha Max" = c(rev(max.vec),complete.alpha), 
								 "Removed Item" = c(colnames(data)[rev(which.max.vec)],"--"))
	
	plot(output[,1],output[,2],t="b",xlab="N.Item", ylab="Alpha Max")
	text(seq(2,k),output[,2],c(colnames(data)[rev(which.max.vec)],""),pos=3,cex=0.6)
	
	return(output)
	}
