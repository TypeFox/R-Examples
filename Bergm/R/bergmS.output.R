bergmS.output <- function(x,
                          ...){
	

# model posterior
par(mfrow = c(1,2),oma=c(0,0,3,0),mar=c(4,3,1.5,1))
     
model.label <-rownames(rep(0,x$nmodels),
                       do.NULL=FALSE,
                       prefix="Model ")                       
     
plot(table(x$M)/length(x$M) , 
     type="h",xlab="",ylab="", 
     axes = FALSE,
     lwd=10)
axis(1,seq(1,x$nmodels),model.label)
axis(2)

plot(x$M, 
     type = "l", 
     xlab = "Iterations",
     ylab="",
     axes = FALSE)
axis(1)     
axis(2,seq(1,x$nmodels),model.label)
title("MCMC output : posterior model probabilities", outer = TRUE)     
    
# parameter posterior (for each model)
mcounts <- abs(sort(-table(x$M)))
m <- rep(0,length(mcounts))
for(i in 1:length(mcounts)){
	m[i] <- as.numeric(names(sort(-table(x$M)))[i])
	Theta <- as.matrix(x$Theta[x$M==m[i],1:x$dims[m[i]]]) 
	if(i==1){ 
		cat("\n", "BEST MODEL","\n", "----------")
	}
	cat("\n", "Model ",m[i],": ",paste(c(x$formula[[m[i]]])),"\n",sep="")
	cat("\n", "Posterior parameter estimate:", "\n")
	overall <- matrix(c(apply(Theta, 2, mean),apply(Theta, 2, sd)),x$dims[m[i]], 2)
	colnames(overall) <- c("Post. mean:", "Post. sd:")
	rownames(overall) <- paste("theta",seq(1,x$dims[m[i]]), 
	                           " (",x$spec[[m[i]]][seq(1,x$dims[m[i]])],")",sep="") 
	overall <- as.table(overall)
	print(overall)
                  
	cat(paste("\n", "Within-model acceptance rate:", 
	          round(x$Waccept[m[i]],2),"\n","\n"))
	          
	if(m[i]!=m[1]){
		cat(paste("BF_",m[1],m[i]," = ",mcounts[1]/mcounts[i],sep="","\n","\n"))
	}
}
T <- as.matrix(x$Theta[x$M==m[1],][,1:x$dims[m[1]]]) 
cat(paste("Between-model acceptance rate:",round(x$Baccept,2)))
	
# parameter posterior plot (for the best model)
G <- mcmc(data=T,start = 1,thin = 1)

dev.new()          
par(mfrow = c(min(4,x$dims[m[1]]), 3),oma=c(0,0,3,0),mar=c(4,3,1.5,1))
for(i in 1:x$dims[m[1]]){
	if(i%in%c(5,9,13)){
		dev.new()
		par(mfrow=c(min(4,x$dims[m[1]]-(i-1)),3),
		            oma=c(0,0,3,0),
		            mar=c(4,3,1.5,1))
	}
	plot(density(T[,i]), 
	     main = "", 
	     axes = FALSE, 
	     xlab=bquote(paste(theta[.(i)],
	     " (",.(x$spec[[m[1]]][i]),")")),
	     ylab = "", lwd = 2)
	axis(1)
	axis(2)
	plot(T[,i], type = "l", xlab = "Iterations",ylab = "")
	autocorr.plot(G[, i], auto.layout = FALSE,...)
	if(i%in%union(x$dims[m[1]],c(4,8,12))) title(paste("MCMC output for Model:",x$formula[m[1]]),outer=TRUE)
}

out=list(Theta=T,
         nchains=1,
         m=m[1],
         dim=x$dims[m[1]],
         model=x$models[[m[1]]],
         formula=x$formulae[[m[1]]],
         specs=x$specs[[m[1]]],
         BF=mcounts[m[1]]/mcounts[-1],
         mcounts=mcounts) 
out               
}

