com.trait.weighted <-
function(my.sample, traits){

mean.output = matrix(NA, nrow=nrow(my.sample), ncol=ncol(traits))
sd.output = matrix(NA, nrow=nrow(my.sample), ncol=ncol(traits))


	for(i in 1:ncol(traits)){
		weight.mean.funk = function(x){
			weighted.mean(traits[names(x[x>0]), i] , x[x > 0])
		}

		weight.sd.funk = function(x){
			wt.sd(traits[names(x[x>0]), i] , x[x>0])	
		}



		mean.output[,i] = apply(my.sample, MARGIN = 1, weight.mean.funk)	
		sd.output[,i] = apply(my.sample, MARGIN = 1, weight.sd.funk)
	}


output = cbind(mean.output, sd.output)
colnames(output)=paste(c(rep("mean",ncol(traits)), rep("sd",ncol(traits)) ), rep(names(traits), 2), sep=".")
rownames(output) = rownames(my.sample)
output

}
