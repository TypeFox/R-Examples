com.trait.moments <-
function(my.sample, traits){

mean.output = matrix(NA, nrow=nrow(my.sample), ncol=ncol(traits))
sd.output = matrix(NA, nrow=nrow(my.sample), ncol=ncol(traits))
skew.output = matrix(NA, nrow=nrow(my.sample), ncol=ncol(traits))
kurt.output = matrix(NA, nrow=nrow(my.sample), ncol=ncol(traits))

skewness = function(x){ 
			m3<-sum((x-mean(x))^3)/length(x)
			s3<-sqrt(var(x))^3
			m3/s3
			}

kurtosis = function(x){
			m4 <- sum((x-mean(x))^4)/length(x)
			s4 <- var(x)^2
			m4/s4 -3


}



	for(i in 1:ncol(traits)){

		kurt.funk = function(x){ 
			kurtosis(traits[names(x[x > 0]), i])
		}

		skew.funk = function(x){ 
			skewness(traits[names(x[x >0]), i])	
		}

		mean.funk = function(x){ 
			mean(traits[names(x[x >0]), i], na.rm = T)
		}

		sd.funk = function(x){ 
			sd(traits[names(x[x >0]), i], na.rm = T)
		}


mean.output[,i] = apply(my.sample, MARGIN = 1, mean.funk)	
sd.output[,i] = apply(my.sample, MARGIN = 1, sd.funk)
skew.output[,i] = apply(my.sample, MARGIN = 1, skew.funk)
kurt.output[,i] = apply(my.sample, MARGIN = 1, kurt.funk)

	}

output = cbind(mean.output, sd.output, skew.output, kurt.output)

colnames(output)=paste(c(rep("mean",ncol(traits)), rep("sd",ncol(traits)), rep("skew",ncol(traits)),rep( "kurtosis",ncol(traits)) ), rep(names(traits), 4), sep=".")
rownames(output) = rownames(my.sample)
output



}
