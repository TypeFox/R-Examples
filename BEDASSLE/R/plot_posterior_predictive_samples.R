plot_posterior_predictive_samples <-
function(posterior.predictive.sample.file,save.figure=NULL,figure.name=NULL){
		posterior.predictive.samples <- load_posterior_predictive_samples(posterior.predictive.sample.file)
		with(posterior.predictive.samples, {	
		distance.array <- array(D,dim=dim(posterior.sample.Fst))
		number.of.populations <- nrow(observed.Fst)
		population.pairs <- combn(1:number.of.populations,2)
		
		if(is.null(save.figure)){
				plot(as.numeric(D[upper.tri(D)]),as.numeric(observed.Fst[upper.tri(observed.Fst)]),col="red",ylim=c(0,max(posterior.sample.Fst)+max(posterior.sample.Fst)/20),pch=19,cex=0.3,main="Posterior Predictive Sampling",xlab=expression(paste("Pairwise Geographic Distance ",(D[ij]),sep="")),ylab=expression(paste("Pairwise ",F[ST],sep="")))
					prog <- txtProgressBar(min=0,ncol(population.pairs),char="|",style=3)
				for(i in 1:ncol(population.pairs)){
					points(as.numeric(distance.array[population.pairs[,i][1],population.pairs[,i][2],]),as.numeric(posterior.sample.Fst[population.pairs[,i][1],population.pairs[,i][2],]),pch=19,col=adjustcolor(1,0.02),cex=0.2)
						setTxtProgressBar(prog,i)
				}
				points(as.numeric(D),as.numeric(observed.Fst),col="red",pch=19,cex=0.4)
		}
		if(!is.null(save.figure)){
			png(figure.name,width=6*200,height=4*200,res=200,pointsize=9)
				plot(as.numeric(D[upper.tri(D)]),as.numeric(observed.Fst[upper.tri(observed.Fst)]),col="red",ylim=c(0,max(posterior.sample.Fst)+max(posterior.sample.Fst)/20),pch=19,cex=0.3,main="Posterior Predictive Sampling",xlab=expression(paste("Pairwise Geographic Distance ",(D[ij]),sep="")),ylab=expression(paste("Pairwise ",F[ST],sep="")))
					prog <- txtProgressBar(min=0,ncol(population.pairs),char="|",style=3)
				for(i in 1:ncol(population.pairs)){
					points(as.numeric(distance.array[population.pairs[,i][1],population.pairs[,i][2],]),as.numeric(posterior.sample.Fst[population.pairs[,i][1],population.pairs[,i][2],]),pch=19,col=adjustcolor(1,0.02),cex=0.2)
						setTxtProgressBar(prog,i)
				}
				points(as.numeric(D),as.numeric(observed.Fst),col="red",pch=19,cex=0.4)
			dev.off()
		}
		})
	}
