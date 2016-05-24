Rud.D <-
function(cutscore, quadrature, sem) 
	 { 
	 	
	 	os <- quadrature[[1]]
		we <- quadrature[[2]]
	 	nn <- length(os)
	 	nc <- length(cutscore)
	
		esacc <- escon <-matrix(NA,length(cutscore), nn, dimnames = list(paste("cut at",round(cutscore,3)), round(os,3)))
		
		for(j in 1:length(cutscore)){
		  
			cuts<-c(-Inf, cutscore[j], Inf)		 	
	 		categ<-cut(os,cuts,labels=FALSE,right=FALSE)
			
		  for(i in 1:nn) {
			  esacc[j,i]<-(pnorm(cuts[categ[i]+1],os[i],sem[i])-pnorm(cuts[categ[i]],os[i],sem[i]))
			  escon[j,i]<-((pnorm(cuts[2], os[i],sem[i]) - pnorm(cuts[1],os[i],sem[i]))^2	+ (pnorm(cuts[3], os[i],sem[i]) - pnorm(cuts[2],os[i],sem[i]))^2	 )
		  }
	 	}
			
			if(nc == 1){
			ans<- (list("Marginal" = cbind("Accuracy" = apply(esacc,1,weighted.mean,we), "Consistency" = apply(escon,1,weighted.mean,we)), "Conditional" = list("Accuracy" =t(esacc), "Consistency" = t(escon))))
			
			return(ans)
			
			} else {
			  
				simul <- matrix(NA, nn, 2, dimnames = list(round(os,3), c("Accuracy", "Consistency")))
				cuts <- c(-Inf, cutscore, Inf)
				categ<-cut(os,cuts,labels=FALSE,right=FALSE)
				
				for(i in 1:nn){
					simul[i,1] <- (pnorm(cuts[categ[i]+1],os[i],sem[i])-pnorm(cuts[categ[i]],os[i],sem[i]))
					
					sha <- 0
					
			    for(j in 1:(nc+1)){
			      sha <- sha + (pnorm(cuts[j+1],os[i],sem[i])-pnorm(cuts[j],os[i],sem[i]))^2
			    }
					
			simul[i,2] <- sha
			}
			
			ans <- (list("Marginal" = rbind(cbind("Accuracy" = apply(esacc,1,weighted.mean,we), "Consistency" = apply(escon,1,weighted.mean,we)), "Simultaneous" = apply(simul,2,weighted.mean,we) ), "Conditional" = list("Accuracy" =cbind(t(esacc), "Simultaneous" =simul[,1]), "Consistency" = cbind(t(escon),"Simultaneous" =simul[,2]))))

				return(ans)
			}
}

