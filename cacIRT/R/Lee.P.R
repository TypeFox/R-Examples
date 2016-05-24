Lee.P <-
function(cutscore, ip, theta, D = 1.7){
	
			ut <- theta
			ni <- dim(ip)[1]
			nn <- length(ut)
			sc <- ni+1
			nc <- length(cutscore)
			
			exp.TS <- rowSums(sapply(1:ni, function(i) ip[i,3] + 
							(1 - ip[i,3])/(1 + exp(-D*ip[i, 1] * 
							(ut-ip[i, 2])))))   
	
      rec.mat <- recursive.raw(ip, ut)
	
	esacc <- escon <-matrix(NA,nc,nn, dimnames = list(paste("cut at",round(cutscore,3)), round(ut,3)))

			for(j in 1:nc){
				cuts<-c(0, cutscore[j], sc)
	   			categ<-cut(exp.TS,cuts,labels=FALSE, right = FALSE)
				bang<-ceiling(cuts)
				rec.s <- list(NA)
					for(i in 1:2){
				rec.s[[i]] <- as.matrix(rec.mat[ , (bang[i]+1):bang[i+1]])}
					
			for(i in 1:nn){
			esacc[j,i]<- sum(rec.s[[categ[i]]][i,])}
			escon[j,]<- rowSums(rec.s[[1]])^2 + rowSums(rec.s[[2]])^2
			}
			
			if(nc > 1){ 
				simul <- matrix(NA,nn, 2, dimnames = list(round(ut,3), c("Accuracy", "Consistency")))
				cuts <- c(0, cutscore, sc)
				categ <- cut(exp.TS,cuts,labels=FALSE, right = FALSE)
				bang <- ceiling(cuts)
				rec.s <- list(NA)
					for(i in 1:(nc+1)){
				rec.s[[i]] <- as.matrix(rec.mat[ , (bang[i]+1):bang[i+1]])}
						
				
				for(i in 1:nn){
			simul[i,1]<- sum(rec.s[[categ[i]]][i,])}
			
			sha <- matrix(0,nn,1)
			for(i in 1:(nc+1)){
				sha <- sha + rowSums(rec.s[[i]])^2}
			simul[,2] <- sha
			
				
			
				ans<- (list("Marginal" = rbind(cbind("Accuracy" = rowMeans(esacc), "Consistency" = rowMeans(escon)), "Simultaneous" = colMeans(simul)), "Conditional" = list("Accuracy" =cbind(t(esacc), "Simultaneous" =simul[,1]), "Consistency" = cbind(t(escon),"Simultaneous" =simul[,2]))))
				return(ans)
				} else {
				
				ans<- (list("Marginal" = cbind("Accuracy" = rowMeans(esacc), "Consistency" = rowMeans(escon)), "Conditional" = list("Accuracy" =t(esacc), "Consistency" = t(escon))))
				
 
        return(ans)
				}
}

