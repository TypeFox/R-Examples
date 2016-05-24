relate.calc <- function(tab.pop, pairs, file.output, value, directory.name, ref.pop)
	{
		

		number.loci <- (length(tab.pop)-2)/2
   	empirical.share.ls <- vector("list",number.loci)
		relate.full.mean <- vector("list",number.loci)
		relate.half.mean <- vector("list",number.loci)
		relate.off.non <- vector("list",number.loci)
		random.pairs.fsib.ls <- vector("list",2)
		random.pairs.hsib1.ls	<- vector("list",2)
		random.pairs.hsib2.ls <- vector("list",2)
    
		# Random pairs for fullsibs overall
		random.pairs.fsib.ls[[1]] <- ref.pop[sample(1:nrow(ref.pop),pairs,replace=TRUE),]
    random.pairs.fsib.ls[[1]][,1]<-paste("FS-",random.pairs.fsib.ls[[1]][,1],"-",seq(1:length(random.pairs.fsib.ls[[1]][,1])),sep="")
		random.pairs.fsib.ls[[2]] <- ref.pop[sample(1:nrow(ref.pop),pairs,replace=TRUE),]
		random.pairs.fsib.ls[[2]][,1]<-paste("FS-",random.pairs.fsib.ls[[2]][,1],"-",seq((length(random.pairs.fsib.ls[[1]][,1])+1):(length(random.pairs.fsib.ls[[2]][,1])+length(random.pairs.fsib.ls[[1]][,1]))),sep="")
		# Random pairs for halfsibs overall 
		random.pairs.hsib1.ls[[1]] <- ref.pop[sample(1:nrow(ref.pop),pairs,replace=TRUE),]
    random.pairs.hsib1.ls[[1]][,1]<-paste("HS-",random.pairs.hsib1.ls[[1]][,1],"-",seq(1:length(random.pairs.hsib1.ls[[1]][,1])),sep="")
		random.pairs.hsib2.ls[[1]] <- ref.pop[sample(1:nrow(ref.pop),pairs,replace=TRUE),]
    random.pairs.hsib2.ls[[1]][,1]<-paste("HS-",random.pairs.hsib2.ls[[1]][,1],"-",seq((length(random.pairs.hsib1.ls[[1]][,1])+1):(length(random.pairs.hsib1.ls[[1]][,1])+length(random.pairs.hsib2.ls[[1]][,1]))),sep="")
		random.pairs.hsib2.ls[[2]] <- ref.pop[sample(1:nrow(ref.pop),pairs,replace=TRUE),]
    random.pairs.hsib2.ls[[2]][,1]<-paste("HS-",random.pairs.fsib.ls[[1]][,1],"-",seq((length(random.pairs.hsib1.ls[[1]][,1])+length(random.pairs.hsib2.ls[[1]][,1])+1):(length(random.pairs.hsib1.ls[[1]][,1])+length(random.pairs.hsib2.ls[[1]][,1])+length(random.pairs.hsib2.ls[[2]][,1]))),sep="")
		
       for (i in 1:number.loci)	
				{
		
		
		message(paste("---","Calculations are performed for Locus",i,"----",Sys.time(),"\n"))
		
			# 1. Empirisches Sharing for each locus calculated
				empirical.share.ls[[i]] <- allele.sharing(tab.pop,tab.pop,(i*2)+1, FALSE, value, ref.pop)
				names(empirical.share.ls)[i] <- names(tab.pop)[(i*2)+1]
	
			
		 	# 4. Random non-related pairs calculated
				  random.pairs.non.ls <- random.pairs(ref.pop,(i*2)+1,pairs)
		  		relate.off.non[[i]] <- allele.sharing(random.pairs.non.ls[[1]],random.pairs.non.ls[[2]],3,TRUE, value, data.frame(ref.pop[,1],ref.pop[,2],ref.pop[,(i*2)+1],ref.pop[,(i*2)+1+1]))
				   

			# Offsprings for reference are calculated
					
					# 1. Random Offspring full
					off.full.ls.1 <- offspring(random.pairs.fsib.ls[[1]],random.pairs.fsib.ls[[2]],(i*2)+1, pairs)
					off.full.ls.2 <- offspring(random.pairs.fsib.ls[[1]],random.pairs.fsib.ls[[2]],(i*2)+1, pairs)
					relate.full.mean[[i]] <- allele.sharing(off.full.ls.1,off.full.ls.2,3,TRUE, value, data.frame(ref.pop[,1],ref.pop[,2],ref.pop[,(i*2)+1],ref.pop[,(i*2)+1+1]))
					
					# 2. Random Offspring half
					off.half.ls.1 <- offspring(random.pairs.hsib1.ls[[1]],random.pairs.hsib2.ls[[1]],(i*2)+1, pairs)
					off.half.ls.2 <- offspring(random.pairs.hsib1.ls[[1]],random.pairs.hsib2.ls[[2]],(i*2)+1, pairs)
					relate.half.mean[[i]] <- allele.sharing(off.half.ls.1,off.half.ls.2,3,TRUE, value, data.frame(ref.pop[,1],ref.pop[,2],ref.pop[,(i*2)+1],ref.pop[,(i*2)+1+1]))
					
		   		}
	
		
			# Empirical

      empirical.Mxy.mean <- do.call("cbind",lapply(empirical.share.ls,array))
			empirical.Mxy.mean <- matrix(rowMeans(empirical.Mxy.mean,na.rm=TRUE),nrow(empirical.share.ls[[1]]))
      row.names(empirical.Mxy.mean) <- paste(tab.pop[,1],tab.pop[,2],sep="-")
		  colnames(empirical.Mxy.mean) <- paste(tab.pop[,1],tab.pop[,2],sep="-")
      
			relate.off.full.Mxy.mean <- do.call("cbind",lapply(relate.full.mean,array))
			relate.off.full.Mxy.mean <- matrix(rowMeans(relate.off.full.Mxy.mean,na.rm=TRUE),nrow(relate.full.mean[[1]]))
						
			relate.off.non.Mxy.mean <- do.call("cbind",lapply(relate.off.non,array))	
			relate.off.non.Mxy.mean <- matrix(rowMeans(relate.off.non.Mxy.mean,na.rm=TRUE),nrow(relate.off.non[[1]]))
			
			relate.off.half.Mxy.mean <- do.call("cbind",lapply(relate.half.mean,array))	
			relate.off.half.Mxy.mean <- matrix(rowMeans(relate.off.half.Mxy.mean,na.rm=TRUE),nrow(relate.half.mean[[1]]))

		# Calculating multiple logistic regression
		return.glm <- glm.prep(empirical.Mxy.mean, relate.off.full.Mxy.mean, relate.off.half.Mxy.mean, relate.off.non.Mxy.mean)
    half <- return.glm[[1]]
    sumlrm <- return.glm[[2]]
		  
    # Threshold calculation
    full <- (sumlrm[[1]][1]-sumlrm[[1]][2])/(sumlrm[[1]][4]-sumlrm[[1]][3])		
		Thres <- data.frame(half,full)
    		    
if (file.output==TRUE)
{
    write.table(file=paste(".","/",directory.name,"/","Random.Fullsib.relatedness.overall.txt",sep=""),x=relate.off.full.Mxy.mean, quote=FALSE, sep=" ")
    write.table(file=paste(".","/",directory.name,"/","Random.Halfsib.relatedness.overall.txt",sep=""),x=relate.off.half.Mxy.mean, quote=FALSE, sep=" ")
    write.table(file=paste(".","/",directory.name,"/","Random.NonRelated.relatedness.overall.txt",sep=""),x=relate.off.non.Mxy.mean, quote=FALSE, sep=" ")

}
		# Assigning Output
		relate.return <- list(empirical.Mxy.mean, relate.off.full.Mxy.mean, relate.off.half.Mxy.mean, relate.off.non.Mxy.mean, Thres)
    names(relate.return) <- c("Relatedness_Empirical"," Randomized_Fullssibs","Randomized_Halfsibs","Ranodmized_Non","Thresholds")
		return(relate.return)
			
	}	

