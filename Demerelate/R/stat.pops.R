stat.pops <- function(Thresholds, tab.pop.pop,  pairs, p.correct, directory.name, out.name, file.output, inputdata, object, value, iteration, ref.pop)

    {

    number.loci <- (ncol(tab.pop.pop)-2)/2
    empirical.share.ls <- vector("list",number.loci)
    empirical.share.ls.pm <- empirical.share.ls
    relate.non.X <- vector("list",number.loci)
    
    							   
    # Calculation of value for each locus in population tab.pop.pop
    for (i in 1:number.loci)
				{
				message(paste("---","Calculations for empirical values are performed for Locus",i,"----",Sys.time()),"\n")

				# Empirical share calculated for each locus
        empirical.share.ls[[i]] <- allele.sharing(tab.pop.pop,tab.pop.pop,(i*2)+1,FALSE, value, ref.pop)
				names(empirical.share.ls)[i] <- paste(names(tab.pop.pop)[(i*2)+1],names(tab.pop.pop)[(i*2)+2],sep="-")
     							
				}
				
    			# Mean empirical
			empirical.list <- do.call("cbind",lapply(empirical.share.ls,array))	
			empirical.list <- matrix(rowMeans(empirical.list,na.rm=TRUE),nrow(empirical.share.ls[[1]]))
			row.names(empirical.list) <- row.names(empirical.share.ls[[1]])
			colnames(empirical.list) <- colnames(empirical.share.ls[[1]])

		

	# Testing on NAs
			if (length(table(is.na(as.dist(empirical.list))))==2)
				{ 

warning(" ############## ERROR FULL STOP ################ ")
warning(" ############## ERROR FULL STOP ################ ")
warning(" ############## ERROR FULL STOP ################ ")
warning(" ## DON'T PANIK JUST READ FURTHER INSTRUCTIONS ## ")
warning("---","\n","\n")
warning("---","\n","\n")
warning(paste("An error occured in calculating allele sharings for population ",tab.pop.pop[1,2],"----",Sys.time(),sep=" "))
warning("---","\n","\n")
warning("Due to too much missing values in inputdata some individual sharings could not be calculated ----")
warning("Check on output for further information")
warning("---","\n","\n") 
warning("NAs mark invalid combinations for individual allele sharings from individuals indicated in row and column names of 'error.individuals' ----")
warning("Please remove at least one individual indicated in 'error.individuals' or NAs from at least one individual from your input.txt and retry analysis")
warning("---","\n","\n")
warning("error.pairings (output=position in emp):")
warning("---","\n","\n")
error.individuals <- which(is.na(as.dist(emp)),arr.ind=TRUE)
error.individuals
warning("---","\n","\n")
warning("If NAs not removable try starting analysis with mode cluster=FALSE; instead of cluster=TRUE (default)")
return(error.individuals)

if (file.output==TRUE)
{ write.table(file="error.individuals.txt",x=error.individuals)
}
stop()
}




          
    # Clusteranalysis
if (file.output==TRUE)
{		
		pdf(paste(".","/",directory.name,"/","Cluster",tab.pop.pop[1,2],out.name,".pdf",sep=""))
        

				par(cex=0.7,font=3)				
				dis.plot <- plot(hclust(as.dist(1-empirical.list),method="complete"),xlab="",sub=paste("Inter-individual Relation in population", tab.pop.pop[1,2],sep=" "), main=paste(value,"normalized dissimilarity"))
				abline(h=1-Thresholds[[1]], col="red", lty="dotdash")
				abline(h=1-Thresholds[[2]], col="blue", lty="dashed")

		for (i in 1:length(empirical.share.ls))
				{ 
					hist(empirical.share.ls[[i]], main = paste("Histogram of", names(empirical.share.ls)[i], "in", tab.pop.pop[1,2], sep=" "), xlab=paste("Empirical relatedness [",value,"]"))
					abline(v=Thresholds[[1]], col="red", lty="dotdash")
					abline(v=Thresholds[[2]], col="blue", lty="dashed")
				}
		
		hist(empirical.list, main = paste("Histogram of normalized mean allele sharing in", tab.pop.pop[1,2], sep=" "), xlab=paste("Empirical relatedness [",value,"]"))
		abline(v=Thresholds[[1]], col="red", lty="dotdash")
		abline(v=Thresholds[[2]], col="blue", lty="dashed")
			
		dev.off()
}
    
    
    if (value=="rxy") 
      {} 
    else {
		
	    
    # length of reference population == empirical population
	
		# Random non-related for X-square
				
				relate.non.X <- lapply(seq(1:number.loci),function(i)
             {
				          random.pairs.non.ls.X <- random.pairs(ref.pop,(i*2)+1,length(empirical.list[!is.nan(as.numeric(empirical.list))]))
                  allele.sharing(random.pairs.non.ls.X[[1]],random.pairs.non.ls.X[[2]],3,TRUE, value)
				      })

		relate.non.X.mean <- do.call("cbind",lapply(relate.non.X,array))	
		relate.non.X.mean <- matrix(rowMeans(relate.non.X.mean,na.rm=TRUE),nrow(relate.non.X[[1]]))
		relate.non.X.mean[relate.non.X.mean>=Thresholds[1,2]] <- "FS"
		relate.non.X.mean[relate.non.X.mean>Thresholds[1,1] & relate.non.X.mean<Thresholds[1,2]] <- "HS"
		relate.non.X.mean[relate.non.X.mean<=Thresholds[1,1]] <- "NON"
    
    empirical.list.values <- empirical.list
    		# Matrix conversion
    		empirical.list[empirical.list>=Thresholds[1,2]] <- "FS"
    		empirical.list[empirical.list>Thresholds[1,1] & empirical.list<Thresholds[1,2]] <- "HS"
    		empirical.list[empirical.list<=Thresholds[1,1]] <- "NON"
		
    
		# P Statistics
		# only FS
		emp <- 0
		non <- 0
		if (!is.na(as.numeric(table(empirical.list)["FS"]))==TRUE){emp <- as.numeric(table(empirical.list)["FS"])}
		if (!is.na(as.numeric(table(relate.non.X.mean)["FS"]))==TRUE){non <- as.numeric(table(relate.non.X.mean)["FS"])}
		stat.p <- prop.test(c(emp,non), c(sum(table(empirical.list)), sum(table(relate.non.X.mean))), conf.level=0.95,correct=p.correct)

		# FS + HS
		emp <- 0
		non <- 0
		if (!is.na(sum(as.numeric(table(empirical.list)["HS"])==TRUE, as.numeric(table(empirical.list)["FS"])==TRUE, na.rm=TRUE))){emp <- sum(as.numeric(table (empirical.list)["FS"]), as.numeric(table(empirical.list)["HS"]),na.rm=TRUE)}
		if (!is.na(sum(as.numeric(table(relate.non.X.mean)["HS"])==TRUE, as.numeric(table(relate.non.X.mean)["FS"])==TRUE, na.rm=TRUE))){non <- sum(as.numeric(table(relate.non.X.mean)["FS"]), as.numeric(table(relate.non.X.mean)["HS"]),na.rm=TRUE)}
		stat.p.HS <- prop.test(c(emp,non),c(sum(table(empirical.list)),sum(table(relate.non.X.mean))),conf.level=0.95,correct=p.correct)
        
		f.p <- as.data.frame(rbind(c(stat.p[[4]][1], stat.p[[4]][2],stat.p[[1]],stat.p[[2]],stat.p[[3]],stat.p[[6]][1],stat.p[[6]][2]),c(stat.p.HS[[4]][1],stat.p.HS[[4]][2],stat.p.HS[[1]],stat.p.HS[[2]],stat.p.HS[[3]],stat.p.HS[[6]][1],stat.p.HS[[6]][2])))
				
		row.names(f.p) <- c("Full Siblings","Full + Half Siblings")
		colnames(f.p) <- c("Prob. of Observed","Prob. of Expected","Chi^2-value","d.f.","p-value","0.95 Lower CI","0.95 Upper CI")

if (file.output==TRUE)
{
  out.file <- file(paste(".","/",directory.name,"/","Relate.mean",tab.pop.pop[1,2],out.name,".txt",sep=""),"w")
  writeLines(
    paste(
      "Demerelate - v.0.8-1", "---","\n","Relatedness outputfile on file:", inputdata,"\n","Analysis had been made using", 
      iteration,"iterations","and",pairs,"pairs","using the",value,"estimator.","\n",

  if (value=="Bxy"){paste("Calculations are based on Li and Horvitz 1953. The values represent an indication on relatedness based on allele sharing.","\n",sep=" ")},
  if (value=="Mxy"){paste("Calculations are based on Bluoin et al. 1996. The values represent relatedness assessment based on genotype sharing.","\n",sep=" ")},
  if (value=="rxy"){paste("Calculations are based on Queller and Goodnight 1989. The values represent relatedness value corrected for total allele diversity.","\n",sep=" ")},
  "For further information mind References at the end of this file.","\n","\n",
  "Calculations had been made for population:", as.character(tab.pop.pop[1,2]),"\n",
  "\n",
  "Relatedness Thresholds","\n",
  "---","\n"),con=out.file)
  write.table(Thresholds, file=out.file, append=T, sep="\t", quote=F, row.names=F) 
              writeLines(paste("\n","---","\n","\n",
  "Relatedness calculations","\n",
  "---","\n",
	"Observed frequencies of full siblings (FS), half siblings (HS) and non related pairs (NON)","\n"),con=out.file)
              write.table(table(empirical.list),file=out.file,append=T,sep="\t", quote=F, row.names=F, col.names=F)
              
  writeLines(paste("\n",
  "---","\n","\n",
  "---","\n",
	"Expected frequencies of full siblings (FS), half siblings (HS) and non related pairs (NON)","\n"),con=out.file)
              
  write.table(table(relate.non.X.mean),file=out.file,append=T,sep="\t", quote=F, row.names=F, col.names=F)
  
  writeLines(paste("---","\n","\n","Chisquare Statistics","\n","---","\n","\n",sep=" "),con=out.file)
    
  write.table(f.p,file=out.file,append=T,sep="\t", quote=F)
              
  writeLines(paste("\n",

  "----------------------------------------------------------------------------------------------------","\n",

      
  "\n","\n","\n","References","\n",
  "Blouin, M.S. et al. (1996) Use of microsatellite loci to classify individuals by relatedness. Molecular Ecology, 5, 393-401.","\n",
  "Li C.C. and Horvitz D.G. (1953) Some methods of estimating the inbreeding coefficient. American Journal of Human Genetics 5, 107-17.","\n",
  "Queller, D.C. and Goodnight, K.F. (1989) Estimating relatedness using genetic markers. Evolution, 43, 258-275.","\n", sep=" "), 
    con=out.file)
    
    close(out.file)

}
		empirical.list <- empirical.list.values
				
    }
  
  if (value=="rxy")
    
  {
    f.p <- "NA"
  }
    
      out.stat <- list(empirical.list, f.p)
      names(out.stat) <- c("Empirical_List", "Chi-square statistics")
      return(out.stat)
    
  
    
  }
