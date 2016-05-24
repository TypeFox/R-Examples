Fis.calc <- function(tab.pop, iteration, number.loci, object, directory.name, out.name)
{
    # input: table of allele frequencies (tab.freq) and table of genotype frequencies (tab.freq.gen)
    empirical.fis <- vector(mode="numeric", number.loci)
    empirical.weir <- vector(mode="numeric", number.loci)
    weir.loci <- vector("list", number.loci)
    p.va <- vector(mode="numeric", number.loci)
    p.weir <- vector(mode="numeric", number.loci)
    b.w <- vector(mode="numeric", number.loci)
    c.w <- vector(mode="numeric", number.loci)
    bootf <- vector("list", number.loci)
   
    
    if(directory.name!="NA" & out.name!="NA") 
    {
    out.file <- file(as.character(paste(".","/",directory.name,"/","Summary",tab.pop[1,2],out.name,".txt",sep="")),"w")
    writeLines(paste(
      "Demerelate - v.0.8-1", "---","\n",
      "Summary outputfile on file:", out.name,"\n",
      "Analysis had been made using", iteration,"iterations.","\n",
      "Populations in inputdata:", tab.pop[1,2],"\n",
      "\n", sep=" "),con=out.file)
    }

    for (i in 1:number.loci)
    {
      allele.column<-i*2+1
      fis.return <- Fis(tab.pop,i*2+1)
      empirical.fis[i] <- fis.return[[5]]
      weir.loci[[i]] <- fis.return[[6]][[1]]
      empirical.weir[i] <- fis.return[[6]][[2]]
      boots.fis <- vector(mode="numeric",iteration)
      boots.weir <- vector(mode="numeric",iteration)
      
      if(directory.name!="NA" & out.name!="NA") 
      {
      writeLines(paste(     
        
        "Locus:", paste(names(tab.pop)[allele.column],names(tab.pop)[allele.column+1]),"\n",
        "Allele diversity:","\n"),con=out.file)
      
      write.table(fis.return[[1]], file=out.file, append=T, quote=F, sep="\t",col.names=F)
      
      writeLines(paste(
        "---","\n","\n",
        "Genotype diversity:","\n"),con=out.file)
      
      write.table(fis.return[[2]],file=out.file, append=T, quote=F, sep="\t",row.names=F,col.names=F)
      
      writeLines(paste(
        "---","\n","\n",
        "Fis per allele (Weir and Cockerham 1984):","\n"),con=out.file)
      
      write.table(fis.return[[6]][[1]], out.file, append=T, quote=F, sep="\t")
      writeLines("\n Fis for locus (Weir and Cockerham 1984):",con=out.file)
      write.table(fis.return[[6]][[2]], out.file, append=T, quote=F, sep="\t", col.names=F, row.names=F)
      
      writeLines(paste(
        "---","\n","\n",
        "Heterozygosity:", 1-fis.return[[3]],"\n","\n",
        "----------------------------------------------------","\n","\n"),con=out.file)
      }
          
        for (j in 1:iteration)
        {
            bootstr1 <- sample(c(tab.pop[,allele.column],tab.pop[,allele.column+1]),length(tab.pop[,1]),replace=TRUE)
            bootstr2 <- sample(c(tab.pop[,allele.column],tab.pop[,allele.column+1]),length(tab.pop[,1]),replace=TRUE)
            bootstr <- data.frame(tab.pop[,1],tab.pop[,2], bootstr1, bootstr2)
            fis.return.boot <- Fis(bootstr,3)
            boots.weir[j] <- fis.return.boot[[6]][[2]]
            c.w[j] <- sum(fis.return.boot[[6]][[1]][3,])
            b.w[j] <- sum(fis.return.boot[[6]][[1]][2,])
            boots.fis[j] <- fis.return.boot[[5]]
	          
        }
      
      bootf[[i]] <- cbind(c.w,b.w)
      
      if (is.nan(empirical.fis[i])==TRUE)
      {p.va[i]  <- NA}
      if (empirical.fis[i]>0)
      {p.va[i]  <- (1+sum(boots.fis >= empirical.fis[i]))/(iteration+1)}
      if (empirical.fis[i]<0)
      {p.va[i]  <- (1+sum(boots.fis <= empirical.fis[i]))/(iteration+1)}
      if (empirical.fis[i]==0)
      {p.va[i]  <- 1}
      
      if (is.nan(empirical.weir[i])==TRUE)
      {p.weir[i] <- NA}  
      if(empirical.weir[i]>0)
      {p.weir[i] <- (1+sum(boots.weir >= empirical.weir[i]))/(iteration+1)}
      if(empirical.weir[i]<0)
      {p.weir[i] <- (1+sum(boots.weir <= empirical.weir[i]))/(iteration+1)}
      if(empirical.weir[i]==0)
      {p.weir[i]  <- 1}
     }
    
     for (k in 2:number.loci)
     {
       bootf[[1]] <- bootf[[1]]+bootf[[k]] 
     }
    
       bootf[[1]] <- 1-(bootf[[1]][,1]/(bootf[[1]][,1]+bootf[[1]][,2]))
    
    # Empirical weighted Weir Fis over loci
    b.weir <- sapply(lapply(weir.loci,function(x){sapply(x[2,],c)}),sum)
    c.weir <- sapply(lapply(weir.loci,function(x){sapply(x[3,],c)}),sum)
    
    weir.overall <- 1-(sum(c.weir)/(sum(c.weir+b.weir)))
    
    if (weir.overall>0)
    {p.weir.overall <- (1+sum(bootf[[1]] >= weir.overall))/(iteration+1)}
    if (weir.overall<0)
    {p.weir.overall <- (1+sum(bootf[[1]] <= weir.overall))/(iteration+1)}
    if (weir.overall==0)
    {p.weir.overall <- 1}
    
    names(empirical.fis) <- names(tab.pop)[3:length(names(tab.pop))][seq(1,length(empirical.fis)*2,2)]
    names(p.va) <- names(tab.pop)[3:length(names(tab.pop))][seq(1,length(empirical.fis)*2,2)]
    names(empirical.weir) <- names(tab.pop)[3:length(names(tab.pop))][seq(1,length(empirical.fis)*2,2)]
    names(p.weir) <- names(tab.pop)[3:length(names(tab.pop))][seq(1,length(empirical.fis)*2,2)]
    
    
 if(directory.name!="NA" & out.name!="NA") 
   {
     
  writeLines(paste("Loci names -- Note -- odd columns are set as loci names for further results","\n"), con=out.file)
  write.table(names(tab.pop)[3:length(names(tab.pop))], out.file, append=T, sep="\t")
  writeLines(paste("\n","\n",
  "Calculations made according to Nei 1983","\n", 
  "Fis values:","\n"), con=out.file)
  write.table(empirical.fis,out.file, append=T, quote=F, sep="\t",col.names=F)
  writeLines(paste("\n",
  "p values:","\n"), con=out.file)
  write.table(p.va,out.file, append=T, quote=F, sep="\t",col.names=F)
  writeLines(paste("\n",
  "Mean fis value:","\n"), con=out.file)
  write.table(mean(empirical.fis, na.rm=TRUE),out.file, append=T, quote=F, sep="\t",col.names=F, row.names=F)
  writeLines(paste("\n",
  "Mean p value:","\n"), con=out.file)
  write.table(mean(p.va, na.rm=TRUE),out.file, append=T, quote=F, sep="\t",col.names=F, row.names=F)
  writeLines(paste("\n","\n",
  "Calculations made according to Weir and Cockerham 1984","\n",
  "Fis values:","\n"), con=out.file)
  write.table(empirical.weir,out.file, append=T, quote=F, sep="\t",col.names=F)
  writeLines(paste("\n",
  "p values:","\n"), con=out.file)
  write.table(p.weir,out.file, append=T, quote=F, sep="\t",col.names=F)
  writeLines(paste("\n",
  "Mean fis value:","\n"),con=out.file)
  write.table(mean(empirical.weir, na.rm=TRUE),out.file, append=T, quote=F, sep="\t",row.names=F,col.names=F)
  writeLines(paste("\n",
  "Mean p value:","\n"),con=out.file)
  write.table(mean(p.weir, na.rm=TRUE),out.file, append=T, quote=F, sep="\t",row.names=F,col.names=F)
  writeLines(paste("\n",
  "Weigthed mean fis value:","\n"),con=out.file)
  write.table(weir.overall,out.file, append=T, quote=F, sep="\t",row.names=F,col.names=F)
  writeLines(paste("\n",
  "Weigthed mean p value:","\n"),con=out.file)
  write.table(p.weir.overall,out.file, append=T, quote=F, sep="\t",row.names=F,col.names=F)
  writeLines(paste("\n",
  "\n","\n","\n","References","\n",
  "Nei, M. (1972) Genetic distance between populations. American Naturalist, 106, 283-292.","\n",
  "Weir, B.S. and Cockerham, C. (1984) Estimating F-Statistics for the Analysis of Population Structure. Evolution, 38, 1358-1370.","\n",
  "Nei, M. and Chesser R.K. (1983) Estimation of fixation indices and gene diversities. Annals of Human Genetics, 47, 253-259.","\n","\n",sep=" "),
              con=out.file)
  close(out.file)

 }
    
    output.fis <- list(empirical.fis, empirical.weir, p.va, p.weir)
    names(output.fis) <- c("Empirical_Fis_Nei", "Empirical_Fis_Weir", "P_values", "P_values_Weir_and_Cockerham")
    return(output.fis)
    
   
}
  
