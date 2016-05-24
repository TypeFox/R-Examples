weir <- function(tab.pop,tab.freq,popsize)

{

# tab.freq = vector of allele p's with names of allele length e.g. 202
# tab.pop = Inputformat, col 3:4 used for calculation
# Weir preparations
  
  y <- data.frame(table(tab.pop[,3],tab.pop[,4])/sum(table(tab.pop[,3],tab.pop[,4])))
  heterozygotes <- sapply(names(tab.freq), function(x){sum(y[which(xor(y[,1]==x, y[,2]==x)),3])})
  
  c.weir <- heterozygotes/2
  
  b.weir<-vector(mode="numeric", length=length(tab.freq))
  fis.weir<-vector(mode="numeric", length=length(tab.freq))
  names(b.weir) <- names(tab.freq)
  names(fis.weir) <- names(tab.freq)
  


    for (i in 1:length(tab.freq))
  {
	if (heterozygotes[which(names(heterozygotes)==names(tab.freq[i]))]!=0)
	{b.weir[i] <- (popsize/(popsize-1))*((tab.freq[i]*(1-tab.freq[i]))-((2*popsize-1)/(4*popsize))*heterozygotes[which(names(heterozygotes)==names(tab.freq[i]))])
   fis.weir[i] <- 1-c.weir[which(names(c.weir)==names(fis.weir[i]))]/sum(c.weir[which(names(c.weir)==names(fis.weir[i]))]+b.weir[i])} 
  else { 
    b.weir[i]<-0
    fis.weir[i]<-0}
  
  }
   
  tab.total.weir <- rbind(tab.freq, b.weir, c.weir, fis.weir)
  rownames(tab.total.weir)<- c("Population.Frequency","b_weir","c_weir","Fis.Weir")
  return(tab.total.weir)
}
