F.stat <- function(tab.pop, object=TRUE, iteration=1000, directory.name="NA", out.name="NA")
{
  # load data
  if (object==FALSE) {tab.pop <- input.txt(tab.pop, "pop")}
  
  number.loci <- (length(tab.pop)-2)/2
  tab.pop.pop <- split(tab.pop,tab.pop[,2])
  output.fis <- vector("list", length(tab.pop.pop))
  
  for (k in 1:length(tab.pop.pop))
              {
                  message(paste("---","Fis calculations are performed for Population",names(tab.pop.pop)[k],"----",Sys.time()),"\n","\n")
                  output.fis[[k]] <- Fis.calc(tab.pop.pop[[k]], iteration, number.loci, object, directory.name, out.name)
                  names(output.fis)[k] <- as.character(tab.pop.pop[[k]][1,2])
              }

  if(directory.name!="NA" & out.name!="NA") {message(paste("--- Statistics on Fis calculations are passed to", directory.name, " ---\n\n"))}
  
  return(output.fis)
}
