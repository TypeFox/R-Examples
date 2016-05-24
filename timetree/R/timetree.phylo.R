#                    timetree.phylo                      ##
##      This code is part of the timetree package        ##
## copyright F.-S. Krah 2015 (last update: 2015-04-15)   ##

timetree.phylo <- function(phy, node.time=c("expert", "mean", "median")){
  nt <- length(phy$tip.label)
  l <- vector("list", nt)
  for (i in 1:(nt-1)){
    cat("Node ", nt+i, "\n")
    d <- unlist(Descendants(phy, nt+i, type="tips"))
    taxa <- phy$tip.label[d][c(1,length(d))]
    div <- timetree(taxa)
    if(length(grep("No molecular data available for this query", div))==0)
    {
      if(length(grep("Expert", div$div[,1]))>0)
      {di <- div$div[grep("Expert", div$div[,1]),][2]
       di <- gsub("\\sMya\\s\\(TimeTree Book\\)","",di)
       expert <- c(nt+i, di)} 
      else{expert <- c(nt+i, "NA")}
      
      if(length(grep("Mean", div$div[,1]))>0)
      {
        di <- div$div[grep("Mean", div$div[,1]),][2]
        di <- gsub("\\sMya","",di)
        mean <- c(di)}
      else{mean <- c("NA")}
      
      if(length(grep("Median", div$div[,1]))>0)
      {
        di <- div$div[grep("Median", div$div[,1]),][2]
        di <- gsub("\\sMya","",di)
        median <- c(di)}
      else{ median <- c("NA") }  
      l[[i]] <- c(expert, mean, median)
    }
    if(length(grep("No molecular data available for this query", div))>0)
    {l[[i]] <- c(nt+i, rep("NA", 3))}
  }
  l <- do.call(rbind,l)
  l <- as.data.frame(l)
  names(l) <- c("node", "expert","mean", "median")
  phy$node.label <- eval(parse(text=paste("l$",node.time,sep="")))
  return(list(ages=l, phy=phy))
}