speciesFigure <- function(xml, species=NULL, type="chr", n=2:11, plot=TRUE){
  if(is.null(species)) stop("Please provide a species object (e.g. the dataset 'species', provided by this package)")
  xmlOne <- xml[[1]]
  for(i in 1:length(xml)){
    xmlOne <- rbind(xmlOne,xml[[i]])
  }
  if(type=="chr") xmlOne <- xmlOne[!is.na(xmlOne$hitChr),]
  xmlOne <- xmlOne[,1]
  
  res <- c()
  for(i in 1:length(species)){
    res[i] <- sum(grepl(species[i],xmlOne))
  }
  names(res) <- species
  
  if(plot){
    plotThis <- sort(res,decreasing=TRUE)[n]
    x <- barplot(plotThis, xaxt="n")
    yOffset <- 0
    text(cex=1, x=x-.25, y=0-yOffset, names(plotThis), xpd=TRUE, srt=45)
  }
  res
}
