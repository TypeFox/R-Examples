##                    timetree                           ##
##      This code is part of the timetree package        ##
## copyright F.-S. Krah 2015 (last update: 2015-04-15)   ##
timetree <- function(taxa){
  url <- paste("http://www.timetree.org/index.php?taxon_a=",
               taxa[1],"&taxon_b=", taxa[2],
               "&submit=Search", sep="")
  parse <- readHTMLTable(url)
  if(length(parse)>0){
  names(parse) <- c("div", "ref")
  parse$div[,1] <- as.character(parse$div[,1])
  parse$div[,2] <- as.character(parse$div[,2])
  parse$div <- rbind(parse$div, names(parse$div))
  names(parse$div) <- c("estimate", "a")
  }
  if(length(parse)==0){parse <- "No molecular data available for this query"}
  return(parse)
}