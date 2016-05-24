addItem <- function(Item, Vec){
#addItem <- function(Item, Vec), output: Vec
#adding an item (Item) to a vector (Vec) sorted in the ascending order

LVec <- length(Vec)
OutList <- findItem(Item, Vec, LVec); Pos <- OutList[[1]]; Status <- OutList[[2]]

if (Status == -1){
  Vec <- c(Item)
}else if (Status == 0){
  Vec <- c(Vec, Item)
}else if (Pos == 1){
  Vec <- c(Item, Vec)
}else{
  Vec <- c(Vec[1:(Pos - 1)], Item, Vec[Pos:LVec])
} #if
return(Vec)
}