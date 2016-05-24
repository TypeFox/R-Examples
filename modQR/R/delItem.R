delItem <- function(Item, Vec){
#delItem <- function(Item, Vec), output: Vec
#deleting an item (Item) from a vector (Vec) sorted in the ascending order

LVec <- length(Vec)
OutList <- findItem(Item, Vec, LVec); Pos <- OutList[[1]]
Vec <- Vec[-Pos]
return(Vec)
}