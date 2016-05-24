findItem <- function(Item, Vec, MaxInd){
#findItem <- function(Item, Vec, MaxInd), output: list(Ind, Status)
#(binary) searching for an item (Item) in the range Vec[1] to Vec[MaxInd],
# determined by a vector (Vec) sorted in the ascending order
#
#Ind         ... the index of Item
#Status      ... the output indicator
#Status = -1 ... Vec is empty
#          0 ... no good index exists (i.e. Item > Vec[MaxInd])
#          1 ... Item <= Vec(Ind) (and Item >= Vec[Ind-1] if Ind > 1)
#
#Vec    ... a vector sorted in the ascending order
#Item   ... the value/item being found
Status <- 0
Ind <- 0

if (1 > MaxInd){
    Status <- -1
    return(list(Ind, Status))
}else{
    L <- 1
    U <- MaxInd
    while (L < U){
        P <- floor((L+U)/2)
        if (Item <= Vec[P]){
             U <- P
        }else{
            L <- P+1
        } #if
    } #while
    #L == U

    if (Item <= Vec[L]){
        Status <- 1
        Ind <- L
    } #if
} #if
return(list(Ind, Status))
}