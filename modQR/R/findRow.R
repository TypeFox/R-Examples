findRow <- function(Row, Mat, NRow, NCol){
#findRow <- function(Row, Mat, NRow, NCol), output: list(Ind, Status)
#(binary) searching for a row (Row) in a matrix (Mat) with
# NCol columns and with its first NRow relevant normalized rows
#  sorted in the lexicographical order
#Ind         ... the row index searched for
#Ind = 0     ... the row index not found
#Status      ... the output indicator
#Status = -1 ... Mat[c(1:NRow),] is empty
#          0 ... no good index exists, i.e. Row > Mat[NRow,]
#          1 ... Row = Mat[Ind,]
#          2 ... Row < Mat[Ind,]

#Mat   ... a matrix with the first NRow normalized rows sorted lexicographically
#Row   ... the row looked for
#NRow    = the number of relevant first sorted rows of Mat
#NCol    = dim(Mat)[2]
Eps <- 2e-11

Status <- 0
Ind <- 0

if (NRow == 0){
  Status <- -1
  return(list(Ind, Status))
} #if

Min     <- 1
Max     <- NRow
for (i in 1:NCol){
  #modified findItem starts here
  StatAux <- 0
  IndAux <- 0

    L <- Min
    U <- Max
    while (L < U){
        P <- floor((L+U)/2)
        if (Row[i] <= Mat[P,i]){
            U <- P
        }else{
            L <- P+1
        } #if
    } #while
    #L == U

    if ((Row[i] > (Mat[L,i]-Eps)) && (Row[i] < (Mat[L,i]+Eps))){
        StatAux <- 1
        IndAux <- L
    }else if (Row[i] >= (Mat[L,i]+Eps)){
        if (Max == NRow){
            Status <- 0
            return(list(Ind, Status))
        }else{
            Status <- 2
            Ind <- Max+1
            return(list(Ind, Status))
        }
    }else if ((L > 1) && (Row[i] < (Mat[L-1,i] + Eps))){
        StatAux <- 1
        IndAux <- L-1
    }else if (Row[i] <= (Mat[L,i]-Eps)){
        Status <- 2
        Ind <- L
        return(list(Ind, Status))
    } #if

    if ((i < NCol) && (StatAux == 1)){
        MinAux <- IndAux
        while ((MinAux >= Min) && (Mat[MinAux,i] > (Row[i]-Eps))){
            MinAux <- MinAux-1
        } #while
        Min <- MinAux+1
        MaxAux <- IndAux
        while ((MaxAux <= Max) && (Mat[MaxAux,i] < (Row[i]+Eps))){
            MaxAux <- MaxAux+1
        } #while
        Max <- MaxAux-1
    } #if
} #for i
if (StatAux == 1){
    Status <- 1
    Ind <- IndAux
} #if
return(list(Ind, Status))
}