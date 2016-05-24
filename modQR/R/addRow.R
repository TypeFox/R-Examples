addRow <- function(Row, Mat, NRow, Pos){
#addRow <- function(Row, Mat, NRow, Pos), output: Mat
#adding a row (Row) to a prespecified position (Pos)
# in a matrix (Mat) with the first NRow recorded rows sorted
#  lexicographically (and with some other rows redundant)
#
#Mat   ... a matrix sorted in the lexicographical order
#Row   ... the row being added
#NRow  ... the number of relevant rows of Mat
#Pos   ... the position where to insert Row
#Pos = 0 means the insertion just after the first NRow rows

if (NRow == 0){
  Mat[1,] <- Row
}else if (Pos == 0){
  Mat[NRow+1,] <- Row
}else{
  Mat[(Pos + 1):(NRow + 1),] <- Mat[Pos:NRow,]
  Mat[Pos,] <- Row
} #if
return(Mat)
}
