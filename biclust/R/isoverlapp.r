isoverlapp<-function(bicResult)

{

rows<-max(rowSums(bicResult@RowxNumber))

cols<-max(colSums(bicResult@NumberxCol))

if(rows>1 & cols>1)
  { cat("There are overlapping rows and cols in the Biclusterresult\n")}
else if(rows>1)
  { cat("There are overlapping rows in the Biclusterresult\n")}
else if(cols>1)
  { cat("There are overlapping cols in the Biclusterresult\n")}
else
  {cat("There are no overlapping rows or cols in the Biclusterresult\n")}
  
return(list(Overlapping= (rows>1 & cols>1), Max.bicluster.Rows=rows, Max.bicluster.cols=cols))


}


  