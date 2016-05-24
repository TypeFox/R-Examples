asym2sym <- function(foodweb, problem, name, single) {
   
  #Check that species names were provided (by checking for duplicates in the first row and column)
  if (any(duplicated(t(foodweb[1,-(1)]))) || any(duplicated(foodweb[-(1),1]))) {
    
    #if duplicates are found, indicate that this is a problematic matrix
    if (single==TRUE) {
      print("Your matrix contains duplicate species names. The first row and column of the matrix must be species names.")} else {
      write.table(cbind(name, "Duplicated species names"), file = problem, append=TRUE, quote=FALSE, sep=",", 
                col.names=FALSE, row.names=FALSE)}      
    } else {
  
  #Check that there are species in the rows that are not in the columns
  #If all species are represented, do nothing. 
  if (length(setdiff(foodweb[1,], foodweb[,1]))+length(setdiff(foodweb[,1], foodweb[1,]))==0){} else {
    
  #Otherwise, proceed to make the matrix symmetrical

  #Assign row and column names

foodweb[1,1] <- 9999
colnames(foodweb) <- foodweb[1,]
row.names(foodweb) <- foodweb[,1]

non.basal.sp <- ncol(foodweb)

#find species mentioned in rows that don't have a corresponding column
no.col <- as.vector(setdiff(row.names(foodweb), colnames(foodweb)))

#Add an all-zero column for each missing species. Give the column the species' name.

for (i in no.col) {
  foodweb <- cbind(foodweb, as.numeric(c(i,rep(0, times=nrow(foodweb)-1))))
  colnames(foodweb)[(which(no.col==i) + non.basal.sp)] <- i
}

#find species mentioned in columns that don't have a corresponding row
no.row <- as.vector(setdiff(colnames(foodweb), row.names(foodweb)))

#Add an all-zero row for each missing species. Give the row the species' name.
for (i in no.row) {
  foodweb <- rbind(foodweb, as.numeric(c(i,rep(0, times=ncol(foodweb)-1))))
  row.names(foodweb)[nrow(foodweb)] <- i
}


#Remove the species names and order the matrix for calculation purposes
foodweb <- foodweb[-(1),-(1)]
foodweb <- foodweb[,order(as.numeric(colnames(foodweb)))]
foodweb <<- foodweb[order(as.numeric(row.names(foodweb))),]

      }
      }
}