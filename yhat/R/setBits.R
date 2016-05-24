setBits <-
function(col, effectBitMap) {
########################################################################################

##DESCRIPTION
##Creates the binary representation of n and stores it in the nth column of the matrix

##REQUIRED ARGUMENTS
##col				Column of matrix to represent in binary image
##effectBitMap	 	Matrix of mean combinations in binary form	

##Initialize variables
row<-1
val<-col
	 ##Create the binary representation of col and store it in its associated column
##One is stored in col 1; Two is stored in col 2; etc.
##While (val >= 1)
##  If the LSB of val is 1; increment the appropriate entry in combo matrix
##  Shift the LSB of val to the right
	while (val!=0){              
	if (odd(val)) {
		effectBitMap[row,col]=1 
	}
	val<-as.integer(val/2)
	row<-row+1
}

##Return matrix
return(effectBitMap)
}

