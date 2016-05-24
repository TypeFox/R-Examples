trans2seg<-function(vect){

# Patrick Giraudoux 16.1.2004
# Convert a transect coordinate file with landmarks 
# into a matrix with segment coordinates
# The argument passed is a matrix or data.frame of two columns
# each row is a transect interval; each column must start (first row)
# and end (last row) with a landmark ; intermediate landmarks must have
# numbers in the two columns of the row. Other rows must be NA values
# The returned value is a matrix of 4 columns to be passed to fonctions
# as segments()

vect<-na.omit(vect)
cbind(vect[1:length(vect[,1])-1,],vect[2:length(vect[,1]),])

}
