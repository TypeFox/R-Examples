prod_array <- function(X,y){

# compute product of each matrix in X by the vector y

Z = 0
for(j in 1:dim(X)[2]) Z = Z+X[,j,]*y[j]
Z

}