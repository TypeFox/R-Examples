rowProds <-
function(matr){

# Call: from din()
# Input: numeric matrix with positive entries
# Output: row products of input matrix

    exp( rowSums( log(matr + 10^(-300) ) ) )
    
}

