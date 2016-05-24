indexDummy <- function(kj, f, p){

# compute the column numbers for the dummy variables for each
# factor (for each factor, the output list contains a vector
# consisting of the column indices)

c <- p - f
JJs <- list()
JJs[[1]] <- c + 1:(kj[1] - 1)

if (f > 1){
for (j in 2:f){
        JJs[[j]] <- max(JJs[[j - 1]]) + 1:(kj[j] - 1)
}}

return(JJs)
}

