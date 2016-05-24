test.rms <- function(observed,expected){
o <- length(observed)
e <- length(expected)
if(o != e){
return(cat("ERROR: Observed and Expected vectors must be same size"))
}

sum <- 0
for(s in 1:o){
sum <- sum + (expected[s]-observed[s])^2
}
rms <- sqrt((1/o)*sum)

return(rms)
}
