extract.val <- 
function(m, samples) {
    i <- cbind(match(samples[,1], rownames(m)), 
               match(samples[,2], colnames(m)))
    m[i]
}
