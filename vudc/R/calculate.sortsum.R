calculate.sortsum <-
function (vector.unsorted) 
{
    vector.sorted <- sort(vector.unsorted, TRUE)
    vector.sorted.sum <- c()
    vector.sorted.sum[1] <- 0
    for (i in 1:length(vector.sorted)) {
        vector.sorted.sum[i + 1] <- sum(vector.sorted[1:i])
    }
    return(vector.sorted.sum)
}
