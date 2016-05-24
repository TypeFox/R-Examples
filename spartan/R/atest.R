atest <-
function(x,y)
{
# Calculate the RankSum required to perform the A-Test. No need to do the Wilcoxon-Mann Test
ranksum <- sum(rank(c(x,y))[1:length(x)])

# Now work out the A-test Value
d <- (ranksum / length(x) - (length(x)+1)/2)/length(y)

}

