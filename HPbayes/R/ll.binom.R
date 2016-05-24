ll.binom <-
function(x, n, p) 
{
    sum(dbinom(x = x, size = n, prob = p, log = TRUE))
}

