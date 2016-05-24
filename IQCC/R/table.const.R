table.const <- function(n)
{
    n <- 2:n
    u <- matrix(c(d2(n), d3(n), c4(n)), max(n) - 1, 3, byrow = FALSE)
    colnames(u) <- c("d2", "d3", "c4")
    rownames(u) <- n
    return(u)    
}