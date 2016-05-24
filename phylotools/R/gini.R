#### Function gini as part of R package phylotools
#### By Jinlong Zhang  <Jinlongzhang01@gmail.com>
#### Institute of Botany, the Chinese Academy of Sciences, Beijing ,China
#### Nov- 01-2010


gini <- function (x) 
{
    n <- length(x)
    x <- sort(x)
    res <- 2 * sum(x * 1:n)/(n * sum(x)) - 1 - (1/n)
    return(res)
}

