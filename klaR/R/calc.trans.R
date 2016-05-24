calc.trans <- function(x)
{
# calculates the transition matrix given the time series of states
# ----------------------------------------------------------------
# x: vector of states
    x <- factor(x)
    if(length(x) > 1){
        tbl <- table(x[1:(length(x)-1)], x[2:length(x)])
        return(tbl / rowSums(tbl))
    }
    else return(table(x[1], x[1]) / sum(table(x[1])))
}
