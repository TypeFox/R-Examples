conditionTable2 <-
function (x, variables, condition) 
{
    out = aperm.default(conditionTable(x, variables, condition), 
        order(c(variables, condition)))
    tmp = patternRepeat(seq_len(prod(dim(x)[c(variables, condition)])), 
        c(variables, condition), dim(x))
    return(array(out[tmp], dim(x)))
}
