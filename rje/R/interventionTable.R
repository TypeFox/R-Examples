interventionTable <-
function (x, variables, condition) 
{
    tmp = conditionTable2(x, variables, condition)
    out = x/tmp
    out
}
