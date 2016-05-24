




orderplot<- function (x, ...) 
UseMethod("orderplot")

orderplot.iBMA.surv<- function(x, ...)
{
    cs<- x 
    nvar<- cs$nVar
    plot(c(0,cs$currIter + 1), c(0, nvar + 1), type='n', main = 'orderplot for iBMA',
        ylab = 'variable', xlab = 'iteration', ...)
    for (i in 1:nvar)
    {
        if (i %in% cs$selected)
        {
            lines(c(cs$first.in.model[i], cs$currIter + 1 ), c(i,i), col = 'blue')
            points(cs$currIter+1, i, col = 'blue', pch= 19)                
        }    
        else 
        {
            if (is.na(cs$first.in.model[i]) | (i %in% cs$new.vars))
                points(0, i, pch = 21, col = 'darkgreen')
            else
            {
                if (is.na(cs$iter.dropped[i]))
                    lines(c(cs$first.in.model[i], cs$currIter+1), c(i,i), col = 'blue')
                else
                {
                    if (cs$first.in.model[i] == (cs$iter.dropped[i] ))
                        points(cs$first.in.model[i], i, pch = '-', col = 'black')
                    else
                        lines(c(cs$first.in.model[i], cs$iter.dropped[i]), c(i,i), col = 'black')
                }
            }
        }
    }
    invisible(x)  
}


orderplot.iBMA.intermediate.glm<- orderplot.iBMA.intermediate.bicreg<- orderplot.iBMA.intermediate.surv<- orderplot.iBMA.glm<- orderplot.iBMA.bicreg<- orderplot.iBMA.surv

