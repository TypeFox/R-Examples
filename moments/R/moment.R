"moment" <-
function(x, order = 1, central = FALSE, absolute = FALSE, na.rm = FALSE) 
{
    if (is.matrix(x)) 
        apply(x, 2, moment, order = order, central = central, absolute = absolute, na.rm = na.rm)
    else if (is.vector(x)) {
          if (na.rm) x = x[!is.na(x)] ;
          if (central) x = x - mean(x)
          if (absolute) x = abs(x)
          sum(x^order)/length(x)
                }
    else if (is.data.frame(x)) 
        sapply(x, moment, order = order, central = central, absolute = absolute, na.rm = na.rm)
    else moment(as.vector(x), order = order, central = central, absolute = absolute, na.rm = na.rm)
}
