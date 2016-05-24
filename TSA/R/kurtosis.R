`kurtosis` <-
function (x, na.rm = FALSE) 
{
var=function(x,...){mean((x-mean(x,...))^2)}
    if (na.rm) 
        x <- x[!is.na(x)]
    sum((x - mean(x))^4)/(length(x) * var(x,na.rm=na.rm)^2) - 3
}

