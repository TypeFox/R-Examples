`skewness` <-
function (x, na.rm = FALSE) 
{
sd=function(x,...){mean((x-mean(x,...))^2)^.5}
    if (na.rm) 
        x <- x[!is.na(x)]
    sum((x - mean(x))^3)/(length(x) * sd(x,na.rm=na.rm)^3)
}

