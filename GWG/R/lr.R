lr <-
function(se,sp)
{
lrpos<- se/(1-sp)
lrneg <- (1-se)/sp
values <- list(lrpos=lrpos, lrneg=lrneg)
return(values)
}
