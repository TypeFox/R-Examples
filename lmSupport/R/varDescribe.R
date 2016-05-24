varDescribe <- function(Data, Detail = 2, Digits = 2)
{
    t3 = psych::describe(Data)
    t3 = data.frame(t3)
    
    if (!is.null(Digits))  t3 = round(t3,Digits)
    
    t2 = t3[,c(1:5,8,9,11,12)]
    t1 = t3[,c(2:4,8,9)]
    t = switch(Detail, t1, t2, t3)
    return(t)
}