cv <-
function (x, norm = FALSE)  
{ v <- sd(x)/mean(x);  
v.norm <- v/sqrt(length(x))

if (norm == FALSE) {
return (v)
}   

if (norm == TRUE) {
return(v.norm)
}
}
