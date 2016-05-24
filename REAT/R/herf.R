herf <-
function (x, norm = FALSE) {
n <- length (x)
a_i <- x/sum(x)
a_i_2 <- a_i^2
H <- sum(a_i_2)   

H.norm <- (H-1/n)/(1-1/n)
if (norm == FALSE) {   
return(H)   
} 
else 
{
return (H.norm)   
}
}
