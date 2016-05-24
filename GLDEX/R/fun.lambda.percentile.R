"fun.lambda.percentile" <-
function(x, u = 0.1)
{
A <- fun.percentile(x, 1 - u)
B <- fun.percentile(x, u)
p1 <- fun.percentile(x, 0.5)
p2 <- A - B
p3 <- (p1 - B)/(A - p1)
p4 <- (fun.percentile(x, 0.75) - fun.percentile(x, 0.25))/p2
return(list("p1"=p1, "p2"=p2, "p3"=p3, "p4"=p4))
}

