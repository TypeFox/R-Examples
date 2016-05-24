"normt3ip" <-
function (mu, sigma) 
{
a <- 1/(sigma^2)
b <- sqrt(2)/sigma
d <- a/b
p <- mu*b

T1 <- ((1-a)*cos(mu*a) + p*b*sin(mu*a)/2)* ic1(p,d)

T2 <- ((1-a)*sin(mu*a) - p*b*cos(mu*a)/2)*is1(p,d)

T3 <- b*exp( - d^2)/2

total <- b*exp(a/2)*(T1+T2+T3)/pi
total
}
