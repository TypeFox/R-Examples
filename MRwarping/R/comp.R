comp <-
function(a,lambda,r1,r2,x)
{

if (lambda==0) {Wx=x}
else
{
r <- min(r1,r2)
Wx <- x
W <- Wx[Wx >(a-r1) & Wx < (a+r2)]
W1 <- W[W < a-3*sqrt(3)/8*lambda*r] 
W2 <- W[W >= a-3*sqrt(3)/8*lambda*r] 
z1 <- (W1-a)/r1
z2 <- (W2-a)/r2

lambda1 <- lambda*r/r1
lambda2 <- lambda*r/r2

c <- 8/(3*sqrt(3)*lambda1)

P <- -4/3-c*z1
Q <- 2/27-2/3*(1+c*z1)-(c^2)/8
R <- Q/2 - sqrt((Q/2)^2+(P/3)^3)
U <- (abs(R)/R)*(abs(R))^(1/3)
U[R==0] <-0
y <- 5/3 - U
y <- y + P/(3*U)
y[R==0] <- 5/3 - U
Ws <- sqrt(-2+2*y)
if (lambda > 0)
	{	u <-  1/2*( Ws - sqrt( -(-6+2*y +(-2*c)/Ws)  )  ) }
else	{	u <-  1/2*(- Ws + sqrt( -(-6+2*y -(-2*c)/Ws)  )  )}
W1 <- a+r1*( u+(1/c)*(1-u^2)^2 )


c2 <- 8/(3*sqrt(3)*lambda2)

P2 <- -1/3-1-c2*z2
Q2 <- 2/27-2/3*(1+c2*z2)-(c2^2)/8
R2 <- Q2/2 - sqrt((Q2/2)^2+(P2/3)^3)
U2 <- (abs(R2)/R2)*(abs(R2))^(1/3)
U2[R2==0] <-0
y2 <- 5/3 - U2
y2 <- y2 + P2/(3*U2)
y2[R2==0] <- 5/3 -U2
Ws2 <- sqrt(-2+2*y2)
if (lambda > 0)
	{u2 <-  1/2*(+ Ws2 - sqrt( -(-6+2*y2 +(-2*c2)/Ws2)  )  ) }
else  {u2 <-  1/2*(- Ws2 + sqrt( -(-6+2*y2 -(-2*c2)/Ws2)  )  ) }
W2 <- a+r2*( u2+(1/c2)*(1-u2^2)^2 )

Wx[Wx>(a-r1) & Wx<(a+r2)] <- c(W1,W2)
}
list(Wx=Wx)
}

