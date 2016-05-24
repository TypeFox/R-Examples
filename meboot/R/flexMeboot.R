
flexMeboot <- function(x, reps = 9, segment = 5, 
  forc = FALSE, myseq = seq(-1, 1, by = 1))
{
  #function to allow for a flexible trend up, flat or down
  #a+bt  replaced by a+Bt, where B=sample(myseq)*b
  if(segment < 5) 
    warning(paste("segment may be too small, segment =", segment))

  n <- length(x)
  newx  <- matrix(NA,nrow=n,ncol=reps) #initialize
  m <- floor(n/segment)

  for (i in 1:m)
  {
    ii <- (i-1) * segment + 1

    if (i < m)
      y <- x[ii:(ii+segment-1)]
    if (i == m)
      y <- x[ii:n]

    reg1 <- lm(y~seq(along=y)) #regress on tim= 1,2,3,4, 5
    b <- coef(reg1)[2]
    B <- b
    if (length(myseq) > 1) 
      B <- sample(myseq)[1]*b

    newy <- coef(reg1)[1]+ B * seq(along=y) + resid(reg1)
    men <- meboot(x=newy,reps=reps,force.clt=forc)$ensemble

    if (i < m) 
      newx [ii:(ii+segment-1),] <- men
    if (i == m) 
      newx [ii:n,] <- men
  }

  newx
}
