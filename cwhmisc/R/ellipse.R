ellipse1 <- function(k, a, b=a, ra=c(-1,361), phi=0 )
{
  ra <- rad(ra)
  alpha <- seq(0, by = (2 * pi)/k, length = k+1)   # define angles, overlapping
  c1 <- alpha >= min(ra)
  c2 <- alpha <= max(ra)
  ra <-  (c1 & c2)
  ran <- alpha[ra]
  Z <- cbind(a*cos(ran), b*sin(ran))
  R <- rotL(phi,1,2)[-3,-3] # 2*2
  res <- Z %*% R
  return( res )
}  # end of ellipse

conf.ellipse <- function(k, a, b, phi, df1, df2, level = 0.95)
{
  ellipse1(k, a, b,, phi=phi) * qf(level, df1, df2)
}
