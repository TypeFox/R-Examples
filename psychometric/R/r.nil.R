"r.nil" <-
function (r, n)
  {
t <- (r*sqrt(n-2))/sqrt(1-r^2)
df <- n-2
p <- pt(t, df)
d <- data.frame("H0:rNot0" = r, t = t, df=df, p=1-p)
return(d)
}

