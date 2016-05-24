size.comparing_probs.mcnemar <-
  function(p1=NULL, p2=NULL, delta=NULL,
           alpha, beta,
           alternative=c("two.sided","one.sided")){
  if ((is.null(p1))&(is.null(p2)))
    {h <- p1 <- 0.5}
  if (is.null(p2)) 
    {p2 <- p1-delta
     h <- p1+p2-2*p1*p2}
  if(alternative == "two.sided") {q1 <- qnorm(1-alpha/2)}
  else {q1 <- qnorm(1-alpha)}
  q2 <- qnorm(1-beta)
  if ((!is.null(p1))&(!is.null(p2)))
    {h <- p1+p2-2*p1*p2
     delta <- p1-p2}
  n <- (q1*h+q2*sqrt(h^2-(3+h)*delta^2/4))/(h*delta^2)
  n <- ceiling(n)
  return(n)
}

