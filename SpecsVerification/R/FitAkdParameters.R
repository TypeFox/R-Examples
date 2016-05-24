FitAkdParameters <- function(ens, obs) {

  #       p(y|x) = 1 / K * sum {dnorm(y, z.i(x), s(x))}
  # where   s(x) = (4/3/K)^0.4 * (s1 + s2 * a^2 * var(x))
  # and   z.i(x) = r1 + r2 * mean(x) + a * x[i]

  # preprocess
  l <- Preprocess(ens=ens, obs=obs)
  ens <- l[["ens"]]
  obs <- as.vector(l[["obs"]])
  stopifnot(nrow(ens) == length(obs))

  # ensemble means
  m.x <- rowMeans(ens, na.rm=TRUE)
  # ensemble variance 
  v.x <-  apply(ens, 1, var, na.rm=TRUE)
  # squared silverman factor 
  K <- max(rowSums(!is.na(ens)))
  sf.2 <- (4 / 3 / K) ^ 0.4

  # initial guesses
  m1 <- lm(obs~m.x)
  m2 <- lm(resid(m1)^2~v.x)
  coef1 <- as.vector(coef(m1))
  coef2 <- as.vector(coef(m2))

  r1 <- coef1[1]
  s1 <- coef2[1] / sf.2
  s2 <- 1
  a <- ifelse(test = coef2[2] <= 0,
               yes = 0,
               no = sqrt(coef2[2] / (1 + sf.2)))
  r2 <- coef1[2] - a

  # the objective function: mean dressing crps with 
  # logarithmic barrier to enforce var.kernel > 0 
  f <- function(parms) {
    parms <- as.list(parms)
    d.ens <- DressEnsemble(ens, "akd", parms)
    sigma <- min(with(parms, sf.2 * (s1 + s2*a*a*v.x)))
    barrier <- ifelse(sigma <= 0, Inf, max(0, -log(sigma)))
    mean(DressCrps(d.ens, obs)) * (1 + 0.01 * barrier)
  }
  parms <- c(a=a, r1=r1, r2=r2, s1=s1, s2=s2)

  # if f cannot be evaluated at initial guesses, 
  # use standard silverman as initial guess
  if (!is.finite(f(parms))) { # catches Inf, NA, NaN
    parms <- c(a=1, r1=0, r2=0, s1=0, s2=1)
  }

  # optimize
  opt <- optim(par=parms, fn=f)
  
  # return
  opt[["par"]]
}

