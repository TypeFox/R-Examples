plot.lqa <-
function (x, y, family, penalty.family, intercept = TRUE, standardize = TRUE, lambdaseq = NULL, offset.values = NULL, show.standardized = FALSE, add.MLE = TRUE, control = lqa.control (), plot.type = "l", ret.true = FALSE, really.plot = TRUE, ...)
{
### show.standardize ... logical, if 'TRUE' the standardized coefficients are plotted, otherwise the unstandardized coefficients are plotted
### add.MLE          ... logical, if 'TRUE' the unrestricted MLE is also plotted. Note this only works for 'n > p' settings. Otherwise this argument is set to 'FALSE' automatically.

  if (is.null (lambdaseq))
    lambdaseq <- exp (seq (-10, 6, length = 60))

  if (!is.function (penalty.family))
    stop ("'penalty.family' not recognized.")

  if (intercept && var (x[,1]) > control$var.eps)    # adds column of ones for the intercept
    x <- cbind (1, x)

  if (!intercept && var (x[,1]) <= control$var.eps)  # deletes column of ones 
    x <- x[,-1]

  if (!is.null (offset.values))
  {
    lambda.to.plot <- which (is.na (offset.values))

    if (length (lambda.to.plot) != 1)
      stop ("'offset.values' must include exactly one 'NA' element.")
  }

  pen.env <- environment (penalty.family)
  ll <- length (lambdaseq)
  nvars <- ncol (x)# + as.integer (intercept)
  beta.mat <- matrix (0, nrow = ll, ncol = nvars)

  if (is.null (colnames (x)))
  {
    colnames (x) <- if (intercept)
                      c ("intercept", paste ("x", 1 : (nvars - 1), sep = ""))
                    else
                      paste ("x", 1 : nvars, sep = "")
  }

  colnames (beta.mat) <- colnames (x)

  for (j in 1 : ll)
  {
    if (is.null (offset.values))
      current.lambda <- lambdaseq[j]
    else
    {
      current.lambda <- offset.values
      current.lambda[lambda.to.plot] <- lambdaseq[j]
    }

    penalty <- eval (penalty.family (current.lambda, ...), pen.env)
    initial.beta <- if (j > 1)
                      beta.mat[j,] 
                    else
                      NULL

    lqa.obj <- lqa.default (x = x, y = y, family = family, penalty = penalty, intercept = intercept, standardize = standardize, control = control, ...)

    beta.mat[j,] <- if (show.standardized)
                      lqa.obj$fit.obj$coefficients          # standardized estimates
                    else
                      lqa.obj$coefficients                  # unstandardized estimates
  }

  s1 <- if (intercept)
          apply (abs (beta.mat[,-1]), 1, sum)
        else 
          apply (abs (beta.mat), 1, sum)

  if (intercept)
    beta.mat <- beta.mat[,-1]


  s1 <- s1 / max (s1)

  main.text <- if (is.null (offset.values))
                 "lqa coefficient built-up"
               else
               {
                 h1.vec <- 1 : length (offset.values)
               }

  y.labtext <- ifelse (show.standardized, "standardized coefficients", "(unstandardized) coefficients")
  ccl <- as.character (current.lambda)

  if (really.plot)
  {

  matplot (s1, beta.mat, type = plot.type, main = bquote (paste ("lqa coefficient built-up (", .(penalty$penalty),")", sep = "")), xlab = "|beta|/max|beta|", ylab = y.labtext)

  if (add.MLE && nrow (x) > nvars)
  {
     glm.obj <- glm (y ~ 0 + x, family = family)

     if (intercept)
       points (x = rep (1.01, nvars - as.integer (intercept)), y = glm.obj$coef[-1])
     else
       points (x = rep (1.01, nvars - as.integer (intercept)), y = glm.obj$coef) 
  }

  axis (4, at = beta.mat[1,], labels = colnames (beta.mat), cex = 0.8, adj = 0, las = 1)
  }

  invisible ()
  ret.obj <- list (s1 = s1, beta.mat = beta.mat)
  
  if (ret.true)
    return (ret.obj)
}

