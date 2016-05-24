get.Amat <-
function (initial.beta = NULL, penalty = NULL, intercept = TRUE, c1 = lqa.control()$c1, x = NULL, ...)
{
  env <- NULL

  if (intercept)
     x <- x[,-1]

  if (is.null (initial.beta))
     stop ("get.Amat: 'initial.beta' is missing.")

  if (is.null (penalty))
     stop ("get.Amat: 'penalty' is missing.")

  nreg <- length (initial.beta) - as.integer (intercept)  

  coefm0 <- if (intercept) 
              initial.beta[-1]   
            else 
              initial.beta


#### Computation of A.lambda:

  if (is.null (penalty$getpenmat))
  {
   a.coefs <- if (is.null (penalty$a.coefs))
            diag (nreg)    # just built for the p regressors! (Intercept not accounted for!!!)
          else
            penalty$a.coefs (coefm0, env = env, x = x, ...)

   A.lambda <- matrix (0, nrow = nreg, ncol = nreg)
   xim0 <- drop (t (a.coefs) %*% coefm0)
   A.lambda2 <- drop (penalty$first.deriv (coefm0, env = env, x = x, ...) * as.integer (xim0 != 0) / sqrt (c1 + xim0^2))

   Jseq <- 1 : ncol (a.coefs)
   l1 <- sapply (Jseq, function (Jseq) {which (a.coefs[,Jseq] != 0)})   # extracts the positions of elements != 0 in the columns of 'a.coefs'
   just.one <- which (sapply (Jseq, function (Jseq) {length (l1[[Jseq]]) == 1}) == TRUE)   # extracts the column indices of 'a.coefs' with just one element != 0
   less.sparse <- setdiff (Jseq, just.one)   # extracts the column indices of 'a.coefs' with more than one element != 0

   if (length (just.one) > 0)  # <FIXME: Does not work if there are more than 'p' columns of 'a.coefs' with just one element != 0!>
   {
      jo2 <- rep (0, nreg)
      jo1 <- A.lambda2[just.one]    ### extracts the elements corresponding with just.one
      sort1 <- sapply (just.one, function (just.one) {l1[[just.one]]})  ### bestimmt die Reihenfolge    
      jo2[sort1] <- jo1
      A.lambda <- A.lambda + diag (jo2)
   }

   for (i in less.sparse)
   {
      ci <- l1[[i]]
      a.ci <- a.coefs[ci, i]
      A.lambda[ci,ci] <- A.lambda[ci,ci] + A.lambda2[i] * outer (a.ci, a.ci)   
   }
  }
  else
    A.lambda <- penalty$getpenmat (beta = coefm0, env = env, x = x, ...)      # if the 'x' argument is not needed (e.g. for penalreg penalty) then it should be accounted for automatically... (hopefully ;-)

   if (intercept)
   {
     A.lambda <- cbind (0, A.lambda)
     A.lambda <- rbind (0, A.lambda)
   }

   A.lambda
}

