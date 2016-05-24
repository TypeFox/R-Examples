### Proposes value from multivariate normal with mean mu, and variance R'R.
### Returns proposal, log importance ratio, and z-values

### Transformation described in: Gentle, J.E. (2009). Computational Statistics. New York: Springer. pp. 315-316.
### http://download.springer.com/static/pdf/677/bok%253A978-0-387-98144-4.pdf?auth66=1426090126_653f471b5793fdd863be1412d1c04822&ext=.pdf

my.proposeB.ID_Norm <- function(mu.prev, mu, R){
  ### Draw from proposal.
  zProp <- rnorm(length(mu))
  prop <- mu + backsolve(R, zProp)

  ### Calculate importance ratio.
  zPrev <- R %*% (mu.prev - mu)
  lir <- -1/2 * (sum(zProp * zProp) - sum(zPrev * zPrev))

  ### Return.
  ret <- list(prop = as.numeric(prop), lir = lir)
  ret
} # End of my.proposeB.ID_Norm().


### In this case, mu.prev == mu since my.drawBConditionalAll.RW_Norm() set this.
my.proposeB.RW_Norm <- function(mu.prev, mu, R,
    b.DrawScale.aa, b.DrawScale.prev.aa){
  ### Draw from proposal.
  zProp <- rnorm(length(mu))
  prop <- mu + backsolve(R / b.DrawScale.aa, zProp)

  ### Check if drawing from the same scale.
  lir <- 0    # no jacobin since no transformation
  if(b.DrawScale.aa != b.DrawScale.prev.aa){
    ### Calculate importance ratio since random walk scale was changed.
    zPrev <- (R / b.DrawScale.prev.aa) %*% (mu.prev - prop)
    lir <- -1/2 * (sum(zProp * zProp) - sum(zPrev * zPrev))
  }

  ### Return.
  ret <- list(prop = as.numeric(prop), lir = lir)
  ret
} # End of my.proposeB.RW_Norm().
