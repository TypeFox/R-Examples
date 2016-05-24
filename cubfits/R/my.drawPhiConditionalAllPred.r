### Function to run MH draw for phi given 
### count y.pred, number of codons n.pred, coefficients b.
###
### Performs vectorized draw.
### Assumes phi.Curr.pred, y.pred, & n.pred are vectors of length # of genes;
### b is 2-dim vector; sigmaWsq is scalar.
###
### Returns list with elements:
###   * phi.New  :   new value for phi resulting from MH draw
###   * accept:   boolean vector indicating if each proposal was accepted

my.drawPhiConditionalAllPred <- function(phi.Curr, y, n, b,
    p.Curr, reu13.df = NULL){
  ### Propose new phi.
  prop <- .cubfitsEnv$my.proposePhiAllPred(phi.Curr)

  ### Calculate acceptance prob.
  lpCurr <- .cubfitsEnv$my.logPosteriorAllPred(
              phi.Curr, y, n, b, p.Curr, reu13.df = reu13.df)
  lpProp <- .cubfitsEnv$my.logPosteriorAllPred(
              prop$phi.Prop, y, n, b, p.Curr, reu13.df = reu13.df)
  logAcceptProb <- lpProp - lpCurr - prop$lir

  ### Run MH acceptance rule.
  u <- runif(length(phi.Curr))
  accept <- u < exp(logAcceptProb)
  phi.New <- phi.Curr
  phi.New[accept] <- prop$phi.Prop[accept]

  ### Extra update and trace.
  my.update.acceptance("phi.pred", accept)
  my.update.adaptive("phi.pred", accept)

  ### Return.
  ret <- phi.New
  ret
} # End of my.drawPhiConditionalAllPred().

