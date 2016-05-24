### Propose new B conditioning on all orhter parameters.
###
### These functions are for all amino acids.

get.my.drawBConditionalAll <- function(type){
  if(!any(type[1] %in% .CF.CT$init.fit)){
    stop("fit is not found.")
  }
  ret <- eval(parse(text = paste("my.drawBConditionalAll.",
                                 type[1], sep = "")))
  assign("my.drawBConditionalAll", ret, envir = .cubfitsEnv)
  ret
} # End of get.my.drawBConditionalAll().


### Draw new B via independent chain. VGAM searching starts from current values.
my.drawBConditionalAll.current <- function(b.Curr, phi.Curr, y, n, reu13.df.obs,
    b.RInitList = NULL){
  ### Note that phi.new = phi.Curr is the E[Phi] rather than phi.Obs.
  b.Fit <- .cubfitsEnv$my.fitMultinomAll(reu13.df.obs, phi.Curr, y, n,
                                         phi.new = phi.Curr, coefstart = b.Curr)

  ### Based on the above new fits of parameters to draw new beta.
  ret <- lapply(1:length(reu13.df.obs),
           function(i.aa){ # i'th amino acid.
             my.drawBConditionalFit.ID_Norm(
               b.Fit[[i.aa]], b.Curr[[i.aa]], phi.Curr, y[[i.aa]],
               n[[i.aa]],
               reu13.df.aa = reu13.df.obs[[i.aa]])
           })

  ### Update beta's acceptance and adaptive.
  accept <- do.call("c", lapply(1:length(ret),
                           function(i.aa){
                             ret[[i.aa]]$accept
                           }))
  my.update.acceptance("b", accept)
  my.update.adaptive("b", accept)

  ret
} # End of my.drawBConditionalAll.current().

### Draw new B via independent chain. VGAM searching starts randomly.
my.drawBConditionalAll.random <- function(b.Curr, phi.Curr, y, n, reu13.df.obs,
    b.RInitList = NULL){
  ### Note that phi.new = phi.Curr is the E[Phi] rather than phi.Obs.
  b.Fit <- .cubfitsEnv$my.fitMultinomAll(reu13.df.obs, phi.Curr, y, n,
                                         phi.new = phi.Curr)

  ### Based on the above new fits of parameters to draw new beta.
  ret <- lapply(1:length(reu13.df.obs),
           function(i.aa){ # i'th amino acid.
             my.drawBConditionalFit.ID_Norm(
               b.Fit[[i.aa]], b.Curr[[i.aa]], phi.Curr, y[[i.aa]],
               n[[i.aa]],
               reu13.df.aa = reu13.df.obs[[i.aa]])
           })

  ### Update beta's acceptance and adaptive.
  accept <- do.call("c", lapply(1:length(ret),
                           function(i.aa){
                             ret[[i.aa]]$accept
                           }))
  my.update.acceptance("b", accept)
  my.update.adaptive("b", accept)

  ret
} # End of my.drawBConditionalAll.random().

### Draw new B via random walk. No VGAM searching.
### b.RInitList is fixed as the initial fits via VGAM without measurement errors.
### Only b.DrawScale and b.DrawScale.prev are changed for adaptive MCMC.
my.drawBConditionalAll.RW_Norm <- function(b.Curr, phi.Curr, y, n, reu13.df.obs,
    b.RInitList){
  ### No VGAM searching.
  b.Fit <- b.Curr

  ### Based on the above new fits of parameters to draw new beta.
  ret <- lapply(1:length(reu13.df.obs),
           function(i.aa){ # i'th amino acid.
             my.drawBConditionalFit.RW_Norm(
               b.Fit[[i.aa]], b.Curr[[i.aa]], phi.Curr, y[[i.aa]],
               n[[i.aa]],
               b.RInitList.aa = b.RInitList[[i.aa]],
               b.DrawScale.aa = .cubfitsEnv$all.DrawScale$b[i.aa],
               b.DrawScale.prev.aa = .cubfitsEnv$all.DrawScale$b.prev[i.aa],
               reu13.df.aa = reu13.df.obs[[i.aa]])
           })

  ### Update beta's acceptance and adaptive.
  accept <- do.call("c", lapply(1:length(ret),
                           function(i.aa){
                             ret[[i.aa]]$accept
                           }))
  my.update.acceptance("b", accept)
  my.update.adaptive("b", accept)

  ret
} # End of my.drawBConditionalAll.RW_Norm().

