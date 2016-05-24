### Get the specific function accroding to the options.
# get.my.drawBConditionalFit <- function(type){
#   if(!any(type[1] %in% .CF.CT$type.B)){
#     stop("type is not found.")
#   }
#   ret <- eval(parse(text = paste("my.drawBConditionalFit.",
#                                  type[1], sep = "")))
#   assign("my.drawBConditionalFit", ret, envir = .cubfitsEnv)
#   ret
# } # End of get.my.drawBConditionalFit().


### Draw new beta from independent chain.
my.drawBConditionalFit.ID_Norm <- function(bFitaa, baa, phi, yaa, naa,
    reu13.df.aa = NULL){
  ### Propose new beta.
  bHat <- bFitaa$coefficients
  R <- bFitaa$R
  proplist <- my.proposeB.ID_Norm(baa, bHat, R)

  ### M-H step.
  ret <- my.drawBConditionalFit.MH(proplist, baa, phi, yaa, naa,
                                   reu13.df.aa = reu13.df.aa)
  ret
} # End of my.drawBConditionalFit.ID_Norm().

### Draw new beta from random walk.
my.drawBConditionalFit.RW_Norm <- function(bFitaa, baa, phi, yaa, naa,
    b.RInitList.aa, b.DrawScale.aa = 1, b.DrawScale.prev.aa = 1,
    reu13.df.aa = NULL){
  ### Propose new beta.
  bHat <- bFitaa
  R <- b.RInitList.aa
  proplist <- my.proposeB.RW_Norm(baa, bHat, R,
                                  b.DrawScale.aa, b.DrawScale.prev.aa)

  ### M-H step.
  ret <- my.drawBConditionalFit.MH(proplist, baa, phi, yaa, naa,
                                   reu13.df.aa = reu13.df.aa)
  ret
} # End of my.drawBConditionalFit.Norm().


### Utility function commonly for all my.drawBConditionalFit.*().
my.drawBConditionalFit.MH <- function(proplist, baa, phi, yaa, naa,
    reu13.df.aa = NULL)
  {
  
  #browser()
  baaProp <- proplist$prop
  lir <- proplist$lir

  
  ### Calculate posterior ratio.
  lpr <- .cubfitsEnv$my.logdmultinomCodOne(baaProp, phi, yaa, naa,
                                           reu13.df.aa = reu13.df.aa) -
         .cubfitsEnv$my.logdmultinomCodOne(baa, phi, yaa, naa,
                                           reu13.df.aa = reu13.df.aa)
  
  ### Calculate prior ratio 
  lprior <- my.drawBPrior(baa, baaProp)
  
  ### log Acceptance probability.
  logAcceptProb <- lpr - lir - lprior
  ### Error handling -- interpreting NaN etc. as ~= 0.
  if(!is.finite(logAcceptProb)){
    warning("log acceptance probability not finite in b draw")
    logAcceptProb <- -Inf
  }
    
  ### Run MH acceptance rule.
  p <- -rexp(1)
  if(p < logAcceptProb){
    bNew <- baaProp
    accept <- 1
  } else{
    bNew <- baa
    accept <- 0
  }

  ### Return.
  ret <- list(bNew = bNew, accept = accept)

  ret
} # End of my.drawBConditionalFit.MH().


## calculates log( (b/b')^-1 )
my.drawBPrior <- function(baa, baaProp)
{
  ncoef <- .cubfitsEnv$my.ncoef #get.my.ncoef(.cubfitsEnv$model, assign.Env = FALSE)
  # on log scale
  priorProp <- 0 # default is uniform 
  
  dmindex <- 1:(length(baa)/ncoef)
  baa <- baa[dmindex]
  baaProp <- baaProp[dmindex]
  if(.CF.CT$prior.dist.M[1] == "normal")
  {
    priorProp <- sum( dnorm(baa, .CF.PARAM$prior.M.a, .CF.PARAM$prior.M.b, log=T)
                  - dnorm(baaProp, .CF.PARAM$prior.M.a, .CF.PARAM$prior.M.b, log=T) )    
  }
  return(priorProp) 
}
