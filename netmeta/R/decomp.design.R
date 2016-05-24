decomp.design <- function(x, tau.preset=x$tau.preset){
  
  
  if (!inherits(x, "netmeta"))
    stop("Argument 'x' must be an object of class \"netmeta\"")
  
  
  tau.within <- tau.within(x)
  decomp.random <- decomp.tau(x, tau.preset=tau.within)
  Q.inc.random <- decomp.random$Q.decomp["Between designs",]
  Q.inc.design.random <- decomp.random$Q.inc.design
  residuals.detach.random <- decomp.random$residuals.detach
  
  
  if (length(tau.preset)==1){
    decomp.random.preset <- decomp.tau(x, tau.preset=tau.preset)
    Q.inc.random.preset <-  decomp.random.preset$Q.decomp["Between designs",]
    Q.inc.design.random.preset <- decomp.random.preset$Q.inc.design
    residuals.inc.detach.random.preset <- decomp.random.preset$residuals.inc.detach
  }
  else{
    tau.preset <- NULL
    Q.inc.random.preset <- NULL
    Q.inc.design.random.preset <- NULL
    residuals.inc.detach.random.preset <- NULL
  }
  
  
  dct <- decomp.tau(x)
  
  
  res <- list(
    Q.decomp=dct$Q.decomp,
    Q.het.design=dct$Q.het.design,
    Q.inc.detach=dct$Q.inc.detach,
    Q.inc.design=dct$Q.inc.design,
    ##
    Q.inc.random=data.frame(Q.inc.random, tau.within),
    Q.inc.random.preset=data.frame(Q.inc.random.preset, tau.preset),
    Q.inc.design.random.preset=Q.inc.design.random.preset,
    ##
    residuals.inc.detach=dct$residuals.inc.detach,
    residuals.inc.detach.random.preset=residuals.inc.detach.random.preset,
    ##
    tau.preset=tau.preset,
    ##
    call=match.call(),
    version=packageDescription("netmeta")$Version
    )
  ##
  class(res) <- "decomp.design"
  ##
  res
}
