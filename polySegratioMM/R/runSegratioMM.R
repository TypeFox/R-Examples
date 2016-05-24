`runSegratioMM` <-
function(seg.ratios, model, priors=setPriors(model),
           inits=setInits(model,priors),
           jags.control= setControl(model, stem, burn.in=burn.in, sample=sample,
             thin=thin), burn.in = 2000, sample = 5000, thin = 1,
                            stem="test", fix.one = TRUE,
           print=TRUE, plots=TRUE,
           print.diagnostics=TRUE, plot.diagnostics=TRUE, run.diagnostics.later=FALSE)
##           print=c(TRUE,FALSE), plots=c(TRUE,FALSE),
##           print.diagnostics=c(TRUE,FALSE), plot.diagnostics=c(TRUE,FALSE))
{
  ## comprehensive wrapper for doing it all in one but with some defaults

  ## check ALL objects passed in args

  #classes
  if (class(seg.ratios) != "segRatio")
    stop("'seg.ratios' must be of class 'segRatio'")
  if (class(model) != "modelSegratioMM")
    stop("'model' must be of class 'modelSegratioMM'")
  if (class(priors) != "priorsSegratioMM")
    stop("'priors' must be of class 'priorsSegratioMM'")
  if (mode(inits) != "list")   #  could do better checking here but ...
    stop("'inits' must be a list")
  if (class(jags.control) != "jagsControl")
    stop("'jags.control' must be of class 'jagsControl'")

  ## check args as match.arg only for characters eg. error for
  ## print <- match.arg(print,c(TRUE,FALSE))
  if (!is.logical(print) | length(print)!=1)
    stop("Error: 'print' must be TRUE or FALSE")
  if (!is.logical(plots) | length(plots)!=1)
    stop("Error: 'plots' must be TRUE or FALSE")
  if (!is.logical(plot.diagnostics) | length(plot.diagnostics)!=1)
    stop("Error: 'plot.diagnostics' must be TRUE or FALSE")
  if (!is.logical(print.diagnostics) | length(print.diagnostics)!=1)
    stop("Error: 'print.diagnostics' must be TRUE or FALSE")


  ## write various files
  dumpData(seg.ratios, model, stem=stem, fix.one=fix.one)   # write data 
  dumpInits(inits, stem=stem)                    # write inits
  writeJagsFile(model, priors, stem=stem)        # jags (.bug) file
  writeControlFile(jags.control)                 #  .cmd file

  if (print) cat("Starting JAGS ...\n")
  run.jags <- runJags(jags.control)  ## just run it
  if (print) print(run.jags)

  read.jags <- readJags(run.jags)  ## read MCMC chain(s)

  ## diagnostics
  if (run.diagnostics.later) {
    ddd <- "Diagnostics to be run later"
  } else {
    ddd <- diagnosticsJagsMix(read.jags, diagnostics=print.diagnostics,
                              plots=plot.diagnostics, return.results=TRUE)
  }
  
  ## summarise parameters
  summary.params <- summary(read.jags, marker.index=NULL)
  if (print) {
    cat("Start:",summary.params$start,", End:", summary.params$end,
        ", Thin:",summary.params$thin,"\n", sep="")
    print(summary.params$statistics)
  }

  DIC <- calculateDIC(read.jags, model, priors, seg.ratios, print.DIC=print)
  ##if (print)
  ##  cat("\nDIC:",DIC,"\n")

  ## summarise/calculate doasages
  doses.jags <- dosagesJagsMix(read.jags, jags.control, seg.ratios) # probs!

  res <- list(seg.ratios=seg.ratios, model=model, priors=priors,
              inits=inits, jags.control=jags.control, stem=stem,
              fix.one = fix.one, 
              run.jags=run.jags, mcmc.mixture=read.jags, diagnostics=ddd,
              summary=summary.params,doses=doses.jags, DIC=DIC)
  
  class(res) <-  "runJagsWrapper"
  return(res)
}

