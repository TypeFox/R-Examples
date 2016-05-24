`setControl` <-
function(model, stem="test",
                       burn.in=2000, sample=5000, thin=1,
                       bugs.file=paste(stem,".bug",sep=""),
                       data.file=paste(stem,"-data.R",sep=""),
                       inits.file=paste(stem,"-inits.R",sep=""),
                       monitor.var=model$monitor.var, seed=1)
{

  if (class(model) != "modelSegratioMM")
    stop("'model' must be of class 'modelSegratioMM'")

  if ((length(thin) != 1) & (length(thin) != length(monitor.var)))
    stop("Error: 'thin' must be a scalar or vector of thinning for each variable set")

  if (length(thin) != 1) { # non-scalar
    if (length(names(thin)) == 0) { # set names
      names(thin) <- monitor.var
      cat("Warning: thin set to\n")
      print(thin)
    } else { # check that names OK
      if (! identical(names(thin), monitor.var)) {# not same in same order
        if (identical(sort(names(thin)), sort(monitor.var))) {# OK but reorder
          thin <- thin[monitor.var]
        } else { # wrong names
          stop("incompatible names for 'thin'")
        }
      }
      # names OK
    }
  }

  ## superseded now that JAGS Version 1.0 required
  ##if (.Platform$OS.type == "windows"){# to be fixed after JAGS 0.90 superseded
  ##  jc <- c(paste("seed", seed),
  ##          paste("model in \"",bugs.file,"\"",sep=""))
  ##} else {
  jc <- c(paste("model in \"",bugs.file,"\"",sep=""))
  ##}

  jc <- c(jc, paste("data in \"",data.file,"\"",sep=""),"compile",
          paste("inits in \"",inits.file,"\"",sep=""),
          "initialize",
          paste("update", burn.in))
  if (length(thin) == 1 & thin[1] == 1) {
    jc <- c(jc, paste("monitor set",monitor.var))
  } else {
    if (length(thin) == 1) {
      thin <- rep(thin,length(monitor.var))
      names(thin) <- monitor.var
    }
    jc <- c(jc, paste("monitor set ",monitor.var,", thin(",thin,")",sep=""))
  } 
  jc <- c(jc, paste("update", sample))

  ## superseded now that JAGS Version 1.0 required
  ##if (.Platform$OS.type == "windows"){# to be fixed after JAGS 0.90 superseded
  ##  jc <- c(jc,      paste("coda *, stem(\"",stem,"\")",sep=""),"exit")
  ##} else {  
  jc <- c(jc,      paste("coda *, stem(\"",stem,"CODA\")",sep=""),"exit")
  ##}

  ## take model out - waste of space??

  res <- list(jags.code=jc, model=model, stem=stem,
              burn.in=burn.in, sample=sample, thin=thin, bugs.file=bugs.file,
              data.file=data.file, inits.file=inits.file,
              monitor.var=monitor.var, call=match.call())
  oldClass(res) <- "jagsControl"
  return(res)

}

