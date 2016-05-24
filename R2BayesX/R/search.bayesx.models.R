search.bayesx.models <-
function(dir)
{
  if(!file.exists(dir))
    stop("directory is not existing!")
  files <- list.files(dir)
  if(length(files) < 1L || is.null(files))
    stop(paste("no files in directory:", dir))
  ok <- FALSE
  model.names <- NULL
  if(length(i <- grep("bayesx.log", files, fixed = TRUE))) {
    logf <- readLines(paste(dir, "/", files[i], sep = ""))
    outfiles <- grep("b.outfile = ", logf, value = TRUE)
    if(length(outfiles) > 0L)
      for(k in 1L:length(outfiles)) {
        of <- splitme(outfiles[k])
        go <- TRUE
        tk <- NULL
        for(j in length(of):1) {
          if(of[j] == "/" || of[j] == "\\")
            go <- FALSE
          if(go)
            tk <- c(tk,of[j])  
        }
        model.names <- c(model.names,resplit(tk[length(tk):1]))
      }
    collect <- NULL
    for(k in 1L:length(model.names))
      if(any(grepl(model.names[k], files)))
        collect <- c(collect, model.names[k])
    if(!is.null(collect)) {
      ok <- TRUE
      model.names <- collect
    }
  }

  ## specifiy possible search endings here
  search.endings <- c("_model_summary.tex", "_graphics.prg", "_stata.do",
    "_r.R", "_FixedEffects.res", "_FixedEffects1.res", "_FixedEffects2.res",
    "_predict.raw", "_modelfit.raw", "_LinearEffects.res", "_LinearEffects1.res",
    "_LinearEffects2.res", ".terms.info", "_predict.res")

  model.names2 <- list()
  n <- 0L
  for(i in 1L:length(search.endings)) {
    if(length(ii <- grep(search.endings[i],files))) {
      n <- n + 1L
      fs <- files[ii]
      k <- length(fs)
      model.names2[[n]] <- rep("NA", k)
      for(j in 1L:length(fs))
        model.names2[[n]][j] <- strsplit(fs[j], search.endings[i], "")[[1L]]
      ok <- TRUE
    }
  }
  model.names <- c(model.names, unlist(model.names2))
  if(!ok)
    stop(paste("no BayesX output files found in:", dir))
  else {
    model.names <- unique(unlist(model.names))
    for(k in 1L:length(model.names)) {
      model.names[k] <- strsplit(model.names[k], "_MAIN_REGRESSION", "")[[1L]][1L]
      model.names[k] <- strsplit(model.names[k], "_RANDOM_EFFECTS", "")[[1L]][1L]
      if(length(grep("hlevel", model.names[k]))) {
        split <- strsplit(model.names[k], "_hlevel")[[1L]]
        if(length(split) < 2L)
          split <- strsplit(model.names[k], ".hlevel")[[1L]]
        mn <- split[1L]
        hlevel <- split[2L]
        if(hlevel < 2L)
          model.names[k] <- paste(mn, "_hlevel", hlevel, "_MAIN_REGRESSION", sep = "")
        else
          model.names[k] <- paste(mn, "_hlevel", hlevel, "_RANDOM_EFFECTS", sep = "")
      }
    }
  }
		
  return(unique(unlist(model.names)))
}

