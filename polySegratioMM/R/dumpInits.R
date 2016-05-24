`dumpInits` <-
function(inits, stem="test",
                      inits.file=paste(stem,"-inits.R",sep=""))
{

  if (mode(inits) != "list")
    stop("'inits' must be a list")


  ##  attach(inits)
  for (i in 1:length(inits)){
    assign(names(inits)[i],inits[[i]])
  }


  ## NB:  control="S_compatible", added for JAGS/R > 2.50 compatability
  dump(names(inits), control="S_compatible", file=inits.file)
  ##dump(c("mu","tau","theta","seed"), file=inits.file)

  ## superseded now that JAGS Version 1.0 required
  ##if (.Platform$OS.type == "windows"){# to be fixed after JAGS 0.90 superseded
  ##  fix.up <- readLines(inits.file)
  ##  fixed <- gsub("`","\"",fix.up)  # replace ` with "
  ##  writeLines(fixed, inits.file)
  ##}
  
}

