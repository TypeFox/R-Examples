# CHNOSZ/examples.R
# run examples from the help files, 
# and a function containing extra examples

examples <- function(do.png=FALSE) {
  # run all the examples in CHNOSZ documentation
  .ptime <- proc.time()
  topics <- c("CHNOSZ-package", "thermo", "sideeffects", "examples",
    "util.args", "util.array", "util.blast", "util.character", 
    "util.data", "util.expression", "util.fasta", "util.formula", "util.matrix", 
    "util.misc", "util.program",
    "util.seq", "util.units", "taxonomy", "info", "protein.info", "hkf", "water", "subcrt",
    "makeup", "basis", "swap.basis", "species", "affinity", "util.affinity", "equil.boltzmann", 
    "diagram", "buffer", "iprotein", "protein", "ionize.aa", "more.aa", "read.expr",
    "objective", "revisit", "transfer", "anim", "EOSregress", "wjd")
  plot.it <- FALSE
  if(is.character(do.png))
    png(paste(do.png,"%d.png",sep=""),width=500,height=500,pointsize=12)
  else if(do.png) plot.it <- TRUE
  for(i in 1:length(topics)) {
    if(plot.it) png(paste(topics[i],"%d.png",sep=""),width=500,height=500,pointsize=12)
    myargs <- list(topic=topics[i],ask=FALSE)
    do.call(example,myargs)
    if(plot.it) dev.off()
  }
  if(is.character(do.png)) dev.off()
  # at the end we attempt to restore the old par() (active as of the first call of thermo.plot.new)
  par(get("thermo")$opar)
  cat("Time elapsed: ", proc.time() - .ptime, "\n")
}

demos <- function(which=c("sources", "NaCl", "density", 
  "nucleobase", "ORP", "revisit", "findit",
  "ionize", "buffer", "yeastgfp", "mosaic",
  "copper", "solubility", "wjd"), do.png=FALSE) {
  # run one or more demos from CHNOSZ with ask=FALSE, and return the value of the last one
  for(i in 1:length(which)) {
    # say something so the user sees where we are
    msgout("------------\n")
    msgout(paste("demos: running '", which[i], "'\n", sep=""))
    if(do.png) png(paste(which[i],"%d.png",sep=""),width=500,height=500,pointsize=12)
    out <- demo(which[i], package="CHNOSZ", character.only=TRUE, echo=FALSE, ask=FALSE)
    if(do.png) dev.off()
  }
  return(invisible(out))
}
