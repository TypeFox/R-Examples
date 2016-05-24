# CHNOSZ/data/thermo.R
# create or restore the 'thermo' data object

# create the CHNOSZ environment if it does not exist
if(!"CHNOSZ" %in% search()) {
  attach(NULL, name="CHNOSZ")
  message("data(thermo): attached environment \"CHNOSZ\"")
}

if(exists("thermo")) {
  # scenario 1: thermo object exists ... restore it to default values
  # (<<- is superassignment to the CHNOSZ environment)
  thermo <<- list(
    # as.is: keep character values as character and not factor
    opt = as.list(read.csv("opt.csv", as.is=TRUE)),
    element = read.csv("element.csv", as.is=1:3),
    obigt = read.csv("OBIGT.csv", as.is=1:7),
    refs = read.csv("refs.csv", as.is=TRUE),
    buffers = read.csv("buffer.csv", as.is=1:3),
    protein = read.csv("protein.csv", as.is=1:4),
    groups = read.csv("groups.csv", row.names=1, check.names=FALSE),
    basis = NULL,
    species = NULL,
    Psat = NULL,
    # all values are restored, except opar (used for plot parameters in examples)
    opar = thermo$opar
  )
} else {
  # scenario 2: thermo object does not exist ... create it
  with(as.environment("CHNOSZ"), 
    thermo <- list(
      # as.is: keep character values as character and not factor
      opt = as.list(read.csv("opt.csv", as.is=TRUE)),
      element = read.csv("element.csv", as.is=1:3),
      obigt = read.csv("OBIGT.csv", as.is=1:7),
      refs = read.csv("refs.csv", as.is=TRUE),
      buffers = read.csv("buffer.csv", as.is=1:3),
      protein = read.csv("protein.csv", as.is=1:4),
      groups = read.csv("groups.csv", row.names=1, check.names=FALSE),
      basis = NULL,
      species = NULL,
      Psat = NULL,
      opar = NULL
    )
  )
}

# give a summary of some of the data
message(paste("thermo$obigt:",
  nrow(thermo$obigt[thermo$obigt$state=="aq",]),
  "aqueous,", nrow(thermo$obigt), "total species"))

# note if there are duplicated species
local({
  idup <- duplicated(paste(thermo$obigt$name, thermo$obigt$state))
  if(any(idup)) warning("thermo$obigt: duplicated species: ", 
    paste(thermo$obigt$name[idup], "(", thermo$obigt$state[idup], ")", sep="", collapse=" "))
})
