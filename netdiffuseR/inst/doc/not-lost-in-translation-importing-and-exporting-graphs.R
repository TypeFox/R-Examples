## ----Loading netdiffuseR, echo=FALSE-------------------------------------
library(netdiffuseR)
knitr::opts_chunk$set(comment = '#')

## ----Ex1: Datasets-------------------------------------------------------
# Loading the datasets
data("fakesurvey")
data("fakeEdgelist")

## ------------------------------------------------------------------------
head(fakesurvey[,c("id", "group")])
head(fakeEdgelist)

## ------------------------------------------------------------------------
# Coercing the edgelist to an adjacency matrix
adjmat  <- edgelist_to_adjmat(
  edgelist   = fakeEdgelist[,1:2], # Should be a two column matrix/data.frame
  w          = fakeEdgelist$value, # An optional vector with weights
  undirected = FALSE,              # In this case, the edgelist is directed
  t          = 5)                  # We use this option to make 5 replicas of it

## ------------------------------------------------------------------------
fakeEdgelist[11,,drop=FALSE]

## ------------------------------------------------------------------------
# Filling the empty data, and checking the outcome
fakeEdgelist[11,"value"] <- 1
fakeEdgelist[11,,drop=FALSE]

# Coercing the edgelist to an adjacency matrix (again)
adjmat  <- edgelist_to_adjmat(
  edgelist   = fakeEdgelist[,1:2], # Should be a two column matrix/data.frame
  w          = fakeEdgelist$value, # An optional vector with weights
  undirected = FALSE,              # In this case, the edgelist is directed
  keep.isolates = TRUE,            # NOTICE THIS NEW ARGUMENT!
  t          = 5)                  # We use this option to make 5 replicas of it

## ------------------------------------------------------------------------
adjmat[[1]]

## ------------------------------------------------------------------------
# Coercing the adjacency matrix and edgelist into a diffnet object
diffnet <- as_diffnet(
  graph               = adjmat,         # Passing a dynamic graph
  toa                 = fakesurvey$toa, # This is required
  vertex.static.attrs = fakesurvey      # Is is optional
  )

# Taking a look at the diffnet object
diffnet

## ------------------------------------------------------------------------
# Before
fakesurvey$id

# Changing the id
fakesurvey$id <- with(fakesurvey, group*100 + id)

# After
fakesurvey$id

## ------------------------------------------------------------------------
diffnet2 <- edgelist_to_diffnet(
  edgelist = fakeEdgelist[,1:2], # Passed to edgelist_to_adjmat
  w        = fakeEdgelist$value, # Passed to edgelist_to_adjmat
  dat      = fakesurvey,         # Data frame with -idvar- and -toavar-
  idvar    = "id",               # Name of the -idvar- in -dat-
  toavar   = "toa",              # Name of the -toavar- in -dat-
  keep.isolates = TRUE           # Passed to edgelist_to_adjmat   
)
diffnet2

## ------------------------------------------------------------------------
# Loading the data
data("fakesurvey")
fakesurvey

## ------------------------------------------------------------------------
fakesurvey[c(4,6),]

## ------------------------------------------------------------------------
# Coercing the survey data into a diffnet object
diffnet_w_unsurveyed <- survey_to_diffnet(
  dat      = fakesurvey,                # The dataset
  idvar    = "id",                      # Name of the idvar (must be integer)
  netvars  = c("net1", "net2", "net3"), # Vector of names of nomination vars
  toavar   = "toa",                     # Name of the time of adoption var
  groupvar = "group",                   # Name of the group var (OPTIONAL)
  no.unsurveyed = FALSE                 # KEEP OR NOT UNSURVEYED
)
diffnet_w_unsurveyed

# Retrieving nodes ids
nodes(diffnet_w_unsurveyed)

## ------------------------------------------------------------------------
# Coercing the survey data into a diffnet object
diffnet_wo_unsurveyed <- survey_to_diffnet(
  dat      = fakesurvey,                # The dataset
  idvar    = "id",                      # Name of the idvar (must be integer)
  netvars  = c("net1", "net2", "net3"), # Vector of names of nomination vars
  toavar   = "toa",                     # Name of the time of adoption var
  groupvar = "group"                    # Name of the group var (OPTIONAL)
)
diffnet_wo_unsurveyed

# Retrieving nodes ids
nodes(diffnet_wo_unsurveyed)

## ------------------------------------------------------------------------
difference <- diffnet_w_unsurveyed - diffnet_wo_unsurveyed
difference

## ------------------------------------------------------------------------
# Taking a look at the data
data("fakeDynEdgelist")
head(fakeDynEdgelist)

## ------------------------------------------------------------------------
data("fakesurveyDyn")
head(fakesurveyDyn)

## ------------------------------------------------------------------------
# Fixing ids
fakesurveyDyn$id <- with(fakesurveyDyn, group*100 + id)

# An individual who is alone
fakeDynEdgelist[11,"value"] <- 1

## ------------------------------------------------------------------------
diffnet <- edgelist_to_diffnet(
  edgelist = fakeDynEdgelist[,1:2], # As usual, a two column dataset
  w        = fakeDynEdgelist$value, # Here we are using weights
  t0       = fakeDynEdgelist$time,  # An integer vector with starting point of spell
  t1       = fakeDynEdgelist$time,  # An integer vector with the endpoint of spell
  dat      = fakesurveyDyn,         # Attributes dataset
  idvar    = "id",                  
  toavar   = "toa",
  timevar  = "time",
  keep.isolates = TRUE              # Keeping isolates (if there's any)
)

diffnet

