################################
##
## Class: Data
##
################################

###produces difficulties in coercing...:

## Initialize methods
#setMethod("initialize", "Dataclass",
#   function(.Object, filename = NULL, Data=.Object@Data, runs=1, obsDim=1) {
#            .Object@filename <- filename
#            .Object@Data <- Data
###changed in 1.8:
#            .Object@version <- "1.8"
#            .Object@runs <- 1
#            .Object@obsDim <- ncol(as.matrix(Data))
#            .Object@samplesize <- nrow(as.matrix(Data))
### old:
#            .Object@runs <- nrow(as.matrix(Data))
#            .Object@samplesize <- ncol(as.matrix(Data))
#            .Object
#          })

### instead generating function ...

################################
##
## Class: Simulation
##
################################

### changed from version 1.8 on:
## ith observation in ith line of datamatrix/array
## jth item/dimension of each observation in jth column of datamatrix/array
## kth run/time of each observation in kth slide of datamatrix/array

## ++old
## +ith run in ith line of datamatrix
## +jth samples of each run in jth column of datamatrix


setMethod("initialize", "Simulation",
          function(.Object,
                   filename = NULL,
                   runs = 100,
                   samplesize = 10,
                   seed = setRNG(),
                   distribution = Norm()) {
            .Object@filename <- filename
            .Object@Data <- NULL
            .Object@runs <- runs
 ##new 03-10-06:
            .Object@version <- "1.9"
            .Object@obsDim <- dim(distribution)
 ##end (new)
            .Object@samplesize <- samplesize
            .Object@seed <- seed
            .Object@distribution <- distribution
            validObject(.Object)
            .Object
          })


################################
##
## Class: Contsimulation
##
################################


setMethod("initialize", "Contsimulation",
          function(.Object,
                   filename = NULL,
                   runs = 100,
                   samplesize = 10,
                   seed = setRNG(),
                   distribution.id = Norm(),
                   distribution.c = Norm(sd = 3),
                   rate = 0.1) {
            .Object@Data <- NULL
            .Object@Data.id <- NULL
            .Object@Data.c <- NULL
### new 031006:
            .Object@version <- "1.9"
            .Object@obsDim <- dim(distribution.id)
###
            .Object@filename <- filename
            .Object@runs <- runs
            .Object@samplesize <- samplesize
            .Object@seed <- seed
            .Object@distribution.id <- distribution.id
            .Object@distribution.c <- distribution.c
            .Object@rate <- rate

            validObject(.Object)
            .Object
          })
