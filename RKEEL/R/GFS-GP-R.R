#Class implementing a Regression Algorithm
  #Implements the GFS-GP-R KEEL regression algorithm

GFS_GP_R <- function(train, test, numLabels=3, numRules=8, popSize=30, numisland=2, steady=1, numIter=100, tourSize=4, mutProb=0.01, aplMut=0.1, probMigra=0.001, probOptimLocal=0.00, numOptimLocal=0, idOptimLocal=0, nichinggap=0, maxindniche=8, probintraniche=0.75, probcrossga=0.5, probmutaga=0.5, lenchaingap=10, maxtreeheight=8, seed=-1){
  alg <- RKEEL::R6_GFS_GP_R$new()
  alg$setParameters(train, test, numLabels, numRules, popSize, numisland, steady, numIter, tourSize, mutProb, aplMut, probMigra, probOptimLocal, numOptimLocal, idOptimLocal, nichinggap, maxindniche, probintraniche, probcrossga, probmutaga, lenchaingap, maxtreeheight, seed)
  return (alg)
}

R6_GFS_GP_R <- R6::R6Class("R6_GFS_GP_R",

  inherit = RegressionAlgorithm,

  public = list(

    #Public properties

    #Number of labels
    numLabels = 3,

    #Number of rules
    numRules = 8,

    #Population size
    popSize = 30,

    #numisland
    numisland = 2,

    #steady
    steady = 1,

    #Number of iterations
    numIter = 100,

    #Tournament size
    tourSize = 4,

    #Mutation probabilitu
    mutProb = 0.01,

    #Amplmut
    aplMut = 0.1,

    #Prob migra
    probMigra = 0.001,

    #Prob local optim
    probOptimLocal = 0.00,

    #Number of optim local
    numOptimLocal = 0,

    #id optim local
    idOptimLocal = 0,

    #nichinggap
    nichinggap = 0,

    #maxindniche
    maxindniche = 8,

    #probintraniche
    probintraniche = 0.75,

    #probcrossga
    probcrossga = 0.5,

    #probmutaga
    probmutaga = 0.5,

    #lenchaingap
    lenchaingap = 10,

    #max tree height
    maxtreeheight = 8,

    #seed
    seed = -1,


    #Public functions

    #Initialize function
    setParameters = function(train, test, numLabels=3, numRules=8, popSize=30,
                          numisland=2, steady=1, numIter=100, tourSize=4,
                          mutProb=0.01, aplMut=0.1, probMigra=0.001,
                          probOptimLocal=0.00, numOptimLocal=0, idOptimLocal=0,
                          nichinggap=0, maxindniche=8, probintraniche=0.75,
                          probcrossga=0.5, probmutaga=0.5, lenchaingap=10,
                          maxtreeheight=8, seed=-1){

      super$setParameters(train, test)

      self$numLabels <- numLabels
      self$numRules <- numRules
      self$popSize <- popSize
      self$numisland <- numisland
      self$steady <- steady
      self$numIter <- numIter
      self$tourSize <- tourSize
      self$mutProb <- mutProb
      self$aplMut <- aplMut
      self$probMigra <- probMigra
      self$probOptimLocal <- probOptimLocal
      self$numOptimLocal <- numOptimLocal
      self$idOptimLocal <- idOptimLocal
      self$nichinggap <- nichinggap
      self$maxindniche <- maxindniche
      self$probintraniche <- probintraniche
      self$probcrossga <- probcrossga
      self$probmutaga <- probmutaga
      self$lenchaingap <- lenchaingap
      self$maxtreeheight <- maxtreeheight


      if(seed == -1) {
        self$seed <- sample(1:1000000, 1)
      }
      else {
        self$seed <- seed
      }

    }

  ),

  private = list(

    #Private properties

    #jar Filename
    jarName = "crispSymRegGAP.jar",

    #algorithm name
    algorithmName = "GFS-GP-R",

    #String with algorithm name
    algorithmString = "Fuzzy Rule Learning, Grammar-based GAP Algorithm",


    #Private functions

    #Get the text with the parameters for the config file
    getParametersText = function(){

      text <- ""
      text <- paste0(text, "seed = ", self$seed, "\n")
      text <- paste0(text, "subAlgorithm = ModelFuzzyGAP", "\n")
      text <- paste0(text, "dataformat = keel", "\n")
      text <- paste0(text, "numlabels = ", self$numLabels, "\n")
      text <- paste0(text, "numrules = ", self$numRules, "\n")
      text <- paste0(text, "outlabel = MFGAP", "\n")
      text <- paste0(text, "popsize = ", self$popSize, "\n")
      text <- paste0(text, "numisland = ", self$numisland, "\n")
      text <- paste0(text, "steady = ", self$steady, "\n")
      text <- paste0(text, "numitera = ", self$numIter, "\n")
      text <- paste0(text, "toursize = ", self$tourSize, "\n")
      text <- paste0(text, "probmuta = ", self$mutProb, "\n")
      text <- paste0(text, "amplmuta = ", self$aplMut, "\n")
      text <- paste0(text, "probmigra = ", self$probMigra, "\n")
      text <- paste0(text, "proboptimlocal = ", self$probOptimLocal, "\n")
      text <- paste0(text, "numoptimlocal = ", self$numOptimLocal, "\n")
      text <- paste0(text, "idoptimlocal = ", self$idOptimLocal, "\n")
      text <- paste0(text, "nichinggap = ", self$nichinggap, "\n")
      text <- paste0(text, "maxindniche = ", self$maxindniche, "\n")
      text <- paste0(text, "probintraniche = ", self$probintraniche, "\n")
      text <- paste0(text, "probcrossga = ", self$probcrossga, "\n")
      text <- paste0(text, "probmutaga = ", self$probmutaga, "\n")
      text <- paste0(text, "lenchaingap = ", self$lenchaingap, "\n")
      text <- paste0(text, "maxtreeheight = ", self$maxtreeheight, "\n")

      return(text)

    }
  )
)
