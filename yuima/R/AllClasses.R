# Class Definitions
# This source MUST be loaded first

# Class 'yuima.pars'

# parameter object included in 'yuima.model'
setClass("model.parameter",representation(all="character",
                                          common="character",
                                          diffusion="character",
                                          drift="character",
                                          jump="character",
                                          measure="character",
# Insert parameters for starting conditions                                           
                                          xinit="character"
                                          )
         )

# Class 'yuima.model'
setClass("yuima.model",representation(drift="expression",
                                      diffusion="list",
                                      hurst="ANY",
                                      jump.coeff="list",
#jump.coeff="expression",
                                      measure="list",
                                      measure.type="character",
                                      parameter="model.parameter",
                                      state.variable="character",
                                      jump.variable="character",
                                      time.variable="character",
                                      noise.number="numeric",
                                      equation.number="numeric",
                                      dimension="numeric",
                                      solve.variable="character",
#                                       xinit="numeric",
                                      xinit="expression",
                                      J.flag="logical"
                                      )
         )

# Class 'carma.info'
setClass("carma.info",
         representation(p="numeric",
                        q="numeric",
                        loc.par="character",
                        scale.par="character",
                        ar.par="character",
                        ma.par="character",
                        lin.par="character",
                        Carma.var="character",
                        Latent.var="character",
                        XinExpr="logical")
         )

# Class 'yuima.carma'

setClass("yuima.carma",
         representation(info="carma.info"),
         contains="yuima.model")

# Class Compound Poisson
setClass("yuima.poisson", contains="yuima.model")


# Class 'yuima.data'

# we want yuimaS4 to use any class of data as input
# the original data will be stored in OrigData
# we convert these objects internally to "zoo" object
# in the future, we may want to use more flexible
# classes

setClass("yuima.data", representation(original.data = "ANY",
                                      zoo.data = "ANY"
                                      )
         )


# Class 'yuima.sampling'

# sampling is now empty, but should give informations on the sampling
# type, rate, deltas, etc.

setClass("yuima.sampling", representation(Initial  = "numeric",
										  Terminal = "numeric",
                                          n = "numeric",
										  delta    = "numeric",
										  grid     = "ANY",
										  random   = "ANY",
										  regular  = "logical",
										  sdelta   = "numeric",
										  sgrid    = "ANY",
										  oindex   = "ANY",
										  interpolation = "character"
                                          )
         )

# Class 'yuima.functional'

# functional model used in 'asymptotic term' procedure

setClass("yuima.functional", representation(F = "ANY",
                                          f = "list",
                                          xinit = "numeric",
                                          e = "numeric"
                                          )
         )
         

# Class 'yuima'

# this is the principal class of yuima project. It may contain up to
# three slots for now: the data, the model and the sampling

setClass("yuima.characteristic", representation(equation.number = "numeric",
                                                time.scale = "numeric"
                                                )
         )


setClass("yuima", representation(data = "yuima.data",
                                 model = "yuima.model",
                                 sampling = "yuima.sampling",
                                 characteristic = "yuima.characteristic",
								 functional = "yuima.functional"
                                 )
         )

# Class yuima.carma.qmle
setClass("yuima.carma.qmle",representation(Incr.Lev = "ANY",
                                           model = "yuima.carma",
                                           logL.Incr = "ANY"
                                           ),
                            contains="mle"
         )

setClass("yuima.qmle",representation(
model = "yuima.model"),
contains="mle"
)

setClass("yuima.CP.qmle",representation(Jump.times = "ANY",
Jump.values = "ANY",
X.values = "ANY",
model = "yuima.model",
threshold="ANY"),
contains="mle"
)

setClass("summary.yuima.carma.qmle",representation(MeanI = "ANY",
                                                   SdI = "ANY",
                                                   logLI = "ANY",
                                                   TypeI = "ANY",
                                                   NumbI = "ANY",
                                                   StatI ="ANY"),
contains="summary.mle"
)

setClass("summary.yuima.CP.qmle",
representation(NJ = "ANY",
MeanJ = "ANY",
SdJ = "ANY",
MeanT = "ANY",
Jump.times = "ANY",
Jump.values = "ANY",
X.values = "ANY",
model = "yuima.model",
threshold = "ANY"),
contains="summary.mle"
)


setClass("summary.yuima.qmle",
representation(
model = "yuima.model",
threshold = "ANY"),
contains="summary.mle"
)

# The yuima.carma.qmle extends the S4 class "mle". It contains three slots: Estimated Levy,
# The description of the carma model and the mle.