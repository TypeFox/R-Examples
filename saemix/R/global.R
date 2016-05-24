####################################################################################
####				Generic functions				####
####################################################################################

setGeneric(name="read.saemixData",
  def=function(object){standardGeneric("read.saemixData")}
)

setGeneric(name="showall",
  def=function(object){standardGeneric("showall")}
)

setGeneric(name="psi",
  def=function(object,indiv.par){standardGeneric("psi")}
)

setGeneric(name="phi",
  def=function(object,indiv.par){standardGeneric("phi")}
)

setGeneric(name="eta",
  def=function(object,indiv.par){standardGeneric("eta")}
)

# setGeneric(name="predict",def=function(object){standardGeneric("predict")})

####################################################################################
####				Library functions				####
####################################################################################

#source("func_main.R")
#source("func_aux.R")
#source("func_plots.R")

####################################################################################
