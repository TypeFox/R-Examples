#################################################################
#
# stmodel.R
#
#######################
# stepp model         #
#######################
setClassUnion("stmodel",
	        c("stmodelGLM","stmodelCOX","stmodelKM","stmodelCI")
		  )  


