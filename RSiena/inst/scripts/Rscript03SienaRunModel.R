################################################################################
###
###
### ---- Rscript03SienaRunModel.R: a script for the introduction to RSiena -----
###
###							version April 16, 2013
################################################################################
#
# The introductory script is divided into the following script files:
# Rscript01DataFormat.R, followed by
# RScriptSNADescriptives.R, code for descriptive analysis of the data, and
# Rscript02SienaVariableFormat.R, which formats data and specifies the model,
# Rscript03SienaRunModel.R, which runs the model and estimates parameters, and
# Rscript04SienaBehaviour.R, which illustrates an example of analysing the
# coevolution of networks and behaviour.
# Written with contributions by Robin Gauthier, Tom Snijders, Ruth Ripley,
# Johan Koskinen, and Paulina Preciado.
#
# This script, Rscript03SienaRunModel.R, runs the estimation in RSiena for the
# model set up and defined in the script Rscript02SienaVariableFormat.R.
#
# A quick version of the model fitting without comments is given at the end
# of this script

########################### ESTIMATION OF PARAMETERS ###########################

# Parameters of the model are estimated by the function siena07.
# This requires the data specification; the effects specification;
# and a number of parameters, or settings, for the estimation algorithm.
# The latter are contained in an object created by the function
# sienaAlgorithmCreate. You can look at the help provided by
# ?sienaAlgorithmCreate
# to find out about options that you may use here;
# for beginning users, only the two options mentioned below are relevant.
#
# Output will be written to a file with name projname.out, where projname is
# whatever name is given; the default (used if no name is given) is Siena.
# This file will be written to your current directory.
# New estimation runs will append to it.
# A new call to print01Report will overwrite it!

		myalgorithm <- sienaAlgorithmCreate(useStdInits = FALSE, projname = 's50_3')

# The useStdInits parameter determines the initial values used for
# the estimation algorithm.
# If useStdInits = TRUE, standard initial values are used;
# if useStdInits = FALSE, the initial values are used that are contained
# in the "initialValue" column of the effects object,
# which are reported by the information request
		myeff
# Below we shall see how these initial values can be altered.

# Let us first redefine the model, to obtain a simpler specification
# that will serve as an illustration here.

		myeff <- getEffects( mydata )
		myeff <- includeEffects( myeff, transTrip, cycle3 )
		myeff <- includeEffects( myeff, egoX, altX, egoXaltX,
								 interaction1 = "alcohol" )
		myeff <- includeEffects( myeff, simX, interaction1 = "smoke1" )
		myeff

# The function siena07 actually fits the specified model to the data
# If you wish the pretty picture of Siena on the screen as information
# about the progress of the algorithm, type

		ans <- siena07( myalgorithm, data = mydata, effects = myeff)

# (ans for "answer").
# If however you do not want the pretty picture, or if this leads to
# difficulties (which may happen e.g. on a Mac), then type

#		ans <- siena07(myalgorithm, data=mydata, effects=myeff, batch=TRUE)

# and intermediate information will be written to the console.

# Function siena07 produces a so-called sienaFit object, here called ans;
# and it fills in a few things in the sienaEffects object myeff,
# if this is the first use of myeff in a siena07 call.
# By using various different effects objects, i.e., with different names,
# you can switch between specifications.

# The batch = FALSE parameters will give a graphical user interface being opened
# which reports on the progress of the estimation algorithm;

# verbose = TRUE leads to extensive diagnostic information being sent
# to the console during the estimation, and results after the estimation
# (these results are also copied to the output file projname.out, see above);
# while batch=TRUE gives only a limited amount of printout sent to the console
# during the estimation (which is seen when clicking in the console,
# or more immediately if the Buffered Output is deselected in the Misc menu)
# which monitors the progress of the estimation algorithm in a different way.

# The call of siena07 leads to output in the file s50_3.out
# (or more generally projname.out,
# where projname is the name given in sienaAlgorithmCreate)
# and to the creation of the object which here is called ans (for "answer").

# To use multiple processors, in the simplest case where your computer has 2
# processors, use

#		ans <- siena07( myalgorithm, data = mydata, effects = myeff,
#					  nbrNodes = 2, useCluster = TRUE)

# Adjust the nbrNodes to the number available.
# If you wish to work on with other programs while running siena07,
# it is advisable to use one node less than the number of available processors.
# If you wish to use other machines as well,
# see the more detailed instructions below.
# You will then need to use the clusterString argument as well.
#
# For more advanced use, it can be helpful to have access to the networks
# simulated in the so-called third phase of the estimation algorithm.
# These networks can be used, e.g., for checking goodness of fit.
# This can be achieved by using the parameter returnDeps=TRUE.
# The fitted object ans will then have a component named "sims"
# which contains a list (each iteration) of lists (each data object)
# of lists (each dependent network or behavior variable) of edgelists for
# networks or vectors for behavior variables.
# See the manual for further explanation.
#
# This option when used with multiple processors would require
# rather a lot of communication between multiple processes,
# slowing down the computations,
# so it might be better to avoid using the two options together.


################### LOOKING AT THE RESULTS ################################

# The file "s50_3.out" will contain the results of the estimation.
# It is contained in the current directory ("getwd()").
# This file can be read by any text editor.
# A summary of the results is obtained on the screen by

		ans

# and a larger summary by

		summary(ans)

# Depending on the random seed and the model specification,
# the results could be something like the following.

# Estimates, standard errors and convergence t-ratios
#
#										Estimate   Standard	  Convergence
#													 Error		t-ratio
#
#Rate parameters:
#  0.1		Rate parameter period 1		 6.6331	 ( 1.1770	)
#  0.2		Rate parameter period 2		 5.2105	 ( 0.8790	)
#
#Other parameters:
#  1.  eval outdegree (density)			-2.7194	 ( 0.1191	)	 0.0024
#  2.  eval reciprocity					 2.4344	 ( 0.2208	)	 0.0013
#  3.  eval transitive triplets			 0.6449	 ( 0.1378	)	 0.0039
#  4.  eval 3-cycles					-0.0881	 ( 0.2917	)	 0.0168
#  5.  eval smoke1 similarity			 0.2239	 ( 0.1991	)	-0.0792
#  6.  eval alcohol alter				-0.0256	 ( 0.0676	)	-0.0048
#  7.  eval alcohol ego					 0.0348	 ( 0.0730	)	 0.0075
#  8.  eval alcohol ego x alcohol alter	 0.1295	 ( 0.0489	)	 0.0126
#

# The results can also be viewed externally in the output file s50_3.out
# It is advisable that you have a look at all three reports and
# understand how information is organized in each of them.

# To understand the table above, note that the "convergence t-ratio"
# is the t-ratio for convergence checking,
# not the t statistic for testing the significance of this effect!
# (See Section 6.1.2 of the manual.)
# In the external output file, these are called
# "t-ratios for deviations from targets".
# The rule of thumb is that all t-ratios for convergence
# should ideally be less than 0.1 in absolute value;
# this signifies good convergence of the algorithm.
# In the example here, this is the case.
# If this would not be the case, the best thing to do would be
# to continue the estimation, using the estimates produced here,
# and contained in ans, as the new initial values.
# This is explained below.
# Because the estimation algorithm is based on random simulations of the
# network evolution, there always will be small differences between
# different runs of the algorithm.
# To obtain "publication grade" estimates, where this variability is minimized,
# choose the value of parameter n3 in sienaAlgorithmCreate()
# ("Number of iterations in phase 3") larger than the default value of 1000;
# e.g., n3=3000.

# With function siena07 we made ans as the object containing
# all the results of the estimation. For example,

		ans$theta

# contains the vector of parameter estimates while

		ans$covtheta

# contains the covariance matrix of the estimates.
# There are several "methods" available for viewing the object
# containing the results of the estimation.
# Above we already mentioned
#		ans
# and
#		summary( ans )
# The command

		xtable( ans )

# will produce a table formatted for inclusion in a LaTeX document
# or formatted in html. Use e.g.

		xtable( ans, type = 'html' )

# to get html, and e.g.

		xtable( ans, file = 'ff.tex' )

# to write the results to a file.
# At http://cran.r-project.org/web/packages/xtable you can find
# a set of vignettes for the xtable package, the xtable gallery,
# which gives more options.
# A function siena.table is available that is specially designed
# for RSiena results and produces html or LaTeX files.
# See

		?print.sienaFit


############## MORE ON INITIALIZING PARAMETERS FOR ESTIMATION ########

# If the estimation algorithm has not produced good estimates
# (it 'has not converged well'),
# as will be indicated by some of the t-ratios for convergence being larger
# than 0.1 (this threshold is not to be taken too precisely, though),
# the best thing to do is continuing the estimation,
# using the estimates produced here,
# and contained in ans, as the new initial values.
# This is done by the option prevAns ('previous ans') as in

		ans <- siena07(myalgorithm, data=mydata, effects=myeff, prevAns=ans)

# the parameter estimates in ans then are extracted and
# used in the new estimation;
# moreover, Phase 1 will be omitted from the algorithm,
# as derivatives and covariance matrix are used from the previous run.
# This should be used only if the model specification in myeff
# has not changed, and if the provisional parameter estimates obtained
# in ans are reasonable; if they are not reasonable,
# omit the prevAns option, use
#		myalgorithm$useStdInits <- TRUE
# to get back on track, and return at the next estimation to
#		myalgorithm$useStdInits <- FALSE
# To understand what happens here, read on:

# Another and more flexible way for determining initial values is by
# using the useStdInits element of the model object,
# and the initial values in the effects object.
# This is done as follows.
# The option useStdInits = TRUE in sienaAlgorithmCreate, will make
# each estimation run start with standard initial values.
# The option useStdInits = FALSE makes the estimation start
# with the initial values in the effects object.
# You can switch between these by commands such as

#		myalgorithm$useStdInits <- FALSE
#		myalgorithm$useStdInits <- TRUE

# Putting the estimates from the results object ans into the
# effects object myeff is done by

		myeff <- updateTheta(myeff, ans)

# A check that the effects object contains the desired initial values is made by

		myeff

# The initial values are in the vector

		myeff$initialValue[myeff$include]

# and this also can be initialised differently, if this is desired.
# Note that this initial vector will be used until you change it again,
# e.g., to the results of a new run, or if you use the prevAns parameter,
# or until you set useStdInits to TRUE.

################################################################################
###
### ---- Testing effects -------------------------------------------------------
###
################################################################################
#
# Two types of tests are available in SIENA.
# 1. t-type tests of single parameters can be carried out by dividing
# the parameter estimate by its standard error.
# Under the null hypothesis that the parameter is 0, these tests have
# approximately a standard normal distribution.
# 2. Score-type tests of single and multiple parameters are described
# in the manual.
# Parameters can be restricted by putting TRUE in the
# include, fix and test columns of the effects object.
# For example, to request a score test for the indegree popularity effect,
# the commands can be as follows.

		myeff <- setEffect(myeff, inPopSqrt, fix=TRUE, test=TRUE,
										  initialValue=0.0)
		ans <- siena07(myalgorithm, data=mydata, effects=myeff)

# After such an operation, again request

		summary(ans)

# to see the results, including those for the score test.


################################################################################
###
### ---- Time test -------------------------------------------------------------
###
################################################################################
#
# An application of the score test is given for the special case of parameter
# heterogeneity by Lospinoso et al. (2010) and implemented in RSiena.
# To apply the test to the results obtained above, request, e.g.,
		tt2 <- sienaTimeTest(ans)
		tt2
# If you wish more information, also
#		summary(tt2)
#		 plot(tt2, effects=3:4)
# If as a consequence of this analysis you wish to add time dummy terms,
# this may be done via
#		myeff <- includeTimeDummy(myeff, transTrip, cycle3)
#		myeff
#		ans3 <- siena07(myalgorithm, data=mydata, effects=myeff)
# and testing again,
#		(tt3 <- sienaTimeTest(ans3))
# and so on.

################################################################################
###
### ---- Summary of model fitted -----------------------------------------------
###
################################################################################

	  friend.data.w1 <- as.matrix(read.table("s50-network1.dat")) # read data
	  friend.data.w2 <- as.matrix(read.table("s50-network2.dat"))
	  friend.data.w3 <- as.matrix(read.table("s50-network3.dat"))
	  drink <- as.matrix(read.table("s50-alcohol.dat"))
	  smoke <- as.matrix(read.table("s50-smoke.dat"))

	  friend.data.w1[ friend.data.w1 %in% c(6,9) ] <- NA # define missing data
	  friend.data.w1[ friend.data.w2 %in% c(6,9) ] <- NA
	  friend.data.w1[ friend.data.w3 %in% c(6,9) ] <- NA

	  friendship <- sienaDependent( array( c( friend.data.w1,
											  friend.data.w2, friend.data.w3 ),
									dim = c( 50, 50, 3 ) ) )

	  drinkingbeh <- sienaDependent( drink, type = "behavior" )
	  smoke1 <- coCovar( smoke[ , 1 ] )
	  alcohol <- varCovar( drink )

	  mydata <- sienaDataCreate( friendship, smoke1, alcohol )

	  myeff <- getEffects( mydata )  # create effects structure

	  print01Report( mydata, myeff, modelname = 's50_3_init' )

	  myeff <- includeEffects( myeff, transTrip, cycle3 )
	  myeff <- includeEffects( myeff, egoX, altX,
							   egoXaltX, interaction1 = "alcohol" )
	  myeff <- includeEffects( myeff, simX, interaction1 = "smoke1" )

	  myalgorithm <- sienaAlgorithmCreate( projname = 's50_3' )
# and finally (commented out because no need to repeat the actual calculations)
#	  ans <- siena07( myalgorithm, data = mydata, effects = myeff)

################################################################################
###
### -- PROCEED TO Rscript04SienaBehaviour.R FOR
###									MODELING NETWORKS AND BEHAVIOUR BY RSIENA --
###
################################################################################
