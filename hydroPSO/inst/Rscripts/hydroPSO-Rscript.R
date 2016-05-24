################################################################################
## File hydroPSO-Rscript.R                                                     #
## Part of the hydroPSO R package:                                             #
##                             http://www.rforge.net/hydroPSO/ ;               #
##                             http://cran.r-project.org/web/packages/hydroPSO #
## Copyright 2011-2012 Mauricio Zambrano-Bigiarini & Rodrigo Rojas             #
## Distributed under GPL 2 or later                                            #
##                                                                             #
## R script created with he PEST2hydroPSO function                             #
##                                                                             #
## Created by Mauricio Zambrano-Bigiarini and Rodrigo Rojas. 08-Nov-2012       #
## Last Update: 08-Nov-2012                                                    #
##              04-Jun-2013                                                    #
################################################################################

###Loading required libraries
library(hydroPSO)
library(hydroGOF)
# library(hydroTSM) # OPTIONAL package

################################################################################
##################   User-defined variables - START  ###########################
################################################################################

###Definition of working directory: input, output and model files paths
model.drty <- "user.model.drty" 
setwd(model.drty)

###Period of analysis. 
###Only required for out.FUN
###To define subperiod of analysis to compute the goodness-of-fit (GoF) measures
Sim.Ini <- "YYYY-MM-DD"
Sim.Fin <- "YYYY-MM-DD"
gof.Ini <- "YYYY-MM-DD"
gof.Fin <- "YYYY-MM-DD"


###Function for reading the simulated equivalents
out.FUN="read_output"   # name of user-defined function for reading model outputs
out.FUN.args=list( ###START out.FUN.args 
             file="filename.out"   # name of the output file
             #,,,                  # additional arguments for 'out.FUN'
                 ) ###END out.FUN.args
                 
                 
###Goodness-of-fit function, either pre-defined from hydroGOF (e.g., ssq) or 
###customized
gof.FUN <- "ssq" # sum of squared residuals. PEST default. 
                 # any other model performance measure could be used
gof.FUN.args <- list()

################################################################################
###################   User-defined variables - END  ############################
################################################################################

###Getting the OBSERVATIONS
obs.fname <- "PEST2hydroPSO_OBS.txt"
obs.fname <- paste(model.drty, "/PSO.in/", obs.fname, sep="")
obs       <- read.table(obs.fname)

###MAIN model function
model.FUN.args=list(
   model.drty=model.drty, 
   param.files=paste(model.drty,"/PSO.in/ParamFiles.txt",sep=""), 
   exe.fname="ModelCommandLine",                                  #TODO
   #stdout="",            
   out.FUN=out.FUN,  
   out.FUN.args=out.FUN.args,
   gof.FUN=gof.FUN,
   gof.FUN.args=gof.FUN.args,
   #gof.Ini=gof.Ini,        # un-comment if 'gof.Ini' is defined
   #gof.Fin=gof.Fin,        # un-comment if 'gof.Fin' is defined
   obs=obs
) ###END model.FUN.args

### MAIN PSO ALGORITHM
### For hydroPSO fine-tuning parameters, see Zambrano-Bigiarini and Rojas, 2012
hydroPSO(
   fn="hydromod",
   model.FUN="hydromod",        
   model.FUN.args=model.FUN.args,
   control=list(     ###START control options
      param.ranges="ParamRanges.txt",                              
      MinMax="min",  # minimisation of ssq
      #maxit=1000,   # case dependent
      #npart=40,     # case dependent
      #,,,           # additional arguments for hydroPSO
      REPORT=5       # frequency of screen messages
   )                 ###END control options
) ###END MAIN hydroPSO ALGORITHM

# Plotting the results
plot_results(MinMax="min", do.png=TRUE)
