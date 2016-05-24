####################################################
### Authors: Simone Padoan and Moreno Bevilacqua.
### Emails: simone.padoan@unibocconi.it,
### moreno.bevilacqua@uv.cl
### Institutions: Department of Decision Sciences,
### University Bocconi of Milan and
### Departamento de Estadistica
### Universidad de Valparaiso
### File name: Utility.r
### Description:
### This file contains a set of procedures
### for the set up of all the package routines.
### Last change: 28/03/2013.
####################################################

### Procedures are in alphabetical order.

CheckCorrModel <- function(corrmodel)
  {
    CheckCorrModel <- NULL
    # Correlation function are in alphabetical order
    CheckCorrModel <- switch(# spatial or temporal correlations
                             corrmodel,
                             cauchy=1,
                             exponential=2,
                             gauss=3,
                             gencauchy=4,
                             spherical=5,
                             stable=6,
                             matern=7,
                             # spatial-temporal non-separable correlations
                             gneiting=21,
                             iacocesare=22,
                             porcu=23,
                             stein=24,
                             porcu2=25,
                             gneiting_GC=26,
                             porcu1=30,
                             # spatial-temporal separable correlations
                             exp_cauchy=41,
                             exp_exp=42,
                             exp_gauss=43,
                             gauss_gauss=44,##not implemented
                             matern_cauchy=45,
                             matern_exp=46,
                             stab_stab=47,
                             #Taper functions
                             Wendland1=15,
                             Wendland2=16,
                             Wendland3=17,
                             Wendland1_Wendland1=100,
                             Wendland1_Wendland2=101,
                             Wendland1_Wendland3=102,
                             Wendland2_Wendland1=103,
                             Wendland2_Wendland2=104,
                             Wendland2_Wendland3=105,
                             Wendland3_Wendland1=106,
                             Wendland3_Wendland2=107,
                             Wendland3_Wendland3=108,
                             qt_time=109,
                             qt_space=110)
    return(CheckCorrModel)
  }

#CheckCorrTap <- function(corrtap)
#  {
#    CheckCorrTap <- NULL
#    # Correlation function are in alphabetical order
#    CheckCorrTap <- switch(corrtap,
#                           Wendland1=15,
#                           Wendland2=16,
#                           Wendland3=17,
#                           Wendland1_Wendland1=100,
#                           Wendland1_Wendland2=101,
#                           Wendland1_Wendland3=102,
#                           Wendland2_Wendland1=103,
#                           Wendland2_Wendland2=104,
#                           Wendland2_Wendland3=105,
#                           Wendland3_Wendland1=106,
#                           Wendland3_Wendland2=107,
#                           Wendland3_Wendland3=108,
#                           qt_time=109,
#                           qt_space=110)
#    return(CheckCorrTap)
#  }


CheckInput <- function(coordx, coordy, coordt, corrmodel, data, distance, fcall, fixed, grid,
                      likelihood, margins, maxdist, maxtime, model, numblock, optimizer, param,
                       replicates, start, taper, tapsep, threshold, type, varest, vartype, weighted)
  {
    error <- NULL
    # START Include internal functions:
    CheckParamRange <- function(param)
    {
        if(!is.na(param['df'])) if(param['df'] <= 0) return(FALSE)
        if(!is.na(param['nugget'])) if(param['nugget'] < 0) return(FALSE)
        if(!is.na(param['power'])) if(param['power'] <=0 || param['power'] > 2) return(FALSE)
        if(!is.na(param['power_s'])) if(param['power_s'] <=0 || param['power_s'] > 2) return(FALSE)
        if(!is.na(param['power_t'])) if(param['power_t'] <=0 || param['power_t'] > 2) return(FALSE)
        if(!is.na(param['power1'])) if(param['power1'] <=0 || param['power1'] > 2) return(FALSE)
        if(!is.na(param['power2'])) if(param['power2'] <= 0) return(FALSE)
        if(!is.na(param['sep'])) if(param['sep'] < 0 || param['sep'] > 1) return(FALSE)
        if(!is.na(param['scale'])) if(param['scale'] <= 0) return(FALSE)
        if(!is.na(param['scale_s'])) if(param['scale_s'] <= 0) return(FALSE)
        if(!is.na(param['scale_t'])) if(param['scale_t'] <= 0) return(FALSE)
        if(!is.na(param['sill'])) if(param['sill'] <= 0) return(FALSE)
        if(!is.na(param['smooth'])) if(param['smooth'] <= 0) return(FALSE)
        if(!is.na(param['smooth_s'])) if(param['smooth_s'] <= 0) return(FALSE)
        if(!is.na(param['smooth_t'])) if(param['smooth_t'] <= 0) return(FALSE)
        return(TRUE)
    }
    # Check if the correlation is spatial or spatial-temporal
    CheckSpaceTime <- function(corrmodel)
    {
        CheckSpaceTime <- NULL
        if(corrmodel >= 1 & corrmodel <= 20) CheckSpaceTime <- FALSE
        else CheckSpaceTime <- TRUE
        return(CheckSpaceTime)
    }
    # Check the type of distances
    CheckDistance<- function(distance)
    {
        CheckDistance <- NULL
        CheckDistance <- switch(distance,
                                eucl=0,
                                Eucl=0,
                                chor=1,
                                Chor=1,
                                geod=2,
                                Geod=2,
                                proj=3,
                                Proj=3)
        return(CheckDistance)
    }
    # END Include internal functions
    if(fcall!="Kriging"){
    ### START Checks inserted input
    # START common check fitting and simulation
    if(missing(coordx) || !is.numeric(coordx)){
        error <- 'insert a suitable set of numeric coordinates\n'
        return(list(error=error))}

    if(!is.null(coordy) & !is.numeric(coordy)){
        error <- 'insert a suitable set of numeric coordinates\n'
        return(list(error=error))}

    if(missing(corrmodel) || !is.character(corrmodel)){
        error <- 'insert the correlation model\n'
        return(list(error=error))}

    if(!is.null(grid) & !is.logical(grid)){
        error <- 'the parameter grid need to be a logic value\n'
        return(list(error=error))}

    if(!is.null(model) & !is.character(model)){
        error <- 'insert the name of the random field\n'
        return(list(error=error))}

    if(is.null(CheckModel(model))){
        error <- 'the model name of the random field is not correct\n'
        return(list(error=error))}

    if(is.null(replicates) || (abs(replicates-round(replicates))>0) || replicates<1){
        error <- 'the parameter replicates need to be a positive integer\n'
        return(list(error=error))}

    if(is.null(CheckDistance(distance))){
        error <- 'the name of the distance is not correct\n'
        return(list(error=error))}

    if(is.null(CheckCorrModel(corrmodel))){
        error <- 'the name of the correlation model is not correct\n'
        return(list(error=error))}

    # END common check fitting and simulation
    # START check fitting
    if(fcall=="Fitting"){
        if(missing(data) || !is.numeric(data)){
            error <- 'insert a numeric vector or matrix of data\n'
            return(list(error=error))}

        if(!is.null(fixed) & !is.list(fixed)){
            error <- 'insert fixed values as a list of parameters\n'
            return(list(error=error))}

        if(!is.null(fixed)){
            namfixed <- names(fixed)
            if(!all(namfixed %in% c(NuisanceParam(model), CorrelationParam(corrmodel)))){
                error <- 'some names of the fixed parameters is/are not correct\n'
                return(list(error=error))}

        if(!CheckParamRange(unlist(fixed))){
            error <- 'some fixed values are out of the range\n'
            return(list(error=error))}}

        if(!is.null(likelihood) & !is.character(likelihood)){
            error <- 'insert the type of likelihood objects\n'
            return(list(error=error))}

        if(!is.null(maxdist)){
            error <- "insert a positive numeric value for the maximum spatial distance\n"
            if(!is.numeric(maxdist)) return(list(error=error))
            else if(maxdist<0) return(list(error=error))}

        if(!is.null(maxtime)){
            error <- "insert a positive numeric value for the maximum time interval\n"
            if(!is.numeric(maxtime)) return(list(error=error))
            else if(maxtime<0) return(list(error=error))}

        if(model=="BinaryGauss"){
            if(length(unique(c(data)))!=2){error <- 'the data are not binary'
            return(list(error=error))}
            if(!is.null(start$sill)) if(start$sill>1){error <- 'some starting values are out of the range\n'
                                                     return(list(error=error))}
            if(!is.null(fixed$sill)) if(fixed$sill>1){error <- 'some starting values are out of the range\n'
                                                     return(list(error=error))}}

        if(model %in% c("BrownResn","ExtGauss","ExtT"))
            if(!margins %in% c("Frechet","Gev")){
                error <- 'insert the correct type of marginal distributions\n'
                return(list(error=error))}

        if(model=="BrowResn")
            if(CheckSpaceTime(CheckCorrModel(corrmodel))){
                if(!corrmodel %in% c("exp_exp","exp_gauss","stable_stable")){
                    error <- "the correlation model is not adequate for the Brown-Resnick process\n"
                    return(list(error=error))}}
            else
                if(!corrmodel %in% c("exponential","gauss","gencauchy","stable")){
                    error <- "the correlation model is not adequate for the Brown-Resnick process\n"
                    return(list(error=error))}

        if(!is.null(optimizer) & !is.character(optimizer)){
            error <- 'insert the type of maximising algorithm\n'
            return(list(error=error))}

        if(!is.null(varest) & !is.logical(varest)){
            error <- 'the parameter std.err need to be a logical value\n'
            return(list(error=error))}

        if(type=="Tapering"){
        if(is.null(taper) || is.null(maxdist)){
          error <- 'tapering need a taper correlation model and/or a compact support\n'
          return(list(error=error))}
        if(!taper %in% c("Wendland1","Wendland2","Wendland3",
                         "Wendland1_Wendland1","Wendland1_Wendland2","Wendland1_Wendland3",
                         "Wendland2_Wendland1","Wendland2_Wendland2","Wendland2_Wendland3",
                         "Wendland3_Wendland1","Wendland3_Wendland2","Wendland3_Wendland3",
                         "qt_time","qt_space")){
            error <- 'insert a correct name for the taper correlation model\n'
            return(list(error=error))}}

        if(!is.null(type) & !is.character(type)){
            error <- 'insert the configuration of the likelihood objects\n'
            return(list(error=error))}

        if(is.null(CheckType(type))){
            error <- 'the type name of the likelihood objects is not correct\n'
            return(list(error=error))}

        if(is.null(CheckLikelihood(likelihood))){
           error <- 'the setting name of the likelihood objects is not correct\n'
           return(list(error=error))}

        if(likelihood == "Full"){
            if(!any(type == c("Restricted", "Standard", "Tapering"))){
                error <- 'insert a type name of the likelihood objects compatible with the full likelihood\n'
                return(list(error=error))}}

        if(likelihood == "Marginal"){
            if(!any(type == c("Difference", "Pairwise"))){
                error <- 'insert a type name of the likelihood objects compatible with the composite-likelihood\n'
                return(list(error=error))}}

        if(varest & (likelihood == "Conditional" || likelihood == "Marginal") & (!is.null(vartype) & !is.character(vartype))){
            error <- 'insert the type of estimation method for the variances\n'
            return(list(error=error))}

        if(varest & is.null(CheckVarType(vartype)) & (likelihood == "Conditional" || likelihood == "Marginal")){
            error <- 'the name of the estimation type for the variances is not correct\n'
            return(list(error=error))}

          if(!is.null(tapsep)){
            error <- "insert a numeric value between 0 and 1 for the separability parameter of the space time taper\n"
            if(!is.numeric(tapsep)) return(list(error=error))
            else if(tapsep<0||tapsep>1) return(list(error=error))}

        if(!is.null(start)){
            if(!is.list(start)){
                error <- 'insert starting values as a list of parameters\n'
                return(list(error=error))}

            namstart <- names(start)

            if(!all(namstart %in% c(NuisanceParam(model), CorrelationParam(corrmodel)))){
                error <- 'some names of the starting parameters is/are not correct\n'
                return(list(error=error))}

            if(any(namstart=='mean') & (type=='Difference' || type=='Restricted')){
                error <- 'the mean parameter is not allow with the difference composite likelihood\n'
                return(list(error=error))}

            if(corrmodel=="gencauchy" && param["power1"]!=2){
                error <- "the parameter power1 need to be equal to 2 for the Brown-Resnick process\n"
                 return(list(error=error))}

            if(any(namstart=='sill') & (model=='BrowResn')){
                error <- 'the sill parameter is not allow with Brown-Renick model\n'
                return(list(error=error))}

            if(!CheckParamRange(unlist(start))){
                error <- 'some starting values are out of the range\n'
                return(list(error=error))}

            if(!is.null(fixed))
                if(any(namstart %in% namfixed)){
                    error <- 'some fixed parameter name/s is/are matching with starting parameter name/s\n'
                    return(list(error=error))}}

        if(!is.null(weighted) & !is.logical(weighted)){
            error <- 'insert if the composite likelihood need to be weighted'
            return(list(error=error))}

        # START - checks the format of the coordinates and dataset
    dimdata <- dim(data) # set the data dimension
    if(is.null(coordt)) # START 1) spatial random field
      {
        if(CheckSpaceTime(CheckCorrModel(corrmodel)))
          {
            error <- 'temporal coordinates are missing\n'
            return(list(error=error))
          }
        if(replicates>1) # a) n iid replicates of a spatial random field
          {
            if(grid) # START regular grid
              {
                if(is.null(dimdata))
                  {
                    error <- c('insert an array d x d of n iid spatial observations\n')
                    return(list(error=error))
                  }
                if(length(dimdata)!=3)
                  {
                    error <- c('the dimension of the data matrix is not correct\n')
                    return(list(error=error))
                  }
                if(length(coordx)!=dimdata[1] || length(coordy)!=dimdata[2])
                  {
                    error <- c('the number of coordinates does not match with the number of spatial observations\n')
                    return(list(error=error))
                  }
                if(dimdata[3]!=replicates)
                  {
                    error <- c('the number of replications does not match with the data observations\n')
                    return(list(error=error))
                  }
              } # END regular grid
            else # START irregular grid
              {
                if(is.null(dimdata))
                  {
                    error <- c('insert a matrix n x d of spatial observations\n')
                    return(list(error=error))
                  }
                if(length(dimdata)!=2)
                  {
                    error <- c('the dimension of the data matrix is not correct\n')
                    return(list(error=error))
                  }
                if(is.null(coordy))
                  {
                    if(is.null(dim(coordx)))
                      {
                        error <- c('insert a matrix d x 2 of spatial coordinates\n')
                        return(list(error=error))
                      }
                    if(dimdata[2]!=nrow(coordx) || ncol(coordx)!=2)
                      {
                        error <- c('the number of coordinates does not match with the number of spatial observations\n')
                        return(list(error=error))
                      }
                  }
                else
                  if(length(coordx)!=length(coordy))
                    {
                      error <- c('the number of the two coordinates does not match\n')
                      return(list(error=error))
                    }
                if(dimdata[1]!=replicates)
                  {
                    error <- c('the number of replications does not match with the data observations\n')
                    return(list(error=error))
                  }
              } # END irregular grid
          }
        else # b) START one realisation of a spatial random field
          {
            if(grid) # START regular grid
              {
                if(is.null(dimdata))
                  {
                    error <- c('insert a matrix d x d of spatial observations\n')
                    return(list(error=error))
                  }
                if(length(dimdata)!=2)
                  {
                    error <- c('the dimension of the data matrix is not correct\n')
                    return(list(error=error))
                  }
                if(length(coordx)!=dimdata[1] || length(coordy)!=dimdata[2])
                  {
                    error <- c('the number of coordinates does not match with the number of spatial observations\n')
                    return(list(error=error))
                  }
              } # END regular grid
            else # START irregular grid
              {
                numsite <- length(data)
                if(is.null(numsite))
                  {
                    error <- c('insert a vector of spatial observations\n')
                    return(list(error=error))
                  }
                if(is.null(coordy))
                  {
                    dimcoord <- dim(coordx)
                    if(is.null(dimcoord))
                      {
                        error <- c('insert a suitable set of coordinates\n')
                        return(list(error=error))
                      }
                    else
                      {
                        if(dimcoord[1]!=numsite || dimcoord[2]!=2)
                          {
                            error <- c('the number of coordinates does not match with the number of spatial observations\n')
                            return(list(error=error))
                          }
                      }
                  }
                else
                  {
                    if(length(coordx)!=length(coordy))
                      {
                        error <- c('the number of the two coordinates does not match\n')
                        return(list(error=error))
                      }
                    if(length(coordx)!=numsite)
                      {
                        error <- c('the number of coordinates does not match with the number of spatial observations\n')
                        return(list(error=error))
                      }
                  }
              }
          } # b) END one realisation of a spatial random field
      } # END 1) spatial random field
    else # 2) case: spatial-temporal random field
      {
        if(!is.numeric(coordt))
          {
            error <- 'insert a numerical vector of temporal coordinates\n'
            return(list(error=error))
          }
        if(length(coordt)<=1)
          {
            error <- 'insert a numerical vector of temporal coordinates\n'
            return(list(error=error))
          }
        if(replicates>1) # START a) n iid replicates of a spatial-temporal random field
          {
            if(grid) # START regular grid
              {
                if(is.null(dimdata))
                  {
                    error <- c('insert an array d x d x t of n iid spatial-temporal observations\n')
                    return(list(error=error))
                  }
                if(length(dimdata)!=4)
                  {
                    error <- c('the dimension of the data matrix is not correct\n')
                    return(list(error=error))
                  }
                if(length(coordx)!=dimdata[1] || length(coordy)!=dimdata[2] )
                  {
                    error <- c('the number of coordinates does not match with the number of spatial observations\n')
                    return(list(error=error))
                  }
                if(length(coordt)!=dimdata[3])
                  {
                    error <- c('the number of the temporal coordinate does not match with the third dimensiona of the data matrix\n')
                    return(list(error=error))
                  }
                if(dimdata[4]!=replicates)
                  {
                    error <- c('the number of replications doen not match with the data observations\n')
                    return(list(error=error))
                  }
              } # END regular grid
            else # START irregular grid
              {
                if(is.null(dimdata))
                  {
                    error <- c('insert an array t x d of n iid spatial-temporal observations\n')
                    return(list(error=error))
                  }
                if(length(dimdata)!=3)
                  {
                    error <- c('the dimension of the data matrix is not correct\n')
                    return(list(error=error))
                  }
                if(is.null(coordy))
                  {
                    if(dimdata[2]!=nrow(coordx) || ncol(coordx)!=2)
                      {
                        error <- c('the number of coordinates does not match with the number of spatial observations\n')
                        return(list(error=error))
                      }
                  }
                else
                  {
                    if(length(coordx)!=length(coordy))
                      {
                        error <- c('the number of the two spatial coordinates does not match\n')
                        return(list(error=error))
                      }
                    if(length(coordx)!=dimdata[2])
                      {
                        error <- c('the number of coordinates does not match with the number of spatial observations\n')
                        return(list(error=error))
                      }
                  }
                if(dimdata[1]!=length(coordt))
                  {
                    error <- c('the time coordinate does not match with the number of rows of the data array\n')
                    return(list(error=error))
                  }
                if(dimdata[3]!=replicates)
                  {
                    error <- c('the number of replications does not match with the data observations\n')
                    return(list(error=error))
                  }
              } # END irregular grid
          }# END a) n iid replicates of a spatial-temporal random field
        else # START b) one realisation of a spatial-temporal random field
          {
            if(grid) # START regular grid
              {
                if(is.null(dimdata))
                  {
                    error <- c('insert an array of d x d x t spatial-temporal observations\n')
                    return(list(error=error))
                  }
                if(length(dimdata)!=3)
                  {
                    error <- c('the dimension of the data matrix is not correct\n')
                    return(list(error=error))
                  }
                if(length(coordx)!=dimdata[1] || length(coordy)!=dimdata[2])
                  {
                    error <- c('the number of coordinates does not match with the number of spatial observations\n')
                    return(list(error=error))
                  }
                if(dimdata[3]!=length(coordt))
                  {
                    error <- c('the time coordinate does not match with the third dimension of the data array\n')
                    return(list(error=error))
                  }
              } # END regular grid
            else # START irregular grid
              {
                if(is.null(dimdata))
                  {
                    error <- c('insert a matrix of t x d spatial-temporal observations\n')
                    return(list(error=error))
                  }
                if(length(dimdata)!=2)
                  {
                    error <- c('the dimension of the data matrix is not correct\n')
                    return(list(error=error))
                  }
                if(is.null(coordy))
                  {
                    if(dimdata[2]!=nrow(coordx) || ncol(coordx)!=2)
                      {
                        error <- c('the number of coordinates does not match with the number of spatial observations\n')
                        return(list(error=error))
                      }
                  }
                else
                  {
                    if(length(coordx)!=length(coordy))
                      {
                        error <- c('the number of the two spatial coordinates does not match\n')
                        return(list(error=error))
                      }
                    if(length(coordx)!=dimdata[2])
                      {
                        error <- c('the number of coordinates does not match with the number of spatial observations\n')
                        return(list(error=error))
                      }
                  }
                if(dimdata[1]!=length(coordt))
                  {
                    error <- c('the time coordinate does not match with the number of the matrix rows\n')
                    return(list(error=error))
                  }
              } # END irregular grid
          } # END b) one realisation of a spatial-temporal random field
      }
      # END - check the format of the inserted coordinates and dataset
    }
    # END check fitting
    # START check simulation
    if(fcall=="Simulation"){
        if(type=="Tapering"){
        #if(is.null(taper) || is.null(maxdist) || is.null(maxtime)) {
        #  error <- 'tapering need a taper correlation model and/or a compact support\n'
        #  return(list(error=error))}
        if(!taper %in% c("Wendland1","Wendland2","Wendland3",
                         "Wendland1_Wendland1","Wendland1_Wendland2","Wendland1_Wendland3",
                         "Wendland2_Wendland1","Wendland2_Wendland2","Wendland2_Wendland3",
                         "Wendland3_Wendland1","Wendland3_Wendland2","Wendland3_Wendland3",
                         "qt_time","qt_space")){
            error <- 'insert a correct name for the taper correlation model\n'
            return(list(error=error))}
            if(!is.null(maxdist)){
            error <- "insert a positive numeric value for the maximum spatial distance\n"
            if(!is.numeric(maxdist)) return(list(error=error))
            else if(maxdist<0) return(list(error=error))}
            if(!is.null(maxtime)){
            error <- "insert a positive numeric value for the maximum temporal distance\n"
            if(!is.numeric(maxtime)) return(list(error=error))
            else if(maxtime<0) return(list(error=error))}
            if(!is.null(tapsep)){
            error <- "separability parameter of spacetime taper must be between 0  and 1\n"
            if(!is.numeric(tapsep)) return(list(error=error))
            else if(tapsep<0||tapsep>1) return(list(error=error))}
            }

        if(is.null(param) || !is.list(param)){
            error <- 'insert the parameters as a list\n'
            return(list(error=error))}

        if(length(param)!=length(c(unique(c(NuisanceParam("Gaussian"),NuisanceParam(model))),CorrelationParam(corrmodel)))){
            error <- "some parameters are missing or does not match with the declared model\n"
            return(list(error=error))}

        if(!all(names(param) %in% c(unique(c(NuisanceParam("Gaussian"),NuisanceParam(model))), CorrelationParam(corrmodel)))){
            error <- 'some names of the parameters are not correct\n'
            return(list(error=error))}

        if(!CheckParamRange(unlist(param))){
            error <- 'some parameters are out of the range\n'
            return(list(error=error))}

        if(model %in% c("BrowResn","ExtT") && !is.null(numblock)){
            error <- "insert the number of observation in each block\n"
            if(!is.numeric(numblock)) return(list(error=error))
            else if(numblock<0) return(list(error=error))}

        if(model=="BrowResn"){
            if(!corrmodel %in% c("exponential","gauss","gencauchy","stable","stable_stable")){
                error <- "the correlation model is not adequate for the Brown-Resnick process\n"
                 return(list(error=error))}
            if(corrmodel=="gencauchy" && param["power1"]!=2){
                error <- "the parameter power1 need to be equal to 2 for the Brown-Resnick process\n"
                 return(list(error=error))}}

    }
    # END check simulation
  }
 else{    ### checking input for kriging
   if(class(corrmodel)!="CovMat")   {
           error <-  "covmatrix is not a object of class covmatrix\n"
           return(list(error=error))}

    if(missing(coordx)) {
    error <- "spatial locations must be a matrix of dimension 2\n"
    return(list(error=error))}
    else     {
    if(is.vector(coordx)&&!length(coordx)==2)    {
           error <- "spatial locations must be a vector of dimension 2\n"
           return(list(error=error))}
    if(is.matrix(coordx)&&!ncol(coordx)==2)       {
           error <- "spatial locations must be  a matrix of dimension 2\n"
           return(list(error=error))}
    }
   if((!is.null(coordt)&&!is.numeric(coordt))  ){
           error <- "time  must be a vector\n"
           return(list(error=error))}
   if(!type %in% c("simple","ordinary","Simple","Ordinary")){
           error <-"kriging type can be  simple or Ordinary\n"
   return(list(error=error))}
   if(missing(data) || !is.numeric(data)){
           error <- "insert a numeric vector of data\n"
           return(list(error=error))}
   if(length(data)!=nrow(corrmodel$covmatrix)) {
           error <-  "data dimension does not correspond to the covmatrix dimension\n"
           return(list(error=error))}
        }
      # END check kriging
  }

CheckLikelihood <- function(likelihood)
  {
    CheckLikelihood <- switch(likelihood,
                              None=0,
                              Conditional=1,
                              Full=2,
                              Marginal=3)
    return(CheckLikelihood)
  }

CheckModel <- function(model)
  {
    CheckModel <- switch(model,
                         None=0,
                         Gaussian=1,
                         BinaryGauss=2,
                         BrowResn=3,
                         ExtGauss=4,
                         ExtT=5)
    return(CheckModel)
  }

CheckVarType <- function(type)
  {
    CheckVarType <- switch(type,
                           Sampling=1,
                           SubSamp=2,
                           Theoretical=3)
    return(CheckVarType)
  }

CheckType <- function(type)
  {
    CheckType <- switch(type,
                        Difference=1,
                        Pairwise=2,
                        Restricted=3,
                        Standard=4,
                        Tapering=5,
                        WLeastSquare=6)
    return(CheckType)
  }

CorrelationParam <- function(corrmodel)
  {
    param <- NULL
    # Cauchy correlation model:
    if(corrmodel=='cauchy'){
      param <- c('power2', 'scale')
      return(param)}
    # Exponential and Gaussian correlation models:
    if(corrmodel=='exponential' || corrmodel=='gauss' || corrmodel=='spherical'){
      param <- c('scale')
      return(param)}
    # Generalised Cauchy correlation model:
    if(corrmodel=='gencauchy'){
      param <- c('power1', 'power2','scale')
      return(param)}
    # Stable correlation model:
    if(corrmodel=='stable'){
      param <- c('power', 'scale')
      return(param)}
    # Whittle-Matern correlation model:
    if(corrmodel=='matern'){
      param <- c('scale', 'smooth')
      return(param)}
    # Non-separable spatial-temporal correlations:
    # Gneiting or Porcu model:
    if(corrmodel=='gneiting'||corrmodel=='porcu'||corrmodel=='porcu1'||corrmodel=='porcu2'||corrmodel=='gneiting_GC'){
      param <- c('power_s', 'power_t','scale_s','scale_t','sep')
      return(param)}
    # Iaco-Cesare model:
    if(corrmodel=='iacocesare'){
      param <- c('power2','power_s', 'power_t','scale_s','scale_t')
      return(param)}
    # Stein model:
    if(corrmodel=='stein'){
      namesparam <- c('power_t','scale_s','scale_t','smooth_s')
      return(param)}
    # Separable spatial-temporal correlations:
    # Exponential-exponential and exponential-Gaussian models:
    if(corrmodel=='exp_exp'||corrmodel=='exp_gauss'){
      param <- c('scale_s','scale_t')
      return(param)}
    # Exponential-Cauchy model:
     if(corrmodel=='exp_cauchy'){
       param <- c('power2','scale_s','scale_t')
       return(param)}
    # Whittle-Matern-exponential model:
     if(corrmodel=='matern_exp'){
       param <- c('scale_s','scale_t','smooth_s')
       return(param)}
    # Whittle-Matern-Cauchy model:
     if(corrmodel=='matern_cauchy'){
       param <- c('power2','scale_s','scale_t','smooth_s')
       return(param)}
    # Stable-Stable model:
     if(corrmodel=='stab_stab'){
       param <- c('power_s','power_t','scale_s','scale_t')
       return(param)}

    return(param)
  }

NuisanceParam <- function(model)
{
  param <- NULL
  # Gaussian random field:
  if(model=='Gaussian' || model=='BinaryGauss'){
    param <- c('mean', 'nugget', 'sill')
    return(param)}
  # Max-stable random field (Extremal Gaussian):
  if(model=='ExtGauss'){
    param <- c('sill')
    return(param)}
  # Max-stable random field (Brown Resnick):
  if(model=='BrowResn'){
    param <- c('sill')
    return(param)}
  # Max-stable random field (Extremal T):
  if(model=='ExtT'){
    param <- c('df', 'sill')
    return(param)}
  return(param)
}


InitParam <- function(coordx, coordy, coordt, corrmodel, data, distance, fcall, fixed, grid,
                      likelihood, margins, maxdist, maxtime, model, numblock, param, parscale,
                      paramrange, replicates, start, taper, tapsep, threshold, type,
                      typereal, varest, vartype, weighted, winconst, winstp)
{
    ### START Includes internal functions:
    # Check if the correlation is spatial or spatial-temporal
    CheckSpaceTime <- function(corrmodel)
    {
        CheckSpaceTime <- NULL
        if(corrmodel >= 1 & corrmodel <= 20) CheckSpaceTime <- FALSE
        else CheckSpaceTime <- TRUE
        return(CheckSpaceTime)
    }
    # Check the type of distances
    CheckDistance<- function(distance)
    {
        CheckDistance <- NULL
        CheckDistance <- switch(distance,
                                eucl=0,
                                Eucl=0,
                                chor=1,
                                Chor=1,
                                geod=2,
                                Geod=2,
                                proj=3,
                                Proj=3)
        return(CheckDistance)
    }
    # Determines the range of the parameters for a given correlation
    SetRangeParam <- function(namesparam, numparam)
    {
        low <- 1e-12
        lower <- NULL
        upper <- NULL
        # Check for the param set:
        for(i in 1:numparam){
            if(namesparam[i]=='mean'){
                lower <- c(lower, -Inf)
                upper <- c(upper, Inf)}
            if(namesparam[i]=='nugget'){
                lower <- c(lower, 0)
                upper <- c(upper, Inf)}
            if(namesparam[i]=='power'){
                lower <- c(lower, low)
                upper <- c(upper, 2)}
            if(namesparam[i]=='power_s'){
                lower <- c(lower, low)
                upper <- c(upper, 2)}
            if(namesparam[i]=='power_t'){
                lower <- c(lower, low)
                upper <- c(upper, 2)}
            if(namesparam[i]=='power1'){
                lower <- c(lower, low)
                upper <- c(upper, 2)}
            if(namesparam[i]=='power2'){
                lower <- c(lower, low)
                upper <- c(upper, Inf)}
            if(namesparam[i]=='scale'){
                lower <- c(lower, low)
                upper <- c(upper, Inf)}
            if(namesparam[i]=='scale_s'){
                lower <- c(lower, low)
                upper <- c(upper, Inf)}
            if(namesparam[i]=='scale_t'){
                lower <- c(lower, low)
                upper <- c(upper, Inf)}
            if(namesparam[i]=='sep'){
                lower <- c(lower, low)
                upper <- c(upper, 1)}
            if(namesparam[i]=='sill'){
                lower <- c(lower, low)
                upper <- c(upper, Inf)}
            if(namesparam[i]=='smooth'){
                lower <- c(lower, low)
                upper <- c(upper, Inf)}}
        return(list(lower=lower, upper=upper))
    }
    ### END Includes internal functions
    ### Set returning variables and initialize the model parameters:
    # Initialises the starting and fixed parameters' names
    error <- NULL
    namesfixed <- namesstart <- namessim <- NULL
    numfixed <- numstart <- 0
    # Set the model, likelihood, correlation and the nuisance parameters:
    namesnuis <- NuisanceParam(model)
    model <- CheckModel(model)
    flagnuis <- NULL
    namescorr <- CorrelationParam(corrmodel)
    numparamcorr <- length(namescorr)
    paramcorr <- rep(1, numparamcorr)
    names(paramcorr) <- namescorr
    flagcorr <- NULL
    #paramrange <- list(lower=NULL, upper=NULL)
    # Set the correlation and if is a space-time random field:
    corrmodel <- CheckCorrModel(corrmodel)
    spacetime <- CheckSpaceTime(corrmodel)
    ### START settings the data structure:
    # set the coordinates sizes:
    if(is.null(coordy)){coordy <- coordx[,2]
                        coordx <- coordx[,1]}
    # checks is the data are on regular grid:
    if(grid) {numcoordx <- length(coordx)
              numcoordy <- length(coordy)
              numcoord <- numcoordx*numcoordy}
    else numcoord <- numcoordx <- numcoordy <- length(coordx)
    # initialize tapering variables:
    tapering<-ia<-idx<-ja<-integer(1)
    nozero<-NULL
    tapmodel=NULL
    cutoff <- FALSE
    distance<-CheckDistance(distance)
    ### END settings the data structure
    # START code for the simulation procedure
    if(fcall=="Fitting"){
        ### Parameters' settings:
        nuisance <- NULL
        likelihood <- CheckLikelihood(likelihood)
        vartype <- CheckVarType(vartype)
        type <- CheckType(type)
        if(model==1){ # Gaussian random field:
           mu <- mean(data)
           if(any(type==c(1, 3, 6)))# Checks the type of likelihood
           if(is.list(fixed)) fixed$mean <- mu# Fixs the mean
           else fixed <- list(mean=mu)
           nuisance <- c(mu, 0, var(c(data)))
           if(likelihood==2 && CheckType(typereal)==5) tapering <- 1}
        if(model==2){ # Binary Gaussian random field:
            # check the threshould:
            if(is.null(threshold) || !is.numeric(threshold))
                threshold <- 0
            p <- mean(data)
            mu <- threshold+qnorm(p)
            nuisance <- c(mu, 0, 1)
            if(!is.null(start$nugget))
                if(length(start)>1) start<-start[!names(start)%in%"nugget"]
                else start<-NULL
            if(is.list(fixed)) fixed$nugget<-0# Fixs the nugget
            else fixed<-list(nugget=0)
            #set the nugget in case the sill is fixed
            if(!is.null(fixed$sill)) fixed$nugget <- 1-fixed$sill}
        if(model>2){ # Max-stable random field:
            if(model==3){# Checks if its the Brown-Resnick model
                if(is.list(fixed)) fixed$sill <- 1 # Fixs the sill
                else fixed <- list(sill=1)
                if(corrmodel==4) fixed$power1 <- 2}
            # Cheks if its the Extremal-t model
            if(model==5) nuisance <- c(nuisance,1)
            nuisance <- c(nuisance,0.5)
            if(margins=="Gev") data <- Dist2Dist(data)}
        # Update the parameter verctor:
        names(nuisance) <- namesnuis
        namesparam <- sort(c(namescorr, namesnuis))
        param <- c(nuisance, paramcorr)
        param <- param[namesparam]
        numparam <- length(param)
        flag <- rep(1, numparam)
        namesflag <- namesparam
        names(flag) <- namesflag
        # Update the parameters with fixed values:
        if(!is.null(fixed)){
            fixed <- unlist(fixed)
            namesfixed <- names(fixed)
            numfixed <- length(namesfixed)
            if(numfixed==numparam){
                error <- 'the are not parameters left to estimate\n'
                return(list(error=error))}
            flag[pmatch(namesfixed, namesflag)] <- 0
            param <- param[-pmatch(namesfixed, namesparam)]
            numparamcorr <- numparamcorr-sum(namesfixed %in% namescorr)
            namesparam <- names(param)
            numparam <- length(param)}
        flagcorr <- flag[namescorr]
        flagnuis <- flag[namesnuis]
        # Update the parameters with starting values:
        if(!is.null(start)){
            start <- unlist(start)
            namesstart <- names(start)
            if(any(type == c(1, 3, 6)))
                if(any(namesstart == 'mean'))
                    start <- start[!namesstart == 'mean']
            namesstart <- names(start)
            numstart <- length(start)
            param[pmatch(namesstart,namesparam)] <- start}
        ### set the scale of the parameters:
        # Insert here!
        # set the range of the parameters if its the case
        if(paramrange) paramrange <- SetRangeParam(namesparam, numparam)
        else paramrange <- list(lower=NULL, upper=NULL)
        ### If the case set the sub-sampling parameters to the default values
        if(missing(winconst) || !is.numeric(winconst)) winconst <- 0
        if(missing(winstp) || !is.numeric(winstp)) winstp <- 0
        ### Set the data format:
        if(spacetime){ # set the number of temporal realisations:
            numtime <- length(coordt)
            if(grid) # if the data are in regular grid:
                if(replicates>1){ # checks if there are iid replicates:
                    dim(data) <- c(numcoord, numtime, replicates)
                    data <- aperm(data, c(2,1,3))}
                else # if there are not iid replicates:
                    data <- matrix(data, ncol=numcoord, nrow=numtime, byrow=TRUE)
                if(typereal=="Tapering"){
                tapering<-1
                idx<-integer((numcoord*numtime)^2)
                ja<-integer((numcoord*numtime)^2)
                ia<-integer(numcoord*numtime+1)
                tapmodel<-CheckCorrModel(taper)
                }}
        else{
            numtime <- 1
            coordt <- 0
            if(grid) data <- matrix(data, ncol=numcoord, nrow=replicates, byrow=TRUE)
            else data <- matrix(data, ncol=numcoord, nrow=replicates)
            if(typereal=="Tapering"){
                tapering<-1
                idx<-integer((numcoord*numtime)^2)
                ja<-integer((numcoord*numtime)^2)
                ia<-integer(numcoord*numtime+1)
                tapmodel<-CheckCorrModel(taper)
                }}
    }
    # END code for the fitting procedure
    # START code for the simulation procedure
    if(fcall=="Simulation"){
        namesnuis <- sort(unique(c(namesnuis,NuisanceParam("Gaussian"))))
        param <- unlist(param)
        numparam <- length(param)
        namesparam <- names(param)
        if(model %in% c(3,5))
            if(is.null(numblock)) numblock <- as.integer(500)
            else numblock <- as.integer(numblock)
        if(model==3){
            param["df"] <- 0
            if(corrmodel==2) param["scale"] <- 2*log(numblock)*param["scale"]
            if(corrmodel %in% c(3,4)) param["scale"] <- 2*sqrt(log(numblock))*param["scale"]
            if(corrmodel==6) param["scale"] <- 2*(log(numblock))^(1/param["power"])*param["scale"]}
        namessim <- c("mean","sill","nugget","scale",namescorr[!namescorr=="scale"])
        if(spacetime) numtime <- length(coordt)
        else {numtime <- 1; coordt <- 0}
        if(typereal=="Tapering"&&type=="Tapering"){
                tapering<-1
                idx<-integer((numcoord*numtime)^2)
                ja<-integer((numcoord*numtime)^2)
                ia<-integer(numcoord*numtime+1)
                tapmodel<-CheckCorrModel(taper)
                }
        if(model==2){ # Binary Gaussian random field:
            # check the threshould:
            if(is.null(threshold) || !is.numeric(threshold))
                threshold <- 0 }
        }
    # END code for the simulation procedure
    ### Compute the spatial and spatial-temporal distances:
    numpairs <- integer(1)
    srange <- double(1)
    trange <- double(1)
    if(is.null(maxdist)) srange<-c(srange,double(1)) else {srange<-c(srange,as.double(maxdist))}                # cutoff<-TRUE
    if(is.null(maxtime)) trange<-c(trange,double(1)) else {trange<-c(trange,as.double(maxtime))}                # cutoff<-TRUE
    isinit <- as.integer(1)
    SG=.C('SetGlobalVar',as.double(coordx),as.double(coordy),as.double(coordt),as.integer(grid),ia=ia,idx=idx,
       isinit=isinit,ja=ja,as.integer(numcoord),as.integer(numcoordx),as.integer(numcoordy),numpairs=numpairs,as.integer(replicates),
       srange=srange, as.double(tapsep), as.integer(numtime),trange=trange,as.integer(tapering),as.integer(tapmodel),
       as.integer(distance),as.integer(weighted),PACKAGE='CompRandFld', DUP=TRUE, NAOK=TRUE)
     srange<-SG$srange
     trange<-SG$trange
     ja<-SG$ja;ia<-SG$ia;idx<-SG$idx
     isinit<-SG$isinit
     numpairs<-SG$numpairs
     if(tapering){ nozero<-numpairs/(numcoord*numtime)^2
                   idx <- idx[1:numpairs]
                   ja <- ja[1:numpairs]}
    ### Returned list of objects:
    return(list(coordx=coordx,coordy=coordy,coordt=coordt,corrmodel=corrmodel,data=data,distance=distance,
                error=error,flagcorr=flagcorr,flagnuis=flagnuis,fixed=fixed,likelihood=likelihood,
                lower=paramrange$lower,model=model,namescorr=namescorr,namesfixed=namesfixed,
                namesnuis=namesnuis,namesparam=namesparam,namessim=namessim,namesstart=namesstart,
                numblock=numblock,numcoord=numcoord,numcoordx=numcoordx,numcoordy=numcoordy,
                numfixed=numfixed,numpairs=numpairs,numparam=numparam,numparamcorr=numparamcorr,
                numrep=replicates,numstart=numstart,numtime=numtime,param=param,setup=list(                ## setup is a list
                ia=ia,idx=idx,ja=ja,nozero=nozero,tapmodel=tapmodel,tapsep=tapsep),                              ## with tapered matrix informations
                spacetime=spacetime,srange=srange,start=start,upper=paramrange$upper,type=type,
                threshold=threshold,trange=trange,vartype=vartype,winconst=winconst,winstp=winstp))
}
