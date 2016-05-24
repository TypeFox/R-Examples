get.mod <-
structure(function # Copy objects between R environments
                   (modelname=ls(FlexParamCurve:::FPCEnv,pattern=".lis"),
                    ### a list of object names
                    from.envir = FlexParamCurve:::FPCEnv,
                    ### R environment currently containing the object(s)
                    to.envir=.GlobalEnv,
                    ### destination R environment to copy the object(s) to
                    write.mod = FALSE
                    ### logical specifying if single models should be assigned or simply returned
                    ) {
    ##description<< Function to copy objects between R environments                
    ##details<< All arguments are optional. With defaults, this function copies any \eqn{nlsList} models
    ## from the FlexParamCurve working environment to the Global Environment. However, user could use 
    ## this function to move any objects between any environments. 
    ##
    ## Default behavior is to assign models to an environment if more than 1 modelname is provided but to
    ## simply return the model from the function if only 1 modelname is given. Notes are printed to the
    ## screen to detail any models moved or any errors encountered.
    ##
if(length(modelname)<1){
}else{
cat("target environment: ")
print(to.envir)
if(length(modelname)>1 | write.mod == TRUE){
  for(i in 1:length(modelname)){
    if(exists(modelname[i] ,envir=from.envir) == F) {
	cat (paste("Object ", modelname[i], " not found in specified environment.
	Note: object name argument should be a character or character vector",sep="\""))
    }else{	
   assign(modelname[i], get(modelname[i],envir = from.envir), envir= to.envir)
   cat( paste(modelname[i],"successfully transfered to this environment", sep = " ")
   , fill=T, labels= as.character(i) )
    }			       }
 }else{
if(exists(modelname[1] ,envir=from.envir) == F) stop ("Object not found in specified environment.
   Note: object name argument should be a character or character vector")
cat(paste("value returned for ",modelname, ":\n ",sep="\""))
return(get(modelname[1] ,envir=from.envir))
 }
}    
}
##value<< If only 1 modelname is provided, the contents of the object is returned. If more
## more than 1 modelname is provided or if write.mod is FALSE then the object(s) will be assigned
## to the environment and no value is returned.
##
##note<< The default function works by detecting the suffix .lis rather than object class, so will
## only return models with this suffix, not necessarily all \eqn{nlsList} models if they have
## different suffixes.
, ex = function(){
   #transfer all \eqn{nlsList} models from the FlexParamCurve working environmment (FPCEnv) 
   #to the Global Environment. Note: unless \code{\link{pn.mod.compare}} or 
   #\code{\link{pn.modselect.step}} have been run, in which case this is default
   #1. subset data object (only 3 individuals) to expediate model selection
   subdata <- subset(posneg.data, as.numeric(row.names (posneg.data) ) < 40)
   #2. run model selection in FPCEnv using \code{\link{pn.mod.compare}}. Only two models (#1 and #5)
   #specified to be run here to reduce processing time. see \code{\link{pn.mod.compare}}
   modseltable <- pn.mod.compare(subdata$age, subdata$mass,
      subdata$id, existing = FALSE, pn.options = "myoptions", mod.subset = c(1,5))
   #3. retrieve models from FlexParamCurve working environmment
   get.mod()
   #transfer an options file called myoptions from FPCEnv to the Global Environment
   #note data are forced to fit a monotonic curve in this example
   modpar(logist.data$age, logist.data$mass, pn.options = "myoptions.1", force4par = TRUE, 
   Envir = FlexParamCurve:::FPCEnv)
   get.mod(modelname = "myoptions.1", write.mod = TRUE)
}
)
