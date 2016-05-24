#
#***************** 0 Objet methodes generiques ********************
#
# A) Affecters 
# B) Allocators
# C) Methode 

# if (!isGeneric("summary")){
#  setGeneric("summary", function(object, ...) standardGeneric("summary"))
# }

##### A) Selecteurs ############################################

# Carto3D MRIaggr
setGeneric(name = "selectContrast", 
           def = function(object, ...){standardGeneric("selectContrast")}
)

# MRIaggr
setGeneric(name = "selectClinic", 
           def = function(object, ...){standardGeneric("selectClinic")}
)

# MRIaggr
setGeneric(name = "selectCoords", 
           def = function(object, ...){standardGeneric("selectCoords")}
)

# Carto3D MRIaggr  
setGeneric(name = "selectDefault_value", 
           def = function(object, ...){standardGeneric("selectDefault_value")}
)

# MRIaggr 
setGeneric(name = "selectDescStats", 
           def = function(object, ...){standardGeneric("selectDescStats")}
)

# Carto3D MRIaggr 
setGeneric(name = "selectFieldDim", 
           def = function(object, ...){standardGeneric("selectFieldDim")}
)

# MRIaggr 
setGeneric(name = "selectHistory", 
           def = function(object, ...){standardGeneric("selectHistory")}
)

# MRIaggr 
setGeneric(name = "selectHemispheres", 
           def = function(object, ...){standardGeneric("selectHemispheres")}
)

# Carto3D MRIaggr 
setGeneric(name = "selectIdentifier", 
           def = function(object, ...){standardGeneric("selectIdentifier")}
)

# MRIaggr 
setGeneric(name = "selectMidplane", 
           def = function(object, ...){standardGeneric("selectMidplane")}
)

# MRIaggr
setGeneric(name = "selectN", 
           def = function(object, ...){standardGeneric("selectN")}
)

# MRIaggr
setGeneric(name = "selectNormalization", 
           def = function(object, ...){standardGeneric("selectNormalization")}
)

# Carto3D MRIaggr 
setGeneric(name = "selectParameter", 
           def = function(object, ...){standardGeneric("selectParameter")}
)

# MRIaggr 
setGeneric(name = "selectTable", 
           def = function(object, ...){standardGeneric("selectTable")}
)

# Carto3D MRIaggr 
setGeneric(name = "selectVoxelDim", 
           def = function(object, ...){standardGeneric("selectVoxelDim")}
)

# MRIaggr 
setGeneric(name = "selectW", 
           def = function(object, ...){standardGeneric("selectW")}
)

##### B) allocants ############################################


# MRIaggr 
setGeneric(name = "allocContrast<-", 
           def = function(object, param = NULL, default_value = NULL, overwrite = FALSE, verbose = TRUE, value){standardGeneric("allocContrast<-")}
)

# MRIaggr
setGeneric(name = "allocClinic<-", 
           def = function(object, add = FALSE, overwrite = FALSE, verbose = TRUE, value){standardGeneric("allocClinic<-")}
)

# MRIaggr 
setGeneric(name = "allocDescStats<-", 
           def = function(object, name, overwrite = FALSE, verbose = TRUE, value){standardGeneric("allocDescStats<-")}
)

# MRIaggr
setGeneric(name = "allocHemisphere<-", 
           def = function(object, overwrite = FALSE, verbose = TRUE, value){standardGeneric("allocHemisphere<-")}
)

# MRIaggr
setGeneric(name = "allocNormalization<-", 
           def = function(object, overwrite = FALSE, verbose = TRUE, value){standardGeneric("allocNormalization<-")}
)

# MRIaggr 
setGeneric(name = "allocTable<-", 
           def = function(object, type, overwrite = FALSE, verbose = TRUE, value){standardGeneric("allocTable<-")}
)

# MRIaggr 
setGeneric(name = "allocW<-", 
           def = function(object, type, overwrite = FALSE, verbose = TRUE, value){standardGeneric("allocW<-")}
)

# MRIaggr
setGeneric(name = "supprContrast<-", 
           def = function(object, verbose = TRUE, value){standardGeneric("supprContrast<-")}
)

# MRIaggr
setGeneric(name = "supprDescStats<-", 
           def = function(object, verbose = TRUE, value){standardGeneric("supprDescStats<-")}
)

##### C) Methodes ############################################

#### calc. ####

# MRIaggr
setGeneric(name  = "calcBrainMask", 
           def = function(object, ...){
             res <-standardGeneric("calcBrainMask")
             
             if(res$update.object == TRUE){
               
               # alloc
               nom_object <- as.character(substitute(object))
               allocContrast(object, param = "mask", overwrite = res$overwrite, verbose = res$verbose) <- res$res$best_group
               
               # update history
               object@history <- c(object@history, 
                                   list(calcBrainMask = list(call = match.call(), date = date(), mask_name = res$res$mask_name))
               )
               
               # update object
               eval(parse(text = paste(
                 "assign(\"", nom_object, "\", value = object, 
                envir = .GlobalEnv)", 
                 sep = "")))             
             }
             
             return(invisible(res$res))
             
           }            
)

# MRIaggr
setGeneric(name = "calcContralateral", 
           def = function(object, ...){
             
             res <- standardGeneric("calcContralateral")
             
             if(res$update.object == TRUE){
               
               # alloc
               nom_object <- as.character(substitute(object))
               newdata <- data.frame(matrix(NA, nrow = selectN(object), ncol = ncol(res$data)))
               names(newdata) <- names(res$data)   
               
               newdata[res$data[,"index"],] <- res$data
               newdata <- newdata[,names(newdata) %in% c("index", "i_hemisphere", "j_hemisphere", "hemisphere") == FALSE]
               
               allocContrast(object, overwrite = res$overwrite, verbose = res$verbose) <- newdata
               
               # update history
               object@history <- c(object@history, 
                                   list(calcContralateral = list(call = match.call(), date = date()))
               )
               
               # update object
               eval(parse(text = paste(
                 "assign(\"", nom_object, "\", value = object, 
                envir = .GlobalEnv)", 
                 sep = "")))         
             }
             
             res$update.object <- NULL
             res$overwrite <- NULL
             res$verbose <- NULL
             
             return(invisible(res))
             
           }
)

# MRIaggr
setGeneric(name  = "calcDistMask", 
           def = function(object, ...){
             res <- standardGeneric("calcDistMask")
             
             if(res$update.object == TRUE){
               
               # alloc
               nom_object <- as.character(substitute(object))
               allocContrast(object, overwrite = res$overwrite, verbose = res$verbose) <- res$res
               
               # update history
               object@history <- c(object@history, 
                                   list(calcDistMask = list(call = match.call(), date = date()))
               )
               
               # update object
               eval(parse(text = paste(
                 "assign(\"", nom_object, "\", value = object, 
                envir = .GlobalEnv)", 
                 sep = "")))             
             }
             
             return(invisible(res$res))
             
           }            
)


# MRIaggr
setGeneric(name = "calcDistTissues", 
           def = function(object, ...){standardGeneric("calcDistTissues")}
)

# MRIaggr
setGeneric(name = "calcFilter", 
           def = function(object, ...){
             res <- standardGeneric("calcFilter")
             
             if(res$update.object == TRUE){
               
               # alloc
               nom_object <- as.character(substitute(object))
               param <- setdiff(names(res$res), c("i", "j", "k"))              
               
               allocContrast(object, overwrite = res$overwrite, verbose = res$verbose) <- res$res[,param, drop = FALSE]
               
               # update history
               object@history <- c(object@history, 
                                   list(calcFilter = list(call = match.call(), date = date()))
               )
               
               # update object
               eval(parse(text = paste(
                 "assign(\"", nom_object, "\", value = object, 
                    envir = .GlobalEnv)", 
                 sep = "")))                 
             }
             
             res$update.object <- NULL
             res$overwrite <- NULL
             res$verbose <- NULL
             
             return(invisible(res))
             
           }
)

# MRIaggr
setGeneric(name = "calcGroupsMask", 
           def = function(object, ...)
           {  res <- standardGeneric("calcGroupsMask")
           
           if(res$update.object == TRUE){
             # alloc
             nom_object <- as.character(substitute(object))
             allocDescStats(object, name = "GroupsLesion", overwrite = res$overwrite, verbose = res$verbose) <- lapply(res$res, function(x){x$group_size})
             
             # update history
             object@history <- c(object@history, 
                                 list(calcGroupsMask = list(call = match.call(), date = date()))
             )
             
             # update object
             eval(parse(text = paste(
               "assign(\"", nom_object, "\", value = object, 
                envir = .GlobalEnv)", 
               sep = "")))  
           }
           
           return(invisible(res$res))
           
           }
)

# MRIaggr
setGeneric(name = "calcHemisphere", 
           def = function(object, ...){
             res <- standardGeneric("calcHemisphere")
             
             if(res$update.object == TRUE){
               
               # alloc
               nom_object <- as.character(substitute(object))
               allocHemisphere(object, overwrite = res$overwrite, verbose = res$verbose) <- list(midplane = res$res$midplane,                                                                                            
                                                                                                 data = res$res$data)
               
               if(!is.null(res$res$hemispheres)){
                 allocHemisphere(object, overwrite = res$overwrite, verbose = res$verbose) <- list(hemispheres = res$res$hemispheres)
               }
               
               # update history
               object@history <- c(object@history, 
                                   list(calcHemisphere = list(call = match.call(), date = date(), optimum = res$res$optimum[,c("position_i", "position_j", "angle_rad")]))
               )
               
               # update object
               eval(parse(text = paste(
                 "assign(\"", nom_object, "\", value = object, 
                envir = .GlobalEnv)", 
                 sep = "")))  
             }  
             
             return(invisible(res$res))
             
           }
)

# MRIaggr
setGeneric(name = "calcROCthreshold", 
           def = function(object, ...){ 
             res <- standardGeneric("calcROCthreshold")
             
             if(res$update.object == TRUE){
               
               # alloc
               nom_object <- as.character(substitute(object))
               allocDescStats(object, name = "Mask_threshold", overwrite = res$overwrite, verbose = res$verbose) <- res$res
               
               eval(parse(text = paste(
                 "assign(\"", nom_object, "\", value = object, 
                 envir = .GlobalEnv)", 
                 sep = "")))
               
               # update history
               object@history <- c(object@history, 
                                   list(calcROCthreshold = list(call = match.call(), date = date()))
               )
               
               # update object
               res$update.object <- NULL
               res$overwrite <- NULL
             }
             return(invisible(res$res))
           }
)

# MRIaggr
setGeneric(name = "calcNormalization", 
           def = function(object, ...){
             
             res <- standardGeneric("calcNormalization")
             
             if(res$update.object == TRUE){
               
               # alloc
               nom_object <- as.character(substitute(object))
               allocNormalization(object, overwrite = res$overwrite, verbose = res$verbose) <- res$res
               
               # update history
               object@history <- c(object@history, 
                                   list(calcNormalization = list(call = match.call(), date = date()))
               )
               
               # update object
               eval(parse(text = paste(
                 "assign(\"", nom_object, "\", value = object, 
                 envir = .GlobalEnv)", 
                 sep = "")))             
             }
             
             return(invisible(res$res))
             
           }
)

# MRIaggr
setGeneric(name = "calcRegionalContrast", 
           def = function(object, ...){
             res <- standardGeneric("calcRegionalContrast")
             
             if(!is.list(res)){
               return(res)
             }else{
               if(res$update.object == TRUE){
                 
                 # alloc
                 nom_object <- as.character(substitute(object))
                 allocContrast(object, overwrite = res$overwrite, verbose = res$verbose) <- res$res
                 
                 # update history
                 object@history <- c(object@history, 
                                     list(calcRegionalContrast = list(call = match.call(), date = date()))
                 )
                 
                 # update object
                 eval(parse(text = paste(
                   "assign(\"", nom_object, "\", value = object, 
                   envir = .GlobalEnv)", 
                   sep = "")))
               }
               
               return(invisible(res$res))
             }
             
           }
)

# MRIaggr
setGeneric(name = "calcSmoothMask", 
           def = function(object, ...){
             res <- standardGeneric("calcSmoothMask")
             
             if(res$update.object == TRUE){
               
               # alloc
               nom_object <- as.character(substitute(object))
               allocContrast(object, param = "mask", overwrite = res$overwrite, verbose = res$verbose) <- res$res$mask
               
               # update history
               object@history <- c(object@history, 
                                   list(calcSmoothMask = list(call = match.call(), date = date()))
               )
               
               # update object
               eval(parse(text = paste(
                 "assign(\"", nom_object, "\", value = object, 
                 envir = .GlobalEnv)", 
                 sep = "")))   
             }
             
             return(invisible(res$res))
             
           }
)

# MRIaggr
setGeneric(name = "calcTableHypoReperf", 
           def = function(object, ...){          
             res <- standardGeneric("calcTableHypoReperf")
             
             if(res$update.object == TRUE){
               
               # alloc
               nom_object <- as.character(substitute(object))
               
               if("volume_hypo" %in% names(res$res)){
                 allocTable(object, type = "hypoperfusion", overwrite = res$overwrite, verbose = res$verbose) <- res$res$volume_hypo
               }
               
               if("volume_reperf" %in% names(res$res)){
                 allocTable(object, type = "reperfusion", overwrite = res$overwrite, verbose = res$verbose) <- res$res$volume_reperf
               }
               
               if("pixel" %in% names(res$res)){
                 nom_param <- names(res$res$pixel)
                 param.reperf_pc <- grep(pattern = "reperf_pc", nom_param, value = TRUE)
                 param.reperf <- grep(pattern = "reperf", nom_param[nom_param %in% param.reperf_pc == FALSE], value = TRUE)
                 param.deperf_pc <- grep(pattern = "deperf_pc", nom_param, value = TRUE)
                 param.deperf <- grep(pattern = "deperf", nom_param[nom_param %in% param.deperf_pc == FALSE], value = TRUE)
                 param.shift <- grep(pattern = "shift", nom_param, value = TRUE)
                 
                 eval(parse(text = paste(
                   "nom_param <- c(\"i\", \"j\", \"k\", ", paste(paste("param.", res$param.update, sep = ""), collapse = ", "), ")", 
                   sep = "")))
                 
                 allocContrast(object, overwrite = res$overwrite, verbose = res$verbose) <- res$res$pixel[,nom_param]
               }
               
               # update history
               object@history <- c(object@history, 
                                   list(calcTableHypoReperf = list(call = match.call(), date = date()))
               )
               
               # update object
               eval(parse(text = paste(
                 "assign(\"", nom_object, "\", value = object, 
                 envir = .GlobalEnv)", 
                 sep = "")))            
             }
             
             return(invisible(res$res))
             
           }
)

# MRIaggr
setGeneric(name = "calcTableLesion", 
           def = function(object, ...){ 
             res <- standardGeneric("calcTableLesion")
             
             if(res$update.object == TRUE){
               
               # alloc
               nom_object <- as.character(substitute(object))
               allocTable(object, type = "lesion", overwrite = res$overwrite, verbose = res$verbose) <- res$res
               
               # update history
               object@history <- c(object@history, 
                                   list(calcTableLesion = list(call = match.call(), date = date()))
               )
               
               # update object
               eval(parse(text = paste(
                 "assign(\"", nom_object, "\", value = object, 
                 envir = .GlobalEnv)", 
                 sep = "")))             
             }
             
             return(invisible(res$res))
             
           }
)

setGeneric(name  = "calcThresholdMRIaggr", 
           def = function(object, ...){
             res <- standardGeneric("calcThresholdMRIaggr")
             
             if(res$update.object == TRUE){
               
               # alloc
               nom_object <- as.character(substitute(object))
               default_value <- data.frame(matrix(TRUE, ncol = length(res$name_newparam)))
               names(default_value) <- res$name_newparam
               
               allocContrast(object, param = res$name_newparam, default_value = default_value, overwrite = res$overwrite, verbose = res$verbose) <- res$res[,res$name_newparam]
               
               # update history
               object@history <- c(object@history, 
                                   list(calcThresholdMRIaggr = list(call = match.call(), date = date()))
               )
               
               # update object
               eval(parse(text = paste(
                 "assign(\"", nom_object, "\", value = object, 
                 envir = .GlobalEnv)", 
                 sep = "")))             
             }
             
             return(invisible(res$res))
             
           }
)

# MRIaggr
setGeneric(name = "calcTissueType", 
           def = function(object, ...){
             res <- standardGeneric("calcTissueType")
             
             if(res$update.object == TRUE){
               
               # alloc
               nom_object <- as.character(substitute(object))
               
               allocContrast(object, param = res$name_newparam, overwrite = res$overwrite, verbose = res$verbose) <- res$res$prob
               
               # update history
               object@history <- c(object@history, 
                                   list(calcTissueType = list(call = match.call(), date = date()))
               )
               
               # update object
               eval(parse(text = paste(
                 "assign(\"", nom_object, "\", value = object, 
                 envir = .GlobalEnv)", 
                 sep = "")))
             }
             
             return(invisible(res$res))
             
           }
)

# MRIaggr
setGeneric(name = "calcW", 
           def = function(object, ...){
             res <- standardGeneric("calcW")
             
             if(res$update.object == TRUE){
               
               # alloc
               nom_object <- as.character(substitute(object))
               if(any(names(res$res) == "W")){names(res$res)[names(res$res) == "W"] <- "Wmatrix"}
               if(any(names(res$res) == "blocks")){names(res$res)[names(res$res) == "blocks"] <- "Wblocks"}

               allocW(object, type = names(res$res), overwrite = res$overwrite, verbose = res$verbose) <- res$res
               
               # update history
               object@history <- c(object@history, 
                                   list(calcW = list(call = match.call(), date = date()))
               )
               
               # update object
               eval(parse(text = paste(
                 "assign(\"", nom_object, "\", value = object, 
                envir = .GlobalEnv)", 
                 sep = "")))      
               
             }
             
             return(invisible(res$res))
             
           }
)

# MRIaggr
setGeneric(name  = "outlineMRIaggr", 
           def = function(object, ...){
             res <- standardGeneric("outlineMRIaggr")
             
             if(res$update.object == TRUE){
               
               # alloc
               nom_object <- as.character(substitute(object))
               default_value <- data.frame(TRUE)
               names(default_value) <- res$name_newparam
               
               allocContrast(object, param = res$name_newparam, default_value = default_value, overwrite = res$overwrite, verbose = res$verbose) <- res$res[,c("i", "j", "k", res$name_newparam)]
               
               # update history
               object@history <- c(object@history, 
                                   list(outlineMRIaggr = list(call = match.call(), date = date()))
               )
               
               # update object
               eval(parse(text = paste(
                 "assign(\"", nom_object, "\", value = object, 
                envir = .GlobalEnv)", 
                 sep = "")))             
             }
             
             return(invisible(res$res))
             
           }
)

#### plot ####

# MRIaggr
setGeneric(name = "boxplotMask", 
           def = function(object, ...){
             standardGeneric("boxplotMask")
           }
)

# MRIaggr
setGeneric(name = "heatmapMRIaggr", 
           def = function(object, ...){
             standardGeneric("heatmapMRIaggr")
           }
)

# Carto3D MRIaggr  
setGeneric(name = "multiplot", 
           def = function(object, ...){
             standardGeneric("multiplot")
           }
)

# MRIaggr
setGeneric(name = "pointsHemisphere", 
           def = function(object, ...){
             standardGeneric("pointsHemisphere")
           }
)

# MRIaggr
setGeneric(name = "plotLesion3D", 
           def = function(object, ...){
             standardGeneric("plotLesion3D")
           }
)

# MRIaggr
setGeneric(name = "plotTableLesion", 
           def = function(object, ...){
             standardGeneric("plotTableLesion")
           }
)

# MRIaggr
setGeneric(name = "plotDistClass", 
           def = function(object, ...){
             standardGeneric("plotDistClass")
           }
)

#### const. ####

# MRIaggr 
setGeneric(name = "constCompressMRIaggr", 
           def = function(object, ...){
             standardGeneric("constCompressMRIaggr")
           }
)

# MRIaggr
setGeneric(name = "constReduceMRIaggr", 
           def = function(object, ...){
             standardGeneric("constReduceMRIaggr")
           }
)

# MRIaggr
setGeneric(name = "writeMRIaggr", 
           def = function(object, ...){
             standardGeneric("writeMRIaggr")
           }
)


#### init.  ####
# Carto3D MRIaggr
setGeneric(name = "initNum", 
           def = function(object, ...){
             standardGeneric("initNum")
           }
)

# MRIaggr
setGeneric(name = "initParameter", 
           def = function(object, ...){
             standardGeneric("initParameter")
           }
)

# MRIaggr
setGeneric(name = "initMask", 
           def = function(object, ...){
             standardGeneric("initMask")
           }
)
