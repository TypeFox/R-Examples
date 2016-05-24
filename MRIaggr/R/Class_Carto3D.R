#**********************************************************************
#**********************************************************************
#*************         0A Classe Carto3D           *******************
#**********************************************************************
#**********************************************************************
#
###### Sommaire #################
# A) Class Definition
# B) Selecters
# C) Allocators
# D) Methods

##### A) Class Definition #############################################################

methods::setClass(
  
  Class = "Carto3D", 
  
  representation(
    identifier  = "character",   # id patient 
    parameter = "character",     # parametre d intensite de la carto
    contrast = "data.frame",      # matrice integrant toutes les cartos d interet de dim (L, nb, l)
    fieldDim = "data.frame",    # vecteur contenant les dimensions des cartos 2D (l, L, nb)
    voxelDim = "data.frame",    # vecteur contenant les dimensions des cartos 2D (l, L, nb)
    default_value  = "character"
  ), 
  
  validity = function(object){
    
    if (optionsMRIaggr("checkArguments")) {
      #### contrast
      validNames(value = object@contrast, name = "@contrast", validLength = 4, validValues = c("i", "j", "k",object@parameter), method = "validity[Carto3D]")
      validDim_matrix(value1 = object@contrast, value2 = c(prod(object@fieldDim), 4), name1 = "@contrast", method = "validity[Carto3D]")
      
      #### fieldDim
      validNames(value = object@fieldDim, name = "@fieldDim", validLength = 3, validValues = c("i", "j", "k"), method = "validity[Carto3D]")
      
      #### voxelDim
      validNames(value = object@voxelDim, name = "@voxelDim", validLength = 4, validValues = c("i", "j", "k", "unit"), method = "validity[Carto3D]")
    }
    
    return(TRUE)
    
  } 
  
)

# Initiateur

methods::setMethod(f = "initialize", 
                   signature = "Carto3D", 
                   definition = function(.Object, identifier , parameter, contrast, fieldDim, voxelDim, default_value){
                     
                     if(missing(identifier)){
                       stop("initialize[Carto3D] : \'identifier\' is missing \n")
                     }
                     .Object@identifier  <- identifier 
                     
                     if(missing(parameter)){
                       if(!missing(contrast) && is.data.frame(contrast) && ncol(contrast) == 4){
                         parameter <- names(contrast)[4]
                       }else{
                         stop("initialize[Carto3D] : \'parameter\' is missing \n")
                       }
                     }
                     .Object@parameter <- parameter
                     
                     if(missing(contrast)){
                       stop("initialize[Carto3D] : \'contrast\' is missing \n")
                     }
                     
                     if(!is.data.frame(contrast)){
                       
                       if(length(dim(contrast)) %in% 3:4 == FALSE){
                         stop("initialize[Carto3D] : wrong specification of \'contrast\'\n", 
                              "required dimension of \'contrast\' : 3 or 4 \n", 
                              "proposed dimension of \'contrast\' : ", length(dim(contrast)), "\n")
                       }
                       
                       if(length(dim(contrast)) == 4 && dim(contrast)[4] == 1){
                         contrast <- contrast[,,,1, drop = TRUE]
                       }
                       
                       if(missing(fieldDim)){
                         fieldDim <- data.frame(i = dim(contrast)[1], j = dim(contrast)[2], k = dim(contrast)[3]) 
                       }
                       
                       if(missing(voxelDim)){
                         voxelDim <- data.frame(i = NA, j = NA, k = NA, unit = NA, stringsAsFactors = FALSE) 
                       }
                       
                       contrast <- data.frame(cbind(which(array(TRUE, dim = dim(contrast)), arr.ind = TRUE), 
                                                    as.vector(contrast)))
                       names(contrast) <- c("i", "j", "k", parameter)
                     }
                     
                     .Object@contrast <- contrast
                     
                     if(missing(fieldDim)){
                       stop("initialize[Carto3D] : \'fieldDim\' is missing \n")
                     }
                     .Object@fieldDim <- fieldDim
                     
                     .Object@voxelDim <- voxelDim
                     
                     if(!missing(default_value)){
                       .Object@default_value <- default_value
                     }else{
                       .Object@default_value <- "NA"
                     }
                     
                     validObject(.Object)
                     return(.Object)
                   }
)

#####  B) Selecters #############################################################

# selectContrast
methods::setMethod(f  = "selectContrast", 
                   signature = "Carto3D", 
                   definition = function(object, num = NULL, na.rm = FALSE, coords = FALSE, 
                                         format = "data.frame")
                   {     
                     #### initialisation ####
                     if(identical(coords, TRUE)){coords <- c("i", "j", "k")}
                     if(identical(coords, FALSE)){coords <- NULL}
                     
                     
                     #### tests ####
                     if (optionsMRIaggr("checkArguments")) {
                       num <- initNum(object = object, num = num, method = "selectContrast[Carto3D]")
                       validLogical(value = na.rm, validLength = 1, method = "selectContrast[Carto3D]")
                       validCharacter(value = format, validLength = 1, validValues = c("vector", "matrix", "data.frame"), method = "selectContrast[Carto3D]")
                       validCharacter(value = coords, validLength = NULL, refuse.NULL = FALSE, validValues = c("i", "j", "k"), method = "selectContrast[Carto3D]")
                       if( (format == "vector") && length(coords)>0 ) {
                         stop("selectContrast[Carto3D] : wrong specification of \'format\' \n", 
                              "vector format is not available when coords are requested \n",
                              "set \'coords\' argument to FALSE \n")
                       }                     
                     }else{
                       num <- initNum(object = object, num = num, test = FALSE, method = "selectContrast[Carto3D]")                     
                     }
                     
                     #### selection ####
                     res <- object@contrast[object@contrast$k %in% num, c(coords,selectParameter(object)),drop = FALSE]
                     
                     if(na.rm == TRUE){
                       res <- res[is.na(res[,selectParameter(object)]) == FALSE,]
                     }
                     
                     switch(format,
                            "matrix" = res <- as.matrix(res),
                            "vector" = res <- res[,,drop = TRUE]
                     )
                     
                     #### export ####
                     return(res)
                     
                   }
)

# selectCoords
methods::setMethod(f  = "selectCoords", 
                   signature  = "Carto3D", 
                   definition = function(object, coords = c("i", "j", "k"), num = NULL, format = "data.frame")
                   { 
                     
                     #### tests ####
                     if (optionsMRIaggr("checkArguments")) {
                       validCharacter(value = coords, validLength = 1:3, validValues = c("i", "j", "k"), method = "selectCoords[Carto3D]")
                       num <- initNum(object = object, num = num, method = "selectContrast[Carto3D]")
                       validCharacter(value = format, validLength = 1, validValues = c("vector", "matrix", "data.frame"), method = "selectCoords[Carto3D]")
                       if( (format == "vector") && length(coords)>1 ) {
                         stop("selectContrast[Carto3D] : wrong specification of \'format\' \n", 
                              "vector format is not available when more than one coordinate is requested \n",
                              "set \'format\' argument to \"matrix\" or \"data.frame\" \n")
                       } 
                     } else {
                       num <- initNum(object = object, test = FALSE, num = num, method = "selectContrast[Carto3D]")
                     }
                     
                     ####  selection ####
                     res <- object@contrast[object@contrast$k %in% num,coords]
                     
                     switch(format,
                            "matrix" = res <- as.matrix(res),
                            "vector" = res <- res[,,drop = TRUE]
                     )
                     
                     #### export ####
                     return(res)
                     
                   }
)

# selectDefault_value
methods::setMethod(f  = "selectDefault_value", 
                   signature  = "Carto3D", 
                   definition = function(object){
                     return(object@default_value) 
                   }
)

# selectIdentifier 
methods::setMethod(f  = "selectIdentifier", 
                   signature  = "Carto3D", 
                   definition = function(object){
                     return(object@identifier )   
                   }
)

# selectN
methods::setMethod(f  = "selectN", 
                   signature  = "Carto3D", 
                   definition = function(object, num = NULL){
                     
                     #### tests ####
                     if (optionsMRIaggr("checkArguments")) {
                       num <- initNum(object = object, num = num, method = "selectContrast[Carto3D]")
                     } else {
                       num <- initNum(object = object, test = FALSE, num = num, method = "selectContrast[Carto3D]")
                     }
                     
                     #### export ####
                     return( sum(object@contrast$k %in% num) )
                   }
) 

# selectParameter
methods::setMethod(f  = "selectParameter", 
                   signature  = "Carto3D", 
                   definition = function(object){
                     return(object@parameter)    
                   }
)

# selectFieldDim
methods::setMethod(f  = "selectFieldDim", 
                   signature  = "Carto3D", 
                   definition = function(object){
                     return(object@fieldDim)  
                   }
)

#####  C) Allocators #############################################################




#####  D) Methods #############################################################


#### plot ####

methods::setMethod(f  = "multiplot", 
                   signature  = "Carto3D", 
                   definition = function(object, num = NULL, col = NULL, pch = NULL, xlim = NULL, ylim = NULL, filename = "multiplot", ...){
                     
                     #### graphical options ####
                     ## get graphical arguments
                     
                     optionsMRIaggr.eg <- optionsMRIaggr()
                     dots.arguments <- list(...)
                     names_dots.arguments <- names(dots.arguments)
                     
                     if("main" %in% names_dots.arguments == FALSE){optionsMRIaggr.eg$main <- paste(selectParameter(object), " : ", selectIdentifier(object), " - slice", sep = "")}
                     if("main.legend" %in% names_dots.arguments == FALSE){optionsMRIaggr.eg$main.legend <- paste("", selectParameter(object), sep = "")}
                     
                     ## tests 
                     if (optionsMRIaggr("checkArguments")) {
                       validCharacter(names_dots.arguments, name = "...", validLength = NULL, validValues = names(optionsMRIaggr.eg), refuse.NULL = FALSE, method = "multiplot[Carto3D]")
                     }
                     
                     ## set specific display
                     if(length(names_dots.arguments) > 0){
                       optionsMRIaggr.eg[names_dots.arguments] <- dots.arguments[names_dots.arguments]
                     }
                     
                     ##                      
                     if(filename == "auto"){
                       filename <- paste("multiplot", selectIdentifier(object), "_", selectParameter(object), sep = "")
                     }

                     #### display ####
                     res <- multiplot(object = selectCoords(object = object, num = num, format = "data.frame"), 
                                      contrast = selectContrast(object = object, num = num, coords = FALSE), 
                                      slice_var = optionsMRIaggr.eg$slice_var, breaks = optionsMRIaggr.eg$breaks, type.breaks = optionsMRIaggr.eg$type.breaks, palette = optionsMRIaggr.eg$palette, col = col, pch = pch, cex = optionsMRIaggr.eg$cex, 
                                      col.NA = optionsMRIaggr.eg$col.NA, pch.NA = optionsMRIaggr.eg$pch.NA, xlim = xlim, ylim = ylim, axes = optionsMRIaggr.eg$axes, 
                                      window = optionsMRIaggr.eg$window, legend = optionsMRIaggr.eg$legend, mfrow = optionsMRIaggr.eg$mfrow, mar = optionsMRIaggr.eg$mar, mgp = optionsMRIaggr.eg$mgp, pty = optionsMRIaggr.eg$pty, asp = optionsMRIaggr.eg$asp, bg = optionsMRIaggr.eg$bg, 
                                      xlab = optionsMRIaggr.eg$xlab, ylab = optionsMRIaggr.eg$ylab, main = optionsMRIaggr.eg$main, num.main = optionsMRIaggr.eg$num.main, cex.main = optionsMRIaggr.eg$cex.main, 
                                      quantiles.legend = optionsMRIaggr.eg$quantiles.legend, digit.legend = optionsMRIaggr.eg$digit.legend, cex.legend = optionsMRIaggr.eg$digit.legend, mar.legend = optionsMRIaggr.eg$mar.legend, 
                                      main.legend = optionsMRIaggr.eg$main.legend,    
                                      filename = filename, width = optionsMRIaggr.eg$width, height = optionsMRIaggr.eg$height, path = optionsMRIaggr.eg$path, unit = optionsMRIaggr.eg$unit, res = optionsMRIaggr.eg$res)
                     
                     #### export ####
                     return(invisible(res))
                   }
)

#### init. ####

methods::setMethod(f  = "initNum", 
                   signature  = "Carto3D", 
                   definition = function(object, num, test = TRUE, init = TRUE, method){
                     
                     if(init == TRUE && is.null(num)){
                       num <- seq(1, object@fieldDim$k)
                     }
                     
                     if(test == TRUE){			
                       validInteger(value = num, validLength = NULL, validValues = seq(1, object@fieldDim$k), 
                                    refuse.NA = TRUE, refuse.NULL = TRUE, refuse.duplicates = TRUE, method = paste(method,"[Carto3D]",sep=""))						 
                     }
                     
                     return(num)
                   }
)

