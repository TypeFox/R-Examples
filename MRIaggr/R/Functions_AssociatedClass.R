#### 1- Construction functions for the objects Carto3D and MRIaggr ####

constCarto3D <- function(array, identifier, param, default_value = NULL, 
                         pos_default_value = c(1, 1, 1), voxelDim = NULL, rm.array = FALSE){
  
  nom.array <- as.character(substitute(array))
  
  
  #### conversion du format ####
  if(class(array) == "anlz"){ # analyse object
    if(is.null(voxelDim)){
      voxelDim <- data.frame(i = NA, j = NA, k = NA, unit = array@vox_units, stringsAsFactors = FALSE)
      voxelDim[,c("i", "j", "k")] <- array@pixdim[2:4] 
    }
    array <- array@.Data
  }else if(class(array) %in% c("nifti", "niftiExtension", "niftiAuditTrail")){ # nifti object
    if(is.null(voxelDim)){
      voxelDim <- data.frame(i = NA, j = NA, k = NA, unit = oro.nifti::convert.units(oro.nifti::xyzt2space(array@xyzt_units)), stringsAsFactors = FALSE)
      voxelDim[,c("i", "j", "k")] <- array@pixdim[2:4]    
    }
    array <- array@.Data
  }else if(is.list(array) && length(array) == 2 && "img" %in% names(array)){ # dicom object
    array <- array$img
    voxelDim <- data.frame(i = NA, j = NA, k = NA, unit = NA, stringsAsFactors = FALSE)
  }else{ # array object
    if(class(array) %in% c("matrix", "array") == FALSE){
      stop("constCarto3D : wrong specification of \'array\' \n", 
           "\'array\' must be either an \"anlz\", \"nifti\" or an array object \n", 
           "class(array) : ", class(array), "\n")
      
    }
    voxelDim <- data.frame(i = NA, j = NA, k = NA, unit = NA, stringsAsFactors = FALSE)
  }
  
  if(length(dim(array)) == 4 && dim(array)[4] == 1){
    array <- array[,,, 1, drop = TRUE]
  }
  if(length(dim(array)) == 2){
    array <- array(array, dim = c(dim(array), 1))
  }
  
  if(length(dim(array)) != 3){
    stop("constCarto3D : wrong specification of \'array\' \n", 
         "\'array\' must have 3 dimensions \n", 
         "proposed \'array\' has ", length(dim(array)), " dimensions : ", paste(dim(array), collapse = " "), " (parameter : ", param, ") \n")
  }
  
  ##### valeurs par defaut ####
  if(is.null(default_value)){
    
    if(is.vector(pos_default_value) == TRUE){pos_default_value <- matrix(pos_default_value, nrow = 1)}
    
    if(ncol(pos_default_value) != length(dim(array))){
      stop("constCarto3D : \'pos_default_value\' and \'array\' dimensions are inconsistent \n", 
           "\'array\' has ", length(dim(array)), " dimensions \n", 
           "\'pos_default_value\' has : ", ncol(pos_default_value), " dimensions \n")          
    }
    
    default_value <- names(which.max(table(array[pos_default_value]))[1, drop = FALSE])
  }
  
  ##### construction de la carto3D ####  
  res <- methods::new(Class = "Carto3D", 
             identifier = identifier, 
             parameter = param, 
             voxelDim = voxelDim, 
             default_value = as.character(default_value), 
             contrast = array)
  
  #### nettoyage ####
  if(rm.array){
    remove(list = nom.array, envir = globalenv())
  }
  
  #### export ####
  return(res)
  
}

constMRIaggr <- function(ls.array, identifier, param, default_value = NULL, 
                         pos_default_value = c(1, 1, 1), 
                         tol = 10^{-10}, voxelDim = NULL, 
                         verbose = optionsMRIaggr("verbose"), rm.ls.array = FALSE){
  
  nom.ls.array <- as.character(substitute(ls.array))
  
  if(class(ls.array) %in% c("anlz", "nifti", "niftiExtension", "niftiAuditTrail") || (is.list(ls.array) && length(ls.array) == 2 && "img" %in% names(ls.array)) ){
    ls.array <- list(ls.array)
  }
  
  if(is.list(ls.array) == FALSE){
    ls.array <- list(ls.array)
  }
  M <- length(ls.array)
  
  if(is.null(param)){
    param <- names(ls.array)
  }
  
  #### test ####
  test.class <- unlist(lapply(ls.array, function(x){class(x) %in% c("anlz", "nifti", "niftiExtension", "niftiAuditTrail", "list", "array", "matrix")}))
  if(any(test.class == FALSE)){
    stop("constMRIaggr : wrong specification of \'ls.array\' \n", 
         "\'ls.array\' must be a list of \"anlz\" or \"nifti\" or \"list\" or \"array\" \n", 
         "elements of \'ls.array\' that are not of type array : ", paste(which(test.class == FALSE), collapse = " "), " \n")
  }
  
  validCharacter(value = identifier, validLength = 1, refuse.NULL = TRUE, method = "constMRIaggr")
  validDim_vector(value1 = param, value2 = list(M), name2 = "ls.array", method  = "constMRIaggr")
  if(!is.null(default_value)){
    validDim_vector(value1 = default_value, value2 = list(M), name2 = "ls.array", method  = "constMRIaggr")
  }
  
  #### transformation en carto3D
  ls.Carto3D <- list()
  for(iter_m in 1:M){
    
    if(is.null(default_value)){
      ls.Carto3D[[iter_m]] <- constCarto3D(ls.array[[iter_m]], 
                                           identifier = identifier, param = param[iter_m], default_value = NULL, 
                                           pos_default_value = pos_default_value, voxelDim = voxelDim)
      
    }else{
      ls.Carto3D[[iter_m]] <- constCarto3D(ls.array[[iter_m]], 
                                           identifier = identifier, param = param[iter_m], default_value = default_value[iter_m], 
                                           pos_default_value = NULL, voxelDim = voxelDim)      
    }
    
  }
  
  #### transformation en MRIaggr
  MRIaggr <- Carto3D2MRIaggr(ls.Carto3D = ls.Carto3D, rm.Carto3D = FALSE, tol = tol, 
                             verbose = verbose)
  
  #### nettoyage 
  if(rm.ls.array){
    rm(list = nom.ls.array, envir = globalenv())
  }
  
  #### export
  return(MRIaggr)
}  

#### 2- Conversion functions ####
array2df <- function(array, coords = NULL, 
                     name_newparam = "res", names_coords = letters[9:(8 + ncol(coords))], na.rm = TRUE){
  
  ### preparation 
  p <- length(dim(array))  
  n <- length(array)  
  index_array <- arrayInd(1:n, .dim = dim(array))
  if(na.rm){
    index_array <- index_array[is.na(array) == FALSE,]
  }
  
  if(is.null(coords)){
    coords <- index_array
  }else{
    if(any(dim(coords) != dim(index_array))){
      stop("array2df : incorrect dimension for \'coords\' : \n", 
           "lenght(array[!is.na(array)]), length(dim(array))  = ", paste(dim(index_array), collapse = " "), "\n", 
           "dim(coords) = ", paste(dim(coords), collapse = " "), "\n")
    }
  }
  
  validDim_vector(value1 = names_coords, value2 = list(p), name1 = "names_coords", name2 = NULL, type = "length", method  = "array2df")
  
  ### integration des donnees 
  data <- data.frame(coords, array[index_array])
  names(data) <- c(names_coords, name_newparam)
  
  ### export
  return(data)
}

Carto3D2MRIaggr <- function(ls.Carto3D, rm.Carto3D = FALSE, tol = 10^{-10}, 
                            num = NULL, verbose = optionsMRIaggr("verbose")){
  
  nom.ls.Carto3D <-  as.character(substitute(ls.Carto3D))[-1]
  
  if(!is.list(ls.Carto3D)){
    ls.Carto3D <- as.list(ls.Carto3D)
  }
  
  #### test preliminaire ####
  if(any(unlist(lapply(ls.Carto3D, function(x){class(x) == "Carto3D"})) == FALSE)){
    stop("Carto3D2MRIaggr : wrong specificition of \'ls.array\' \n", 
         "\'ls.array\' must be a list of \"Carto3D\" objects \n", 
         "elements of \'ls.array\' that are not \"Carto3D\" objects : ", paste(which(unlist(lapply(ls.Carto3D, function(x){class(x) == "Carto3D"})) == FALSE), collapse = " "), " \n")    
  }
  
  #### preparation ####
  
  if(is.null(num)){
    num <- seq(1, selectFieldDim(ls.Carto3D[[1]])$k, by = 1)
  }
  fieldDim <- data.frame(i = selectFieldDim(ls.Carto3D[[1]])$i, 
                         j = selectFieldDim(ls.Carto3D[[1]])$j, 
                         k = length(num), stringsAsFactors = FALSE)
  
  default_value <- data.frame(matrix("NA", ncol = length(ls.Carto3D), nrow = 1), stringsAsFactors = FALSE)
  identifiant <- selectIdentifier(ls.Carto3D[[1]])
  size <- ls.Carto3D[[1]]@voxelDim
  
  data_global <- data.frame(matrix(NA, 
                                   nrow = selectN(ls.Carto3D[[1]], num = num), 
                                   ncol = length(ls.Carto3D)))
  data_global[,1:3] <- selectCoords(ls.Carto3D[[1]], num = num, format = "data.frame")
  names(data_global)[1:3] <- c("i", "j", "k")
  
  #### assemblage ####
  if(verbose){cat("Merging : ")}
  for(iter in 1:length(ls.Carto3D)){
    
    if(verbose){cat("(", iter, ")", sep = "")}

    data_tempo <- selectContrast(ls.Carto3D[[iter]], num = num, coords = TRUE)
    nom_param <- selectParameter(ls.Carto3D[[iter]])
    if(verbose){cat(" ", nom_param, " ", sep = "")}
    
    default_value[iter] <- selectDefault_value(ls.Carto3D[[iter]])
    names(default_value)[iter] <- nom_param
    
    if(sum(abs(data_global[,1:3]-data_tempo[,1:3])) > tol){
      stop("Carto3D2MRIaggr - some \"Carto3D\" objects of \'ls.Carto3D\' have different coordinates \n", 
           "Difference between \'ls.Carto3D\'[[1]] and \'ls.Carto3D\'[[", iter, "]]  \n")
    }
    
    if(selectIdentifier(ls.Carto3D[[iter]]) != identifiant){
      stop("Carto3D2MRIaggr -  some \"Carto3D\" objects of \'ls.Carto3D objects correspond to different patients \n", 
           "patient of \'ls.Carto3D\'[[1]] : ", selectIdentifier(ls.Carto3D[[1]]), "\n", 
           "patient of \'ls.Carto3D\'[[", iter, "]] : ", selectIdentifier(ls.Carto3D[[iter]]), "\n")
    }
    
    if(nom_param %in% names(data_global)){
      warning("Carto3D2MRIaggr - some \"Carto3D\" objects of \'ls.Carto3D\' correspond to the same parameter \n", 
              "common parameter : ", nom_param, "\n", 
              "index of the elements in \'ls.Carto3D\' : ", paste(which(names(data_global) == nom_param) - 2, collapse = " "), ") \n")
    }
    
    data_global[,iter + 3] <- data_tempo[,4]
    names(data_global)[iter + 3] <- nom_param
  } 
  if(verbose){cat("\n")}
  
  #### nettoyage
  if(rm.Carto3D){
    rm(list = nom.ls.Carto3D, envir = globalenv())
  }
  
  #### export ####
  res <- methods::new(Class = "MRIaggr", 
                      identifier = identifiant, 
                      contrast = data_global, 
                      default_value = default_value, 
                      fieldDim = fieldDim, 
                      voxelDim = size)
  
  return(res)
}

df2array <- function(contrast, coords, format = "any", default_value = NA, 
                     range.coords = NULL){
  
  ### preparation
  validClass(value = coords, validClass = c("matrix", "data.frame"), superClasses = TRUE, method = "df2array")
  
  contrast <- as.data.frame(contrast)
  M <- ncol(contrast)
  p <- ncol(coords)
  n <- nrow(coords)
  
  validDim_matrix(value1 = nrow(contrast), value2 = n, name1 = "contrast", name2 = "coords", type = "nrow", method = "df2array")
  
  if(p > 3){
    stop("df2array : wrong specification of \'coords\' \n", 
         "\'coords\' can have up to 3 columns \n", 
         "number of columns of \'coords\' : ", p, "\n")
  }
  
  validCharacter(value = format, validLength = 1, validValues = c("any", "matrix", "data.frame", "list"), 
                 refuse.NULL = TRUE, method = "df2array")
  
  if(!is.null(range.coords) && p != length(range.coords)){
    stop("df2array : wrong specification of \'range.coords\' \n", 
         "if not NULL, \'range.coords\' must be a vector of length ", p, " \n", 
         "is(range.coords) : ", paste(is(range.coords), collapse = " "), " \n", 
         "length(range.coords) : ", length(range.coords), "\n")
  }  
  
  #### definition des coordonnees des points dans le repere de la matrice
  scale <- rep(0, p)
  if(is.null(range.coords)){ 
    Mcoords <- apply(coords, 2, 
                     function(x){scale_tempo = min(x) - 1;
                     return(c(scale_tempo, x-scale_tempo))}
    )
    scale <- Mcoords[1,]
    Mcoords <- Mcoords[-1,, drop = FALSE]    
  }else{
    Mcoords <- coords
  }
  
  if(any(Mcoords %% 1 != 0)){    
    Mcoords <- apply(Mcoords, 2, function(x){as.numeric(as.factor(x))})
  }
  
  dimnames <- names(Mcoords)
  Mcoords <- as.matrix(Mcoords)
  if(is.null(range.coords)){
    range.coords <- unlist(apply(Mcoords, 2, function(x){max(x)}))
  }
  if(any(range.coords < apply(Mcoords, 2, max))){
    stop("df2array : wrong specification of \'range.coords\' \n", 
         "\'range.coords\' must be at least ", paste(apply(Mcoords, 2, max), collapse = " "), " \n", 
         "requested \'range.coords\' : ", paste(range.coords, collapse = " "), "\n")    
  }
  
  #### index correspondant aux points dans la matrice
  
  Mindex <- eval(parse(text = paste("1+", paste( c(1, if(ncol(Mcoords) > 1){range.coords[1]}, if(ncol(Mcoords) > 2){range.coords[1]*range.coords[2]}), 
                                                 " * (Mcoords[,", 1:ncol(Mcoords), "] - 1)", sep = "", collapse = "+"), sep = "")))
  
  ### integration des donnees
  dataM <- list()
  
  for(iter_m in 1:M){
    
    dataM[[iter_m]] <- array(default_value, dim = range.coords)
    dataM[[iter_m]][Mindex] <- contrast[,iter_m]
    
    if(format == "matrix" && p == 2){
      dataM[[iter_m]] <- as.matrix(dataM[[iter_m]])
    }
    if(format == "data.frame" && p == 2){
      dataM[[iter_m]] <- as.data.frame(dataM[[iter_m]])
    }
  }
  names(dataM) <- names(contrast)
  
  ### mise en forme
  if(format %in% c("matrix", "data.frame")  && M == 1){
    dataM <- dataM[[1]]
  }
  
  unique_coords <- lapply(1:ncol(Mcoords), function(x){scale[x]+seq(min(Mcoords[,x]), max(Mcoords[,x]), by = 1)})
  names(unique_coords) <- names(coords)
  
  return(list(contrast = dataM, 
              coords = coords, 
              unique_coords = unique_coords))
  
  
  
}

readMRI <- function (filename, format, na.value = 0, 
                     what = "numeric", size = "NA_integer_", dimensions = NULL, 
                     SPM = FALSE, reorient = FALSE, flipud = FALSE){
  
  validCharacter(value = format, validLength = 1, validValues = c("rawb.gz", "analyze", "nifti", "dicom"),
                 refuse.NULL = TRUE, method = "readMRI")
  
  if (format == "rawb.gz") {
    
    validInteger(value = dimensions, validLength = 3, 
                 refuse.NA = TRUE, refuse.NULL = TRUE, refuse.duplicates = FALSE, method = "readMRI")
    
    f <- gzfile(filename, open = "rb")
    on.exit(close(f))
    data <- readBin(con = f, what = what, n = prod(dimensions), size = size)
    res <- array(data, dimensions)
	
  } else if (format == "analyze") {
    
    # initPackage(package = "oro.nifti", argument = "format = \"oro.nifti\"", method = "readMRI")
    res <- oro.nifti::readANALYZE(filename, SPM = SPM)
	
  } else if (format == "nifti") {
  
    # initPackage(package = "oro.nifti", argument = "format = \"oro.nifti\"", method = "readMRI")
    res <- oro.nifti::readNIfTI(filename, reorient = reorient)
	
  } else if (format == "dicom")  {
    # initPackage(package = "oro.dicom", argument = "format = \"oro.dicom\"", method = "readMRI")
    res <- oro.dicom::readDICOMFile(filename, flipud = flipud)
  }
  
  if(!is.na(na.value)){
    
    if(format == "dicom"){
      test.na <- is.na(res$img)
      if(any(test.na)){res$img[test.na] <- na.value}        
    }else{
      test.na <- is.na(res)
      if(any(test.na)){res[test.na] <- na.value}
    }
    
  }
  
  return(res)
  
}

writeMRI <- function (data, filename, format, gzipped = TRUE, verbose = optionsMRIaggr("verbose"), size = "NA_integer_"){
  
  validCharacter(value = format, validLength = 1, validValues = c("rawb.gz", "analyze", "nifti", "dicom"), 
                 refuse.NULL = TRUE, method = "writeMRI")
  
  objClass <- class(data)
  if (objClass != "array" || length(dim(data)) != 3){
    stop("writeMRI : incorrect \'data\' : \n", 
         "data has to be a 3 dimensional array  \n", 
         "proposed data (dimension) : ", paste(is(data), collapse = " "), " (", paste(dim(data), collapse = " "), ") \n")
  }
  
  if (format == "rawb.gz") {
  
    f <- gzfile(filename, open = "wb")
    writeBin(con = f, object = as.vector(data), size = size)
    on.exit(close(f))
	
  } else if (format == "analyze") {
  
    # initPackage(package = "oro.nifti", argument = "format = \"oro.nifti\"", method = "writeMRI")
    oro.nifti::writeANALYZE(as(data, "anlz"), filename = filename, gzipped = gzipped, verbose = verbose)
	
  } else if (format == "nifti") {
  
    # initPackage(package = "oro.nifti", argument = "format = \"oro.nifti\"", method = "writeMRI")
    oro.nifti::writeNIfTI(as(data, "nifti"), filename = filename, gzipped = gzipped, verbose = verbose)
	
  } else if (format == "dicom")  {
  
    # initPackage(package = "oro.dicom", argument = "format = \"oro.dicom\"", method = "writeMRI")
    cat("writeMRI for dicom files is not implemented \n")
    invisible(return(FALSE))    
	
  }
  
}

#### 3- Calculation functions ####

calcBlockW <- function(W, site_order = NULL, dist.center = NULL, dist.max = Inf, verbose = optionsMRIaggr("verbose")){
  
  ## test
  # initPackage(package = "spam", method = "calcBlockW")
  # initPackage(package = "Matrix", method = "calcBlockW")
  
  
  if(!is.null(dist.center)){
    site_order <- order(dist.center) - 1
  }else{
    dist.center <- rep(0, length(W@p) - 1)
    if(dist.max < Inf){
      stop("calcBlockW : \'dist.max\' cannot be specified if argument \'dist.center\' is NULL \n")
    }
  }
  
  if(is.null(site_order)){
    site_order <- -1
  }else{
    if(length(site_order) != (length(W@p) - 1)){
      stop("calcBlockW : wrong specification of \'site_order\' \n", 
           "length(site_order) : ", length(site_order), "\n", 
           "length(W@p) - 1 : ", length(W@p) - 1, "\n") 
    }
    
    if(min(site_order) != 0 || max(site_order) != (length(W@p) - 2) ){
      stop("calcBlockW : wrong specification of \'site_order\' \n", 
           "must contains values between 0 and ", (length(W@p) - 2), " \n", 
           "range(site_order) = ", paste(range(site_order), collapse = " "), "\n") 
    }
  }
  
  ## conversion des indices depuic C a R
  res <- calcBlockW_cpp(W_i = W@i, W_p = W@p, site_order = site_order, dist_center = dist.center, dist_max = dist.max, verbose = verbose)
  res$ls_groups <- lapply(res$ls_groups, function(x){x + 1})
  
  ## export
  return(res)
  
}

calcGroupsCoords <- function(coords, array = NULL, Neighborhood, max_groups = 10000, verbose = optionsMRIaggr("verbose")){
  
  if(is.null(array)){
    if(any(coords < 0)){
      stop("calcGroupsCoords : wrong specification of  \'coords\' \n", 
           "\'coords must\' be positve \n", 
           "min(coords) : ", min(coords), "\n")
    }
    
    if(nrow(coords) == 0){
      return(list(indexArray.group = NULL,
                  indexCoords.group = NULL, 
                  df.group = NULL, 
                  group_size = 0)
      )
    }
    
    if(nrow(coords) == 1){
      return(list(indexArray.group = NULL,
                  indexCoords.group = list(1),  
                  df.group = data.frame(coords, index = 1, group = 1), 
                  group_size = 1)
      )
    }
    
    coords_NNA <- apply(coords, 2, function(x){x-min(x)})
    dim_coordsNNA <- apply(coords_NNA, 2, max) + 1
    p <- length(dim_coordsNNA)
    index_NNA <- rowSums(sweep(coords_NNA, 2, c(1, cumprod(dim_coordsNNA)[-p]), FUN = "*"))
  }else{
    coords_NNA <- which(!is.na(array), arr.ind = TRUE) - 1
    
    if(nrow(coords_NNA) == 0){
      return(list(indexArray.group = NULL,
                  indexCoords.group = NULL, 
                  df.group = NULL, 
                  group_size = 0)
      )
    }
    
    if(nrow(coords_NNA) == 1){
      return(list(ls.group = list(1), 
                  df.group = data.frame(coords_NNA, index = 1, group = 1), 
                  group_size = 1)
      )
    }
    
    dim_coordsNNA <- dim(array)
    p <- length(dim_coordsNNA)
    index_NNA <- which(!is.na(array)) - 1
  }
  
  #### definition de Neighborhood
  if(length(Neighborhood) == 1 && is.character(Neighborhood)){
    Neighborhood <- initNeighborhood(Neighborhood, method = "calcGroupsCoords")
  }
  
  validDim_matrix(value1 = Neighborhood, value2 = p, name1 = "Neighborhood", type = "ncol", method = "calcGroupsCoords")
  
  #### Rcpp
  res_cpp <-  calcGroupsCoords_cpp(coords_NNA = coords_NNA, 
                                   index_NNA = index_NNA, 
                                   min_index_NNA = 0,#index_NNA[1],
                                   max_index_NNA = index_NNA[length(index_NNA)],
                                   Neighborhood = Neighborhood, 
                                   coords_max = dim_coordsNNA, 
                                   max_groups = max_groups, 
                                   verbose = verbose)
  
  if(res_cpp$cv == FALSE){
    warning("calcGroupsCoords : maximum number of groups reached \n", 
            "the identification of the spatial groups may not be complet \n", 
            "set \'max_groups\' higher to allow more groups \n")
  }
  
  #### export  
  indexArray.group <- list()
  for(iter_group in 1:length(res_cpp$group_size)){
    indexArray.group[[iter_group]] <- index_NNA[res_cpp$group == iter_group]
  }
  
  indexCoords.group <- list()
  for(iter_group in 1:length(res_cpp$group_size)){
    indexCoords.group[[iter_group]] <- which(res_cpp$group == iter_group)
  }
 
  if(is.null(array)){
    df.group <- data.frame(coords, group = res_cpp$group)
  }else{
    df.group <- NULL
  }
  
  return(list(indexArray.group = indexArray.group,
              indexCoords.group = indexCoords.group, 
              df.group = df.group, 
              group_size = res_cpp$group_size)
  )
}

calcThreshold <- function(contrast, param, hemisphere = NULL, rm.CSF = FALSE, threshold = 1:10, decreasing = FALSE, 
                          GRalgo = FALSE, W = NULL, seed = NULL, numeric2logical = FALSE, verbose = optionsMRIaggr("verbose")){
  

  #### pre test
  if( optionsMRIaggr("checkArguments")) {
    
    validClass(value = contrast, validClass = c("matrix","data.frame"), superClasses = FALSE, method = "calcThreshold")
    validLogical(value = numeric2logical, validLength = 1, method = "calcThreshold")
    
    if (is.character(seed) == TRUE) {  
    
      validCharacter(seed, validLength = NULL, validValues = names(contrast), method = "calcThreshold")
      sapply(seed,function(x){
        if(numeric2logical == TRUE){
          validNumeric(contrast[,x], name = "seed", validLength = NULL, method = "calcThreshold")
        }else{
          validLogical(contrast[,x], name = "seed", validLength = NULL, method = "calcThreshold")
        }
      })
      
    }
    
  }
  
  #### pre initialization
  #if(GRalgo == TRUE){
    # initPackage(package = "spam", argument = "GRalgo = \"TRUE\"", method = "calcThreshold")
    # initPackage(package = "Matrix", argument = "GRalgo = \"TRUE\"", method = "calcThreshold")
  #}
  
  p <- length(param)
  n <- nrow(contrast)
  if(decreasing == TRUE){
    threshold <- - threshold
    contrast[,param] <- - contrast[,param]
  }
  threshold <- sort(threshold)
  n.threshold <- length(threshold)
  
  if(is.character(rm.CSF) == TRUE){
    param.CSF <- rm.CSF
    rm.CSF <- TRUE
  }else{ 
    param.CSF <- "CSF"
  }
  
  if (is.character(seed) == TRUE) {
    seed <- contrast[,seed, drop = FALSE]
    if(numeric2logical == TRUE){ seed <- apply(seed, 2, as.logical)}
    seed <-  which(rowSums(seed) > 0)
  }
 
   #### test
  if( optionsMRIaggr("checkArguments")) {
    
    validNames(value = param, validValues = names(contrast), method = "calcThreshold")
    
    if (rm.CSF == TRUE && (param.CSF %in% names(contrast) == FALSE)) {
      stop("calcThreshold : wrong specification of \'contrast\' \n", 
           "contrast must be contains a column \"", param.CSF, "\" if rm.CSF = TRUE \n", 
           "column names of \'contrast\' : ", paste(names(contrast), collapse = " "), " \n")
    }
    
    if (!is.null(hemisphere)) {
      
      if ("hemisphere" %in% names(contrast) == FALSE) {
        stop("calcThreshold : wrong specification of \'contrast\' \n", 
             "contrast must be contains a columns \"hemisphere\" if \'hemisphere\' is not NULL \n", 
             "column names of \'contrast\' : ", paste(names(contrast), collapse = " "), " \n")
      }
      
      if (hemisphere %in% unique(contrast$hemisphere) == FALSE) {
        stop("calcThreshold : wrong specification of \'hemisphere\' \n", 
             "\'hemisphere\' does not match element of the hemisphere column in contrast \n", 
             "unique(contrast$hemisphere) : ", unique(contrast$hemisphere), "\n", 
             "proposed \'hemisphere\' : ", hemisphere, " \n")
      }
      
    }
    
    validLogical(value = rm.CSF, validLength = 1, method = "calcThreshold")
    validNumeric(value = threshold, validLength = NULL, refuse.duplicates = TRUE, method = "calcThreshold")
    validLogical(value = decreasing, validLength = 1, method = "calcThreshold")
    validLogical(value = GRalgo, validLength = 1, method = "calcThreshold")
    
    if(GRalgo == TRUE){
      
      validDim_matrix(value1 = W, value2 = c(n,n), name1 = "W", name2 = "contrast", type = "both", method = "calcThreshold")
      validInteger(value = seed, validLength = NULL, min = 1 , max = n, refuse.duplicates = TRUE, method = "calcThreshold")
      
    }
    
    validLogical(value = verbose, validLength = 1, method = "calcThreshold")
    
  }
  
  #### step 1 - CSF et hemi
  if(verbose == TRUE){cat("Step 1 : ")}
  index.perf <- 1:n
  
  if(!is.null(hemisphere)){
    if(verbose == TRUE){cat("keep only the \"", hemisphere, "\" hemisphere", sep = "")}
    index.perf <- intersect(index.perf, which(contrast$hemisphere %in% hemisphere))    
  } else{
    if(verbose == TRUE){cat("keep both hemipheres")}
  }
  
  if(rm.CSF == TRUE){
    if(verbose == TRUE){cat(", remove CSF ")}
    index.perf <- intersect(index.perf, which(contrast[,param.CSF] < 0.5))
  } else{
    if(verbose == TRUE){cat(", keep CSF ")}
  }
  
  contrast <- contrast[,param, drop = FALSE]
  min_perf <- threshold[1] - 1
  contrast[-index.perf,] <- min_perf
  if(verbose == TRUE){cat("\n")}
  
  #### etape 2 - seuillage
  if(verbose == TRUE){cat("Step 2 : thresholding ")}
  tempo <- matrix(min_perf, nrow = n, ncol = p)
  colnames(tempo) <- param
  
  
  for(iter_threshold in threshold){
    if(verbose == TRUE){cat("*")}
    for(iter_param in 1:p){
      tempo[contrast[,param[iter_param]] >= iter_threshold, iter_param] <- iter_threshold      
    }
  }
  
  contrast <- tempo
  rm(tempo)
  gc()
  if(verbose == TRUE){cat("\n")}
  
  #### etape 3 - GR
  if(GRalgo == TRUE){
    if(verbose == TRUE){cat("Step 3 : Growing Region \n", sep = "")}
    
    for(iter_param in param){
      if(verbose == TRUE){cat(iter_param, " ", sep = "")}
      
      for(iter_threshold in 1:n.threshold){
        if(verbose == TRUE){cat("*")}
        
        index_threshold <- which(contrast[,iter_param] >= threshold[iter_threshold])
        
        contrastBis <- rep(0, n)
        contrastBis[index_threshold] <- 100
        
        resGR <- GRalgo(contrastBis, W = W, seed = seed[seed >= threshold[iter_threshold]], 
                        sigma_max = 0.0001, range = c(threshold[iter_threshold], Inf), breaks = seq(-1, 109, by = 10), step = 10, 
                        operator = "sd", iter_max = 1000, keep.lower = FALSE, keep.upper = FALSE,
                        history.sigma = FALSE, history.step = FALSE, history.front = FALSE)
        
        if(length(setdiff(index_threshold, resGR$GR)) == 0){next}
        
        if(iter_threshold == 1){
          contrast[setdiff(index_threshold, resGR$GR), iter_param] <- min_perf
        }else{
          contrast[setdiff(index_threshold, resGR$GR), iter_param] <- threshold[iter_threshold - 1]
        }
     
      }
      if(verbose == TRUE){cat("\n")}
    }
    if(verbose == TRUE){cat("\n")}
  }
  
  if(decreasing == TRUE){
    contrast <- - contrast  
  }
  names(contrast) <- param
  
  return(as.data.frame(contrast))
}

calcGroupsW <- function(W, subset = NULL, max_groups = 10000, verbose = optionsMRIaggr("verbose")){ 
  
  #### tests
  # initPackage(package = "spam", method = "calcGroupsW")
  # initPackage(package = "Matrix", method = "calcGroupsW")
  
  validClass(value = W, validClass = "dgCMatrix", superClasses = TRUE, method = "calcGroupsW")
  
  
  #### initialization
  if(is.null(subset)){
    subset <- -1
  }else{
    subset <- subset - 1
  }
  
  #### call C++ function
  resCpp <- calcGroupsW_cpp(W_i = W@i, W_p = W@p, subset = subset, max_groups = max_groups, verbose = verbose)

  if(resCpp$cv == FALSE){
    warning("calcGroupsW : maximum number of groups reached \n", 
            "the identification of the spatial groups may not be complet \n", 
            "set \'max_groups\' higher to allow more groups \n")
  }
  
  #### export
  res <- list(group = resCpp$group, 
              group_subset = resCpp$group_subset, 
              group_size = resCpp$group_size, 
              group_number = resCpp$nb_groups, 
              group_max = which.max(resCpp$group_size)              
  )
  
  return(res)
}

methods::setMethod(f  = "calcW", 
                   signature  = "data.frame", 
                   definition = function(object, range, method = "euclidean", upper = NULL, format = "dgCMatrix", row.norm = FALSE, 
                                         spatial_res = rep(1, ncol(object)), calcBlockW = FALSE)
                   { 
                     p <- ncol(object)
                     
                     #### tests
                     # initPackage(package = "spam", method = "calcW")
                     # initPackage(package = "Matrix", method = "calcW")
                     
                     validCharacter(value = format, validLength = 1, validValues = c("spam", "dgCMatrix"), refuse.NULL = TRUE, method = "calcW")
                     validNumeric(value = range, validLength = 1, min = 0, refuse.NA = TRUE, refuse.NULL = TRUE, method = "calcW")
                     
                     validDim_vector(value1 = spatial_res, value2 = list(p), name2 = NULL, method = "calcW")
                     
                     #### initialisation
                     for(iter_c in 1:p){
                       object[,iter_c] <-  object[,iter_c]*spatial_res[iter_c]
                     }
                     
                     #### computing            
                     W <- spam::nearest.dist(object, method = "euclidean", delta = range, upper = upper)
                     
                     if(format == "dgCMatrix"){
                       W <- spam::as.dgCMatrix.spam(W)
                       if(row.norm){
                         pSum <- spam::rowSums(W)
                         pSum[pSum == 0] <- -1
                         W <- W / pSum
                       }
                       W <- Matrix::drop0(W, is.Csparse = TRUE)
                     }      
                     
                     
                     if(calcBlockW == TRUE){
                       dist.center <- sqrt(rowSums(sweep(object, MARGIN = 2, STATS = apply(object, 2, median), FUN = "-")^2))
                       blocks <- calcBlockW(W = W, site_order = NULL, dist.center = dist.center, dist.max = range, verbose = FALSE)
                     }else{
                       blocks <- NULL
                     }
                     
                     #### export
                     return(list(res = list(W = W, blocks = blocks, upper = upper), 
                                 update.object = FALSE, 
                                 overwrite = FALSE)
                     )
                   }
)

#### 3- Display functions ####

legendMRI <- function(breaks, palette, mar, cex, cex.main, main, quantiles, digit){
  
  ### definition de scale 
  if(length(breaks) > 8){
    at_legend <- seq(min(breaks), max(breaks), length.out = 8)
  }else{
    at_legend <- breaks
  }
  at_legend <- signif(at_legend, digit)
  n.legend <- length(at_legend)
  n.breaks <- length(breaks)
  seq_scale <- 1:(n.breaks - 1)
  
  if(!is.null(mar)){graphics::par(mar = mar)}
  graphics::image(1, (breaks[-n.breaks] + breaks[-1]) / 2, rbind(seq_scale), col = palette, 
                  axes = FALSE, ylab = "", xlab = "")
  
  graphics::title(main, cex.main = cex.main)
  graphics::axis(2, cex.axis = cex, at = c(at_legend[-n.legend], 0.99 * at_legend[n.legend]), labels = at_legend, las = 2)            
  
  if(!identical(quantiles, FALSE)){
    graphics::abline(h = quantiles, col = c("red", "blue", "green", "blue", "red"), lwd = 3, lty = c(1, 2, 3, 2, 1))
  }
  
  # export
  return(invisible(TRUE))
}

legendMRI2 <- function(param, palette, mar, cex, cex.main){
  
  ### definition de scale 
  if(!is.null(mar)){graphics::par(mar = mar)}
  n.param <- length(param)
  seq_FALSE <- 0
  seq_TRUE <- seq(0, 1, 0.1)
  
  for(iter_param in 1:n.param){
    eval(parse(text = paste(
      "graphics::image(1, seq(0, 1, 0.1), rbind(1:10), col = ", palette, "(seq_", iter_param == 1, ", seq_", iter_param == 2, ", seq_", iter_param == 3, "), 
          axes = FALSE, ylab = \"\", xlab = \"\")", 
      sep = "")))
    graphics::title(param[iter_param], cex.main = cex.main)
    graphics::axis(2, cex.axis = cex, at = c(0, 0.25, 0.5, 0.75, 1), labels = c(0, 0.25, 0.5, 0.75, 1), las = 2)            
  }
  
  # export
  return(invisible(TRUE))
}

methods::setMethod(f  = "multiplot", 
                   signature  = "data.frame", 
                   definition = function(object, contrast = NULL, num = NULL, index1 = NULL, index2 = NULL, index3 = NULL, 
                                         col = NULL, pch = NULL, xlim = NULL, ylim = NULL, filename = "multiplot", ...){
                   
                     #### graphical options
                     ## get graphical arguments
                     optionsMRIaggr.eg <- optionsMRIaggr() 
                     dots.arguments <- list(...)
                     names_dots.arguments <- names(dots.arguments)
                     
                     validCharacter(names_dots.arguments, name = "...", validLength = NULL, 
                                    validValues = setdiff(names(optionsMRIaggr.eg), "hemisphere"),
                                    refuse.NULL = FALSE, method = "multiplot[data.frame]")
                     
                     ## set specific display
                     if(length(names_dots.arguments) > 0){
                       optionsMRIaggr.eg[names_dots.arguments] <- dots.arguments[names_dots.arguments]
                     }
                     
                     ## create all variables
                     slice_var <- optionsMRIaggr.eg$slice_var
                     norm_mu <- optionsMRIaggr.eg$norm_mu
                     norm_sigma <- optionsMRIaggr.eg$norm_sigma
                     numeric2logical <- optionsMRIaggr.eg$numeric2logical
                     breaks <- optionsMRIaggr.eg$breaks
                     type.breaks <- optionsMRIaggr.eg$type.breaks
                     palette <- optionsMRIaggr.eg$palette
                     cex <- optionsMRIaggr.eg$cex
                     col.NA <- optionsMRIaggr.eg$col.NA
                     pch.NA <- optionsMRIaggr.eg$pch.NA
                     axes <- optionsMRIaggr.eg$axes
                     window <- optionsMRIaggr.eg$window
                     legend <- optionsMRIaggr.eg$legend
                     mfrow <- optionsMRIaggr.eg$mfrow
                     mar <- optionsMRIaggr.eg$mar
                     mgp <- optionsMRIaggr.eg$mgp
                     pty <- optionsMRIaggr.eg$pty
                     asp <- optionsMRIaggr.eg$asp
                     bg <- optionsMRIaggr.eg$bg
                     xlab <- optionsMRIaggr.eg$xlab
                     ylab <- optionsMRIaggr.eg$ylab
                     main <- optionsMRIaggr.eg$main
                     num.main <- optionsMRIaggr.eg$num.main
                     cex.main <- optionsMRIaggr.eg$cex.main
                     quantiles.legend <- optionsMRIaggr.eg$quantiles.legend
                     digit.legend <- optionsMRIaggr.eg$digit.legend
                     cex.legend <- optionsMRIaggr.eg$cex.legend
                     mar.legend <- optionsMRIaggr.eg$mar.legend
                     main.legend <- optionsMRIaggr.eg$main.legend
                     outline.index <- optionsMRIaggr.eg$outline.index
                     cex.index <- optionsMRIaggr.eg$cex.index
                     pch.index <- optionsMRIaggr.eg$pch.index
                     col.index <- optionsMRIaggr.eg$col.index
                     filter.index <- optionsMRIaggr.eg$filter.index
                     width <- optionsMRIaggr.eg$width
                     height <- optionsMRIaggr.eg$height
                     path <- optionsMRIaggr.eg$path
                     unit <- optionsMRIaggr.eg$unit
                     res <- optionsMRIaggr.eg$res
                     
                     #### test ####
                     validCharacter(legend, validLength = NULL, validValues = c(FALSE, TRUE, "only"), refuse.NULL = FALSE, method = "multiplot")
                     
                     if(is.null(main.legend)){
                       main.legend <- paste("", names(contrast), sep = "")
                     }   
                     
                     if(is.data.frame(contrast)){
                       param <- names(contrast)
                       contrast <- as.matrix(contrast)
                     }else if(is.matrix(contrast)){
                       param <- NULL
                     }else if(is.vector(contrast)){
                       contrast <- as.matrix(contrast)
                       param <- NULL
                     }else if(is.null(contrast)){
                       contrast <- matrix(1, nrow = nrow(object), ncol = 1)
                       param <- NULL
                     }else{
                       stop("multiplot[data.frame] : wrong specification of \'contrast\' \n", 
                            "must be either NULL or data.frame or matrix or vector \n", 
                            "is(contrast) : ", paste(is(contrast), collapse = " "), "\n")
                     }
                     
                     if(ncol(object) == 2){
                       object <- cbind(object, 1)
                     }
                     
                     if(ncol(object) != 3){
                       stop("multiplot[data.frame] : \'object\' must be a 3 column matrix with the observation coordinates \n", 
                            "ncol(object) : ", ncol(object), "\n")
                     }
                   
                     validDim_matrix(value1 = contrast, value2 = object, name1 = "contrast", name2 = "object", type = "nrow", method ="multiplot[data.frame]")
                     
                     if(is.null(num)){num <- unique(object[,3])}
                     
                     if(any(num %in% object[,3]) == FALSE){
                       stop("multiplot[data.frame] : \'object\' inconsistent with \'num\' \n", 
                            "requested slices : ", paste(num, collapse = " "), " \n", 
                            "slices in \'object\' : ", paste(unique(object[,3]), collapse = " "), " \n")
                     }
                     
                     n.num <- length(num)
                     if(!is.null(legend)){
                       if(legend == FALSE){n.plot <- n.num}
                       if(legend == TRUE){n.plot <- n.num + 1}
                       if(legend == "only"){n.plot <- 1}
                     }else{
                       n.plot <- n.num
                     }
                     
                     contrast <- contrast[num %in% object[,3],, drop = FALSE]
                     object <- object[num %in% object[,3],]
                     n.px <- nrow(contrast)
                     n.contrast <- ncol(contrast)
                     
                     #### initialization and tests ####
                     
                     ## windows
                     par.init <- graphics::par()
                     
                     res.init <- initWindow(window = window, filename = filename, path = path, width = width, height = height, unit = unit, res = res, 
                                            n.plot = n.plot, mfrow = mfrow, xlim = xlim, ylim = ylim, 
                                            method = "multiplot[data.frame]")
                     scale <- res.init$scale
                     mfrow <- res.init$mfrow
                     n.graph_par_window <- res.init$n.graph_par_window
                     xlim.plot <- res.init$xlim.plot
                     ylim.plot <- res.init$ylim.plot
                     
                     ## color and breaks   
                     res.init <- initCol(contrast = contrast, coords = object, param = param, pch = pch, col = col, palette = palette, breaks = breaks, legend = legend, type.breaks = type.breaks, 
                                         method = "multiplot[data.frame]")
                     
                     contrast <- res.init$contrast 
                     palette <- res.init$palette 
                     breaks <- res.init$breaks
                     
                     col <- res.init$col
                     pch <- res.init$pch
                     palette_sauve <- res.init$palette_sauve
                     breaks_sauve <- res.init$breaks_sauve
                     index_duplicated <- res.init$index_duplicated
                     index_order <- res.init$index_order
                     
                     ## index            
                     if(!is.null(index1)){
                       res.init <- initIndex(object = NULL, index = index1, num = num, hemisphere = "both", numeric2logical = numeric2logical, 
                                             method = "multiplot[data.frame]", indexNum = 1, outline.default = outline.index, 
                                             cex.default = cex.index[1], pch.default = pch.index[1], col.default = col.index[1], filter.default = filter.index)
                       
                       index1 <- res.init$coords
                       indexindex1 <- res.init$index
                       pch_index1 <- res.init$pch
                       cex_index1 <- res.init$cex
                       col_index1 <- res.init$col
                     }
                     
                     if(!is.null(index2)){
                       res.init <- initIndex(object = NULL, index = index2, num = num, hemisphere = "both", numeric2logical = numeric2logical, 
                                             method = "multiplot[data.frame]", indexNum = 2, outline.default = outline.index, 
                                             cex.default = cex.index[2], pch.default = pch.index[2], col.default = col.index[2], filter.default = filter.index)
                       
                       index2 <- res.init$coords
                       indexindex2 <- res.init$index
                       pch_index2 <- res.init$pch
                       cex_index2 <- res.init$cex
                       col_index2 <- res.init$col
                     }
                     
                     if(!is.null(index3)){
                       test.intersection1 <- length(index3) == 1 && index3 == "intersection"
                       test.intersection2 <- is.list(index3) && "coords" %in% names(index3) && length(index3$coords) == 1 && index3$coords == "intersection"
                       
                       if(test.intersection1 || test.intersection2){
                         
                         if(test.intersection1){index3 <- list(coords = "intersection")}
                         if(is.null(index1)){
                           stop("multiplot[data.frame] : \'index3\' cannot request \"interaction\" parameter if \'index1\' is missing  \n")
                         }
                         if(is.null(index2)){
                           stop("multiplot[data.frame] : \'index3\' cannot request \"interaction\" parameter if \'index2\' is missing  \n")
                         }
                         
                         test1.intersect <- indexindex1 %in% indexindex2
                         test2.intersect <- indexindex2 %in% indexindex1
                         
                         index1_sauve <- index1
                         index3$coords <- index1_sauve$coords[test1.intersect == TRUE, c("i", "j", "k")]
                         index1 <- index2$coords[test1.intersect == FALSE, c("i", "j", "k")]
                         index2 <- index1_sauve$coords[test2.intersect == FALSE, c("i", "j", "k")]
                         
                       }
                       
                       res.init <- initIndex(object = NULL, index = index3, num = num, hemisphere = "both", numeric2logical = numeric2logical, 
                                             method = "multiplot[data.frame]", indexNum = 3, outline.default = outline.index, 
                                             cex.default = cex.index[3], pch.default = pch.index[3], col.default = col.index[3], filter.default = filter.index)
                       
                       index3 <- res.init$coords   
                       pch_index3 <- res.init$pch
                       cex_index3 <- res.init$cex
                       col_index3 <- res.init$col     
                       
                     }
                     
                     
                     #### display plot ####
                     compteur <- 1
                     
                     if(is.null(legend) || legend != "only"){
                       for(iter_num in 1:n.num){
                         
                         # device
                         if(!is.null(window) && compteur == 1){
                           
                           if(window %in% c("png", "eps", "svg", "pdf") && iter_num > 1){grDevices::dev.off()}
                           filename_all <- paste(filename, "_", paste(param, collapse = "-"), "(slice", num[iter_num], "-", min(max(num), num[iter_num + n.graph_par_window - 1], na.rm = TRUE), ")", sep = "")                  
                           
                           initDisplayWindow(window = window, filename = filename_all, path = path, width = width, height = height, scale = scale, res = res, 
                                             mfrow = mfrow, bg = bg, pty = pty, mar = mar, mgp = mgp, 
                                             n.contrast = if(identical(legend, TRUE) && ((n.num-iter_num) < prod(mfrow) ) ){n.contrast}else{1})
                         }
                         
                         # data
                         index_k <- which(object[,3] == num[iter_num])
                         contrastK <- contrast[index_k, , drop = FALSE]
                         coordsK <- object[index_k,, drop = FALSE]
                         if(is.null(col)){colK <- NULL}else{colK <- col[index_k]}
                         mainK <- paste(main, 
                                        if(num.main == TRUE){paste(" ", num[iter_num], sep = "")}, 
                                        if(!is.null(param)){paste(" (", paste(param, collapse = "-"), ")", sep = "")}, sep = "")
                         
                         # xlim-ylim
                         if(is.null(xlim)){                
                           if(is.null(coordsK) || nrow(coordsK) == 0){
                             xlim.plot <- NULL
                           }else{
                             xlim.plot <- c(min(coordsK[,1]) - 0.5, max(coordsK[,1]) + 0.5)
                           }                 
                         }
                         
                         if(is.null(ylim)){
                           if(is.null(coordsK) || nrow(coordsK) == 0){
                             ylim.plot <- NULL
                           }else{
                             ylim.plot <- c(min(coordsK[,2]) - 0.5, max(coordsK[,2]) + 0.5)
                           }
                         }
                         
                         plot.test <- plotMRI(contrast = contrastK, coords = coordsK[,1:2], breaks = breaks, col = colK, palette = palette, 
                                              asp = asp, 
                                              xlim = xlim.plot, ylim = ylim.plot, pch = pch, cex = cex, axes = axes, col.NA = col.NA, pch.NA = pch.NA, 
                                              xlab = xlab, ylab = ylab, main = mainK, cex.main = cex.main)
                         
                         if(!is.null(index1) && any(index1[,3] == num[iter_num])){
                           graphics::points(index1[index1[,3] == num[iter_num],1:2, drop = FALSE], 
                                            pch = pch_index1, cex = cex_index1, col = col_index1)
                         }
                         
                         if(!is.null(index2) && any(index2[,3] == num[iter_num])){
                           graphics::points(index2[index2[,3] == num[iter_num],1:2, drop = FALSE], 
                                            pch = pch_index2, cex = cex_index2, col = col_index2)
                         }
                         
                         if(!is.null(index3) && any(index3[,3] == num[iter_num])){
                           graphics::points(index3[index3[,3] == num[iter_num],1:2, drop = FALSE], 
                                            pch = pch_index3, cex = cex_index3, col = col_index3)
                         }
                         
                         compteur <- compteur + 1
                         if(compteur > n.graph_par_window)
                         {compteur <- 1}
                         
                       }            
                     }
                     
                     #### legend ####
                     if(is.null(legend)){
                       compteur <- 1
                       mfrow <- c(1, 1)
                       legend <- TRUE
                     }
                     
                     if(legend %in% c(TRUE, "only")){
                       
                       if(!is.null(window) && compteur == 1){
                         
                         if(window %in% c("png", "eps", "svg", "pdf") && iter_num > 1){grDevices::dev.off()}
                         filename_all <- paste(filename, "_", paste(param, collapse = "-"), "(slice", min(num), "-", max(num), ") - legend", sep = "")
                         
                         initDisplayWindow(window = window, filename = filename_all, path = path, width = width, height = height, scale = scale, res = res, 
                                           mfrow = mfrow, bg = bg, pty = pty, mar = NULL, mgp = mgp, n.contrast = n.contrast)
                         
                       }
                       
                       
                       if(n.contrast == 1){
                         if(quantiles.legend == TRUE){
                           quantiles.legend <- stats::quantile(contrast[,1], na.rm = TRUE)
                         }else{
                           quantiles.legend <- quantiles.legend
                         }   
                         
                         if(length(unique(breaks_sauve)) == 1){
                           
                           if(is.null(col)){
                             breaks_sauve <- c(breaks_sauve - 10^{-12}, breaks_sauve, breaks_sauve + 10^{-12})
                           }else{
                             breaks_sauve <- 1:length(unique(col))
                           }
                           
                         }
                         
                         plot.test <- legendMRI(breaks = breaks_sauve, palette = palette_sauve, mar = mar.legend, 
                                                cex = cex.legend, main = paste(main.legend, collapse = " "), cex.main = cex.main, quantiles = quantiles.legend, digit = digit.legend)
                         
                       }else{
                         if(is.null(param)){param <- 1:n.contrast}
                         plot.test <- legendMRI2(param = param, palette = palette, mar = mar, cex = cex.legend, cex.main = cex.main)
                         
                       }              
                     }
                     
                     
                     
                     if(!is.null(window) && window %in% c("eps", "svg", "png", "pdf")){
                       grDevices::dev.off()
                     }
                     graphics::par(bg = par.init$bg, pty = par.init$pty, mar = par.init$mar, mgp = par.init$mgp)
                     
                     #### export 
                     return(invisible(list(breaks.plot = breaks, 
                                           palette.plot = palette, 
                                           breaks.legend = breaks_sauve, 
                                           palette.legend = palette_sauve, 
                                           quantiles.legend = quantiles.legend)))
                   }
) 

outline <- function(n = 50, sequential = TRUE, min_dist = 1, 
                    col = c("blue", "red", "grey"), pch = 20, cex = c(0.75, 1, 0.75)){
  
  # initPackage(package = "RANN", method = "outline")
  
  #### one shot outline 
  if(sequential == FALSE){
    res_locator <- graphics::locator(n = n)
    res_locator$x <- c(res_locator$x, res_locator$x[1])
    res_locator$y <- c(res_locator$y, res_locator$y[1])
    n <- length(res_locator$x)
    
    # round
    res_locator$x <- round(res_locator$x)
    res_locator$y <- round(res_locator$y)  
    
    
    df_points <- matrix(NA, nrow = 0, ncol = 4)
    
    for(iter_point in 1:(n - 1)){
      n_tempo <- max(1 + abs(res_locator$x[iter_point]-res_locator$x[iter_point + 1]), 
                     1 + abs(res_locator$y[iter_point]-res_locator$y[iter_point + 1]))
      n_tempo <- ceiling(n_tempo)
      i_tempo <- seq(res_locator$x[iter_point], res_locator$x[iter_point + 1], length.out = n_tempo)
      j_tempo <- seq(res_locator$y[iter_point], res_locator$y[iter_point + 1], length.out = n_tempo)
      
      i_tempo <- round(i_tempo)
      j_tempo <- round(j_tempo)    
      
      if(n_tempo > 1){
        df_points <- rbind(df_points, 
                           cbind(i = i_tempo[-n_tempo], j = j_tempo[-n_tempo], edge = iter_point, points = c(1, rep(0, n_tempo-2))))
      }
    }
    
  }
  
  
  #### sequential outline
  if(sequential == TRUE){
    
    iter <- 1
    dist <- min_dist + 1
    df_points <- matrix(NA, nrow = 0, ncol = 4)
    
    while(iter <= n && dist > min_dist){
      res_locator <- graphics::locator(n = 1)
      
      # if user enter echap directly
      if(is.null(res_locator)){
        return(list(edge = NULL, 
                    surface = NULL)
        )
      }
      
      if(iter > 1){
        dist <- sqrt( (df_points[1, "i"] - res_locator$x)^2 + (df_points[1, "j"] - res_locator$y)^2)
        
        res_locator$x <- c(df_points[nrow(df_points), "i"], round(res_locator$x))
        res_locator$y <- c(df_points[nrow(df_points), "j"], round(res_locator$y))
        
        n_tempo <- max(1 + abs(res_locator$x[1] - res_locator$x[2]), 
                       1 + abs(res_locator$y[1] - res_locator$y[2]))
        n_tempo <- ceiling(n_tempo)
        if(n_tempo == 1){i_tempo <- c() ; next} # cas ou l on selectionne le meme points
        
        i_tempo <- round(seq(res_locator$x[1], res_locator$x[2], length.out = n_tempo)[-1])
        j_tempo <- round(seq(res_locator$y[1], res_locator$y[2], length.out = n_tempo)[-1])
        
        points <- c(rep(0, n_tempo - 2), 1)
      }else{
        
        # round
        res_locator$x <- round(res_locator$x)
        res_locator$y <- round(res_locator$y)  
        
        i_tempo <- res_locator$x
        j_tempo <- res_locator$y
        points <- 1  
      }
      
      if(dist < 1){
        if(length(i_tempo) > 1){
          df_points <- rbind(df_points, 
                             cbind(i = i_tempo[ - (n_tempo - 1)], j = j_tempo[ - (n_tempo - 1)], edge = iter, points = 0))
        }
      }else{
        df_points <- rbind(df_points, 
                           cbind(i = i_tempo, j = j_tempo, edge = iter, points = points))
      }
      
      graphics::points(i_tempo, j_tempo, 
                       col = col[points + 1], pch = pch, cex = cex[points + 1])
      
      
      
      iter <- iter + 1
    }
    
    # relier au premier point   
    n_tempo <- max(1 + abs(df_points[nrow(df_points),"i"] - df_points[1,"i"]), 
                   1 + abs(df_points[nrow(df_points),"j"] - df_points[1,"j"]))
    n_tempo <- ceiling(n_tempo)
    
    if(n_tempo > 2){
      i_tempo <- round(seq(df_points[nrow(df_points),"i"], df_points[1,"i"], length.out = n_tempo))
      j_tempo <- round(seq(df_points[nrow(df_points),"j"], df_points[1,"j"], length.out = n_tempo))
      
      df_points <- rbind(df_points, 
                         cbind(i = i_tempo[c(-1, -n_tempo)], j = j_tempo[c(-1, -n_tempo)], edge = iter, points = 0))
      
      graphics::points(i_tempo[c(-1, -n_tempo)], j_tempo[c(-1, -n_tempo)], 
                       col = col[points + 1], pch = pch, cex = cex[points + 1])
    }
  }
  df_points <- as.data.frame(df_points)
  
  # enlever d eventuels dupliques
  df_points <- df_points[duplicated(df_points[,c("i", "j")]) == FALSE,]
  
  
  #### filling the outline
  df_points <- cbind(df_points, valid = TRUE)
  j_level <- unique(df_points[,"j"])
  i_tempo <- c()
  j_tempo <- c()  
  
  for(iter_j in 1:length(j_level)){
    
    index_j <- which(df_points[,"j"] == j_level[iter_j])
    
    if(length(index_j) == 1){
      i_tempo <-  c(i_tempo, df_points[index_j, "i"])
      j_tempo <-  c(j_tempo, df_points[index_j, "j"])
    }else{
      index_j <- index_j[order(df_points[index_j, "i"])]
      
      # enlever les doublons d un meme trait
      diff_tempo <- c(FALSE, abs(diff(index_j)) %in% c(1, nrow(df_points) - 1))
      
      while(TRUE %in% diff_tempo && length(diff_tempo) > 2){   
        pos_FALSE <- which(diff_tempo == TRUE)[1]
        if(pos_FALSE %% 2 == 0){
          df_points$valid[index_j[pos_FALSE]] <- FALSE
          index_j <- index_j[-pos_FALSE]          
        }else{          
          df_points$valid[index_j[pos_FALSE - 1]] <- FALSE
          index_j <- index_j[ - (pos_FALSE - 1)]
        }
        diff_tempo <- diff_tempo[-pos_FALSE]
      }
      
      # enlever les V inverses
      df_Vpoints <- cbind(index = 1:nrow(df_points), df_points)[df_points$valid == TRUE,]
      index_Vj <- sapply(index_j, function(x){which(df_Vpoints$index == x)})
      index_m1 <- c(nrow(df_Vpoints), 1:nrow(df_Vpoints) - 1)
      index_p1 <- c(2:nrow(df_Vpoints), 1)
      
      diff_tempo <- rbind(df_Vpoints[,"i"]-df_Vpoints[index_m1,"i"], 
                          df_Vpoints[,"j"]-df_Vpoints[index_m1,"j"], 
                          df_Vpoints[,"i"]-df_Vpoints[index_p1,"i"], 
                          df_Vpoints[,"j"]-df_Vpoints[index_p1,"j"])[,index_Vj, drop = FALSE]
      
      test_V <- apply(diff_tempo, 2, function(x){
        test.H <- all(x == c(1,-1,-1,-1)) + all(x == c(0,-1,-1,-1)) + all(x == c(0,-1,0,-1)) + all(x == c(1,-1,0,-1))
        test.antiH <- all(x == c(-1,-1,1,-1)) + all(x == c(-1,-1,0,-1)) + all(x == c(0,-1,1,-1))
        return(test.H + test.antiH)       
      }) 
      
      test_Vinv <- apply(diff_tempo, 2, function(x){
        test.H <- all(x == c(1,1,-1,1)) + all(x == c(0,1,-1,1)) + all(x == c(0,1,0,1)) + all(x == c(1,1,0,1))
        test.antiH <- all(x == c(-1,1,1,1)) + all(x == c(-1,1,0,1)) + all(x == c(0,1,1,1))
        return(test.H + test.antiH)        
      }) 
      
      i_tempo <- c(i_tempo, df_points[index_j[test_V == 1], "i"])
      j_tempo <- c(j_tempo, df_points[index_j[test_V == 1], "j"])
      index_j <- index_j[c(test_Vinv + test_V) == 0]
      
      n.i <- length(index_j)
      
      # remplissage
      if(n.i > 0){
        for(iter_i in 1:floor(n.i / 2)){
          seq_i <- seq(df_points[index_j[2 * (iter_i - 1) + 1], "i"], df_points[index_j[2 * iter_i], "i"], by = 1)
          i_tempo <-  c(i_tempo, seq_i)
          j_tempo <-  c(j_tempo, rep(j_level[iter_j], length(seq_i)))
        }
      }
      
    }
    
  }
  
  df_fill <- data.frame(i = i_tempo, j = j_tempo)
  
  if(sequential == TRUE){
    test <- RANN::nn2(data = df_points[,c("i", "j")], query = df_fill[,c("i", "j")], k = 1)
    graphics::points(i_tempo[test$nn.dists > 0], j_tempo[test$nn.dists > 0], 
                     col = col[3], pch = pch, cex = cex[3])
  }
  
  
  return(list(edge = df_points, 
              surface = df_fill)
  )
}

plotMRI <- function(contrast, coords, breaks, palette, col, asp, 
                    xlim, ylim, pch, cex, axes, col.NA, pch.NA, xlab, ylab, main, cex.main){
  
  #### gestion des coupes ####
  if(is.null(main) && !is.null(names(contrast)))
  {main <- names(contrast)}
  
  if(nrow(contrast) == 0 || sum(is.na(contrast) == FALSE) == 0 || is.null(contrast)){
    
    if(is.logical(xlim) || is.logical(ylim))
    {graphics::plot(1, 1, col = col, xlab = "", ylab = "", axes = axes)}else
    {graphics::plot(1, 1, col = col, xlim = xlim, ylim = ylim, xlab = "", ylab = "", axes = axes)}
    graphics::legend("center", c("Missing", "data"), cex = 2, pch = 62, bty = "n")
    
    graphics::title(main = main)
    graphics::axis(1)
    graphics::axis(2)
    
    return(invisible(FALSE))
  }
  
  #### mis en forme des donnees ####
  if(is.null(col)){
    array <- df2array(contrast = contrast, coords = coords, format = "matrix", default_value = NA)
    spdf <- array$contrast
    coords_x <- array$unique_coords[[1]]
    coords_y <- array$unique_coords[[2]]
  }
  
  #### taille du pixel ####
  if(!is.null(col)){
    
    unique_x <- sort(unique(coords[,1]))
    n_x <- (unique_x[length(unique_x)] - unique_x[1]) # /  min(diff(unique_x))
    
    unique_y <- sort(unique(coords[,2]))
    n_y <- (unique_y[length(unique_y)] - unique_y[1]) # /  min(diff(unique_y))
    
    cex <- max(3 * graphics::par("pin") / (graphics::par("cin") * max(n_y,n_x))) * cex
    #     cex = 100 * min(graphics::par()$plt[2] - graphics::par()$plt[1], graphics::par()$plt[4] - graphics::par()$plt[3]) / max(xlim[2] - xlim[1], ylim[2] - ylim[1]) 
  }
  
  #### affichage #### 
  if(is.null(col)){
    eval(parse(text = paste(
      "graphics::image(x = coords_x, y = coords_y, z = spdf, 
      breaks = breaks, col = palette, 
      ", if(!is.null(xlim)){"xlim = xlim, "}, if(!is.null(ylim)){"ylim = ylim, "}, 
      "xlab = xlab, ylab = ylab, axes = axes, asp = asp)", 
      sep = "")))
    
  }else{
    
    graphics::plot(coords, 
                   pch = pch, cex = cex, col = col, xlim = xlim, ylim = ylim, 
                   axes = axes, asp = asp, xlab = xlab, ylab = ylab)
    
  }
  
  graphics::title(main = main, cex.main = cex.main)
  # affichage des NA
  if(!is.null(col.NA)){
    test_NA <- which(is.na(contrast))
    
    if(sum(test_NA) > 0){
      graphics::points(coords[test_NA,], col = col.NA, pch = pch.NA, cex = cex[1] / 5)
    }
    
  }
  
  # export
  return(invisible(TRUE))
  
}

pointsOutline <- function(coords, array = NULL, filter = "2D_N4"){

  if(is.null(array)){
    p <- ncol(coords)
    
    # order voxels by coordinates
    names_coords <- names(coords)
    if(is.null(names_coords)){names_coords <- letters[9:(8 + p)]}
    eval(parse(text = paste("coords <- coords[with(coords, order(", paste(names_coords[length(names_coords):1], collapse = ", "), ")),]", sep = "")))
    
    n.neighbors <- as.numeric(paste(strsplit(filter, split = "", fixed = TRUE)[[1]][ - (1:4)], collapse = ""))
    if(nrow(coords) <= n.neighbors){
      return(coords)
    }      

    array <- df2array(contrast = rep(1, nrow(coords)), 
                      coords = coords)$contrast[[1]]    
  }else{
    coords <- array2df(array)[,1:length(dim(array))]
  }
 
  res <- calcFilter(array, filter = filter, norm.filter = FALSE)
  index_array <- which(!is.na(array))
  index_outline <- which(res$res < max(res$res, na.rm = TRUE))
  
  coords <- coords[index_array %in% index_outline,]
  
  return(coords)
  
}


#### 4- Filtering functions ####

methods::setMethod(f  = "calcFilter", 
                   signature  = "array", 
                   definition = function(object, filter, norm.filter = TRUE, bilateral = FALSE, na.rm = FALSE){
                     
                     #### test
                     validLogical(value = norm.filter, validLength = 1, refuse.NULL = TRUE, refuse.NA = TRUE, method = "calcFilter")
                     validLogical(value = bilateral, validLength = 1, refuse.NULL = TRUE, refuse.NA = TRUE, method = "calcFilter")
                     validLogical(value = na.rm, validLength = 1, refuse.NULL = TRUE, refuse.NA = TRUE, method = "calcFilter")
                     
                     if(length(dim(object)) %in% c(2, 3) == FALSE){
                       stop("calcFilter[array] : wrong specification of \'object\' \n", 
                            "\'object\' must be an array of dimension 2 or 3 \n", 
                            "dim(object) : ", paste(dim(object), collapse = " "), "\n")
                     }
                     
                     #### Filtres disponibles 
                     if(length(filter) == 1 && is.character(filter)){
                       
                       filter_split <- strsplit(filter, split = "")[[1]]
                       
                       if(filter_split[[4]] == "N"){
                         
                         res <- initNeighborhood(filter, method = "calcFilter")
                         filter_split <- list(ncol(res), "D", "_", filter_split[[4]], nrow(res))
                         
                         if(filter_split[[1]] == 2){
                           filter <- matrix(NA, nrow = 3, ncol = 3)
                           for(iter in 1:filter_split[[5]]){
                             filter[res[iter, 1] + 2, res[iter, 2] + 2] <- 1  
                           }
                         }else{
                           filter <- array(NA, dim = c(3, 3, 3))
                           for(iter in 1:filter_split[[5]]){
                             filter[res[iter, 1] + 2, res[iter, 2] + 2, res[iter, 3] + 2] <- 1  
                           }
                         }
                         
                       }else{
                         
                         res <- initFilter(filter, method = "calcFilter")
                         filter <- res$filter
                         filter_split <- res$filter_split
                         
                       }
                     }else{ # Perso
                       filter_split <- c(length(dim(filter)), "D", "_", "P", dim(filter)[1])
                     }
                     
                     if(filter_split[[4]] == "S" && na.rm == FALSE){
                       warning("calcFilter[array] : gradient values may be incorrect at the edges \n", 
                               "set \'na.rm\' to TRUE to remove these values \n")
                     }
                     if(filter_split[[4]] == "M" && bilateral == TRUE){
                       warning("calcFilter[array] : there is no edge preserving median filtering \n", 
                               "\'bilateral\' will be ignored \n")
                     }
                     
                     ### preparation
                     p <- dim(filter)
                     p.dim <- length(p)
                     p_ref <- sapply(1:p.dim, function(x){stats::median(1:p[x])})
                     if(p.dim > 3){
                       stop("calcFilter[array] : wrong specification of \'filter\' \n", 
                            "can only handle 1D 2D and 3D filters \n", 
                            "dimension of the proposed filter : ", p.dim, "\n")
                     }
                     
                     M.dim <- dim(object)
                     Mres <- array(NA, dim = dim(object))
                     
                     Ind.operateur_compr <- which(as.vector(filter) != 0)
                     Vec.operateur_compr <- as.vector(filter)[Ind.operateur_compr]
                    
                     # filtrage 2D
                     if(filter_split[[4]] %in% c("N", "G", "I", "P", "S") && length(M.dim) == 2){
                       
                       resCpp <- filtrage2D_cpp(M_data = object, 
                                                M_operateur = filter, 
                                                index_data = which(!is.na(object), arr.ind = TRUE) - 1, 
                                                bilateral = bilateral, 
                                                na_rm = na.rm)
                       
                       if(norm.filter == TRUE){
                         Mres <- resCpp$Mres / resCpp$Wres
                       }else{
                         Mres <- resCpp$Mres
                       }
                       
                     }
                     
                     
                     # filtrage 3D
                     if(filter_split[[4]] %in% c("N", "G", "I", "P", "S") && length(M.dim) == 3){
                       if(filter_split[[1]] == 2){
                         for(iter_k in 1:(M.dim[3])){
                    
                           resCpp <- filtrage2D_cpp(M_data = as.matrix(object[,,iter_k, drop = TRUE]), 
                                                    M_operateur = filter, 
                                                    index_data = which(!is.na(object[,,iter_k, drop = FALSE]), arr.ind = TRUE) - 1, 
                                                    bilateral = bilateral, 
                                                    na_rm = na.rm)
                           
                           if(norm.filter == TRUE){
                             Mres[,, iter_k] <- resCpp$Mres / resCpp$Wres
                           }else{
                             Mres[,, iter_k] <- resCpp$Mres
                           }      
                         }
                       }else{  
                         
                         resCpp <- filtrage3D_cpp(Vec_data = as.vector(object), p_data = dim(object), 
                                                  Vec_operateur = as.vector(filter), p_operateur = dim(filter), 
                                                  index_data = which(!is.na(object), arr.ind = TRUE) - 1, 
                                                  bilateral = bilateral, 
                                                  na_rm = na.rm)                
                         
                         if(norm.filter == TRUE){
                           Mres <- resCpp$Mres / resCpp$Wres
                         }else{
                           Mres <- resCpp$Mres
                         }
                       }
                       
                     }           
                     
                     # filtrage median 2D
                     if(filter_split[[4]] == "M" && length(M.dim) == 2){
                       Mres <- filtrage2Dmed_cpp(M_data = object, 
                                                 M_operateur = filter, 
                                                 index_data = which(!is.na(object), arr.ind = TRUE) - 1, 
                                                 na_rm = na.rm)
                     }
                     
                     # filtrage median 3D
                     if(filter_split[[4]] == "M" && length(M.dim) == 3){
                       if(filter_split[[1]] == 2){
                         for(iter_k in 1:(M.dim[3])){
                           Mres[,, iter_k] <- filtrage2Dmed_cpp(M_data = as.matrix(object[,,iter_k, drop = TRUE]), 
                                                                M_operateur = filter, 
                                                                index_data = which(!is.na(object[,,iter_k, drop = FALSE]), arr.ind = TRUE) - 1, 
                                                                na_rm = na.rm)
                         }
                       }else{
                         
                         Mres <- filtrage3Dmed_cpp(Vec_data = as.vector(object), p_data = dim(object), 
                                                   Vec_operateur = as.vector(filter), p_operateur = dim(filter), 
                                                   index_data = which(!is.na(object), arr.ind = TRUE) - 1, 
                                                   na_rm = na.rm) 
                       }
                     }
                     
                     
                     ### export
                     return(list(res = Mres, 
                                 filter = filter, 
                                 update.object = FALSE, 
                                 overwrite = FALSE)
                     )
                   }
)

#### 5- init functions ####

initCol <- function(contrast, coords, param = NULL, pch, col, palette, breaks, legend, type.breaks, method){
  
  # > image is used in all but three cases
  # col has been specifed by the user
  # there are more than one parameter (color have to be defined with a multiparametric palette)
  # user has specified a pch value
  
  p <- ncol(contrast)
  n.px <- nrow(contrast)
  if(is.null(param)){param <- names(contrast)}
  
  #### tests ####
  validCharacter(value = type.breaks, validLength = NULL, validValues = c("range", "quantile", "range_center"), method = method)  
  validCharacter(value = legend, validLength = 1, validValues = c(TRUE,FALSE,"only"), refuse.NULL = FALSE, method = method)
  
  if(!is.numeric(breaks)){
    stop(method, " : wrong specificaiton of \'breaks\' \n", 
         "\'breaks\' must be either a numeric indicating the number of breaks or a numeric vector containing the break values \n", 
         "proposed \'break\' is(break) : ", paste(is(breaks), collapse = " "), "\n")
  }
  
  if((p %in% 1:3) == FALSE){
    if(method == "multiplot[MRIaggr]"){
      stop(method, " : wrong specification of \'param\' \n", 
           "\'param\' must contains between 1 and 3 parameters \n", 
           "proposed parameters (nb) : ", paste(names(contrast), collapse = " "), " (", ncol(contrast), ") \n")  
    }else{
      stop(method, " : wrong specification of \'contrast\' \n", 
           "\'contrast\' must have between 1 and 3 columns \n", 
           "columns of \'data\' (ncol) : ", paste(names(contrast), collapse = " "), " (", ncol(contrast), ") \n")    
    }
  }
  
  
  #### dealing with multiple parameters ####
  index_duplicated <- NULL
  index_order <- NULL
  
  if(p %in% c(2:3)){
    
    validCharacter(value = palette, validLength = 1, validValues = c("rgb","hsv"), refuse.NULL = TRUE, method = method)
    
    if(any(contrast > 1) || any(contrast < 0)){
      if(method == "multiplot[MRIaggr]"){
        stop(method, " : wrong specification of \'param\' \n", 
             "parameters must take value in [0;1] if several parameters are intended to be displayed on the same map \n",                   
             "current range : ", paste(range(contrast), collapse = " "), "\n")
      }else{
        stop(method, " : wrong specification of \'contrast\' \n", 
             "contrast must contains values be in [0;1] if several parameters are intended to be displayed on the same map \n",                   
             "current range : ", paste(range(contrast), collapse = " "), "\n")
      }
    }
    
    if(ncol(contrast) == 2){      
      col <- eval(parse(text = paste("grDevices::", palette, "(contrast[,1], contrast[,2], 1)", sep = "")))      
      contrast <- contrast[,1, drop = FALSE]
    }else{
      
      col <- eval(parse(text = paste("grDevices::", palette, "(contrast[,1], contrast[,2], contrast[,3])", sep = "")))
      
      # ordering dataset by Tier
      order1 <- order(contrast[,1], decreasing = TRUE)[seq(1, n.px, length.out = round(n.px / 3))]            
      order2 <- setdiff(order(contrast[,2], decreasing = TRUE), order1)[seq(1, n.px - round(n.px / 3), length.out = round(n.px / 3))]
      order3 <- setdiff(order(contrast[,3], decreasing = TRUE), c(order1, order2))[seq(1, n.px - 2 * round(n.px / 3), length.out = n.px - 2 * round(n.px / 3))]
      
      contrast[order2, 1] <- contrast[order2, 2] + 1
      contrast[order3, 1] <- contrast[order3, 3] + 2 * 1
      contrast <- contrast[,1, drop = FALSE]            
      
      index_duplicated <- which(duplicated(contrast[,1]) == FALSE)
      index_order <- order(contrast[index_duplicated,1])
    }
    
  }
  
  #### come back to one parameter
  if(is.null(col)){
    
    ## test palette
    available_palette <- c("grey.colors", "gray.colors", "rainbow", "heat.colors", "terrain.colors", "topo.colors", "cm.colors", "tim.colors")
    
    if(length(palette) == 1 && palette %in% available_palette == FALSE){
      unique_tempo <-  unique(contrast[,1])
      if(length(unique_tempo) == 1){
        palette <- rep("red", 2)
        breaks <- c(unique_tempo - 10^{-12}, unique_tempo, unique_tempo + 10^{-12})
      }else{
        stop(method, " : wrong specification of \'palette\' \n", 
             "available palette for 1 parameter : \"", paste(available_palette, collapse = "\" \""), "\" \n", 
             "proposd palette : ", palette, "\n")
      }
    }
    
    ## infinte values
    index.infinite <- which(is.infinite(contrast[,1]))
    
    if(length(index.infinite) > 0){
      maxInf <-  max(c(contrast[-index.infinite,1], 99999))
      if(method == "multiplot[MRIaggr]"){
        warning(method, " : \'param\' values in \'object\' contains Inf values \n", 
                "they are set to ", maxInf, "\n")     
      }else{
        warning(method, " : \'contrast\' values contains Inf values \n", 
                "they are set to ", maxInf, "\n")          
      }
      contrast[index.infinite,1] <- maxInf   
    }
    
    ## breaks 
    if(length(breaks) == 1){
      range_data <- range(stats::na.omit(contrast[,1]))
      
      if(length(palette) > 1){breaks <- length(palette) + 1}
      
      breaks_sauve <- switch(type.breaks, 
                             "range" = seq(range_data[1], range_data[2], length.out = min(500, breaks)), 
                             "range_center" = seq(-max(abs(range_data)), max(abs(range_data)), length.out = min(500, breaks)), 
                             "quantile" = stats::quantile(contrast[,1], probs = seq(0, 1, length.out = min(500, breaks)))
      )
      breaks_sauve <- unique(breaks_sauve)
      #if(length(breaks_sauve) == 1){breaks_sauve <- c(breaks_sauve - 1 * 10^{-12}, breaks_sauve, breaks_sauve + 1 * 10^{-12})}
      
      breaks <- breaks_sauve
      if(length(breaks) == 1){breaks <- c(breaks - 10^{-12}, breaks, breaks + 10^{-12})}
      breaks[1] <- breaks[1] - 10^{-12}
      breaks[length(breaks)] <- breaks[length(breaks)] + 10^{-12}
      
    }else{
      breaks <- sort(breaks)
      breaks_sauve <- breaks
    }
    
    ## palette
    if(length(palette) == 1 && palette %in% available_palette){
      if(palette %in% c("grey.colors", "gray.colors", "rainbow", "heat.colors", "terrain.colors", "topo.colors", "cm.colors")){
        palette <- eval(parse(text = paste("grDevices::", palette, "(length(breaks) - 1)", sep = "")))
      }else if(palette == "tim.colors"){
        palette <- eval(parse(text = paste("fields::", palette, "(length(breaks) - 1)", sep = "")))
      }
    }
    palette_sauve <- palette
    
    ## integrate extreme contrast in range
    
    if(max(breaks) <= max(contrast[,1], na.rm = TRUE)){
      breaks[length(breaks)] <-max(contrast[,1], na.rm = TRUE) + 10^{-12}
    }
    if(min(breaks) >= min(contrast[,1], na.rm = TRUE)){
      breaks[1] <- min(contrast[,1], na.rm = TRUE) - 10^{-12}           
    }
    
    ## check breaks - palette
    if(length(breaks) != (length(palette) + 1)){
      stop(method, " : \'breaks\' and \'palette\' are incompatible \n", 
           "length(palette) must be equal to length(breaks) + 1 \n", 
           "length(palette) : ", length(palette), " \n", 
           "length(breaks) : ", length(breaks), " \n")
    }
    
    if( (length(unique(breaks)) != length(breaks)) || sum( abs(sort(breaks) - breaks) > 10^{-16}) > 0){
      stop(method, " : the elements of \'breaks\' must be distinct and in ascending order \n", 
           "length(breaks) : ", length(breaks), "\n", 
           "length(unique(breaks)) : ", length(unique(breaks)), "\n", 
           "rank(breaks) : ", paste(rank(breaks), collapse = " "), "\n")
    }
    
    if(!is.null(pch) || any(coords %% 1 > 0)){
      col <- palette[findInterval(contrast[,1], breaks, all.inside = TRUE)]
      col[is.na(contrast[,1])] <- NA
      if(is.null(pch)){pch <- 15}
    }
    
  }else{
    palette_sauve <- unique(col[order(contrast[,1])])
    breaks_sauve <- seq(min(contrast[,1], na.rm = TRUE), max(contrast[,1], na.rm = TRUE), length.out = length(palette_sauve) - 1)
    
    #     breaks_sauve <- c(min(contrast, na.rm = TRUE) - 10^{-12}, 
    #                       seq(min(contrast[,1], na.rm = TRUE), max(contrast[,1], na.rm = TRUE), length.out = length(palette_sauve) - 1), 
    #                       max(contrast, na.rm = TRUE) + 10^{-12})
    
    if(any(is.na(contrast[,1]))){
      col[is.na(contrast[,1])] <- NA
    }
    
    if(is.null(pch)){pch <- 15}
  }
  
  #### export ####
  res <- list()
  res$contrast <- contrast
  res$palette <- palette
  res$breaks <- breaks
  res$col <- col
  res$pch <- pch
  res$palette_sauve <- palette_sauve
  res$breaks_sauve <- breaks_sauve
  res$index_duplicated <- index_duplicated
  res$index_order <- index_order
  return(res)
  
}

initDisplayWindow <- function(window, filename, path, width, height, scale, res, 
                              mfrow, bg, pty, mar, mgp, n.contrast = 1){
  
  if(window %in% c("png", "eps", "svg", "pdf")){
    switch(window, 
           "eps" = grDevices::postscript(file = paste(path, filename, ".eps", sep = ""), width = width * scale / 90, height = height * scale / 90, horizontal = FALSE, onefile = FALSE, paper = "special"), 
           "svg" = grDevices::svg(filename = paste(path, filename, ".svg", sep = ""), width = width * scale / 90, height = height * scale / 90, onefile = FALSE), 
           "png" = grDevices::png(filename = paste(path, filename, ".png", sep = ""), width = width * scale, height = height * scale, res = res), 
           "pdf" = grDevices::pdf(file = paste(path, filename, ".pdf", sep = ""), width = width * scale / 90, height = height * scale / 90, onefile = FALSE, paper = "special")
    )
  }
  
  if(window == TRUE){grDevices::dev.new()}
  
  if(!is.null(mfrow)){
    
    if(n.contrast == 1){ # cas uniparametrique
      graphics::par(mfrow = mfrow)
      
    }else{ # cas multiparametrique avec legende
      M.layout <- matrix(NA, nrow = mfrow[1], ncol = mfrow[2]+n.contrast - 1)
      M.layout[,1:mfrow[2]] <- 1:prod(mfrow)
      M.layout[1:(mfrow[1] - 1),  - (1:mfrow[2])] <-  M.layout[1:(mfrow[1] - 1), mfrow[2]]    
      M.layout[mfrow[1],  - (1:mfrow[2])] <-  seq(prod(mfrow) + 1, prod(mfrow) + n.contrast - 1)  
      widths.layout <- c(rep(1 / mfrow[2], mfrow[2] - 1), 
                         rep(1 / (mfrow[2] * n.contrast), n.contrast)
      )
      
      graphics::layout(M.layout, widths = widths.layout)
    }
  }
  
  if(!is.null(bg)){graphics::par(bg = bg)}
  if(!is.null(pty)){graphics::par(pty = pty)}
  if(!is.null(mar)){graphics::par(mar = mar)}   
  if(!is.null(mgp)){graphics::par(mgp = mgp)}     
  
}

initFilter <- function(filter, method){
  
  filter_name <- as.character(filter)
  filter_split <- as.list(strsplit(filter_name, split = "")[[1]])
  
  #### tests ####
  # 1
  if(filter_split[[1]] %in% c("2", "3") == FALSE){
    stop(method, " : wrong specification of \'filter\' \n", 
         "valid filter[1] : \"2\" \"3\" \n", 
         "proposed filter[1] :  ", filter_split[[1]], "\n")
  }
  filter_split[[1]] <- as.numeric(filter_split[[1]])
  
  # 2-3
  if(filter_split[[2]] != "D" || filter_split[[3]] != "_" ){
    stop(method, " : wrong specification of \'filter\' \n", 
         "valid filter[2:3] : \"D_\" \n", 
         "proposed filter[2:3] :  ", filter_split[[2]], " ", filter_split[[3]], "\n")
  }
  
  # 4
  if(filter_split[[4]] %in% c("G", "M", "S", "I") == FALSE){
    stop(method, " : wrong specification of \'filter\' \n", 
         "valid filter[4] : \"G\", \"M\", \"S\" or \"I\" \n", 
         "proposed filter[4] :  ", filter_split[[4]], "\n")
  }
  
  # 5-6
  if(filter_split[[4]] == "S"){  
    
    if(filter_split[[1]] == 2 && filter_split[[5]] %in% c("x", "y") == FALSE){
      stop(method, " : wrong specification of \'filter\' \n", 
           "if the fourth letter of \'filter\' is \"S\" then the fifth must be \"x\" or \"y\" \n", 
           "proposed letter: ", filter_split[[5]], "\n")
    }
    if(filter_split[[1]] == 3 && filter_split[[5]] %in% c("x", "y", "z") == FALSE){
      stop(method, " : wrong specification of \'filter\' \n", 
           "if the fourth letter of \'filter\' is \"S\" then the fifth must be \"x\", \"y\" or \"z\" \n", 
           "proposed letter: ", filter_split[[5]], "\n")
    }
    
    if(length(filter_split) != 5){
      stop(method, " : wrong specification of \'filter\' \n", 
           "if the fourth letter of \'filter\' is \"S\" then it must contains only five letters \n", 
           "proposed nomber of letter: ", length(filter_split), "\n")
    }
    
  }else{
    
    if(length(filter_split) == 5){
      if(filter_split[[5]] %in% c("3", "5", "7", "9") == FALSE){
        stop(method, " : wrong specification of \'filter\' \n", 
             "the fifth letter of \'filter\' must correspond to \"3\", \"5\", \"7\", or \"9\" \n", 
             "proposed letter: ", filter_split[[5]], "\n")
      }
      filter_split[[5]] <- as.numeric(filter_split[[5]])
    }else if(length(filter_split) == 6){
      if(filter_split[[5]] %in% c("1", "3", "5", "7", "9") == FALSE){
        stop(method, " : wrong specification of \'filter\' \n", 
             "the fifth letter of \'filter\' be odd \n", 
             "proposed letter: ", filter_split[[5]], "\n")
      }
      
      if(filter_split[[6]] %in% as.character(1:9) == FALSE){
        stop(method, " : wrong specification of \'filter\' \n",            
             "the six letter of \'filter\' must correspond to \"1\", \"2\" ... \"9\" \n", 
             "proposed letter: ", filter_split[[6]], "\n")
      }
      
      filter_split[[5]] <- as.numeric(filter_split[[6]]) + 10 * as.numeric(filter_split[[5]])
      filter_split[[6]] <- NULL
    }else{
      stop(method, " : wrong specification of \'filter\' \n", 
           "\'filter\' must contains 5 or 6 letters \n", 
           "nb of letters proposed : ", length(filter_split), "\n")
    }
    
  }
  
  
  #### creation of the filters ####  
  if(filter_split[[1]] == 2){
    if(filter_split[[4]] == "I"){ # immediate neighborhood
      filter <- matrix(0, nrow = filter_split[[5]], ncol = filter_split[[5]])
      index_1 <- rowSums((which(filter == 0, arr.ind = TRUE) - stats::median(1:filter_split[[5]]))^2) <= filter_split[[5]]
      filter[index_1] <- 1
    }
    if(filter_split[[4]] == "G"){ # gaussian
      val_Filter <- stats::dbinom(0:(filter_split[[5]] - 1), filter_split[[5]] - 1, 0.5)
      filter <- matrix(val_Filter, nrow = filter_split[[5]]) %*% matrix(val_Filter, ncol = filter_split[[5]])
    }
    if(filter_split[[4]] == "S"){ # sobel
      if(filter_split[[5]] == "x"){filter <- c(1,2,1) %*% t(c(1,0,-1))}
      if(filter_split[[5]] == "y"){filter <- c(1,0,-1) %*% t(c(1,2,1))}      
    }
    if(filter_split[[4]] == "M"){ # median
      filter <- matrix(1, nrow = filter_split[[5]], ncol = filter_split[[5]])
    }     
  }
  if(filter_split[[1]] == 3){
    if(filter_split[[4]] == "I"){
      filter <- array(0, dim = rep(filter_split[[5]], 3))     
      index_1 <- rowSums((which(filter == 0, arr.ind = TRUE) - stats::median(1:filter_split[[5]]))^2) <= filter_split[[5]]
      filter[index_1] <- 1
    }
    if(filter_split[[4]] == "G"){
      val_Filter <- stats::dbinom(0:(filter_split[[5]] - 1), filter_split[[5]] - 1, 0.5)
      filter <- array(val_Filter, dim = filter_split[[5]]) %o% array(val_Filter, dim = filter_split[[5]]) %o% array(val_Filter, dim = filter_split[[5]])
    }
    if(filter_split[[4]] == "S"){      
      if(filter_split[[5]] == "x"){filter <- c(1,2,1) %o% c(1,0,-1) %o% c(0,1,0)}
      if(filter_split[[5]] == "y"){filter <- c(1,0,-1) %o% c(1,2,1) %o% c(0,1,0)}
      if(filter_split[[5]] == "z"){filter <- c(0,1,0) %o% c(1,2,1) %o% c(1,0,-1)}
    }  
    if(filter_split[[4]] == "M"){
      filter <- array(1, dim = rep(filter_split[[5]], 3))  
    }      
  }
  
  return(list(filter = filter, 
              filter_split = filter_split)
  )
}

initIndex <- function(object, index, num, hemisphere = "both", numeric2logical, indexNum = NULL, 
                      outline.default, cex.default, pch.default, col.default, filter.default, method){
  
  #### initialization ####
  if(length(index) == 1 && is.character(index)){
    index <- list(coords = index)
  }
  if(is.matrix(index)){
    index <- as.data.frame(index)
  }
  if(is.data.frame(index)){
    index <- list(coords = index)
  }
  
  #### parameter ####
  if(is.list(index) && length(index$coords) == 1 && is.character(index$coords)){
    
    param_index <- index$coords
    
    if(class(object) == "MRIaggr"){
      initParameter(object, param = param_index, test = TRUE, init = FALSE, accept.coords = FALSE, accept.mask = TRUE, accept.index = FALSE, 
                    arg_name = "param", long_name = "parameters", method = "multiplot")        
      index$coords <- selectContrast(object, param = param_index, num = num, hemisphere = hemisphere, coords = TRUE)
    }else{        
      stop(method, " : wrong specification of \'index", indexNum, "\' \n", 
           "\'index", indexNum, "\' it must be a list containing an element named \"coords\" containing the coordinates of the points to display \n", 
           "index", indexNum, "$coords : ", index$coords, "\n")
    }
    
    if(numeric2logical == TRUE){index$coords[,param_index] <- as.logical(index$coords[,param_index])}
    
    if(is.logical(index$coords[,param_index]) == FALSE){
      stop(method, " : \'index", indexNum, "$coords\'is not of type logical \n", 
           "requested parameter : ", param_index, "\n", 
           "type of ", param_index, " values : ", paste(is(index$coords[,param_index]), collapse = " "), "\n", 
           "to force the conversion to logical set \'numeric2logical\'= TRUE \n")
    }
    
    index$index <- selectContrast(object, param = "index", num = num, hemisphere = hemisphere, format = "vector")[index$coords[,param_index] == TRUE]
    index$coords <- index$coords[index$coords[,param_index] == TRUE, c("i", "j", "k")]             
  }
  
  if(!is.list(index) || "coords" %in% names(index) == FALSE ){
    stop(method, " : wrong specification of \'index", indexNum, "\' \n", 
         "\'index", indexNum, "\' must be a list containing an element named \"coords\" containing the coordinates of the points to display \n", 
         "names(index", indexNum, ") : ", paste(names(index), collapse = " "), "\n")
  }
  
  if(class(object) == "MRIaggr"){
    if(length(index$coords) == 1 || ncol(index$coords) != 3 || any(names(index$coords) %in% c("i", "j", "k") == FALSE)){
      stop(method, " : wrong specification of \'index", indexNum, "\' \n", 
           "\'index", indexNum, "$coords\' must have 3 columns named \"i\" \"j\" \"k\" \n", 
           "names(index", indexNum, "$coords) : ", paste(names(index), collapse = " "), "\n")
    }
  }else{  
    if(length(index$coords) == 1 || ncol(index$coords) != 3){
      stop(method, " : wrong specification of \'index", indexNum, "\' \n", 
           "\'index", indexNum, "$coords\' must have 3 columns \n", 
           "names(index", indexNum, "$coords) : ", paste(names(index), collapse = " "), "\n")
    }
  }  
  
  if("cex" %in% names(index) == FALSE){index$cex <- cex.default}
  if("pch" %in% names(index) == FALSE){index$pch <- pch.default}
  if("col" %in% names(index) == FALSE){index$col <- col.default}
  
  test.IndexOutlineT <- ("outline" %in% names(index) == TRUE) && (index$outline == TRUE)
  test.IndexOutlineF <- ("outline" %in% names(index) == TRUE) && (index$outline == FALSE)

  if(test.IndexOutlineT || (outline.default == TRUE && !test.IndexOutlineF)){
    if("filter" %in% names(index) == TRUE){filter <- index$filter}else{filter <- filter.default}
    index$coords <-  pointsOutline(index$coords, filter = filter)
  }
  
  #### export ####
  return(index)
  
}

initNeighborhood <- function(Neighborhood, method){
  
  filter_name <- as.character(Neighborhood)
  valid_names <- c("2D_N4", "2D_N8", "3D_N4", "3D_N6", "3D_N8", "3D_N10", "3D_N18", "3D_N26")
  
  validCharacter(value = Neighborhood, validLength = 1, 
                 validValues = c("2D_N4", "2D_N8", "3D_N4", "3D_N6", "3D_N8", "3D_N10", "3D_N18", "3D_N26"), 
                 refuse.NULL = TRUE, method = method)
  
  Neighborhood_split <- unlist(strsplit(Neighborhood, split = ""))
  
  p.Neighborhood <- as.numeric(Neighborhood_split[[1]])
  n.Neighborhood <- if(length(Neighborhood_split) == 5){
    as.numeric(Neighborhood_split[5])
  }else{
    sum(as.numeric(Neighborhood_split[5:6]) * c(10, 1))
  }
  
  Neighborhood <- matrix(0, nrow = n.Neighborhood, ncol = p.Neighborhood)
  
  Neighborhood[1:4,1:2] <- rbind(c(-1,0), 
                                 c(0,-1), 
                                 c(1,0), 
                                 c(0,1))
  
  if(n.Neighborhood %in% c(8,10,18,26)){
    Neighborhood[5:8,1:2] <- rbind(c(-1,-1), 
                                   c(1,1), 
                                   c(-1,1), 
                                   c(1,-1)
    )
  }
  
  if(n.Neighborhood %in% c(6,10,18,26)){
    row_tempo <-  min(which(rowSums(abs(Neighborhood)) == 0))
    Neighborhood[seq(row_tempo, row_tempo + 1),] <- rbind(c(0,0,1), 
                                                          c(0,0,-1)
    )
  }
  
  if(n.Neighborhood %in% c(18,26)){
    row_tempo <-  min(which(rowSums(abs(Neighborhood)) == 0))
    Neighborhood[seq(row_tempo, row_tempo + 7),] <- rbind(c(1,0,1), 
                                                          c(0,1,1), 
                                                          c(-1,0,1), 
                                                          c(0,-1,1), 
                                                          c(1,0,-1), 
                                                          c(0,1,-1), 
                                                          c(-1,0,-1), 
                                                          c(0,-1,-1)
    )
  }
  
  if(n.Neighborhood == 26){
    row_tempo <-  min(which(rowSums(abs(Neighborhood)) == 0))
    Neighborhood[seq(row_tempo, row_tempo + 7),] <- rbind(c(1,1,1), 
                                                          c(-1,1,1), 
                                                          c(-1,-1,1), 
                                                          c(1,-1,1), 
                                                          c(1,1,-1), 
                                                          c(-1,1,-1), 
                                                          c(-1,-1,-1), 
                                                          c(1,-1,-1)
    )
  }
  
  return(Neighborhood)
}

initPackage <- function(package, argument = NULL, tryAttach = FALSE, method){
  
  test.package <- requireNamespace(package, quietly = TRUE)
  if(test.package == FALSE){
    stop(method, " : this function ", if(!is.null(argument)){paste("with argument ", argument, " ", sep = "")}, "requires to have installed the ", package, " package to work \n")
  }
  if(tryAttach == TRUE && (paste("package:", package, sep = "") %in% search() == FALSE) ){
    try(attachNamespace(package))
  }
}

initSpace <- function(character){
  
  nchar.character <- lapply(character, function(x){nchar(x)})
  nchar_max.character <- max(unlist(nchar.character))
  
  space.character <- sapply(character, function(x){do.call(paste, c(as.list(rep(" ", nchar_max.character-nchar(x))), sep = ""))})
  space.character <- unlist(lapply(space.character, function(x){if(length(x) == 0){""}else{x}}))
  return(space.character)
}

initWindow <- function(window, filename, path, width, height, unit, res, 
                       n.plot, mfrow, xlim, ylim, method){
  
  #### tests ####
  validCharacter(window, validLength = 1, validValues = c(TRUE,FALSE,"png","eps","svg","pdf"), refuse.NULL = FALSE, method = method)
  
  validCharacter(filename, validLength = 1, method = method)
  
  validNumeric(width, validLength = 1, min = 0, max = NULL, refuse.NA = TRUE, method = method)
  
  validNumeric(height, validLength = 1, min = 0, max = NULL, refuse.NA = TRUE, method = method)
  
  validCharacter(path, validLength = 1, refuse.NULL = FALSE, method = method)

  validPath(path, method = method)
  
  validCharacter(unit, validLength = 1, validValues = c("px","in","cm","mm"), refuse.NULL = TRUE, method = method)
  
  #### initialization ####
  scale <- switch(unit, 
                  "px" = 1, 
                  "in" = 90, 
                  "cm" = 35.43, 
                  "mm" = 3.543)
  
  if(is.null(mfrow)){
    if(n.plot <= 9){
      mfrow <- c(ceiling(n.plot / ceiling(sqrt(n.plot))), ceiling(sqrt(n.plot)))
    }else{
      mfrow <- c(3, 3)
    }              
  }else{
    mfrow <- mfrow
  }
  
  n.graph_par_window <- prod(mfrow)
  
  if(!is.null(xlim)){xlim.plot <- xlim}else{xlim.plot <- NULL}
  if(!is.null(ylim)){ylim.plot <- ylim}else{ylim.plot <- NULL}
  
  #### export ####
  res <- list()
  res$scale <- scale
  res$mfrow <- mfrow
  res$n.graph_par_window <- n.graph_par_window
  res$xlim.plot <- xlim.plot
  res$ylim.plot <- ylim.plot
  
  return(res)  
  
}
