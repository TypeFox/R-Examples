#' @title Extract HAZUS damage functions for specified function type
#'
#' @param func_type Flood damage or depreciation function type. 
#' Choose one of depth (depth-damage functions), velocity (velocity-depth-damage 
#' functions), ag (damage functions for agriculture, based on duration of 
#' flooding), bridge (damage function for bridges based on the severity of the 
#' flood) or deprec (depreciation with age).
#' @param long_format Logical flag to indicate whether raw data is desired or in
#' a format suited for plotting using ggplot2. Damage function data from HAZUS
#' are typically in the wide format.
#' 
#' @return data frame, the number of rows and columns depend on the first 
#' argument of the function.
#' 
#' @author Gopi Goteti
#' 
#' @export
#' 
#' @examples
#' # depth-damage functions
#' fl_dept <- extract_hazus_functions()
#' # depth-damage functions, raw data only
#' fl_dept <- extract_hazus_functions(long_format = FALSE)
#' # velocity-depth-damage functions
#' fl_velo <- extract_hazus_functions(func_type = "velocity")
#' # agriculture damage functions
#' fl_agri <- extract_hazus_functions(func_type = "ag")
#' # bridge damage functions
#' fl_bridge <- extract_hazus_functions(func_type = "bridge")
#' # depreciation functions
#' fl_depr <- extract_hazus_functions(func_type = "deprec") 
#' # columns names of all flood damage functions
#' lapply(ls(pattern = "fl_"), FUN = function(x) colnames(get(x)))
#' # flood occupancy types and description
#' data(haz_fl_occ)
#' head(haz_fl_occ)
extract_hazus_functions <- function(func_type = "depth", long_format = TRUE) {
  
  func_opts <- c("depth", "velocity", "ag", "bridge", "deprec")
  
  if (!(func_type %in% func_opts)) {
    stop("argument func_type has to be one of - ", paste(func_opts, collapse = ","))
  }
  
  if (!is.logical(long_format)) {
    stop("argument long_format has to be either TRUE or FALSE!")
  }
  
  if (func_type == "depth") {  
    haz_fl_dept <- NULL
    data(haz_fl_dept, envir = environment())
    
    if (long_format) {
      fx <- melt(data = haz_fl_dept, 
                 id.vars = c("Occupancy", "DmgFnId", "Source", "Description", 
                             "Source_Table", "Occupy_Class", "Cover_Class", "Comment"))
      
      fx$depth <- ifelse(grepl("\\_", fx$variable), 
                         0.1, 
                         ifelse(grepl("m", fx$variable), 
                                -1.0, 
                                1.0))
      fx$variable <- gsub("ft", "", fx$variable)
      fx$variable <- gsub("m", "", fx$variable)
      fx$variable <- gsub("\\_", "", fx$variable)
      fx$depth <- as.double(fx$variable) * fx$depth
      fx$variable <- NULL
      colnames(fx)[9] <- "damage"
      
      # remove missing damage percentages
      fx <- fx[!is.na(fx$damage), ]
      
      return (fx)
    } 
    else {
      return (haz_fl_dept)
    }
  } 
  else if (func_type == "velocity") {  
    haz_fl_velo <- NULL
    data(haz_fl_velo, envir = environment())
    
    if (long_format) {
      fx <- NULL
      for (ro in 1:nrow(haz_fl_velo)) {
        if (haz_fl_velo$struc_type[ro] == "Wood") {
          fl_tmp <- data.frame(struc_type = haz_fl_velo$struc_type[ro],
                               num_story = haz_fl_velo$num_story[ro],
                               vel = c(0, haz_fl_velo$thresh_v[ro], 
                                       seq(haz_fl_velo$thresh_v[ro] + 0.05, 10, 0.1)),
                               dep = NA,
                               stringsAsFactors = FALSE)
          fl_tmp$dep <- ifelse(fl_tmp$vel <= haz_fl_velo$thresh_v[ro], 
                               haz_fl_velo$thresh_d[ro], 
                               haz_fl_velo$coef1[ro] * (fl_tmp$vel ^ haz_fl_velo$expo1[ro]))
        } 
        else {
          fl_tmp <- data.frame(struc_type = haz_fl_velo$struc_type[ro],
                               num_story = haz_fl_velo$num_story[ro],
                               vel = seq(haz_fl_velo$thresh_v[ro], 10, 0.1),                                       
                               dep = NA,
                               stringsAsFactors = FALSE)
          fl_tmp$dep <- haz_fl_velo$coef1[ro] * (fl_tmp$vel ^ haz_fl_velo$expo1[ro]) + 
            haz_fl_velo$coef2[ro] * fl_tmp$vel +  haz_fl_velo$coef3[ro]
          if (max(fl_tmp$dep) < 35.0) {
            fl_tmp2 <- data.frame(struc_type = haz_fl_velo$struc_type[ro],
                                  num_story = haz_fl_velo$num_story[ro],
                                  vel = haz_fl_velo$thresh_v[ro],                                       
                                  dep = 35,
                                  stringsAsFactors = FALSE)
            fl_tmp <- rbind(fl_tmp2, fl_tmp)
          }
        }
        fx <- rbind(fx, fl_tmp)
      }
      
      fx$struc_type <- factor(fx$struc_type)
      fx$num_story <- factor(fx$num_story)
      
      return (fx)
    } 
    else {
      return (haz_fl_velo)  
    }
  } 
  else if (func_type == "ag") {
    haz_fl_agri <- NULL
    data(haz_fl_agri, envir = environment())
    
    if (long_format) {
      fx <- melt(data = haz_fl_agri, 
                 id.vars = c("Crop", "FunctionSource", "JulianDay"))
      # rename select columns
      colnames(fx)[4:5] <- c("loss_type", "damage")
      
      return (fx)
    } 
    else {
      return (haz_fl_agri)
    }    
  } 
  else if (func_type == "bridge") {    
    haz_fl_bridge <- NULL
    data(haz_fl_bridge, envir = environment())
    
    if (long_format) {
      fx <- melt(data = haz_fl_bridge, 
                 id.vars = c("BridgeDmgFnId", "Occupancy", "Source", "Description"))
      # change RPs to integer
      fx$variable <- gsub("RP", "", fx$variable)
      fx$variable <- as.integer(fx$variable)
      # rename select columns
      colnames(fx)[5:6] <- c("ret_period", "damage")
      
      return (fx)
    } 
    else {
      return (haz_fl_bridge)
    }    
  } 
  else {    
    haz_fl_depr <- NULL
    data(haz_fl_depr, envir = environment())
    
    if (long_format) {
      meas_vars <- colnames(haz_fl_depr)[!(colnames(haz_fl_depr) %in% c("Age"))]
      fx <- melt(data = haz_fl_depr, measure.vars = meas_vars)
      # rename select columns
      colnames(fx) <- c("Age", "Occupancy", "deprec")
      
      return (fx)
    } 
    else {
      return (haz_fl_depr)
    }
  }
}
