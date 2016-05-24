##' This function avoids double counting of China.
##'
##' This function should only be used when performing aggregations.
##'
##' We decide to use the smaller subsets in the regional level because
##' weighting variable may not exist for other variables for the
##' larger subsets.
##'
##' The function only work for FAOST_CODE, if the country coding
##' system is not in FAOST_CODE then use the translateCountryCode
##' function to translate it.
##'
##' @param var The variables that require to be sanitized.
##' @param data The data frame which contains the data
##' @param year The column which correspond to the year.
##' @export


CHMT = function(var, data, year = "Year"){
    for(i in 1:length(var)){
      if(length(unique(data[data$FAOST_CODE %in% c(41, 351, 357) &
                            !is.na(data[, var[i]]), "FAOST_CODE"])) > 1) {
        cat(paste0("\nNOTE: Multiple China detected in '", var[i],
                   "' sanitization is performed\n"))
        for(j in sort(unique(data[, year]))){
          ## If the China, mainland exist then we will not use
          ## (China + Taiwan) nor the (China + Hong Kong + Macau
          ## + Taiwan).
          if(NROW(data[data$FAOST_CODE == 41 &
                         data[, year] == j, var[i]]) == 1){
            if(!is.na(data[data$FAOST_CODE == 41 &
                             data[, year] == j, var[i]])){
              data[data$FAOST_CODE %in% c(351, 357), var[i]] = NA
              ## If China mainland does not exist then we
              ## will use (China + Taiwan).
            } else if(NROW(data[data$FAOST_CODE == 357 &
                                  data[, year] == j, var[i]]) == 1){
              if(!is.na(data[data$FAOST_CODE == 357 &
                               data[, year] == j, var[i]])){
                data[data$FAOST_CODE %in% c(41, 214, 351),
                     var[i]] = NA
                ## If both fails then we will use (China +
                ## Hong Kong + Macau + Taiwan)
              } else if(NROW(data[data$FAOST_CODE == 351 &
                                    data[, year] == j, var[i]]) == 1){
                if(!is.na(data[data$FAOST_CODE == 351 &
                                 data[, year] == j, var[i]])){
                  data[data$FAOST_CODE %in%
                         c(41, 96, 128, 214, 357), var[i]] = NA
                }
              }
            }
          }
        }
      } 
    }
    data
}
