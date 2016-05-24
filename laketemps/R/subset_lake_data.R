#'@importFrom dplyr filter select %>%
#'@importFrom reshape2 acast
subset_lake_data = function(lake_name, types){
  
  check_lake(lake_name)  
  
  # -- fix for R CMD check 'no visible binding for global variable'
  siteID <- "_private"
  variable <- "_private"
  year <- "_private"
  value <- "_private"
  # -- fix for R CMD check 'no visible binding for global variable'
  IDs <- get_site_ID(lake_name) # can have 1+ sites w/ same lake name
  df <- data.frame()
  df = tryCatch({
    for (i in 1:length(IDs)){
      vals <- filter(gltc_values, tolower(variable) %in% tolower(types), siteID == IDs[i]) %>%
        select(variable, year, value)
      df <- rbind(vals, df)
    }
      
    df <- acast(df, year ~ variable)
    df <- cbind(data.frame(year = as.numeric(row.names(df))), df)
    rownames(df) <- NULL
    df
  }, error = function(e) {
    return(df)
  })
  if (nrow(df) == 0) df = data.frame()
  return(df)
}