`RODM_store_settings` <- function(
#
# Store model settings in the ODM settings table
# This is an internal auxilliary function that should not be
# called directly by the end-user.
#
  database,
  rodmset, 
  auto_data_prep,
  sql.log.file = NULL,
  bias_frame = NULL,
  bias_weights = FALSE
)
{
   # Populate a table that holds bias information (e.g., weights)
   if (!is.null(bias_frame)) {

     if (ncol(bias_frame) == 2) {
       if (bias_weights == TRUE) {
         bias_set_name <- "CLAS_WEIGHTS_TABLE_NAME"
       } else {
         bias_set_name <- "CLAS_PRIORS_TABLE_NAME"
       }
       if (class(bias_frame[,1]) == "numeric") {
         bias_set_value <- "RODM_NUM_PRIORS"
         quot <- ""
       } else {
         bias_set_value <- "RODM_CAT_PRIORS"
         quot <- "'"
       }
       bias_col_list <- " (TARGET_VALUE, PRIOR_PROBABILITY)"
     } else if (ncol(bias_frame) == 3) {
       bias_set_name <- "CLAS_COST_TABLE_NAME"
       if (class(bias_frame[,1]) == "numeric") {
         bias_set_value <- "RODM_NUM_COSTS"
         quot <- ""
       } else {
         bias_set_value <- "RODM_CAT_COSTS"
         quot <- "'"
       }
       bias_col_list <- " (ACTUAL_TARGET_VALUE, PREDICTED_TARGET_VALUE, COST)"
     } else {
       stop("Invalid priors, weights, or costs data frame")
     }  

     # truncate previous bias information
     query.string <- paste("TRUNCATE TABLE ", bias_set_value, sep="")
     if (!is.null(sql.log.file)) write(query.string, file = sql.log.file, append = TRUE, ncolumns = 1000)
     sqlQuery(database, query = query.string)

     for (i in 1:length(bias_frame[,1])) {
       # Insert by hand (instead of using sqlSave) for performance reasons
       query.string <- paste(
           "INSERT INTO ", bias_set_value, bias_col_list, " VALUES (",
           quot, bias_frame[i,1], quot, ",", 
           quot, bias_frame[i,2], quot, sep="")
       if (ncol(bias_frame) == 3) {
         query.string <- paste(query.string,
           ",", quot, bias_frame[i,3], quot, sep="")
       }
       query.string <- paste(query.string, ")", sep="")
       if (!is.null(sql.log.file)) write(query.string, file = sql.log.file, append = TRUE, ncolumns = 1000)   
       sqlQuery(database, query = query.string)
     }
     D <- data.frame(matrix(c(bias_set_name, bias_set_value), nrow=1, ncol=2, byrow=TRUE,
                          dimnames = list(NULL,c("SETTING_NAME", "SETTING_VALUE"))))
     rodmset <- rbind(rodmset, D)
   }

   query.string <- "TRUNCATE TABLE RODM_SETTINGS_TABLE"
   if (!is.null(sql.log.file)) write(query.string, file = sql.log.file, append = TRUE, ncolumns = 1000)
   sqlQuery(database, query = query.string)

   if (auto_data_prep == TRUE) {
     D <- data.frame(matrix(c("PREP_AUTO", "ON"), nrow=1, ncol=2, byrow=TRUE,
                            dimnames = list(NULL,c("SETTING_NAME", "SETTING_VALUE"))))
     rodmset <- rbind(rodmset, D)
   }

   # Insert by hand (instead of using sqlSave) for performance reasons
   for (i in 1:length(rodmset[,1])) {
     query.string <- paste(
         "INSERT INTO RODM_SETTINGS_TABLE (SETTING_NAME, SETTING_VALUE) ",
         "VALUES ('", rodmset[i,1],
         "','", rodmset[i,2], "')", sep="")
     if (!is.null(sql.log.file)) write(query.string, file = sql.log.file, append = TRUE, ncolumns = 1000)   
     sqlQuery(database, query = query.string)
   }

   return()
} # end of RODM_store_settings
