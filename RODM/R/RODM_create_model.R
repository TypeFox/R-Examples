`RODM_create_model` <- function(
#                                    
# This function creates an Oracle Data Mining model, and returns
# basic information if desired (settings and signature).
# If a model with the same name already exists, the existing model
# is first dropped.
# This is an internal auxilliary function that should not be
# called directly by the end-user.
#
   database,                          # Database ODBC channel identifier
   model_name,                        # ODM Model name
   mining_function_type,              # ODM mining function
   data_table_name,                   # Database table containing the input dataset
   case_id_column_name = "",          # Row unique case identifier in data_table_name				
   target_column_name = "",           # Target column name in data_table_name					
   retrieve_outputs_to_R = TRUE,      # Flag controlling if the output results are moved to the R environment 
   sql.log.file = NULL)               # File where to append the log of all the SQL calls made by this function
{
   # Drop the model if it already exists
   if (sqlQuery(database, 
        paste("select count(*) from user_mining_models where model_name = '",
              toupper(model_name), "'", sep="", collapse="")) > 0) {
     RODM_drop_model(database, model_name, sql.log.file)
   }

   # Create the model
   query.string <- paste("BEGIN dbms_data_mining.create_model(",
      "model_name => '", model_name, "'",
      ",mining_function => ", mining_function_type,
      ",data_table_name => '", data_table_name, "'",
      ",case_id_column_name => '", case_id_column_name, "'",
      ",target_column_name => '", target_column_name, "'",
      ",settings_table_name => 'RODM_SETTINGS_TABLE'); END;",
     sep="", collapse="")
   if (!is.null(sql.log.file)) write(query.string, file = sql.log.file, append = TRUE, ncolumns = 1000)   
   res <- sqlQuery(database, query = query.string)
   if (!identical(res, character(0))) stop(res)

   model.list <- NULL
   if (retrieve_outputs_to_R == TRUE) { 

     # Extract model settings
     query.string <- 
       paste("select setting_name, setting_value, setting_type ",
             "from user_mining_model_settings ",
             "where model_name = '",
             toupper(model_name), "' order by 1", sep="")
     if (!is.null(sql.log.file)) write(query.string, file = sql.log.file, append = TRUE, ncolumns = 1000)   
     model.model_settings <- data.frame(sqlQuery(database, query.string))

     # Extract model attributes
     query.string <- 
       paste("select attribute_name, attribute_type, data_type, data_length, ",
             "data_precision, data_scale, usage_type, target ",
             "from user_mining_model_attributes ",
             "where model_name = '",
             toupper(model_name), "' order by 1", sep="")
     model.model_attributes <- data.frame(sqlQuery(database, query.string))
     if (!is.null(sql.log.file)) write(query.string, file = sql.log.file, append = TRUE, ncolumns = 1000)   

     model.list <- c(model.list, list(
                 "model.model_settings" = model.model_settings,
                 "model.model_attributes" = model.model_attributes
                  ))
   }

   return(model.list)
} # End of RODM_create_model
