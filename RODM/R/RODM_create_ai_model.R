`RODM_create_ai_model` <- function(
#
# This function creates an ODM Attribute Importance model. 
#
   database,                          # Database ODBC channel identifier
   data_table_name,                   # Database table/view containing the training dataset
   case_id_column_name = NULL,        # Row unique case identifier in data_table_name				
   target_column_name,                # Target column name in data_table_name					
   model_name = "AI_MODEL",           # ODM Model name				  
   auto_data_prep = TRUE,             # Setting to perform automatic data preparation
   retrieve_outputs_to_R = TRUE,      # Flag controlling if the output results are moved to the R environment 
   leave_model_in_dbms = TRUE,        # Flag controlling if the model is deleted or left in RDBMS               
   sql.log.file = NULL)               # File where to append the log of all the SQL calls made by this function
{
   if (!is.null(sql.log.file)) write(paste("--- SQL calls by ODM function: RODM_create_ai_model ", 
              date(), "---"), file = sql.log.file, append = TRUE, ncolumns = 1000)

   # Store settings in the RDBMS RODM settings table
   AI.settings.table <- data.frame(matrix(c(
       "ALGO_NAME", "ALGO_AI_MDL"),
       nrow = 1, ncol=2, byrow=TRUE))
   names(AI.settings.table) <- c("SETTING_NAME", "SETTING_VALUE")
   RODM_store_settings(database, AI.settings.table, auto_data_prep, 
                       sql.log.file)

   # Create the ODM Attribute Importance model, retrieving
   # basic details (settings and attributes) if desired
   ai.list <- RODM_create_model(
     database, model_name, "dbms_data_mining.attribute_importance",
     data_table_name, case_id_column_name, target_column_name, 
     retrieve_outputs_to_R, sql.log.file)

   # Retrieve AI-specific details if desired
   if (retrieve_outputs_to_R == TRUE) { 
      query.string <- paste("SELECT ", 
              "attribute_name, ",
              "attribute_subname, ",
              "importance_value, ",
              "rank ", 
              "FROM TABLE(DBMS_DATA_MINING.GET_MODEL_DETAILS_AI('",
                           model_name, "')) ORDER BY RANK;", sep="")
      importance <- sqlQuery(database, query = query.string)
      if (!is.null(sql.log.file)) write(query.string, file = sql.log.file, append = TRUE, ncolumns = 1000)
      ai.list <- c(ai.list, list("ai.importance" = importance))
   } 

   # Clean up as requested
   if (leave_model_in_dbms == FALSE) RODM_drop_model(database, model_name, sql.log.file)
   
   return(ai.list)

} # End of RODM_create_ai_model
