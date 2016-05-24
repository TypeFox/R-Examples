`RODM_create_nb_model` <- function(
#                                    
# This function creates an ODM Naive Bayes model. 
#
   database,                          # Database ODBC channel identifier
   data_table_name,                   # Database table/view containing the training dataset
   case_id_column_name = NULL,        # Row unique case identifier in data_table_name				
   target_column_name,                # Target column name in data_table_name					
   model_name = "NB_MODEL",           # ODM Model name
   auto_data_prep = TRUE,             # Setting to perform automatic data preparation
   class_priors = NULL,               # Data frame containing target class priors
   retrieve_outputs_to_R = TRUE,      # Flag controlling if the output results are moved to the R environment 
   leave_model_in_dbms = TRUE,        # Flag controlling if the model is deleted or left in RDBMS               
   sql.log.file = NULL)               # File where to append the log of all the SQL calls made by this function
{
   if (!is.null(sql.log.file)) write(paste("--- SQL calls by ODM function: RODM_create_nb_model ", 
              date(), "---"), file = sql.log.file, append = TRUE, ncolumns = 1000)

   # Store settings in the RDBMS RODM settings table
   NB.settings.table <- data.frame(matrix(c(
       "ALGO_NAME", "ALGO_NAIVE_BAYES"),
       nrow = 1, ncol=2, byrow=TRUE))
   names(NB.settings.table) <- c("SETTING_NAME", "SETTING_VALUE")
   RODM_store_settings(database, NB.settings.table, auto_data_prep, 
                       sql.log.file, class_priors)

   # Create the ODM Naive Bayes classification model, retrieving
   # basic details (settings and attributes) if desired
   nb.list <- RODM_create_model(
     database, model_name, "dbms_data_mining.classification",
     data_table_name, case_id_column_name, target_column_name, 
     retrieve_outputs_to_R, sql.log.file)

   # Retrieve NB-specific details if desired
   if (retrieve_outputs_to_R == TRUE) { 
      query.string <- paste("SELECT ",
              "TARGET_ATTRIBUTE_NAME, ",
              "TARGET_ATTRIBUTE_STR_VALUE, ",
              "TARGET_ATTRIBUTE_NUM_VALUE, ",
              "PRIOR_PROBABILITY, ",
              "ATTRIBUTE_NAME, ",
              "ATTRIBUTE_SUBNAME, ",
              "ATTRIBUTE_STR_VALUE, ",
              "ATTRIBUTE_NUM_VALUE, ",
              "CONDITIONAL_PROBABILITY ",
              "FROM table(dbms_data_mining.get_model_details_nb('",
              model_name, "')) t, table(t.conditionals) s ",
              "order by 1,2,3,5,6,7,8", sep="");
      conditionals <- sqlQuery(database, query = query.string)
      if (!is.null(sql.log.file)) write(query.string, file = sql.log.file, append = TRUE, ncolumns = 1000)
      nb.list <- c(nb.list, list("nb.conditionals" = conditionals))
   } 

   # Clean up as requested
   if (leave_model_in_dbms == FALSE) RODM_drop_model(database, model_name, sql.log.file)
   
   return(nb.list)
} # End of RODM_create_nb_model
