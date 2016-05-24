`RODM_create_glm_model` <- function(
#                                        
# This function creates an ODM Generalized Linear Model. 
#
   database,                          # Database ODBC channel identifier
   data_table_name,                   # Database table/view containing the training dataset
   case_id_column_name = NULL,        # Row unique case identifier in data_table_name				
   target_column_name,                # Target column name in data_table_name					
   model_name = "GLM_MODEL",          # ODM Model name				  
   mining_function = "classification",# Type of GLM model: "classification" or "regression"
   auto_data_prep = TRUE,             # Setting to perform automatic data preparation
   class_weights = NULL,              # Data frame containing target class weights
   weight_column_name = NULL,         # Column in data_table_name to weight rows differently
   conf_level = NULL,
   reference_class_name = NULL,
   missing_value_treatment = NULL,    # Setting for specifying missing value treatment
   ridge_regression = NULL,
   ridge_value = NULL,
   vif_for_ridge = NULL,
   diagnostics_table_name = NULL,     # Table to hold per-row diagnostics
   retrieve_outputs_to_R = TRUE,      # Flag controlling if the output results are moved to the R environment 
   leave_model_in_dbms = TRUE,        # Flag controlling if the model is deleted or left in RDBMS               
   sql.log.file = NULL)               # File where to append the log of all the SQL calls made by this function
{
   if (!is.null(sql.log.file)) write(paste("--- SQL calls by ODM function: RODM_create_glm_model ", 
                       date(), "---"), file = sql.log.file, append = TRUE, ncolumns = 1000)

   # Validate mining function
   if (mining_function == "classification") {
       mining_function_value <- "dbms_data_mining.classification"    
   } else if (mining_function == "regression") {
       mining_function_value <- "dbms_data_mining.regression"    
   } else {
     stop("Invalid mining_function specified for RODM_create_glm_model")
     return()
   }

   # Store settings in the RDBMS RODM settings table
   GLM.settings.table <- data.frame(matrix(c(
       "ALGO_NAME", "ALGO_GENERALIZED_LINEAR_MODEL"),
       nrow = 1, ncol=2, byrow=TRUE))
   names(GLM.settings.table) <- c("SETTING_NAME", "SETTING_VALUE")
   if (!is.null(weight_column_name)) {
     GLM.settings.table <- rbind(GLM.settings.table, 
         data.frame(matrix(c("ODMS_ROW_WEIGHT_COLUMN_NAME", weight_column_name),
           nrow=1, ncol=2, byrow=TRUE,
           dimnames = list(NULL,c("SETTING_NAME", "SETTING_VALUE")))))
   }
   if (!is.null(conf_level)) {
     GLM.settings.table <- rbind(GLM.settings.table, 
         data.frame(matrix(c("GLMS_CONF_LEVEL", conf_level),
           nrow=1, ncol=2, byrow=TRUE,
           dimnames = list(NULL,c("SETTING_NAME", "SETTING_VALUE")))))
   }
   if (!is.null(reference_class_name)) {
     GLM.settings.table <- rbind(GLM.settings.table, 
         data.frame(matrix(c("GLMS_REFERENCE_CLASS_NAME", reference_class_name),
           nrow=1, ncol=2, byrow=TRUE,
           dimnames = list(NULL,c("SETTING_NAME", "SETTING_VALUE")))))
   }
   if (!is.null(missing_value_treatment)) {
     GLM.settings.table <- rbind(GLM.settings.table, 
         data.frame(matrix(c("ODMS_MISSING_VALUE_TREATMENT", missing_value_treatment),
           nrow=1, ncol=2, byrow=TRUE,
           dimnames = list(NULL,c("SETTING_NAME", "SETTING_VALUE")))))
   }
   if (!is.null(ridge_regression)) {
     GLM.settings.table <- rbind(GLM.settings.table, 
         data.frame(matrix(c("GLMS_RIDGE_REGRESSION", ridge_regression),
           nrow=1, ncol=2, byrow=TRUE,
           dimnames = list(NULL,c("SETTING_NAME", "SETTING_VALUE")))))
     if (ridge_regression == "GLMS_RIDGE_REG_ENABLE") {
       if (!is.null(ridge_value)) {
         GLM.settings.table <- rbind(GLM.settings.table, 
             data.frame(matrix(c("GLMS_RIDGE_VALUE", ridge_value),
               nrow=1, ncol=2, byrow=TRUE,
               dimnames = list(NULL,c("SETTING_NAME", "SETTING_VALUE")))))
       }
       if (!is.null(vif_for_ridge)) {
         GLM.settings.table <- rbind(GLM.settings.table, 
             data.frame(matrix(c("GLMS_VIF_FOR_RIDGE", vif_for_ridge),
               nrow=1, ncol=2, byrow=TRUE,
               dimnames = list(NULL,c("SETTING_NAME", "SETTING_VALUE")))))
       }
     }
   }
   if (!is.null(diagnostics_table_name)) {
     GLM.settings.table <- rbind(GLM.settings.table, 
         data.frame(matrix(c("GLMS_DIAGNOSTICS_TABLE_NAME", diagnostics_table_name),
           nrow=1, ncol=2, byrow=TRUE,
           dimnames = list(NULL,c("SETTING_NAME", "SETTING_VALUE")))))
   }
   RODM_store_settings(database, GLM.settings.table, auto_data_prep, 
                       sql.log.file, class_weights, TRUE)

   # Create the ODM Generalized Linear model, retrieving
   # basic details (settings and attributes) if desired
   glm.list <- RODM_create_model(
     database, model_name, mining_function_value,
     data_table_name, case_id_column_name, target_column_name, 
     retrieve_outputs_to_R, sql.log.file)

   # Retrieve GLM-specific details if desired
   if (retrieve_outputs_to_R == TRUE) { 
     query.string <- paste("SELECT ",
            "GLOBAL_DETAIL_NAME, ",
            "GLOBAL_DETAIL_VALUE ",
            "FROM table(dbms_data_mining.get_model_details_global('",
            model_name, "')) t ",
            "order by 1", sep="");
     globals <- sqlQuery(database, query = query.string)
     if (!is.null(sql.log.file)) write(query.string, file = sql.log.file, append = TRUE, ncolumns = 1000)
     glm.list <- c(glm.list, list("glm.globals" = globals))
     query.string <- paste("SELECT ",
            "CLASS, ",
            "ATTRIBUTE_NAME, ",
            "ATTRIBUTE_SUBNAME, ",
            "ATTRIBUTE_VALUE, ",
            "COEFFICIENT, ",
            "STD_ERROR, ",
            "TEST_STATISTIC, ",
            "P_VALUE, ",
            "VIF, ",
            "STD_COEFFICIENT, ",
            "LOWER_COEFF_LIMIT, ",
            "UPPER_COEFF_LIMIT, ",
            "EXP_COEFFICIENT, ",
            "EXP_LOWER_COEFF_LIMIT, ",
            "EXP_UPPER_COEFF_LIMIT ",
            "FROM table(dbms_data_mining.get_model_details_glm('",
            model_name, "')) t ",
            "order by 1,2,3,4", sep="");
     coefficients <- sqlQuery(database, query = query.string)
     if (!is.null(sql.log.file)) write(query.string, file = sql.log.file, append = TRUE, ncolumns = 1000)
     glm.list <- c(glm.list, list("glm.coefficients" = coefficients))
   }

   # Clean up as requested
   if (leave_model_in_dbms == FALSE) RODM_drop_model(database, model_name, sql.log.file)

   return(glm.list)
} # End of RODM_create_glm_model

