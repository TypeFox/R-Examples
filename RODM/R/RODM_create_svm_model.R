`RODM_create_svm_model` <- function(
#                      
# This function creates an ODM SVM model. 
#
   database,                          # Database ODBC channel identifier
   data_table_name,                   # Database table/view containing the training dataset
   case_id_column_name = NULL,        # Row unique case identifier in data_table_name				
   target_column_name = NULL,         # Target column name in data_table_name					
   model_name = "SVM_MODEL",          # ODM Model name				  
   mining_function = "classification",# Type of SVM model: "classification", "regression" or "anomaly_detection"
   auto_data_prep = TRUE,             # Setting to perform automatic data preparation
   class_weights = NULL,              # Data frame containing target class weights
   active_learning = TRUE,            # Setting for enabling active learning
   complexity_factor = NULL,          # Setting that specifies the complexity factor for SVM.
   conv_tolerance = NULL,             # Setting that specifies tolerance for SVM.
   epsilon = NULL,                    # Setting that specifies epsilon for SVM Regression.
   kernel_cache_size = NULL,          # Setting that specifies the Gaussian kernel cache size (bytes) for SVM.
   kernel_function = NULL,            # Setting for speifying the kernel function for SVM (SVMS_GAUSSIAN or SVMS_LINEAR)
   outlier_rate = NULL,               # Setting specifying the desired rate of outliers in the training data (one-class SVM)
   std_dev = NULL,                    # Setting that specifies standard deviation for SVM Gaussian kernel.
   retrieve_outputs_to_R = TRUE,      # Flag controlling if the output results are moved to the R environment 
   leave_model_in_dbms = TRUE,        # Flag controlling if the model is deleted or left in RDBMS               
   sql.log.file = NULL)               # File where to append the log of all the SQL calls made by this function
{
   if (!is.null(sql.log.file)) write(paste("--- SQL calls by ODM function: RODM_create_svm_model ", 
                       date(), "---"), file = sql.log.file, append = TRUE, ncolumns = 1000)

   # Validate mining function
   if (mining_function == "classification") {
       mining_function_value <- "dbms_data_mining.classification"    
   } else if (mining_function == "regression") {
       mining_function_value <- "dbms_data_mining.regression"    
   } else if (mining_function == "anomaly_detection") {
       if (!is.null(target_column_name)) {
         stop("Anomaly Detection requires NULL target column name")
         return()
       }
       mining_function_value <- "dbms_data_mining.classification"  
   } else {
     stop("Invalid mining_function specified for RODM_create_svm_model")
     return()
   }

   # Store settings in the RDBMS RODM settings table
   SVM.settings.table <- data.frame(matrix(c(
       "ALGO_NAME", "ALGO_SUPPORT_VECTOR_MACHINES"),
       nrow = 1, ncol=2, byrow=TRUE))
   names(SVM.settings.table) <- c("SETTING_NAME", "SETTING_VALUE")
   if (active_learning == FALSE) {
       SVM.settings.table <- rbind(SVM.settings.table, 
           data.frame(matrix(c("SVMS_ACTIVE_LEARNING", "SVMS_AL_DISABLE"), 
             nrow=1, ncol=2, byrow=TRUE,
             dimnames = list(NULL,c("SETTING_NAME", "SETTING_VALUE")))))
   }
   if (!is.null(complexity_factor)) {
       SVM.settings.table <- rbind(SVM.settings.table, 
           data.frame(matrix(c("SVMS_COMPLEXITY_FACTOR", complexity_factor),
             nrow=1, ncol=2, byrow=TRUE,
             dimnames = list(NULL,c("SETTING_NAME", "SETTING_VALUE")))))
   }
   if (!is.null(conv_tolerance)) {
       SVM.settings.table <- rbind(SVM.settings.table, 
           data.frame(matrix(c("SVMS_CONV_TOLERANCE", conv_tolerance),
             nrow=1, ncol=2, byrow=TRUE,
             dimnames = list(NULL,c("SETTING_NAME", "SETTING_VALUE")))))
   }
   if (!is.null(epsilon)) {
       SVM.settings.table <- rbind(SVM.settings.table, 
           data.frame(matrix(c("SVMS_EPSILON", epsilon),
             nrow=1, ncol=2, byrow=TRUE,
             dimnames = list(NULL,c("SETTING_NAME", "SETTING_VALUE")))))
   }
   if (!is.null(kernel_function)) {
     if (kernel_function == "SVMS_GAUSSIAN") {
       if (!is.null(kernel_cache_size)) {
           SVM.settings.table <- rbind(SVM.settings.table, 
               data.frame(matrix(c("SVMS_KERNEL_CACHE_SIZE", kernel_cache_size),
                 nrow=1, ncol=2, byrow=TRUE,
                 dimnames = list(NULL,c("SETTING_NAME", "SETTING_VALUE")))))
       }
       if (!is.null(std_dev)) {
           SVM.settings.table <- rbind(SVM.settings.table, 
               data.frame(matrix(c("SVMS_STD_DEV", std_dev),
                 nrow=1, ncol=2, byrow=TRUE,
                 dimnames = list(NULL,c("SETTING_NAME", "SETTING_VALUE")))))
       }
     }
     SVM.settings.table <- rbind(SVM.settings.table, 
         data.frame(matrix(c("SVMS_KERNEL_FUNCTION", kernel_function), 
           nrow=1, ncol=2, byrow=TRUE,
           dimnames = list(NULL,c("SETTING_NAME", "SETTING_VALUE")))))
   }
   if (mining_function == "anomaly_detection") {
     if (!is.null(outlier_rate)) {
         SVM.settings.table <- rbind(SVM.settings.table, 
             data.frame(matrix(c("SVMS_OUTLIER_RATE", outlier_rate),
               nrow=1, ncol=2, byrow=TRUE,
               dimnames = list(NULL,c("SETTING_NAME", "SETTING_VALUE")))))
     }
   }
   RODM_store_settings(database, SVM.settings.table, auto_data_prep, 
                       sql.log.file, class_weights)

   # Create the ODM Support Vector Machine model, retrieving
   # basic details (settings and attributes) if desired
   svm.list <- RODM_create_model(
     database, model_name, mining_function_value,
     data_table_name, case_id_column_name, target_column_name, 
     retrieve_outputs_to_R, sql.log.file)

   # Retrieve SVM-specific details if desired
   if (retrieve_outputs_to_R == TRUE) { 
     if (1 == length(svm.list$model.model_settings[svm.list$model.model_settings$SETTING_VALUE == "SVMS_LINEAR",]$SETTING_VALUE)) {
       query.string <- paste("SELECT ",
              "CLASS, ",
              "ATTRIBUTE_NAME, ",
              "ATTRIBUTE_SUBNAME, ",
              "ATTRIBUTE_VALUE, ",
              "COEFFICIENT ",
              "FROM table(dbms_data_mining.get_model_details_svm('",
              model_name, "')) t, table(t.attribute_set) s ",
              "order by 1,2,3,4", sep="");
       coefficients <- sqlQuery(database, query = query.string)
       if (!is.null(sql.log.file)) write(query.string, file = sql.log.file, append = TRUE, ncolumns = 1000)
       svm.list <- c(svm.list, list("svm.coefficients" = coefficients))
     }
   } 

   # Clean up as requested
   if (leave_model_in_dbms == FALSE) RODM_drop_model(database, model_name, sql.log.file)
   
   return(svm.list)
} # End of RODM_create_svm_model
