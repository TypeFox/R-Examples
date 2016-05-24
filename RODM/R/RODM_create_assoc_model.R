RODM_create_assoc_model <- function(
   database,                    # Database ODBC channel identifier
   data_table_name,             # Database table/view containing the training dataset
   case_id_column_name,         # Row unique case identifier in data_table_name		
   model_name = "AR_MODEL",     # Model name
   min_support = NULL,          # Min. support
   min_confidence = NULL,       # Min. confidence
   max_rule_length = NULL,      # Max. rule length
   retrieve_outputs_to_R = TRUE,# Flag controlling if the output results are moved to the R environment (optional)
   leave_model_in_dbms = TRUE,  # Flag controlling if the model is deleted or left in RDBMS
   sql.log.file = NULL)         # File where to append the log of all the SQL calls made by this function (optional)
{
   if (!is.null(sql.log.file)) write(paste("--- SQL calls by ODM function: RODM_create_assoc_model ", 
                       date(), "---"), file = sql.log.file, append = TRUE, ncolumns = 1000)

   # Store settings in the RDBMS RODM settings table
   AR.settings.table <- data.frame(matrix(c(
       "ALGO_NAME", "ALGO_APRIORI_ASSOCIATION_RULES"),
       nrow = 1, ncol=2, byrow=TRUE))
   names(AR.settings.table) <- c("SETTING_NAME", "SETTING_VALUE")
   if (!is.null(min_support)) {
       AR.settings.table <- rbind(AR.settings.table, 
           data.frame(matrix(c("ASSO_MIN_SUPPORT", min_support),
             nrow=1, ncol=2, byrow=TRUE,
             dimnames = list(NULL,c("SETTING_NAME", "SETTING_VALUE")))))
   }
   if (!is.null(min_confidence)) {
       AR.settings.table <- rbind(AR.settings.table, 
           data.frame(matrix(c("ASSO_MIN_CONFIDENCE", min_confidence),
             nrow=1, ncol=2, byrow=TRUE,
             dimnames = list(NULL,c("SETTING_NAME", "SETTING_VALUE")))))
   }
   if (!is.null(max_rule_length)) {
       AR.settings.table <- rbind(AR.settings.table, 
           data.frame(matrix(c("ASSO_MAX_RULE_LENGTH", max_rule_length),
             nrow=1, ncol=2, byrow=TRUE,
             dimnames = list(NULL,c("SETTING_NAME", "SETTING_VALUE")))))
   }
   RODM_store_settings(database, AR.settings.table, FALSE, sql.log.file)

   # Create the ODM Association Rules model, retrieving
   # basic details (settings and attributes) if desired
   ar.list <- RODM_create_model(
     database, model_name, "dbms_data_mining.association",
     data_table_name, case_id_column_name, NULL, 
     retrieve_outputs_to_R, sql.log.file)

   # Retrieve AR-specific details if desired
   if (retrieve_outputs_to_R == TRUE) { 
    query.string <- paste("SELECT ",
            "T.RULE_ID, ",
            "T.RULE_SUPPORT, ",
            "T.RULE_CONFIDENCE, ",
            "T.RULE_LIFT, ",
            "T.NUMBER_OF_ITEMS, ",
            "A.ATTRIBUTE_NAME ANT_NAME, ",
            "A.ATTRIBUTE_SUBNAME ANT_SUBNAME, ",
            "A.ATTRIBUTE_NUM_VALUE ANT_NUM_VALUE, ",
            "A.ATTRIBUTE_STR_VALUE ANT_STR_VALUE, ",
            "T.ANTECEDENT_SUPPORT ANT_SUPPORT, ",
            "C.ATTRIBUTE_NAME CONS_NAME, ",
            "C.ATTRIBUTE_SUBNAME CONS_SUBNAME, ",
            "C.ATTRIBUTE_NUM_VALUE CONS_NUM_VALUE, ",
            "C.ATTRIBUTE_STR_VALUE CONS_STR_VALUE, ",
            "T.CONSEQUENT_SUPPORT CONS_SUPPORT ",
            "FROM table(dbms_data_mining.get_association_rules('",
            model_name, "')) T, ",
            "table(T.consequent) C, ",
            "table(T.antecedent) A ",
            "order by 3 desc, 2 desc,1,8,9,10,11", sep="");
     rules <- sqlQuery(database, query = query.string)
     if (!is.null(sql.log.file)) write(query.string, file = sql.log.file, append = TRUE, ncolumns = 1000)
     query.string <- paste("SELECT ",
            "T.ITEMSET_ID, ",
            "T.SUPPORT, ",
            "T.NUMBER_OF_ITEMS, ",
            "A.ATTRIBUTE_NAME, ",
            "A.ATTRIBUTE_SUBNAME, ",
            "A.ATTRIBUTE_NUM_VALUE, ",
            "A.ATTRIBUTE_STR_VALUE ",
            "FROM table(dbms_data_mining.get_frequent_itemsets('",
            model_name, "')) T, ",
            "table(T.items) A ",
            "order by 3 desc, 2 desc,1,4,5,6,7", sep="");
     itemsets <- sqlQuery(database, query = query.string)
     if (!is.null(sql.log.file)) write(query.string, file = sql.log.file, append = TRUE, ncolumns = 1000)
     ar.list <- c(ar.list, list(
                 "ar.rules" = rules,
                 "ar.itemsets" = itemsets
                 ))
   } 

   # Clean up as requested
   if (leave_model_in_dbms == FALSE) RODM_drop_model(database, model_name, sql.log.file)
   
   return(ar.list)
} # End of RODM_create_assoc_model
