`RODM_create_oc_model` <- function(
#
# Create an ODM O-Cluster model
#
   database,                   # Database ODBC channel identifier
   data_table_name,            # Database table/view containing the training dataset
   case_id_column_name,        # name of the column of data_table_frame containing the case id
   model_name = "OC_MODEL",    # model name
   auto_data_prep = TRUE,      # Setting to perform automatic data preparation
   num_clusters = NULL,        # Setting that specifies the number of clusters for a clustering model.
   max_buffer = NULL,          # Buffer size for O-Cluster. Default is 50,000.
   sensitivity = NULL,         # A fraction that specifies the peak density required for separating a new cluster. 
   retrieve_outputs_to_R = TRUE, # Flag controlling if the output results are moved to the R environment (optional)
   leave_model_in_dbms = TRUE,   # Flag controlling if the model is deleted or left in RDBMS
   sql.log.file = NULL)          # File where to append the log of all the SQL calls made by this function (optional)

{
   if (!is.null(sql.log.file)) write(paste("--- SQL calls by ODM function: RODM_create_oc_model ", 
                       date(), "---"), file = sql.log.file, append = TRUE, ncolumns = 1000)

   # Store settings in the RDBMS RODM settings table
   OC.settings.table <- data.frame(matrix(c(
       "ALGO_NAME", "ALGO_O_CLUSTER"),
       nrow = 1, ncol=2, byrow=TRUE))
   names(OC.settings.table) <- c("SETTING_NAME", "SETTING_VALUE")
   if (!is.null(num_clusters)) {
       OC.settings.table <- rbind(OC.settings.table, 
           data.frame(matrix(c("CLUS_NUM_CLUSTERS", num_clusters),
             nrow=1, ncol=2, byrow=TRUE,
             dimnames = list(NULL,c("SETTING_NAME", "SETTING_VALUE")))))
   }
   if (!is.null(max_buffer)) {
       OC.settings.table <- rbind(OC.settings.table, 
           data.frame(matrix(c("OCLT_MAX_BUFFER", max_buffer),
             nrow=1, ncol=2, byrow=TRUE,
             dimnames = list(NULL,c("SETTING_NAME", "SETTING_VALUE")))))
   }
   if (!is.null(sensitivity)) {
       OC.settings.table <- rbind(OC.settings.table, 
           data.frame(matrix(c("OCLT_SENSITIVITY", sensitivity),
             nrow=1, ncol=2, byrow=TRUE,
             dimnames = list(NULL,c("SETTING_NAME", "SETTING_VALUE")))))
   }
   RODM_store_settings(database, OC.settings.table, auto_data_prep, sql.log.file)

   # Create the ODM O-Cluster model, retrieving
   # basic details (settings and attributes) if desired
   oc.list <- RODM_create_model(
     database, model_name, "dbms_data_mining.clustering",
     data_table_name, case_id_column_name, NULL, 
     retrieve_outputs_to_R, sql.log.file)

   # Retrieve OC-specific details if desired
   if (retrieve_outputs_to_R == TRUE) { 
     query.string <- paste("SELECT ",
            "ID CLU_ID, ",
            "CLUSTER_ID CLU_NAME, ",
            "RECORD_COUNT REC_CNT, ",
            "PARENT PARENT_CLU_ID, ",
            "TREE_LEVEL ",
            "FROM table(dbms_data_mining.get_model_details_oc('",
            model_name, "', null, null, 0, 0, 0)) ",
            "order by 1", sep="");
     clusters <- sqlQuery(database, query = query.string)
     if (!is.null(sql.log.file)) write(query.string, file = sql.log.file, append = TRUE, ncolumns = 1000)
     query.string <- paste("SELECT ",
            "ID CLU_ID, ",
            "ATTRIBUTE_NAME, ",
            "CONDITIONAL_OPERATOR, ",
            "ATTRIBUTE_NUM_VALUE, ",
            "ATTRIBUTE_STR_VALUE ",
            "FROM table(dbms_data_mining.get_model_details_oc('",
            model_name, "', null, null, 0, 0, 0)) t, ","
            table(t.split_predicate) ",
            "order by 1,2,3,4,5", sep="");
     split_predicate <- sqlQuery(database, query = query.string)
     if (!is.null(sql.log.file)) write(query.string, file = sql.log.file, append = TRUE, ncolumns = 1000)
     query.string <- paste("SELECT ",
            "T.ID PARENT_CLU_ID, ",
            "S.ID CHILD_CLU_ID ",
            "FROM table(dbms_data_mining.get_model_details_oc('",
            model_name, "', null, null, 0, 0, 0)) T, ","
            table(t.child) S ",
            "order by 1,2", sep="");
     taxonomy <- sqlQuery(database, query = query.string)
     if (!is.null(sql.log.file)) write(query.string, file = sql.log.file, append = TRUE, ncolumns = 1000)
     query.string <- paste("SELECT ",
            "T.ID CLU_ID, ",
            "ATTRIBUTE_NAME, ",
            "MEAN, ",
            "MODE_VALUE ",
            "FROM table(dbms_data_mining.get_model_details_oc('",
            model_name, "', null, null, 1, 0, 0)) T, ","
            table(t.centroid) S ",
            "order by 1,2,3,4", sep="");
     centroid <- sqlQuery(database, query = query.string)
     if (!is.null(sql.log.file)) write(query.string, file = sql.log.file, append = TRUE, ncolumns = 1000)
     query.string <- paste("SELECT ",
            "T.ID CLU_ID, ",
            "ATTRIBUTE_NAME, ",
            "BIN_ID, ",
            "LABEL, ",
            "COUNT ",
            "FROM table(dbms_data_mining.get_model_details_oc('",
            model_name, "', null, null, 0, 1, 0)) T, ","
            table(t.histogram) S ",
            "order by 1,2,3", sep="");
     histogram <- sqlQuery(database, query = query.string)
     if (!is.null(sql.log.file)) write(query.string, file = sql.log.file, append = TRUE, ncolumns = 1000)
     query.string <- paste("SELECT ",
            "S.ATTRIBUTE_NUM_VALUE CONSEQUENT_CLU_ID, ",
            "S.ATTRIBUTE_SUPPORT CONSEQUENT_SUPPORT, ",
            "S.ATTRIBUTE_CONFIDENCE CONSEQUENT_CONFIDENCE, ",
            "T.RULE.RULE_SUPPORT ANTECEDENT_SUPPORT, ",
            "T.RULE.RULE_CONFIDENCE ANTECEDENT_CONFIDENCE, ",
            "R.ATTRIBUTE_NAME ANT_ATTRIBUTE_NAME, ",
            "R.CONDITIONAL_OPERATOR ANT_CONDITIONAL_OPERATOR, ",
            "R.ATTRIBUTE_NUM_VALUE ANT_ATTRIBUTE_NUM_VALUE, ",
            "R.ATTRIBUTE_STR_VALUE ANT_ATTRIBUTE_STR_VALUE, ",
            "R.ATTRIBUTE_SUPPORT ANT_ATTRIBUTE_SUPPORT, ",
            "R.ATTRIBUTE_CONFIDENCE ANT_ATTRIBUTE_CONFIDENCE ",
            "FROM table(dbms_data_mining.get_model_details_oc('",
            model_name, "', null, null, 0, 0, 2)) T, ",
            "table(T.rule.consequent) S, ",
            "table(T.rule.antecedent) R ",
            "order by 1,6", sep="");
     rule <- sqlQuery(database, query = query.string)
     if (!is.null(sql.log.file)) write(query.string, file = sql.log.file, append = TRUE, ncolumns = 1000)
     query.string <- paste(
            "SELECT CLU_ID, COUNT(*) CNT FROM (",
            "SELECT ",
            "CLUSTER_ID(", model_name, " USING *) CLU_ID ",
            "FROM ", data_table_name,
            ") GROUP BY CLU_ID ORDER BY CLU_ID", sep="");
     leaf_cluster_count <- sqlQuery(database, query = query.string)
     if (!is.null(sql.log.file)) write(query.string, file = sql.log.file, append = TRUE, ncolumns = 1000)
     query.string <- paste(
            "SELECT ", case_id_column_name, ", ",
            "CLUSTER_ID(", model_name, " USING *) CLU_ID, ",
            "CLUSTER_PROBABILITY(", model_name, " USING *) CLU_PROBABILITY ",
            "FROM ", data_table_name, sep="");
     assignment <- sqlQuery(database, query = query.string)
     if (!is.null(sql.log.file)) write(query.string, file = sql.log.file, append = TRUE, ncolumns = 1000)
     oc.list <- c(oc.list, list(
                 "oc.clusters" = clusters,
                 "oc.split_predicate" = split_predicate,
                 "oc.taxonomy" = taxonomy,
                 "oc.centroid" = centroid,
                 "oc.histogram" = histogram,
                 "oc.rule" = rule,
                 "oc.leaf_cluster_count" = leaf_cluster_count,
                 "oc.assignment" = assignment
                 ))
   } 

   # Clean up as requested
   if (leave_model_in_dbms == FALSE) RODM_drop_model(database, model_name, sql.log.file)
   
   return(oc.list)
} # End of RODM_create_oc_model

