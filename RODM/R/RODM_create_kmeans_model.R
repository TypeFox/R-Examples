`RODM_create_kmeans_model` <- function(
#
# Create an ODM Hierarchical k-Means model
#
   database,                           # Database ODBC channel identifier
   data_table_name,                    # Database table/view containing the training dataset
   case_id_column_name = NULL,         # name of the column of data_table_frame containing the case id
   model_name = "KM_MODEL",            # model name
   auto_data_prep = TRUE,              # Setting to perform automatic data preparation
   num_clusters = NULL,                # Setting that specifies the number of clusters for a clustering model.
   block_growth = NULL,                # Setting that specifies the growth factor for memory to hold cluster data for k-Means.
   conv_tolerance = NULL,              # Setting that specifies the convergence tolerance for k-Means.
   euclidean_distance = TRUE,          # Distance function set to Euclidean (if false, use Cosine)
   iterations = NULL,                  # Setting that specifies the number of iterations for k-Means.
   min_pct_attr_support = NULL,        # Setting that specifies the minimum % required for attributes in rules 
   num_bins = NULL,                    # Setting that specifies the number of histogram bins k-Means.
   variance_split = TRUE,              # Use variance for splits (if false, use size)
   retrieve_outputs_to_R = TRUE,       # Flag controlling if the output results are moved to the R environment (optional)
   leave_model_in_dbms = TRUE,         # Flag controlling if the model is deleted or left in RDBMS
   sql.log.file = NULL)                # File where to append the log of all the SQL calls made by this function (optional)

{
   if (!is.null(sql.log.file)) write(paste("--- SQL calls by ODM function: RODM_create_kmeans_model ", 
                       date(), "---"), file = sql.log.file, append = TRUE, ncolumns = 1000)

   # Store settings in the RDBMS RODM settings table
   KM.settings.table <- data.frame(matrix(c(
       "ALGO_NAME", "ALGO_KMEANS"),
       nrow = 1, ncol=2, byrow=TRUE))
   names(KM.settings.table) <- c("SETTING_NAME", "SETTING_VALUE")
   if (!is.null(num_clusters)) {
       KM.settings.table <- rbind(KM.settings.table, 
           data.frame(matrix(c("CLUS_NUM_CLUSTERS", num_clusters),
             nrow=1, ncol=2, byrow=TRUE,
             dimnames = list(NULL,c("SETTING_NAME", "SETTING_VALUE")))))
   }
   if (!is.null(block_growth)) {
       KM.settings.table <- rbind(KM.settings.table, 
           data.frame(matrix(c("KMNS_BLOCK_GROWTH", block_growth),
             nrow=1, ncol=2, byrow=TRUE,
             dimnames = list(NULL,c("SETTING_NAME", "SETTING_VALUE")))))
   }
   if (!is.null(conv_tolerance)) {
       KM.settings.table <- rbind(KM.settings.table, 
           data.frame(matrix(c("KMNS_CONV_TOLERANCE", conv_tolerance),
             nrow=1, ncol=2, byrow=TRUE,
             dimnames = list(NULL,c("SETTING_NAME", "SETTING_VALUE")))))
   }
   if (euclidean_distance == FALSE) {
       KM.settings.table <- rbind(KM.settings.table, 
           data.frame(matrix(c("KMNS_DISTANCE", "KMNS_COSINE"), 
             nrow=1, ncol=2, byrow=TRUE,
             dimnames = list(NULL,c("SETTING_NAME", "SETTING_VALUE")))))
   }
   if (!is.null(iterations)) {
       KM.settings.table <- rbind(KM.settings.table, 
           data.frame(matrix(c("KMNS_ITERATIONS", iterations),
             nrow=1, ncol=2, byrow=TRUE,
             dimnames = list(NULL,c("SETTING_NAME", "SETTING_VALUE")))))
   }
   if (!is.null(min_pct_attr_support)) {
       KM.settings.table <- rbind(KM.settings.table, 
           data.frame(matrix(c("KMNS_MIN_PCT_ATTR_SUPPORT", min_pct_attr_support),
             nrow=1, ncol=2, byrow=TRUE,
             dimnames = list(NULL,c("SETTING_NAME", "SETTING_VALUE")))))
   }
   if (!is.null(num_bins)) {
       KM.settings.table <- rbind(KM.settings.table, 
           data.frame(matrix(c("KMNS_NUM_BINS", num_bins),
             nrow=1, ncol=2, byrow=TRUE,
             dimnames = list(NULL,c("SETTING_NAME", "SETTING_VALUE")))))
   }
   if (variance_split == FALSE) {
       KM.settings.table <- rbind(KM.settings.table, 
           data.frame(matrix(c("KMNS_SPLIT_CRITERION", "KMNS_SIZE"), 
             nrow=1, ncol=2, byrow=TRUE,
             dimnames = list(NULL,c("SETTING_NAME", "SETTING_VALUE")))))
   }
   RODM_store_settings(database, KM.settings.table, auto_data_prep, sql.log.file)

   # Create the ODM Heirarchical K-Means model, retrieving
   # basic details (settings and attributes) if desired
   km.list <- RODM_create_model(
     database, model_name, "dbms_data_mining.clustering",
     data_table_name, case_id_column_name, NULL, 
     retrieve_outputs_to_R, sql.log.file)

   # Retrieve KM-specific details if desired
   if (retrieve_outputs_to_R == TRUE) { 
     query.string <- paste("SELECT ",
            "ID CLU_ID, ",
            "CLUSTER_ID CLU_NAME, ",
            "RECORD_COUNT REC_CNT, ",
            "PARENT PARENT_CLU_ID, ",
            "TREE_LEVEL, ",
            "DISPERSION ",
            "FROM table(dbms_data_mining.get_model_details_km('",
            model_name, "', null, null, 0, 0, 0)) ",
            "order by 1", sep="");
     clusters <- sqlQuery(database, query = query.string)
     if (!is.null(sql.log.file)) write(query.string, file = sql.log.file, append = TRUE, ncolumns = 1000)
     query.string <- paste("SELECT ",
            "T.ID PARENT_CLU_ID, ",
            "S.ID CHILD_CLU_ID ",
            "FROM table(dbms_data_mining.get_model_details_km('",
            model_name, "', null, null, 0, 0, 0)) T, ","
            table(t.child) S ",
            "order by 1,2", sep="");
     taxonomy <- sqlQuery(database, query = query.string)
     if (!is.null(sql.log.file)) write(query.string, file = sql.log.file, append = TRUE, ncolumns = 1000)
     query.string <- paste("SELECT ",
            "T.ID CLU_ID, ",
            "ATTRIBUTE_NAME, ",
            "ATTRIBUTE_SUBNAME, ",
            "MEAN, ",
            "MODE_VALUE, ",
            "VARIANCE ",
            "FROM table(dbms_data_mining.get_model_details_km('",
            model_name, "', null, null, 1, 0, 0)) T, ","
            table(t.centroid) S ",
            "order by 1,2,3,4,5", sep="");
     centroid <- sqlQuery(database, query = query.string)
     if (!is.null(sql.log.file)) write(query.string, file = sql.log.file, append = TRUE, ncolumns = 1000)
     query.string <- paste("SELECT ",
            "T.ID CLU_ID, ",
            "ATTRIBUTE_NAME, ",
            "ATTRIBUTE_SUBNAME, ",
            "BIN_ID, ",
            "LOWER_BOUND, ",
            "UPPER_BOUND, ",
            "LABEL, ",
            "COUNT ",
            "FROM table(dbms_data_mining.get_model_details_km('",
            model_name, "', null, null, 0, 1, 0)) T, ","
            table(t.histogram) S ",
            "order by 1,2,3,4", sep="");
     histogram <- sqlQuery(database, query = query.string)
     if (!is.null(sql.log.file)) write(query.string, file = sql.log.file, append = TRUE, ncolumns = 1000)
     query.string <- paste("SELECT ",
            "S.ATTRIBUTE_NUM_VALUE CONSEQUENT_CLU_ID, ",
            "S.ATTRIBUTE_SUPPORT CONSEQUENT_SUPPORT, ",
            "S.ATTRIBUTE_CONFIDENCE CONSEQUENT_CONFIDENCE, ",
            "T.RULE.RULE_SUPPORT ANTECEDENT_SUPPORT, ",
            "T.RULE.RULE_CONFIDENCE ANTECEDENT_CONFIDENCE, ",
            "R.ATTRIBUTE_NAME ANT_ATTRIBUTE_NAME, ",
            "R.ATTRIBUTE_SUBNAME ANT_ATTRIBUTE_SUBNAME, ",
            "R.CONDITIONAL_OPERATOR ANT_CONDITIONAL_OPERATOR, ",
            "R.ATTRIBUTE_NUM_VALUE ANT_ATTRIBUTE_NUM_VALUE, ",
            "R.ATTRIBUTE_STR_VALUE ANT_ATTRIBUTE_STR_VALUE, ",
            "R.ATTRIBUTE_SUPPORT ANT_ATTRIBUTE_SUPPORT, ",
            "R.ATTRIBUTE_CONFIDENCE ANT_ATTRIBUTE_CONFIDENCE ",
            "FROM table(dbms_data_mining.get_model_details_km('",
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
     if (!is.null(case_id_column_name)) {
       query.string <- paste(
            "SELECT ", case_id_column_name, ", ",
            "CLUSTER_ID(", model_name, " USING *) CLU_ID, ",
            "CLUSTER_PROBABILITY(", model_name, " USING *) CLU_PROBABILITY ",
            "FROM ", data_table_name, sep="");
       assignment <- sqlQuery(database, query = query.string)
       if (!is.null(sql.log.file)) write(query.string, file = sql.log.file, append = TRUE, ncolumns = 1000)
     } else {
       assignment <- NULL
     }
     km.list <- c(km.list, list(
                 "km.clusters" = clusters,
                 "km.taxonomy" = taxonomy,
                 "km.centroid" = centroid,
                 "km.histogram" = histogram,
                 "km.rule"= rule,
                 "km.leaf_cluster_count" = leaf_cluster_count,
                 "km.assignment" = assignment
                 ))
   } 

   # Clean up as requested
   if (leave_model_in_dbms == FALSE) RODM_drop_model(database, model_name, sql.log.file)
   
   return(km.list)
} # End of RODM_create_kmeans_model

