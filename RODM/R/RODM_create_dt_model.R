`RODM_create_dt_model` <- function(
#
# This function creates an ODM Decision Tree model. 
#
   database,                     # Database ODBC channel identifier
   data_table_name,              # Database table/view containing the training dataset
   case_id_column_name = NULL,   # Row unique case identifier in data_table_name		
   target_column_name,           # Target column name in data_table_name		
   model_name = "DT_MODEL",      # ODM Model name				  
   auto_data_prep = TRUE,        # Setting to perform automatic data preparation
   cost_matrix = NULL,           # Data frame containing target class cost matrix
   gini_impurity_metric = TRUE,  # Use gini impurity metric (false for entropy)
   max_depth = NULL,             # Tree max depth
   minrec_split = NULL,          # Minimum records in a node to split the node
   minpct_split = NULL,          # Minimum % of records in a node to split the node
   minrec_node = NULL,           # Minimum records to form a node
   minpct_node = NULL,           # Minimum % of records to form a node
   retrieve_outputs_to_R = TRUE, # Flag controlling if the outpout results are moved to the R environment 
   leave_model_in_dbms = TRUE,   # Flag controlling if the model is deleted or left in RDBMS               
   sql.log.file = NULL)          # File where to append the log of all the SQL calls made by this function
{
   if (!is.null(sql.log.file)) write(paste("--- SQL calls by ODM function: RODM_create_dt_model ", 
              date(), "---"), file = sql.log.file, append = TRUE, ncolumns = 1000)

   # Store settings in the RDBMS RODM settings table
   DT.settings.table <- data.frame(matrix(c(
       "ALGO_NAME", "ALGO_DECISION_TREE"),
       nrow = 1, ncol=2, byrow=TRUE))
   names(DT.settings.table) <- c("SETTING_NAME", "SETTING_VALUE")
   if (gini_impurity_metric == FALSE) {
       DT.settings.table <- rbind(DT.settings.table, 
           data.frame(matrix(c("TREE_IMPURITY_METRIC", "TREE_IMPURITY_ENTROPY"), 
             nrow=1, ncol=2, byrow=TRUE,
             dimnames = list(NULL,c("SETTING_NAME", "SETTING_VALUE")))))
   }
   if (!is.null(max_depth)) {
       DT.settings.table <- rbind(DT.settings.table, 
           data.frame(matrix(c("TREE_TERM_MAX_DEPTH", max_depth),
             nrow=1, ncol=2, byrow=TRUE,
             dimnames = list(NULL,c("SETTING_NAME", "SETTING_VALUE")))))
   }
   if (!is.null(minrec_split)) {
       DT.settings.table <- rbind(DT.settings.table, 
           data.frame(matrix(c("TREE_TERM_MINREC_SPLIT", minrec_split),
             nrow=1, ncol=2, byrow=TRUE,
             dimnames = list(NULL,c("SETTING_NAME", "SETTING_VALUE")))))
   }
   if (!is.null(minpct_split)) {
       DT.settings.table <- rbind(DT.settings.table, 
           data.frame(matrix(c("TREE_TERM_MINPCT_SPLIT", minpct_split),
             nrow=1, ncol=2, byrow=TRUE,
             dimnames = list(NULL,c("SETTING_NAME", "SETTING_VALUE")))))
   }
   if (!is.null(minrec_node)) {
       DT.settings.table <- rbind(DT.settings.table, 
           data.frame(matrix(c("TREE_TERM_MINREC_NODE", minrec_node),
             nrow=1, ncol=2, byrow=TRUE,
             dimnames = list(NULL,c("SETTING_NAME", "SETTING_VALUE")))))
   }
   if (!is.null(minpct_node)) {
       DT.settings.table <- rbind(DT.settings.table, 
           data.frame(matrix(c("TREE_TERM_MINPCT_NODE", minpct_node),
             nrow=1, ncol=2, byrow=TRUE,
             dimnames = list(NULL,c("SETTING_NAME", "SETTING_VALUE")))))
   }
   RODM_store_settings(database, DT.settings.table, auto_data_prep, 
                       sql.log.file, cost_matrix)

   # Create the ODM Decision Tree classification model, retrieving
   # basic details (settings and attributes) if desired
   dt.list <- RODM_create_model(
     database, model_name, "dbms_data_mining.classification",
     data_table_name, case_id_column_name, target_column_name, 
     retrieve_outputs_to_R, sql.log.file)

   # Retrieve DT-specific details if desired
   if (retrieve_outputs_to_R == TRUE) { 

      query.string <- paste(
        "SELECT * FROM ",
        "  XMLTable('for $s in /PMML/TreeModel//ScoreDistribution ",
        "       return ",
        "         <scores id=\"{$s/../@id}\" ",
        "                 tvalue=\"{$s/@value}\" ",
        "                 tcount=\"{$s/@recordCount}\" ",
        "         />' ",
        "  passing dbms_data_mining.get_model_details_xml('", model_name, "') ",
        "     COLUMNS  ",
        "       node_id      NUMBER PATH '/scores/@id', ",
        "       target_value VARCHAR2(4000) PATH '/scores/@tvalue', ",
        "       target_count NUMBER PATH '/scores/@tcount') ", sep="");
      distributions <- sqlQuery(database, query = query.string)
      if (!is.null(sql.log.file)) write(query.string, file = sql.log.file, append = TRUE, ncolumns = 1000)
      dt.list <- c(dt.list, list("dt.distributions" = distributions))

      query.string <- paste(
       "WITH X as ",
       "(SELECT * FROM  ",
       "  XMLTable('for $n in /PMML/TreeModel//Node ",
       "      let $rf :=  ",
       "        if (count($n/CompoundPredicate) > 0) then ",
       "          $n/CompoundPredicate/*[1]/@field ",
       "        else ",
       "          if (count($n/SimplePredicate) > 0) then ",
       "            $n/SimplePredicate/@field ",
       "          else ",
       "            $n/SimpleSetPredicate/@field ",
       "      let $ro :=  ",
       "        if (count($n/CompoundPredicate) > 0) then ",
       "          if ($n/CompoundPredicate/*[1] instance of  ",
       "              element(SimplePredicate)) then ",
       "            $n/CompoundPredicate/*[1]/@operator ",
       "          else if ($n/CompoundPredicate/*[1] instance of  ",
       "              element(SimpleSetPredicate)) then ",
       "            (\"in\") ",
       "          else () ",
       "        else ",
       "          if (count($n/SimplePredicate) > 0) then ",
       "            $n/SimplePredicate/@operator ",
       "          else if (count($n/SimpleSetPredicate) > 0) then ",
       "            (\"in\") ",
       "          else () ",
       "      let $rv :=  ",
       "        if (count($n/CompoundPredicate) > 0) then ",
       "          if ($n/CompoundPredicate/*[1] instance of  ",
       "              element(SimplePredicate)) then ",
       "            $n/CompoundPredicate/*[1]/@value ",
       "          else  ",
       "            $n/CompoundPredicate/*[1]/Array/text() ",
       "        else ",
       "          if (count($n/SimplePredicate) > 0) then ",
       "            $n/SimplePredicate/@value ",
       "          else ",
       "            $n/SimpleSetPredicate/Array/text() ",
       "      let $sf :=  ",
       "        if (count($n/CompoundPredicate) > 0) then ",
       "          $n/CompoundPredicate/*[2]/@field ",
       "        else () ",
       "      let $so :=  ",
       "        if (count($n/CompoundPredicate) > 0) then ",
       "          if ($n/CompoundPredicate/*[2] instance of  ",
       "              element(SimplePredicate)) then ",
       "            $n/CompoundPredicate/*[2]/@operator ",
       "          else if ($n/CompoundPredicate/*[2] instance of  ",
       "              element(SimpleSetPredicate)) then ",
       "            (\"in\") ",
       "          else () ",
       "        else () ",
       "      let $sv :=  ",
       "        if (count($n/CompoundPredicate) > 0) then ",
       "          if ($n/CompoundPredicate/*[2] instance of  ",
       "              element(SimplePredicate)) then ",
       "            $n/CompoundPredicate/*[2]/@value ",
       "          else ",
       "            $n/CompoundPredicate/*[2]/Array/text() ",
       "        else () ",
       "      return ",
       "        <pred id=\"{$n/../@id}\" ",
       "              score=\"{$n/@score}\" ",
       "              rec=\"{$n/@recordCount}\" ",
       "              cid=\"{$n/@id}\" ",
       "              rf=\"{$rf}\" ",
       "              ro=\"{$ro}\" ",
       "              rv=\"{$rv}\" ",
       "              sf=\"{$sf}\" ",
       "              so=\"{$so}\" ",
       "              sv=\"{$sv}\" ",
       "        />' ",
       "passing dbms_data_mining.get_model_details_xml('", model_name, "') ",
       "      COLUMNS  ",
       "        parent_node_id   NUMBER PATH '/pred/@id', ",
       "        child_node_id    NUMBER PATH '/pred/@cid', ",
       "        rec              NUMBER PATH '/pred/@rec', ",
       "        score            VARCHAR2(4000) PATH '/pred/@score', ",
       "        rule_field       VARCHAR2(4000) PATH '/pred/@rf', ",
       "        rule_op          VARCHAR2(20) PATH '/pred/@ro', ",
       "        rule_value       VARCHAR2(4000) PATH '/pred/@rv', ",
       "        surr_field       VARCHAR2(4000) PATH '/pred/@sf', ",
       "        surr_op          VARCHAR2(20) PATH '/pred/@so', ",
       "        surr_value       VARCHAR2(4000) PATH '/pred/@sv')) ",
       "select pid parent_node, nid node, rec record_count,  ",
       "       score prediction, rule_pred local_rule, surr_pred local_surrogate,  ",
       "       rtrim(replace(full_rule,'$O$D$M$'),' AND') full_simple_rule from ( ",
       " select row_number() over (partition by nid order by rn desc) rn, ",
       "  pid, nid, rec, score, rule_pred, surr_pred, full_rule from ( ",
       "  select rn, pid, nid, rec, score, rule_pred, surr_pred,  ",
       "    sys_connect_by_path(pred, '$O$D$M$') full_rule from ( ",
       "   select row_number() over (partition by nid order by rid) rn, ",
       "     pid, nid, rec, score, rule_pred, surr_pred,  ",
       "     nvl2(pred,pred || ' AND ',null) pred from( ",
       "    select rid, pid, nid, rec, score, rule_pred, surr_pred, ",
       "      decode(rn, 1, pred, null) pred from ( ",
       "     select rid, nid, rec, score, pid, rule_pred, surr_pred, ",
       "      nvl2(root_op, '(' || root_field || ' ' || root_op || ' ' || root_value || ')', null) pred, ",
       "      row_number() over (partition by nid, root_field, root_op order by rid desc) rn from ( ",
       "      SELECT  ",
       "        connect_by_root(parent_node_id) rid, ",
       "        child_node_id nid, ",
       "        rec, score, ",
       "        connect_by_root(rule_field) root_field, ",
       "        connect_by_root(rule_op) root_op, ",
       "        connect_by_root(rule_value) root_value, ",
       "        nvl2(rule_op, '(' || rule_field || ' ' || rule_op || ' ' || rule_value || ')',  null) rule_pred, ",
       "        nvl2(surr_op, '(' || surr_field || ' ' || surr_op || ' ' || surr_value || ')',  null) surr_pred, ",
       "        parent_node_id pid ",
       "        FROM ( ",
       "         SELECT parent_node_id, child_node_id, rec, score, rule_field, surr_field, rule_op, surr_op,  ",
       "                replace(replace(rule_value,'&quot; &quot;', ''', '''),'&quot;', '''') rule_value, ",
       "                replace(replace(surr_value,'&quot; &quot;', ''', '''),'&quot;', '''') surr_value ",
       "         FROM ( ",
       "           SELECT parent_node_id, child_node_id, rec, score, rule_field, surr_field, ",
       "                  decode(rule_op,'lessOrEqual','<=','greaterThan','>',rule_op) rule_op, ",
       "                  decode(rule_op,'in','('||rule_value||')',rule_value) rule_value, ",
       "                  decode(surr_op,'lessOrEqual','<=','greaterThan','>',surr_op) surr_op, ",
       "                  decode(surr_op,'in','('||surr_value||')',surr_value) surr_value ",
       "           FROM X) ",
       "        ) ",
       "        CONNECT BY PRIOR child_node_id = parent_node_id ",
       "      ) ",
       "     ) ",
       "    ) ",
       "   ) ",
       "   CONNECT BY PRIOR rn = rn - 1 ",
       "          AND PRIOR nid = nid ",
       "   START WITH rn = 1 ",
       " ) ",
       ") ",
       "where rn = 1", sep="");
      nodes <- sqlQuery(database, query = query.string)
      if (!is.null(sql.log.file)) write(query.string, file = sql.log.file, append = TRUE, ncolumns = 1000)
      dt.list <- c(dt.list, list("dt.nodes" = nodes))

   } 

   # Clean up as requested
   if (leave_model_in_dbms == FALSE) RODM_drop_model(database, model_name, sql.log.file)
   
   return(dt.list)
} # End of RODM_create_dt_model

