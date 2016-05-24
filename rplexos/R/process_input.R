#' @rdname process_folder
#' 
#' @useDynLib rplexos
#' @export
process_input <- function(file) {
  # Check that inputs are valid
  stopifnot(is.character(file), length(file) == 1L)
  
  # Check that file exists
  if (!file.exists(file)) {
    warning(file, " does not exist and was ignored.", call. = FALSE, immediate. = TRUE)
    return(invisible(""))
  }

  # Database name will match that of the zip file
  db.name <- gsub(".xml|.XML", "-input.db", file)
  
  # Delete old file, if possible
  if (file.exists(db.name)) {
    stop_ifnot_delete(db.name)
  }
  
  # Read content from the XML file
  read.con <- file(file, open = "r")
  xml.content.temp <- NULL
  try(xml.content.temp <- readLines(read.con, warn = FALSE))
  if (is.null(xml.content.temp)) {
    stop("Error reading XML file into memory", call. = FALSE)
  }
  xml.content <- paste(xml.content.temp, collapse = " ")
  close(read.con)
  
  # Check that XML is a valid PLEXOS file
  plexos.check <- grep("MasterDataSet", xml.content)
  if (length(plexos.check) == 0L) {
    rplexos_message("Invalid XML content in ", file)
    warning(file, " is not a PLEXOS input file and was ignored.", call. = FALSE, immediate. = TRUE)
    return(invisible(""))
  }
  
  # Create an empty database and add the XML information
  rplexos_message("  - Input: '", file, "'", sep = "")
  
  # Open connection to SQLite for R
  dbf <- src_sqlite(db.name, create = TRUE)
  
  # Add basic XML structure and delete cached XML file
  new_database(dbf, xml.content, is.solution = FALSE)
  rm(xml.content)
  
  # Add a few tables that will be useful later on
  add_extra_tables_input(dbf)
  
  # Close database connections
  DBI::dbDisconnect(dbf$con)
  
  # Message that file processing is done
  rplexos_message("Finished processing file ", file, "\n")
  
  # Return the name of the database that was created
  invisible(db.name)
}

# Add a few tables that will be useful later on
add_extra_tables_input <- function(db) {
  rplexos_message("Adding extra tables to the database")
  
  # View with class and class_group
  sql <- "CREATE VIEW [class] AS
          SELECT c.class_id class_id,
                 c.name,
                 g.name class_group,
                 c.is_enabled,
                 c.description
          FROM t_class c
          INNER JOIN t_class_group g
            ON c.class_group_id = g.class_group_id"
  DBI::dbGetQuery(db$con, sql)
  
  # View with object, category, class, class_group
  sql <- "CREATE VIEW [object] AS
          SELECT o.object_id,
                 o.name name,
                 o.description,
                 cat.name category,
                 c.class_group,
                 c.name class,
                 c.class_id
          FROM t_object o
          JOIN [class] c
            ON o.class_id = c.class_id
          JOIN t_category cat
            ON o.category_id = cat.category_id"
  DBI::dbGetQuery(db$con, sql)
  
  # View with property
  sql <- "CREATE VIEW [property] AS
          SELECT p.property_id,
                 p.name,
                 p.default_value,
                 p.period_type_id,
                 p.is_enabled,
                 p.is_dynamic,
                 p.is_multi_band,
                 p.max_band_id,
                 pg.name property_group,
                 c.name collection,
                 u.value unit,
                 p.description
          FROM t_property p
          JOIN t_property_group pg
            ON p.property_group_id = pg.property_group_id
          JOIN t_collection c
            ON c.collection_id = p.collection_id
         JOIN t_unit u
           ON u.unit_id = p.unit_id"
  DBI::dbGetQuery(db$con, sql)
  
  # View for attribute
  sql <- "CREATE VIEW [attribute] AS
          SELECT a.attribute_id,
                 a.class_id,
                 a.name,
                 u.value unit,
                 a.default_value,
                 a.is_enabled,
                 a.description
          FROM t_attribute a
          JOIN t_unit u
            ON a.unit_id = u.unit_id"
  DBI::dbGetQuery(db$con, sql)

  # View for attribute data
  sql <- "CREATE VIEW [attribute_data] AS
          SELECT o.class_group,
                 o.class,
                 o.name,
                 a.name attribute,
                 a.unit,
                 a.default_value,
                 d.value given_value,
                 ifnull(d.value, a.default_value) value
          FROM object o
          INNER JOIN attribute a
            ON a.class_id = o.class_id
          LEFT JOIN t_attribute_data d
            ON d.object_id = o.object_id 
            AND d.attribute_id = a.attribute_id"
  DBI::dbGetQuery(db$con, sql)
  
  # View with memberships, collection, parent and child objects
  col.table <- tbl(db, "t_collection")
  if ("complement_name" %in% col.table$select) {
    txt.comp <- "c.complement_name"
  } else {
    txt.comp <- ""
  }
  
  sql <- sprintf("CREATE VIEW [membership] AS
          SELECT m.membership_id,
                 m.parent_object_id parent_object_id,
                 m.child_object_id child_object_id,
                 c.name collection,
                 %s comp_collection,
                 p.name parent_name,
                 p.class parent_class,
                 p.class_group parent_group,
                 p.category parent_category,
                 ch.name child_name,
                 ch.class child_class,
                 ch.class_group child_group,
                 ch.category child_category
          FROM t_membership m
          JOIN t_collection c
            ON c.collection_id = m.collection_id
          JOIN [object] p
            ON p.object_id = m.parent_object_id
          JOIN [object] ch
            ON ch.object_id = m.child_object_id", txt.comp)
  DBI::dbGetQuery(db$con, sql)
  
  # Create table with all the tags
  rplexos_message("Creating tag table")
  if (db_has_table(db$con, "t_tag")) {
    # Query and reformat the data
    sql <- "CREATE VIEW [temp_tag] AS
            SELECT t.data_id, o.category, o.class, o.name
            FROM t_tag t
            JOIN object o
            ON t.object_id = o.object_id"
    DBI::dbGetQuery(db$con, sql)
    
    tag.table <- tbl(db, "temp_tag") %>%
      collect %>%
      tidyr::spread(class, name) %>%
      as.data.frame
    
    DBI::dbGetQuery(db$con, "DROP VIEW [temp_tag]")
    
    # Avoid space in the table name
    names(tag.table) <- gsub("Data File", "DataFile", names(tag.table))
    
    # Add missing names
    for (col in c("DataFile", "Escalator", "Scenario")) {
      if (!col %in% names(tag.table))
        tag.table[[col]] <- rep(NA, nrow(tag.table))
    }
  } else {
    tag.table <- data.frame(data_id = integer(0),
                            DataFile = character(0),
                            Escalator = character(0),
                            Scenario = character(0))
  }
  
  DBI::dbWriteTable(db$con, "table_tag", tag.table, row.names = FALSE)
  
  # Create table with data text entries
  rplexos_message("Creating text table")
  if (db_has_table(db$con, "t_text")) {
    # Query and reformat the data
    sql <- "CREATE VIEW [temp_text] AS
            SELECT t.data_id, c.name class, t.value
            FROM t_text t
            JOIN t_class c
            ON t.class_id = c.class_id"
    DBI::dbGetQuery(db$con, sql)
      
    text.table <- tbl(db, "temp_text") %>%
      collect %>%
      tidyr::spread(class, value) %>%
      as.data.frame
    
    DBI::dbGetQuery(db$con, "DROP VIEW [temp_text]")
    
    # Avoid space in the table name
    names(text.table) <- gsub("Data File", "DataFile", names(text.table))
    
    # Add missing names
    for (col in c("DataFile", "Scenario", "Timeslice")) {
      if (!col %in% names(text.table))
        text.table[[col]] <- rep(NA, nrow(text.table))
    }
  } else {
    text.table <- data.frame(data_id = integer(0),
                             DataFile = character(0),
                             Scenario = character(0),
                             Timeslice = character(0))
  }
  
  DBI::dbWriteTable(db$con, "table_text", text.table, row.names = FALSE)
  
  # Add t_memo_data if it doesn't exist
  if (!db_has_table(db$con, "t_memo_data")) {
    rplexos_message("Adding t_memo_data")
    memo.table <- data.frame(data_id = integer(0),
                             value = character(0))
    DBI::dbWriteTable(db$con, "t_memo_data", memo.table, row.names = FALSE)
  }
  
  # Add t_band if it doesn't exist
  if (!db_has_table(db$con, "t_band")) {
    rplexos_message("Adding t_band")
    band.table <- data.frame(data_id = integer(0),
                             band_id = integer(0))
    DBI::dbWriteTable(db$con, "t_band", band.table, row.names = FALSE)
  }
  
  # Add t_date_from if it doesn't exist
  if (!db_has_table(db$con, "t_date_from")) {
    rplexos_message("Adding t_date_from")
    datefrom.table <- data.frame(data_id = integer(0),
                                 date = character(0))
    DBI::dbWriteTable(db$con, "t_date_from", datefrom.table, row.names = FALSE)
  }
  
  # Add t_date_to if it doesn't exist
  if (!db_has_table(db$con, "t_date_to")) {
    rplexos_message("Adding t_date_to")
    dateto.table <- data.frame(data_id = integer(0),
                               date = character(0))
    DBI::dbWriteTable(db$con, "t_date_to", dateto.table, row.names = FALSE)
  }
  
  # Add data view
  sql <- "CREATE VIEW data AS
          SELECT m.collection, m.parent_class, m.parent_group, m.parent_category, m.parent_name,
                 m.child_class, m.child_group, m.child_category, m.child_name,
                 p.name property,
                 d.value,
                 p.default_value,
                 p.period_type_id,
                 p.unit,
                 df.date date_from, dt.date date_to,
                 table_text.Timeslice Timeslice,
                 table_tag.Escalator Escalator,
                 ifnull( table_tag.DataFile, table_text.DataFile )  DataFile,
                 ifnull( table_tag.Scenario, table_text.Scenario )  Scenario,
                 ifnull( b.band_id, 1 )  band,
                 memo.value memo
          FROM t_data d
          INNER JOIN membership m
               ON m.membership_id = d.membership_id
          INNER JOIN property p
               ON p.property_id = d.property_id
          LEFT JOIN t_date_from df
               ON d.data_id = df.data_id
          LEFT JOIN t_date_to dt
               ON d.data_id = dt.data_id
          LEFT JOIN table_tag
               ON d.data_id = table_tag.data_id
          LEFT JOIN table_text
               ON d.data_id = table_text.data_id
          LEFT JOIN t_band b
               ON d.data_id = b.data_id
          LEFT JOIN t_memo_data memo
               ON d.data_id = memo.data_id"
  DBI::dbGetQuery(db$con, sql)
  
  0
}
