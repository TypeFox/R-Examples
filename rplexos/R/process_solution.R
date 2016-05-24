#' @rdname process_folder
#' 
#' @useDynLib rplexos
#' @importFrom utils packageVersion unzip
#' @export
process_solution <- function(file, keep.temp = FALSE) {
  # Check that inputs are valid
  stopifnot(is.character(file), length(file) == 1L)
  
  # Check that file exists
  if (!file.exists(file)) {
    warning(file, " does not exist and was ignored.", call. = FALSE, immediate. = TRUE)
    return(invisible(""))
  }

  # Database name will match that of the zip file
  db.temp <- gsub(".zip", "-temp.db", file)
  db.name <- gsub(".zip", "-rplexos.db", file)
  
  # Delete old files, if possible
  if (file.exists(db.temp)) {
    stop_ifnot_delete(db.temp)
  }
  if (file.exists(db.name)) {
    stop_ifnot_delete(db.name)
  }
  
  # Read list of files in the zip file
  zip.content <- unzip(file, list = TRUE)
  
  # Check that zip file has valid XML, Log and BIN files
  xml.pos <- grep("^Model.*xml$", zip.content$Name)
  bin.pos <- grep("^t_data_[0-4].BIN$", zip.content$Name)
  log.pos <- grep("^Model.*Log.*.txt$", zip.content$Name)
  if ((length(xml.pos) == 0L) | (length(bin.pos) == 0L) | (length(log.pos) == 0L)) {
    # If in debug mode, give some more explanation
    if (length(xml.pos) == 0L)
      rplexos_message("No XML file found in file ", file)
    if (length(bin.pos) == 0L)
      rplexos_message("No BIN file found in file ", file)
    if (length(log.pos) == 0L)
      rplexos_message("No LOG file found in file ", file)
    
    warning(file, " is not a PLEXOS solution file and was ignored.", call. = FALSE, immediate. = TRUE)
    return(invisible(""))
  }
  
  # Check if the interval binary file is not being read correctly within the zip file
  #   This seems to happen when R is using 32-bit versions of the zip libraries
  #   No clear solution as to how to solve this yet
  bin.int.pos <- grep("^t_data_0", zip.content$Name)
  if (length(bin.int.pos) == 1L) {
    if (zip.content$Length[bin.int.pos] == (2^32 - 1)) {
      warning("Skipped file ", file, "\n",
              "  Interval data is too large and cannot be processed in Mac/Linux.\n",
              "  Reduce number of outputs in PLEXOS or process with the 64-bit Windows version of R.\n",
              "  No need to report this message; we are working on a solution.",
              call. = FALSE, immediate. = TRUE)
      
      return(invisible(""))
    }
  }
  
  # Read content from the XML file
  xml.content <- try(read_file_in_zip(file, xml.pos), silent = !is_debug_rplexos())
  if (inherits(xml.content, "try-error")) {
    stop("Error reading XML file into memory", call. = FALSE)
  }
  
  # Check that XML is a valid PLEXOS file
  plexos.check <- grep("SolutionDataset", xml.content[1])
  if (length(plexos.check) == 0L) {
    rplexos_message("Invalid XML content in ", file)
    warning(file, " is not a PLEXOS database and was ignored.", call. = FALSE, immediate. = TRUE)
    return(invisible(""))
  }
  
  # Create an empty database and add the XML information
  rplexos_message("  - Solution: '", file, "'", sep = "")
  
  # Open connection to SQLite for R
  dbt <- src_sqlite(db.temp, create = TRUE)
  
  # Add basic XML structure and delete cached XML file
  new_database(dbt, xml.content)
  rm(xml.content)
  
  # Add a few tables that will be useful later on
  add_extra_tables(dbt)
  
  # Create SQLite database to store final results
  rplexos_message("Creating final database and adding basic data")
  dbf <- src_sqlite(db.name, create = TRUE)
  
  # Store time stamps
  sql <- "CREATE TABLE data_time (phase_id INT, interval INT, time TEXT)"
  DBI::dbGetQuery(dbf$con, sql)
  sql <- "CREATE VIEW time AS
          SELECT phase_id, interval, datetime(time) time
          FROM data_time"
  DBI::dbGetQuery(dbf$con, sql)
  
  # Turn PRAGMA OFF
  DBI::dbGetQuery(dbf$con, "PRAGMA synchronous = OFF")
  DBI::dbGetQuery(dbf$con, "PRAGMA journal_mode = MEMORY")
  DBI::dbGetQuery(dbf$con, "PRAGMA temp_store = MEMORY")
  
  # Attach final database to temporary database
  DBI::dbGetQuery(dbt$con, sprintf("ATTACH '%s' AS new", db.name))
  
  # Add config table
  DBI::dbGetQuery(dbt$con, "CREATE TABLE new.config AS SELECT * FROM t_config")
  sql <- sprintf("INSERT INTO new.config VALUES ('rplexos', '%s')", packageVersion("rplexos"))
  DBI::dbGetQuery(dbt$con, sql)
  
  # Add time data
  sql <- "INSERT INTO new.data_time
          SELECT phase_id, interval_id, time
          FROM temp_period_0"
  DBI::dbGetQuery(dbt$con, sql)
  
  # Collate information to key (first period data, then summary data)
  sql <- "CREATE TABLE new.key (key INT PRIMARY KEY,
                                table_name TEXT,
                                collection TEXT,
                                property TEXT,
                                unit TEXT,
                                name TEXT,
                                parent TEXT,
                                category TEXT,
                                region TEXT,
                                zone TEXT,
                                class TEXT,
                                class_group TEXT,
                                phase_id INT,
                                period_type_id INT,
                                timeslice TEXT,
                                band INT,
                                sample TEXT)"
  DBI::dbGetQuery(dbt$con, sql)
  
  sql <- "INSERT INTO new.key
          SELECT *
          FROM temp_key"
  DBI::dbGetQuery(dbt$con, sql)
  
  # Detach database
  DBI::dbGetQuery(dbt$con, "DETACH new");
  
  # Define columns from the key table that go into the views
  view.k <- paste0("k.", c("key", "collection", "property", "unit", "name", "parent", "category",
                            "region", "zone", "phase_id", "period_type_id", "timeslice",
                            "band", "sample"))
  view.k2 <- paste(view.k, collapse = ", ")
  
  # For each summary time, create a table and a view
  rplexos_message("Creating data tables and views")
  times <- c("day", "week", "month", "year")
  for (i in times) {
    sql <- sprintf("CREATE TABLE data_%s (key integer, time real, value double)", i);
    DBI::dbGetQuery(dbf$con, sql)
    
    sql <- sprintf("CREATE VIEW %s AS    
                    SELECT %s, datetime(d.time) AS time, d.value 		
                    FROM data_%s d NATURAL LEFT JOIN key k ", i, view.k2, i);		
    DBI::dbGetQuery(dbf$con, sql)
  }
  
  # Create interval data tables and views
  sql <- "SELECT DISTINCT table_name
          FROM key
          WHERE period_type_id = 0"
  props <- DBI::dbGetQuery(dbf$con, sql)
  
  for (p in props$table_name) {
    sql <- sprintf("CREATE TABLE '%s' (key INT, time_from INT, time_to INT, value DOUBLE)", p)
    DBI::dbGetQuery(dbf$con, sql)
    
    view.name <- gsub("data_interval_", "", p)
    sql <- sprintf("CREATE VIEW %s AS
                    SELECT %s, t1.time time_from, t2.time time_to, d.value 
                    FROM %s d
                    NATURAL JOIN key k
                    JOIN time t1
                      ON t1.interval = d.time_from
                     AND t1.phase_id = k.phase_id
                    JOIN time t2
                      ON t2.interval = d.time_to
                     AND t2.phase_id = k.phase_id
                    WHERE k.table_name = '%s'", view.name, view.k2, p, p);
    DBI::dbGetQuery(dbf$con, sql)
  }
  
  # Create table for list of properties
  sql <- "CREATE TABLE property AS
          SELECT DISTINCT class_group,
                          class,
                          collection,
                          property,
                          unit,
                          phase_id,
                          period_type_id AS is_summary,
                          table_name,
                          COUNT(DISTINCT band) AS count_band,
                          COUNT(DISTINCT sample) AS count_sample,
                          COUNT(DISTINCT timeslice) AS count_timeslice
          FROM key
          GROUP BY class_group, class, collection, property, unit, phase_id, period_type_id, table_name
          ORDER BY phase_id, period_type_id, class_group, class, collection, property"
  DBI::dbGetQuery(dbf$con, sql)
  
  # Add binary data
  for (period in 0:4) {
    # Check if binary file exists, otherwise, skip this period
    period.name <- ifelse(period == 0, "interval", times[period])
    bin.name <- sprintf("t_data_%s.BIN", period)
    if(!bin.name %in% zip.content$Name)
      next
    bin.con <- unz(file, bin.name, open = "rb")
    
    # Print debug message
    rplexos_message("Reading ", period.name, " binary data")
    
    # Check if length in t_key_index is correct
    correct.length <- correct_length(dbt, period)
    
    # Read time data
    t.time <- tbl(dbt, sprintf("temp_period_%s", period)) %>%
      arrange(phase_id, period_id) %>%
      select(phase_id, period_id, time, interval_id) %>%
      collect()
    
    # Read t_key_index entries for period data
    sql <- sprintf("SELECT nk.[key], nk.phase_id, nk.table_name, tki.period_offset, tki.length
                    FROM t_key_index tki
                    JOIN temp_key nk
                    ON tki.key_id = nk.[key]
                    WHERE tki.period_type_id = %s
                    ORDER BY tki.position", period)
    tki <- DBI::dbSendQuery(dbt$con, sql)
    
    # All the data is inserted in one transaction
    DBI::dbBegin(dbf$con)
      
    # Read one row from the query
    num.rows <- ifelse(period == 0, 1, 1000)
    trow <- DBI::dbFetch(tki, num.rows)
    num.read <- 0
    
    # Iterate through the query results
    while (nrow(trow) > 0) {
      # Fix length if necessary
      if (!correct.length)
        trow <- trow %>% mutate(length = length - period_offset)
      
      # Expand data
      tdata <- trow %>%
        select(key_id = key, phase_id, period_offset, length) %>%
        expand_tkey
      
      # Query data
      value.data <- readBin(bin.con,
                            "double",
                            n = nrow(tdata),
                            size = 8L,
                            endian = "little")
      num.read <- num.read + length(value.data)
      
      # Check the size of data (they won't match if there is a problem)
      if (length(value.data) < nrow(tdata)) {
        rplexos_message("   ", num.read, " values read")
        stop("Problem reading ", period.name, " binary data (reached end of file).\n",
             "  ", nrow(tdata), " values requested, ", length(value.data), " returned.\n",
             "  This is likely a bug in rplexos. Please report it.", call. = FALSE)
      }
      
      # Copy data
      tdata$value <- value.data
      
      # Join with time
      tdata2 <- tdata %>%
        inner_join(t.time, by = c("phase_id", "period_id"))
      
      # Add data to SQLite
      if (period > 0) {
        tdata3 <- tdata2 %>% select(key, time, value)
        
        RSQLite::dbGetPreparedQuery(dbf$con,
          sprintf("INSERT INTO data_%s VALUES(?, ?, ?)", times[period]),
          bind.data = tdata3 %>% as.data.frame)
      } else {
        # Eliminate consecutive repeats
        default.interval.to.id <- max(tdata2$interval_id)
        tdata3 <- tdata2 %>%
          arrange(interval_id) %>%
          filter(value != lag(value, default = Inf)) %>%
          mutate(interval_to_id = lead(interval_id - 1, default = default.interval.to.id)) %>%
          select(key, time_from = interval_id, time_to = interval_to_id, value)
        
        RSQLite::dbGetPreparedQuery(dbf$con,
          sprintf("INSERT INTO %s (key, time_from, time_to, value)
                  VALUES(?, ?, ?, ?)", trow$table_name),
          bind.data = tdata3 %>% as.data.frame)
      }
      
      # Read next row from the query
      trow <- DBI::dbFetch(tki, num.rows)
    }
    
    # Finish transaction
    rplexos_message("   ", num.read, " values read")
    DBI::dbClearResult(tki)
    DBI::dbCommit(dbf$con)
    
    # Close binary file connection
    close(bin.con)
  }
  
  # Read Log file into memory
  rplexos_message("Reading and processing log file")
  log.content <- try(read_file_in_zip(file, log.pos), silent = !is_debug_rplexos())
  if (inherits(log.content, "try-error")) {
    # Error reading log file, throw a warning
    warning("Could not read Log in solution '", file, "'\n",
            "    Data parsed correctly if no other errors were found.",
            call. = FALSE)
  } else {
    # Success reading file, try to parse it
    log.result <- plexos_log_parser(log.content)
    
    if (length(log.result) < 2L) {
      warning("Log in solution '", file, "' did not parse correctly.\n",
              "    Data parsed correctly if no other errors were found.",
              call. = FALSE)
    }
    
    for (i in names(log.result)) {
      DBI::dbWriteTable(dbf$con, i, log.result[[i]] %>% as.data.frame, row.names = FALSE)
    }
  }
  
  # Close database connections
  DBI::dbDisconnect(dbt$con)
  DBI::dbDisconnect(dbf$con)
  
  # Message that file processing is done
  rplexos_message("Finished processing file ", file, "\n")
  
  # Delete temporary database
  if (!keep.temp) {
    stop_ifnot_delete(db.temp)
  }
  
  # Return the name of the database that was created
  invisible(db.name)
}

# Read a file in a zip file onto memory
#' @importFrom utils unzip
read_file_in_zip <- function(zip.file, position) {
  zip.content <- unzip(zip.file, list = TRUE)
  read.file <- zip.content[position, ]
  read.con <- unz(zip.file, read.file$Name, "rb")
  .nBytes <- 2^30
  
  # readChar cannot read files that are very large
  if (read.file$Length > .nBytes) {
    rplexos_message("File '", read.file$Name, "' is large (", read.file$Length, " bytes)")
  }
  
  nchunks <- ceiling(read.file$Length / .nBytes)
  read.content <- character(nchunks)
  
  for (i in seq_len(nchunks)) {
    nread <- min(.nBytes, read.file$Length - (i-1) * .nBytes)
    read.content[i] <- readChar(read.con, nread, TRUE)
  }
  
  close(read.con)
  read.content
}

# Populate new database with XML contents
new_database <- function(db, xml, is.solution = TRUE) {
  rplexos_message("Reading XML file and saving content")
  
  # Read XML and convert to a list
  xml.list <- process_xml(xml)
  
  # Print PLEXOS version when debuging
  ver.pos <- grep("Version|version", xml.list$t_config[, 1])
  if (length(ver.pos) == 1L) {
    rplexos_message("   PLEXOS version:  '", xml.list$t_config[ver.pos, 2], "'")
    rplexos_message("   rplexos version: '", packageVersion("rplexos"), "'")
  } else {
    rplexos_message("   PLEXOS version: N/A")
    rplexos_message("   rplexos version: '", packageVersion("rplexos"), "'")
  }
  
  # Turn PRAGMA OFF
  DBI::dbGetQuery(db$con, "PRAGMA synchronous = OFF");
  DBI::dbGetQuery(db$con, "PRAGMA journal_mode = MEMORY");
  DBI::dbGetQuery(db$con, "PRAGMA temp_store = MEMORY");
  
  # Do the following for solution files
  if (is.solution) {
    # Create data tables
    for (i in 0:4) {
      sql <- sprintf("CREATE TABLE t_data_%s (key_id integer, period_id integer, value double)", i)
      DBI::dbGetQuery(db$con, sql)
    }
    
    # Create phase tables
    for (i in 0:4) {
      sql <- sprintf("CREATE TABLE t_phase_%s (interval_id integer, period_id integer)", i)
      DBI::dbGetQuery(db$con, sql)
    }
    
    # Create t_key_index table
    DBI::dbGetQuery(db$con, "CREATE TABLE t_key_index (key_id integer, period_type_id integer, position long, length integer, period_offset integer)");
  }
    
  # Write tables from XML file
  for (t in names(xml.list))
    DBI::dbWriteTable(db$con, t, xml.list[[t]], append = TRUE, row.names = FALSE)
  
  0
}

# Add a few tables that will be useful later on
add_extra_tables <- function(db) {
  rplexos_message("Adding extra tables to the database")
  
  # View with class and class_group
  sql <- "CREATE VIEW temp_class AS
          SELECT tc.class_id class_id,
                 tc.name class,
                 tcg.name class_group
          FROM t_class tc
          LEFT JOIN t_class_group tcg
            ON tc.class_group_id = tcg.class_group_id"
  DBI::dbGetQuery(db$con, sql)
  
  # View with object, category, class, class_group
  sql <- "CREATE VIEW temp_object AS
          SELECT o.object_id object_id,
                 o.name name,
                 cat.name category,
                 c.class_group class_group,
                 c.class class
          FROM t_object o
          JOIN temp_class c
            ON o.class_id = c.class_id
          JOIN t_category cat
            ON o.category_id = cat.category_id"
  DBI::dbGetQuery(db$con, sql)
  
  # Create a new table making long version of property
  sql <- "CREATE TABLE temp_property AS
          SELECT p.property_id property_id,
                 '0' period_type_id,
                 p.name property,
                 c.name collection,
                 u.value unit
          FROM t_property p
          JOIN t_collection c
            ON c.collection_id = p.collection_id
         JOIN t_unit u
           ON u.unit_id = p.unit_id"
  DBI::dbGetQuery(db$con, sql)
  sql <- "INSERT INTO temp_property
          SELECT p.property_id property_id,
                 '1' period_type_id,
                 p.summary_name property,
                 c.name collection,
                 u.value unit
          FROM t_property p
          JOIN t_collection c
            ON c.collection_id = p.collection_id
         JOIN t_unit u
           ON u.unit_id = p.summary_unit_id"
  DBI::dbGetQuery(db$con, sql)
  
  # View with memberships, collection, parent and child objects
  sql <- "CREATE VIEW temp_membership AS
          SELECT m.membership_id membership_id,
                 m.parent_object_id parent_object_id,
                 m.child_object_id child_object_id,
                 c.name collection,
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
          JOIN temp_object p
            ON p.object_id = m.parent_object_id
          JOIN temp_object ch
            ON ch.object_id = m.child_object_id"
  DBI::dbGetQuery(db$con, sql)
  
  # Views to list zones
  sql <- "CREATE VIEW temp_zones_id AS
          SELECT child_object_id,
                 min(parent_object_id) parent_object_id
          FROM temp_membership
          WHERE collection = 'Generators' 
                AND parent_class = 'Region'
          GROUP BY child_object_id"
  DBI::dbGetQuery(db$con, sql)
  sql <- "CREATE VIEW temp_zones AS
          SELECT a.child_object_id,
                 b.name region,
                 b.category zone
          FROM temp_zones_id a
          JOIN temp_object b
          WHERE a.parent_object_id = b.object_id"
  DBI::dbGetQuery(db$con, sql)
  
  # Read key data and transform it
  sql <- "SELECT k.key_id key,
                 m.child_name child_name,
                 m.parent_name parent_name,
                 m.child_class child_class,
                 m.parent_class parent_class,
                 m.child_group child_group,
                 m.parent_group parent_group,
                 m.child_category child_category,
                 m.parent_category parent_category,
                 m.collection child_collection,
                 p.property property,
                 p.unit unit,
                 ts.name timeslice,
                 k.band_id band,
                 k.sample_id sample,
                 k.period_type_id period_type_id,
                 k.phase_id phase_id,
                 z.region region,
                 z.zone zone
          FROM t_key k
          JOIN temp_membership m
            ON m.membership_id = k.membership_id
          JOIN t_timeslice ts
            ON ts.timeslice_id = k.timeslice_id
          JOIN temp_property p
            ON p.property_id = k.property_id
               AND
               p.period_type_id = k.period_type_id
          LEFT OUTER JOIN temp_zones z
               on z.child_object_id = m.child_object_id"
  key <- DBI::dbGetQuery(db$con, sql) %>%
    mutate(zone = ifelse(is.na(zone), "", zone),
           region = ifelse(is.na(region), "", region))
  
  # Add collection (to match PLEXOS GUI) and table name
  key2 <- key %>%
    mutate(collection = ifelse(parent_class == "System",
                               child_class,
                               paste0(parent_class, ".", child_collection)),
           table_name = paste0("data_interval_", clean_string(collection), "_", clean_string(property)))
  
  # Change sample to identify stats
  key3 <- key2 %>%
    mutate(sample = gsub("^0$",  "Mean",  sample)) %>%
    mutate(sample = gsub("^-1$", "StDev", sample)) %>%
    mutate(sample = gsub("^-2$", "Min",   sample)) %>%
    mutate(sample = gsub("^-3$", "Max",   sample))
  
  # Write results into a new table
  key3 %>%
    select(key, table_name, collection, property, unit, name = child_name,
           parent = parent_name, category = child_category, region, zone,
           class = child_class, class_group = child_group, phase_id, period_type_id,
           timeslice, band, sample) %>%
    DBI::dbWriteTable(db$con, "temp_key", ., row.names = FALSE)
  
  # Check that t_key and temp_key have the same number of rows
  t.key.length <- tbl(db, "t_key") %>% summarize(n = n()) %>% collect
  temp.key.length <- tbl(db, "temp_key") %>% summarize(n = n()) %>% collect
  rplexos_message("   t_key has    ", t.key.length$n,    " rows")
  rplexos_message("   temp_key has ", temp.key.length$n, " rows")
  
  # Create tables to hold interval, day, week, month, and yearly timestamps
  for (i in 0:4) {
    sql <- sprintf("CREATE TABLE temp_period_%s (phase_id INT, period_id INT, interval_id INT, time real)", i)
    DBI::dbGetQuery(db$con, sql)
  }
  
  # For each phase add corresponding values to the time tables
  column.times <- c("day_id", "week_id", "month_id", "fiscal_year_id")
  for (phase in 1:4) {
    # Join t_period_0 and t_phase
    sql <- sprintf("CREATE VIEW temp_phase_%s AS
                    SELECT p.*, ph.period_id, julianday(p.year || '-' || substr(0 || p.month_of_year, -2)
                    || '-' || substr(0 || p.day_of_month, -2) || 'T' || substr(p.datetime, -8)) AS correct_time
                    FROM t_period_0 p
                    INNER JOIN t_phase_%s ph
                    ON p.interval_id = ph.interval_id", phase, phase)
    DBI::dbGetQuery(db$con, sql)
    
    # Fix time stamps in t_period_0 (interval)
    sql <- sprintf("INSERT INTO temp_period_0
                    SELECT %s, period_id, interval_id, correct_time time
                    FROM temp_phase_%s", phase, phase)
    DBI::dbGetQuery(db$con, sql)
    
    # Fix time stamps in t_period_X (summary data)
    for (i in seq_along(column.times)) {
      sql <- sprintf("INSERT INTO temp_period_%s
                      SELECT %s, %s, min(interval_id), min(correct_time) time
                      FROM temp_phase_%s
                      GROUP BY %s", i, phase, column.times[i], phase, column.times[i])
      DBI::dbGetQuery(db$con, sql)
    }
  }
  
  0
}

# Does t_key_index have the correct length?
correct_length <- function(db, p) {
  res <- tbl(db, "t_key_index") %>%
    filter(period_type_id == p) %>%
    summarize(JustLength            = max(position / 8 + length),
              JustLengthMinusOffset = max(position / 8 + length - period_offset),
              SumLength             = sum(length),
              SumLengthMinusOffset  = sum(length - period_offset)) %>%
    collect
  
  # Debug output
  rplexos_message("   Max position:           ", res$JustLength)
  rplexos_message("   Adjusted max position:  ", res$JustLengthMinusOffset)
  rplexos_message("   Sum of length:          ", res$SumLength)
  rplexos_message("   Sum of adjusted length: ", res$SumLengthMinusOffset)
  
  if (res$JustLength == res$SumLength) {
    rplexos_message("   ", res$JustLength, " entries expected in t_data_", p, ".BIN")
    return(TRUE)
  } else if (res$JustLengthMinusOffset == res$SumLengthMinusOffset) {
    rplexos_message("   Length correction is needed")
    rplexos_message("   ", res$JustLengthMinusOffset, " entries expected in t_data_", p, ".BIN")
    return(FALSE)
  }
  
  # This case is 
  warning("Problem with length of 't_key_index' for period ", p, "\n",
          "in file '", db$path, "'",
          call. = FALSE, immediate. = TRUE)
  
  TRUE
}
