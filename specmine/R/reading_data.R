# reads a dataset from CSV files: one for data and (optionally) one for metadata

"read_dataset_csv" = function(filename.data, filename.meta= NULL, type = "undefined", description = "", 
                              label.x = NULL, label.values = NULL, sample.names = NULL,
                              format = "row", header.col = TRUE, header.row = TRUE, sep = ",", 
                              header.col.meta = TRUE, header.row.meta = TRUE, sep.meta = ",")
{
  if (!is.null(filename.meta))
    metadata = read_metadata(filename.meta, header.col = header.col.meta, header.row = header.row.meta, sep = sep.meta)
  else metadata = NULL
  
  data = read_data_csv(filename.data, format = format, header.col = header.col, header.row = header.row, sep = sep)
  if (!is.null(sample.names))
      dataset = create_dataset(data, type = type, metadata = metadata, description = description, 
                             sample.names = sample.names, label.x = label.x, label.values = label.values)
  else
      dataset = create_dataset(data, type = type, metadata = metadata, description = description, 
                           label.x = label.x, label.values = label.values)
  dataset
}


# reads a CSV file creating a data matrix 
# format: "row" -> samples are in rows (if header.row is T, sample names are in the first column;
#               if header.col is T, value labels are in the first row)
#   	    "column" -> samples in columns (if header.col is T, sample names are in the first row;
#               if header.row is T, value labels are in the first column)
# returns dataset with the default structure

"read_data_csv" = function(filename, format = "row", header.col = TRUE, header.row = TRUE, sep = ",")
{
  if (header.row) rownames = 1
  else rownames = NULL
  df = read.table(filename, header = header.col, row.names = rownames, sep = sep, check.names = FALSE)
  
  if( sum(apply(df,c(1,2), is.numeric)) != nrow(df) * ncol(df) )
    stop("There are non-numerical values in data")
  
  if (format == "row"){
    data = as.matrix(t(df))
  } 
  else if (format == "col") {
    data = as.matrix(df)
  }
  data
}

"read_metadata" = function(filename, header.col = T, header.row = T, sep = ",")
{
  if (header.row) rownames = 1
  else rownames = NULL
  metadata = read.table(filename, header = header.col, row.names = rownames, sep = sep, check.names = FALSE)
  metadata
}