docs_bulk_file <- function(filename) {
  conn <- es_get_auth()
  url <- paste0(conn$base, ":", conn$port, '/_bulk')
  POST(url, body=upload_file(filename), encode = "json")
}

docs_bulk_data <- function(df) {
  conn <- es_get_auth()
  url <- paste0(conn$base, ":", conn$port, '/_bulk')
  POST(url, body=data, encode = "json")
}



make_bulk <- function(df, filename = "hello_world.json") {
  for (i in seq_len(NROW(df))) {
    dat <- list(index = list(`_index` = 'hello', `_type` = 'world', `_id` = i-1))
    cat(proc_doc(dat), sep = "\n", file = filename, append = TRUE)
    cat(proc_doc(df[i,]), sep = "\n", file = filename, append = TRUE)
  }
  message("file written to ", filename)
}

make_bulk2 <- function(df, filename = "hello_world.json") {
  metadata_fmt <- '{"index":{"_index":"hello","_type":"world","_id":%d}}'
  metadata <- sprintf(metadata_fmt, seq_len(nrow(df)) - 1L)
  data <- jsonlite::toJSON(df, collapse=FALSE)
  writeLines(paste(metadata, data, sep="\n"), filename)
}

make_bulk2_cat <- function(df, filename = "hello_world.json") {
  metadata_fmt <- '{"index":{"_index":"hello","_type":"world","_id":%d}}'
  metadata <- sprintf(metadata_fmt, seq_len(nrow(df)) - 1L)
  data <- jsonlite::toJSON(df, collapse=FALSE)
  cat(paste(metadata, data, sep="\n"), file = filename)
}

make_bulk2(iris, "~/hello_world.json")

library(ggplot2)
make_bulk2(diamonds[1:20000L,], "~/hello_world_writelines.json")
make_bulk2_cat(diamonds[1:20000L,], "~/hello_world_cat.json")
