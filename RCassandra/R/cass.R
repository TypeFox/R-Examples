.keyspace.cache <- new.env(TRUE, emptyenv())

RC.connect <- function(host = NULL, port = 9160L) .Call("RC_connect", host, port, PACKAGE="RCassandra")

RC.close <- function(conn) .Call("RC_close", conn, PACKAGE="RCassandra")

RC.use <- function(conn, keyspace, cache.def = TRUE) {
  res <- .Call("RC_use", conn, keyspace, PACKAGE="RCassandra")
  if (isTRUE(cache.def)) .keyspace.cache[[keyspace]] <- RC.describe.keyspace(conn, keyspace)
  .keyspace.cache[["*current keyspace*"]] <- keyspace
  invisible(res)
}    

RC.get <- function(conn, c.family, key, c.names, comparator=NULL, validator=NULL) .Call("RC_get_list", conn, key, c.family, c.names, length(c.names), 0, comparator, validator, PACKAGE="RCassandra")

.current.comparator <- function(c.family) {
  if(!is.null(ks <- .keyspace.cache[["*current keyspace*"]])) {
    ct <- .keyspace.cache[[ks]][["cf_defs"]][[c.family]][["comparator_type"]]
    if (is.character(ct)) return(gsub(".*\\.","",ct[1]))
  }
  NULL
}

RC.get.range <- function(conn, c.family, key, first="", last="", reverse=FALSE, limit=1e7, comparator=NULL, validator=NULL) {
  if (is.null(comparator)) comparator <- .current.comparator(c.family)
  .Call("RC_get_range", conn, key, c.family, first, last, limit, reverse, comparator, validator, PACKAGE="RCassandra")
}

RC.mget.range <- function(conn, c.family, keys, first="", last="", reverse=FALSE, limit=1e7, comparator=NULL, validator=NULL) {
  if (is.null(comparator)) comparator <- .current.comparator(c.family)
  .Call("RC_mget_range", conn, keys, c.family, first, last, limit, reverse, comparator, validator, PACKAGE="RCassandra")
}

RC.get.range.slices <- function(conn, c.family, k.start="", k.end="", first="", last="", reverse=FALSE, limit=1e7, k.limit=1e7, tokens=FALSE, fixed=FALSE, comparator=NULL, validator=NULL) {
  if (is.null(comparator)) comparator <- .current.comparator(c.family)
  .Call("RC_get_range_slices", conn, k.start, k.end, c.family, first, last, limit, reverse, k.limit, tokens, fixed, comparator, validator, PACKAGE="RCassandra")
}

RC.insert <- function(conn, c.family, key, column, value=NULL, comparator=NULL, validator=NULL) {
  if (is.null(comparator)) comparator <- .current.comparator(c.family)
  .Call("RC_insert", conn, key, c.family, column, value, comparator, validator, PACKAGE="RCassandra")
}

RC.mutate <- function(conn, mutation) .Call("RC_mutate", conn, mutation)

RC.cluster.name <- function(conn) .Call("RC_cluster_name", conn, PACKAGE="RCassandra")

RC.version <- function(conn) .Call("RC_version", conn, PACKAGE="RCassandra")

RC.login <- function(conn, username="default", password="") .Call("RC_login", conn, username, password, PACKAGE="RCassandra")

RC.write.table <- function(conn, c.family, df) .Call("RC_write_table", conn, c.family, df, row.names(df), names(df), PACKAGE="RCassandra")

RC.read.table <- function(conn, c.family, convert = TRUE, na.strings = "NA", as.is = FALSE, dec = ".") {
  df <- RC.get.range.slices(conn, c.family, fixed=TRUE)
  if (convert) for (i in seq.int(df)) df[[i]] <- type.convert(df[[i]], na.strings, as.is, dec)
  df
}

RC.consistency <- function(conn, level = c("one", "quorum", "local.quorum", "each.quorum", "all", "any", "two", "three")) {
  level <- match.arg(level)
  level <- match(level, c("one", "quorum", "local.quorum", "each.quorum", "all", "any", "two", "three"))
  invisible(.Call("RC_set_cl", conn, level))
}

.map.KsDef <- function(res) {
  sn <- c("name", "strategy_class", "strategy_options", "replication_factor", "cf_defs", "durable_writes")
  names(res) <- sn[as.integer(names(res))]
  sn <- c("keyspace", "name", "column_type", "4", "comparator_type", "subcomparator_type", "7",
          "9", "10", "11", "comment", "read_repair_chance", "column_metadata", "gc_grace_seconds",
          "default_validation_class", "id", "min_compaction_threshold", "max_compaction_threshold",
          "19", "20", "21", "22", "23", "replicate_on_write", "merge_shards_chance", "key_validation_class",
          "27", "key_alias", "comparaction_strategy", "comparaction_strategy_options", "31",
          "compression_options", "bloom_filter_fp_chance", 34:50)
  cm <- c("name", "validation_class", "index_type", "index_name", "index_options", 6:10)
  res$cf_defs <- lapply(res$cf_defs, function(x) {
    names(x) <- sn[as.integer(names(x))]
    if (length(x$column_metadata))
      names(x$column_metadata) <- cm[as.integer(names(x$column_metadata))]
    x
  })
  names(res$cf_defs) <- sapply(res$cf_defs, function(x) x$name)
  res
}

RC.describe.keyspace <- function(conn, keyspace) {
  res <- .Call("RC_call_ks", conn, "describe_keyspace", keyspace)
  if (is.list(res)) .map.KsDef(res) else res  
}

RC.describe.keyspaces <- function(conn) {
  res <- .Call("RC_call_void", conn, "describe_keyspaces")
  if (is.list(res)) {
    res <- lapply(res, .map.KsDef)
    names(res) <- sapply(res, function(x) x[["name"]])
    res
  } else res
}
