check.NameChanges <- function(conn, order.by = "spec"){
  x <- paste("SELECT spec, spec_ncbi FROM taxonomy", 
             "WHERE spec != spec_ncbi AND spec_ncbi IS NOT NULL AND spec_ncbi != 'NA'",
             "ORDER BY", order.by)
  dbGetQuery(conn, x)
}