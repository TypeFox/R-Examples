addNbAcc <- function(tip, dbname = "fungi", pw = ""){
  conn <- dbConnect(PostgreSQL(), dbname = dbname, user = "postgres", password = pw)
  tip <- dbGetQuery(conn, paste("SELECT * FROM taxonomy WHERE spec='", 
                                tip, "'", sep = ""))
  dbDisconnect(conn)
  tip <- tip[, grep("spec|_sel", names(tip))]
  names(tip) <- gsub("_sel", "", names(tip))
  tip[is.na(tip)] <- 0
  tip <- tip[, tip[1, ] > 0]
  tip <- paste(names(tip), tip, sep = "=")
  tip <- gsub("^(spec=|_)", "", tip)
  tip[1] <- gsub("_", " ", tip[1])
  tip[1] <- paste("italic(\"", tip[1], "\")", sep = "")
  tip[-1] <- paste("plain(\"", tip[-1], "\")", sep = "")
  tip <- paste(tip, collapse = "*\"   \"*")
  parse(text = tip)
}