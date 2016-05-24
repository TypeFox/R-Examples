
setClass("H2Driver", contains = "JDBCDriver")
setClass("H2Connection", contains = "JDBCConnection")
setClass("H2Result", contains = "JDBCResult")

# H2("RH2", jar = NULL, morePaths = "C:/tmp2/h2-1.3.155.jar") can be used
# to specify a different jar file.
H2 <- function(driverClass='org.h2.Driver', 
               identifier.quote="\"", jars = getOption("RH2.jars"), ...) {
  #identifier.quote="\"", parameters = "-Dh2.identifiersToUpper=false", ...) {
  # options(java.parameters=parameters)
  if (is.null(jars)) jars <- "*"
  .jpackage("RH2", jars = jars, ...)
  if (nchar(driverClass) && is.jnull(.jfindClass(as.character(driverClass)[1])))
    stop("Cannot find H2 driver class ",driverClass)
  jdrv <- .jnew(driverClass, check=FALSE)
  .jcheck(TRUE)
  if (is.jnull(jdrv)) jdrv <- .jnull()
  new("H2Driver", identifier.quote=as.character(identifier.quote), jdrv=jdrv)
}

setMethod("dbConnect", "H2Driver", def=function(drv, url = "jdbc:h2:mem:", 
                                                user='sa', password='', DATABASE_TO_UPPER = getOption("RH2.DATABASE_TO_UPPER"), ...) {
  if (is.null(DATABASE_TO_UPPER) ||
        (is.logical(DATABASE_TO_UPPER) && !DATABASE_TO_UPPER) ||
        (is.character(DATABASE_TO_UPPER) && DATABASE_TO_UPPER == "FALSE"))
    url <- paste(url, "DATABASE_TO_UPPER=FALSE", sep = ";")
  jc <- .jcall("java/sql/DriverManager","Ljava/sql/Connection;","getConnection", as.character(url)[1], as.character(user)[1], as.character(password)[1], check=FALSE)
  if (is.jnull(jc) || !is.jnull(drv@jdrv)) {
    # ok one reason for this to fail is its interaction with rJava's
    # class loader. In that case we try to load the driver directly.
    oex <- .jgetEx(TRUE)
    p <- .jnew("java/util/Properties")
    if (length(user)==1 && nchar(user)) .jcall(p,"Ljava/lang/Object;","setProperty","user",user)
    if (length(password)==1 && nchar(password)) .jcall(p,"Ljava/lang/Object;","setProperty","password",password)
    jc <- .jcall(drv@jdrv, "Ljava/sql/Connection;", "connect", as.character(url)[1], p)
  }
  .verify.JDBC.result(jc, "Unable to connect JDBC to ",url)
  new("H2Connection", jc=jc, identifier.quote=drv@identifier.quote)},
  valueClass="H2Connection")

setMethod("dbWriteTable", "H2Connection", def=function(conn, name, value, overwrite=TRUE, ...) {
  dots <- list(...)
  temporary <- "temporary" %in% names(dots) && dots$temporary
  ac <- .jcall(conn@jc, "Z", "getAutoCommit")
  if (is.vector(value) && !is.list(value)) value <- data.frame(x=value)
  if (length(value)<1) stop("value must have at least one column")
  if (is.null(names(value))) names(value) <- paste("V",1:length(value),sep='')
  if (length(value[[1]])>0) {
    if (!is.data.frame(value)) value <- as.data.frame(value, row.names=1:length(value[[1]]))
  } else {
    if (!is.data.frame(value)) value <- as.data.frame(value)
  }
  fts <- sapply(value, dbDataType, dbObj=conn)
  if (dbExistsTable(conn, name)) {
    if (overwrite) dbRemoveTable(conn, name)
    else stop("Table `",name,"' already exists")
  }
  fdef <- paste(.sql.qescape(names(value), FALSE, conn@identifier.quote),fts,collapse=',')
  # qname <- .sql.qescape(name, TRUE, conn@identifier.quote)
  # cat("conn@identifier.quote:", conn@identifier.quote, "\n")
  qname <- .sql.qescape(name, FALSE, conn@identifier.quote)
  ct <- paste(if (temporary) "CREATE TEMPORARY TABLE" else "CREATE TABLE ",
              qname," (",fdef,")",sep= '')
  # cat("ct:", ct, "\n")
  if (ac) {
    .jcall(conn@jc, "V", "setAutoCommit", FALSE)
    on.exit(.jcall(conn@jc, "V", "setAutoCommit", ac))
  }
  dbSendUpdate(conn, ct)
  if (length(value[[1]])) {
    inss <- paste("INSERT INTO ",qname," VALUES(", paste(rep("?",length(value)),collapse=','),")",sep='')
    # send Date variables as character strings
    is.Date <- sapply(value, inherits, what = "Date")
    for(i in which(is.Date)) {
      value[[i]] <- as.character(value[[i]])
    }
    # send times variables as character strings
    is.times <- sapply(value, function(x) identical(class(x), "times"))
    for(i in which(is.times)) {
      value[[i]] <- as.character(value[[i]])
    }
    if (NCOL(value) > 0) {
      for (j in seq_along(value[[1]]))
        dbSendUpdate(conn, inss, list=as.list(value[j,]))
    }
  }
  if (ac) dbCommit(conn)            
})

setMethod("dbDataType", signature(dbObj="H2Connection", obj = "ANY"),
          def = function(dbObj, obj, ...) {
            if (is.integer(obj)) "INTEGER"
            else if (inherits(obj, "Date")) "DATE"
            else if (identical(class(obj), "times")) "TIME"
            else if (inherits(obj, "POSIXct")) "TIMESTAMP"
            else if (is.numeric(obj)) "DOUBLE PRECISION"
            else "VARCHAR(255)"
          }, valueClass = "character")

setMethod("fetch", signature(res="H2Result", n="numeric"), def=function(res, n, ...) {
  # Column count
  cols <- .jcall(res@md, "I", "getColumnCount")
  
  # Row count  
  nc <- .jcall(res@jr, "I", "getRow")
  .jcall(res@jr, "Z", "absolute",as.integer(-1)) # Last row
  nl <- .jcall(res@jr, "I", "getRow")
  .jcall(res@jr, "Z", "absolute",as.integer(nc)) # Return to previous row
  nn <- nl - nc # Remaining rows
  if (n<0 | nn<n) n <- nn
  
  # Get column names
  cns <- sapply( 1:cols,function(i) .jcall(res@md, "S", "getColumnName", i) )
  
  # Get column types
  cts <- sapply( 1:cols,function(i) .jcall(res@md, "I", "getColumnType", i) )
  
  # Corresponding R data types
  rts <- rep("character",cols)
  I = which(cts == -5 | cts ==-6 | (cts >= 2 & cts <= 8))
  rts[I] <- "numeric"
  I = which(cts == 91)
  rts[I] <- "Date"
  I = which(cts == 92)
  rts[I] <- "times"
  I = which(cts == 93)
  rts[I] <- "POSIXct"
  
  jts <- c(numeric="getDouble",character="getString",Date="getString",
           times="getString",POSIXct="getString")[rts]
  jnis <- c(numeric="D",character="S",Date="S",times="S",POSIXct="S")[rts] 
  pps <- c(numeric=as.numeric,character=as.character,Date=as.Date,times=times,
           POSIXct=as.POSIXct)[rts]
  
  # Pre-create result list
  l <- list()
  for (i in 1:cols) l[[i]] = (pps[[i]])(rep(NA,n))
  
  names(l) <- cns
  
  if (n==0) return(as.data.frame(l,stringsAsFactors=FALSE))
  
  for (j in 1:n) {
#     if (j %% 1000==0) print(paste("Retrieving row",j))
    if (!(.jcall(res@jr, "Z", "next")))
      stop("Row not found") # Should never happen
    
    for (i in 1:cols) {
      val <- .jcall(res@jr, jnis[[i]], jts[[i]], i)
      if (!(.jcall(res@jr, "Z", "wasNull"))) l[[i]][j] <- (pps[[i]])(val)
    }
  }
  
  # Just changes attributes to avoid large copy
  attr(l, 'class') <- 'data.frame'
  attr(l, "row.names") <- c(NA_integer_,n)
  return(l)
})

setMethod("dbSendQuery", signature(conn="H2Connection", statement="character"),  def=function(conn, statement, ..., list=NULL) {
  s <- .jcall(conn@jc, "Ljava/sql/PreparedStatement;", "prepareStatement", as.character(statement)[1],
              as.integer(.jfield("java/sql/ResultSet","I","TYPE_SCROLL_INSENSITIVE")), 
              as.integer(.jfield("java/sql/ResultSet","I","CONCUR_READ_ONLY")), 
              check=FALSE)
  .verify.JDBC.result(s, "Unable to execute JDBC statement ",statement)
  if (length(list(...))) .fillStatementParameters(s, list(...))
  if (!is.null(list)) .fillStatementParameters(s, list)
  r <- .jcall(s, "Ljava/sql/ResultSet;", "executeQuery", check=FALSE)
  .verify.JDBC.result(r, "Unable to retrieve JDBC result set for ",statement)
  md <- .jcall(r, "Ljava/sql/ResultSetMetaData;", "getMetaData", check=FALSE)
  .verify.JDBC.result(md, "Unable to retrieve JDBC result set meta data for ",statement, " in dbSendQuery")
  new("H2Result", jr=r, md=md)
})

setMethod("dbGetQuery", signature(conn="H2Connection", statement="character"),  def=function(conn, statement, ...) {
  r <- dbSendQuery(conn, statement, ...)
  fetch(r, -1)
})

.verify.JDBC.result <- RJDBC:::.verify.JDBC.result
.fillStatementParameters <- RJDBC:::.fillStatementParameters
.sql.qescape <- RJDBC:::.sql.qescape

