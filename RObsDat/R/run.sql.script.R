run.sql.script <- function(con, script.file){

	script  <- readLines(script.file)
	script <- gsub("--.*", "", script)
	script <- gsub("\t", "", script)
	script <- paste(script, collapse=" ")
	scriptparts <- strsplit(script, ";")[[1]]


	for(i in seq(along=scriptparts)){

	  statement <- gsub("COMMENT.'[^']*'", "",scriptparts[i] )

	  tablename <- gsub("CREATE TABLE *([[:alpha:]]+).* ","\\1",statement)
	  tablename <- gsub("[[:space:]]*", "", tablename)

	  if(class(con)=="PostgreSQLConnection"){
		statement <- gsub("ENGINE.*", "",statement )
		statement <- gsub("int\\(11\\)  NOT NULL auto_increment", "SERIAL",statement)
		statement <- gsub("int\\(11\\)", "integer",statement)
		statement <- gsub("double", "double precision",statement)
		statement <- gsub("datetime", "timestamp",statement)
	  }
	  if(class(con)=="SQLiteConnection"){
		statement <- gsub("ENGINE.*", "",statement )
		statementa <- gsub("int\\(11\\)  NOT NULL auto_increment", "INTEGER PRIMARY KEY",statement)
		if(statementa != statement){
			statement <- gsub(", *PRIMARY KEY *\\([^)]*\\)","",statementa)
		}
	  }
	  if(!mdbExistsTable(con, tablename)){
		  if(getOption("verbose.queries", default=FALSE)) print(statement)
		  rs1 <- dbSendQuery(con, statement)
	  }
	}

}
