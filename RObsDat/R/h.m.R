h.m <- function(object){ #handler.message
	requireNamespace("RSQLite")
	m <- dbDriver("SQLite")
	con <- dbConnect(m, dbname = "RODM.db")

	#dbGetQuery(con, "SELECT * FROM Versions")
	sqhandler <-  new("odm1_1Ver", con=con)
	options(odm.handler=sqhandler)
        stop("Setting default database scheme (odm1_1 with version management) and file RODM.db\nExecute command again")

}

