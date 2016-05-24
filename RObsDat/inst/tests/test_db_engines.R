#options(testPostgreSQL = TRUE, testMySQL = TRUE)
#options(testMySQL = TRUE)
#options(testPostgreSQL = TRUE)
context("db_engines")

test_that("postgreSQL works", {
	#Setup
	if(!as.logical(getOption("testPostgreSQL", FALSE))){
		cat("Not testing example on PostgreSQL. Use 'options(testPostgreSQL = TRUE)' to do so\n")
		return()
	}
		require("RObsDat")
		on.exit( detach("package:RPostgreSQL", unload=TRUE))
		require("RPostgreSQL")
		postDrive <- dbDriver("PostgreSQL")
		conpost <- dbConnect(postDrive, user="reusser", password="4nsp", dbname="obsdat_test", port="5433", host="localhost")
		sqhandler <-  new("odm1_1Ver", con=conpost)
		options(odm.handler=sqhandler)
		on.exit({
			cat("Exiting and cleaning PostgreSQL database\n")
			dbGetQuery(conpost, "DROP SCHEMA public CASCADE;")
			dbGetQuery(conpost, "CREATE SCHEMA public ;")
			odm.close()
			rm(postDrive, conpost)
			detach("package:RPostgreSQL", unload=TRUE)
			cat("Exiting and cleaning finished\n")
				})

		exampleCommands()
		})

test_that("MySQL works", {
	#Setup
	if(!as.logical(getOption("testMySQL", FALSE))){
		cat("Not testing example on MySQL. Use 'options(testMySQL = TRUE)' to do so\n")
		return()
	}
		require("RObsDat")
		on.exit( detach("package:RMySQL", unload=TRUE))
		require("RMySQL")
		sqDrive <- dbDriver("MySQL")
		mycon <- dbConnect(sqDrive, user="a_user", password="secret", dbname="obsdat_test")
		sqhandler <-  new("odm1_1Ver", con=mycon)
		options(odm.handler=sqhandler)
		on.exit({
			cat("Exiting and cleaning MySQL database\n")
			cleanupMySQL(mycon)
			odm.close()
			rm(sqDrive, mycon)
			detach("package:RMySQL", unload=TRUE)
			cat("Exiting and cleaning finished\n")
				})


		exampleCommands()


		})

