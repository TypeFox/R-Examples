#options(testLongExample = TRUE, testMySQL = TRUE)
#options(testLongExample = TRUE, testPostgreSQL = TRUE, testMySQL = TRUE)
#options(testLongExample = TRUE, testPostgreSQL = TRUE)
context("Example with multiple files")

test_that("LongExample Data works", {
	cat("SQLite\n")
	on.exit( { cat("Exiting and deleting RODM.db\n")
				unlink("RODM.db")})
	unlink("RODM.db")
	getDefaultDB()
	cat("SQLite file initiated\n")
	longExample()
	cat("SQLite example processed\n")
		})

test_that("postgreSQL works with LongExample Data", {
	#Setup
	if(!(as.logical(getOption("testPostgreSQL", FALSE)))){
		cat("Not testing long example on PostgreSQL. Use 'options(testLongExample = TRUE, testPostgreSQL = TRUE)' to do so\n")
		return()
	}
		require("RObsDat")
		on.exit( detach("package:RPostgreSQL", unload=TRUE))
		require("RPostgreSQL")
		postDriver <- dbDriver("PostgreSQL")
		con <- dbConnect(postDriver, user="reusser", password="4nsp", dbname="obsdat_test", port="5433", host="localhost")
		sqhandler <-  new("odm1_1Ver", con=con)
		options(odm.handler=sqhandler)
		on.exit({
			cat("Exiting and cleaning PostgreSQL database\n")
			dbGetQuery(con, "DROP SCHEMA public CASCADE;")
			dbGetQuery(con, "CREATE SCHEMA public ;")
			odm.close()
			rm(postDriver, con)
			detach("package:RPostgreSQL", unload=TRUE)
			cat("Exiting and cleaning finished\n")
			})

	        try(getMetadata("Site"), silent=TRUE)
		longExample()



		})
test_that("MySQL works with LongExampleData", {
	#Setup
	if(!(as.logical(getOption("testMySQL", FALSE)))){
		cat("Not testing long example on MySQL. Use 'options(testLongExample = TRUE, testMySQL = TRUE)' to do so\n")
		return()
	}
		require("RObsDat")
		on.exit( detach("package:RMySQL", unload=TRUE))
		require("RMySQL")
		sqDriver <- dbDriver("MySQL")
		con <- dbConnect(sqDriver, user="a_user", password="secret", dbname="obsdat_test")
		sqhandler <-  new("odm1_1Ver", con=con)
		options(odm.handler=sqhandler)
		on.exit({
			cat("Cleaning MySQL database\n")
			cleanupMySQL(con)
			odm.close()
			rm(sqDriver, con)
			detach("package:RMySQL", unload=TRUE)
			cat("Exiting and cleaning finished\n")
			}
				)


	        try(getMetadata("Site"), silent=TRUE)
		longExample()




		})

