.onLoad <- function(libname,pkgname){

	# Environments for some reason lose their attributes
	# prior to .onload 
	class(Mappers) <- c('dataMappers',class(Mappers))

	if(.Platform$OS.type != "windows") {

		# A simple example of a read-only mapper for the UNIX command 'top'.
		newMapper(
			type="UNIX:top",
			init= function(map){
				# Nothing realy to do but tell datamap
				# the name of our objects
				install('top',map)
				return(TRUE)
			},
			get = function(x){
				con <- textConnection(system("top -b -n 1", intern=TRUE)[-1:-6])
				on.exit(close(con))
				read.table(con,header=TRUE)
			}
		)
	}

	# Again R CMD check, be quiet.
	con <- NULL
	newMapper(
		type="DBI:tables",
		init= function(map,...){

			# Various initialization code
			con <- dbConnect(...)

			# Assign persistent state into our map
			map$con <- con

			# Install symbols within our map. Note that 'install' is
			# different than assigning state. Installed symbols become
			# the objects visible to R
			install(dbListTables(con),map)

			# Returning FALSE means something has failed
			return(TRUE)
		},
		get = function(x) {
			if (dbExistsTable(con,x))
				return(dbReadTable(con,x))
			else
				return(NULL)
		},
		assign = function(x,val) dbWriteTable(con,x,val),
		finalize = function(map) dbDisconnect(con),
		delete = function(x) dbRemoveTable(con,x)
	)
}
