get_binaryIDs<-function(path, net = 1){
	# connect to binaryID.db to import vector of binaryIDs
	driver<-dbDriver("SQLite")
  # open a connection to the binary ID table
	connect<-try(dbConnect(SQLite(), paste(path, "/binaryID.db", sep = "")), silent = TRUE)
	if(class(connect) == "try-error"){
	  dbUnloadDriver(driver)
	  stop("Unable to open a connection to the binary table - directory could be incorrect")
	}
	# List all database tables - there should be one for each network
	networks<-dbListTables(connect)
	IDTable <- try(dbGetQuery(connect, statement = paste("select * from net", net, sep = "")),
                 silent = TRUE)
	if(class(IDTable) == "try-error"){
	  dbDisconnect(connect)
	  dbUnloadDriver(driver)
    nnets<-length(networks)
    stop(paste("Error finding the network number specified, number should be between 1 and ", nnets, sep = ""))
	}
	rids<-IDTable[,1]
	IDTable<-IDTable[order(rids),]
	rids.sort<-rids[order(rids)]

	### Must close SQLite Driver
  dbDisconnect(connect)
  dbUnloadDriver(driver)
	IDTable
	}
