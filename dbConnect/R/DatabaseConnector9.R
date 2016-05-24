
#' Convert data types from SQL to R
#'
#' Takes a vector of SQL data types and and returns the corresponding
#' data types in R.
#'
#' @param data The vector of SQL data types to convert
#' @return A vector of R data types.
#' @author Dason Kurkiewicz \email{dasonk@@iastate.edu}
sqlToR <- function(data){
	data[data %in% c("varchar", "char", "date", "text")] <- "character"
	data[data %in% c("int", "decimal", "bigint", "double")] <- "numeric"
	return(data)
}

#' Convert times to a more useful form
#'
#' Converts times to a more useful form
#'
#' @param time A single time
#' @param type Either 'second' or 'millisecond' specifying whether the data is seconds
#'   or milliseconds from the epoch
#' @return A more useful time
unix2POSIXct <- function(time, type = "second"){
	if(type == "second"){
		return(structure(time, class = c("POSIXt", "POSIXct")))
	}else if(type == "millisecond"){
		return(structure(time/1000, class = c("POSIXt", "POSIXct")))
	}else{
		#Add error message
	}
}


#~ insertFunction <- function(variables, func){
	#~ output <- paste(func,"(", variables[1],")",sep="")
	#~ for(i in variables[-1]){
		#~ output <- paste(output,", ",func,"(",i,")",sep="")
	#~ }
	#~ return(output)
#~ }




#' Connect to a database and view it
#'
#' Provides a connection dialog and a viewer to explore a MySQL database
#'
#' @param toolkit Which gWidgets toolkit to use when creating the GUI
#'
#' @author Dason Kurkiewicz \email{dasonk@@iastate.edu}
DatabaseConnect <- function(toolkit = "RGtk2"){
.dbhistory <- ""

dlgConnect <- function(toolkit = "RGtk2") {
	library(gWidgets)
	options("guiToolkit"=toolkit)
	#library(gWidgetstcltk)

	# Save options so we can reset them later
	op.p <- options()$prompt
	op.c <- options()$continue
	options(prompt=" ")
	options(continue=" ")

	win = gwindow("DB Connection Dialog")
	g = ggroup(horizontal = FALSE, cont=win)
	addSpring(g)

	connectionGroup = ggroup(vertical=TRUE, horizontal=FALSE, cont=g)
	serverGroup = ggroup(horizontal=TRUE, cont=connectionGroup)
	addSpring(serverGroup)
	glabel(text="Server Host: ", cont=serverGroup)
	server = gedit(text="headnode.stat.iastate.edu", editable=T, cont=serverGroup, width=16)
	addSpace(serverGroup, 20)
	glabel(text="Port: ", cont=serverGroup)
	port = gedit(text="3306", editable=T, cont=serverGroup, width=8)

##########
	userGroup = ggroup(horizontal=TRUE, cont=connectionGroup)
	addSpring(userGroup)
	glabel(text="Username: ", cont=userGroup)
	user = gedit(text="2009Expo", editable=T, width=34, cont=userGroup)

##########
	pwdGroup = ggroup(horizontal=TRUE, cont=connectionGroup)
	addSpring(pwdGroup)
	glabel(text="Password: ", cont=pwdGroup)
	pwd = gedit(text="R R0cks", editable=T, width=34, cont=pwdGroup)

##########
	dbGroup = ggroup(horizontal=TRUE, cont=connectionGroup)
	addSpring(dbGroup)
	glabel(text="Database: ", cont=dbGroup)
	dbname = gedit(text="data_expo_2009", editable=T, width=34, cont=dbGroup)
	size(connectionGroup) <- c(200,125)

##########
	buttonGroup = ggroup(horizontal = TRUE, cont=g)
	#helpButton = gbutton("help", cont=buttonGroup)
	addSpring(buttonGroup)
	m <- NULL
	con <- NULL
	okButton = gbutton("Connect", cont=buttonGroup,#icon = getStockIcons()$ok,
	    handler=function(h,...) {
		library(RMySQL)
		    .dbhistory <<- paste(.dbhistory, "library(RMySQL)", sep = "\n")
		    drv <- dbDriver("MySQL")
		    .dbhistory <<- paste(.dbhistory, "drv <- dbDriver(\"MySQL\")", sep = "\n")
		    con <<- dbConnect(
		        drv,
		        user = svalue(user),
	            password = svalue(pwd),
	            port = as.numeric(svalue(port)),
	            dbname = svalue(dbname),
	            host = svalue(server)
	            ) 
		    .dbhistory <<- paste(.dbhistory, "\ncon <- dbConnect(drv, user = \"",
		       svalue(user), "\", password = \"",svalue(pwd), "\", port = ",
		       as.numeric(svalue(port)),", dbname = \"", svalue(dbname),
		       "\", host = \"",svalue(server),"\")", sep = ""
		    )

		svalue(okButton) <- "Connected"
		    dispose(win)
		    databaseViewer(con)
		}
	)

	addSpace(buttonGroup, 10)
	cancelButton = gbutton("Close",
	    cont=buttonGroup,
		handler=function(h,...){dispose(win)}
	)
	
	#Restoring options	
	#options(optionsSave)	
	options(prompt = op.p)
	options(continue = op.c)
	
} # End of dlgConnect


databaseViewer <- function(con){
	#Store options and then change them.  Reset later.
	op.p <- options()$prompt
	op.c <- options()$continue
	options(prompt=" ")
	options(continue=" ")
	
	#Getting wanted values
	Tables = list()
	Tables$names <- dbListTables(con)
	Tables$count <- integer(length(Tables$names))
	names(Tables$count) <- Tables$names
	Tables$totalNumber <- length(Tables$names)
	Fields <- list()
	for(i in Tables$names){
		Fields[[i]] <- dbListFields(con, i)
		Tables$count[i] <- as.numeric(dbGetQuery(con, paste("SELECT COUNT(*) FROM",i)))
	}
	
	#Creating Main window

	saveHist <- function(dbhist){
		win <- gwindow("Save History?")
		gg <- ggroup(hor = F, cont = win)
		g1 <- ggroup(hor = T, cont = gg)
		glabel("Save history to file?", cont = g1)
		g2 <- ggroup(hor = T, cont = gg)
		dirlab <- glabel("Save in directory:", cont = g2)
		dirline <- gedit(paste(getwd()), expand = T, cont = g2)
		g3 <- ggroup(hor = T, cont = gg)
		dirlab <- glabel("Save as:", cont = g3)
		saveline <- gedit("dbhistory.R", expand = T, cont = g3)
		g4 <- ggroup(hor = T, cont = gg)
		savebutton <- gbutton("Save", handler = function(h, ...){
	  	    cat(dbhist,
			  file = paste(svalue(dirline),"/",svalue(saveline), sep = "")
		    )
		    dispose(win)
		}, cont = g4)
		cancelbutton <- gbutton("Cancel", handler = function(h,...){dispose(win)}, cont = g4)
	}
        
	destroyHandler <- function(h,...){
		options(prompt = op.p)
		options(continue = op.c)
		saveHist(.dbhistory)
		dbDisconnect(con)
	}
	mainWindow <- gwindow("Main Browser", handler = destroyHandler)
	mainGroup <- ggroup(hor=F, cont=mainWindow, expand = T)
	mainNotebook <- gnotebook(cont = mainGroup, expand = T)
	
	.db <- list(tables = Tables,
	            fields = Fields,
	            main = list(
	                   group = mainGroup,
	                   notebook = mainNotebook,
	                   nbList = list(
					            line = list(),
	                            variable = list(),
	                            subsample = list(),
	                            query = list()
						)
				)
	)
	
	## Create main tabs
	.db$main$nbList$line$group <- ggroup(
	       hor = F, cont = .db$main$notebook, label = "Data Viewer", expand = T
	)

	.db$main$nbList$variable$group <- ggroup(
	       hor = F, cont = .db$main$notebook, label = "Variable Browser" , expand = T
	)

	.db$main$nbList$subsample$group <- ggroup(
	       hor = T, cont = .db$main$notebook, label = "Subsample", expand = T
	)

	.db$main$nbList$query$group <- ggroup(
	       hor = F, cont = .db$main$notebook, label = "Query", expand = T
	)

	svalue(.db$main$notebook) <- 1


	#######################################
	#### Make Line Viewer Tab
	#######################################
	

	#Create Notebook for the line viewer
	.db$main$nbList$line$notebook <- gnotebook(cont = .db$main$nbList$line$group, expand = T)

	.db$main$nbList$line$nbList <- list()

	.dbhistory <<- paste(.dbhistory, "\n\n# Preview of first 5 lines of each nonempty table", sep = "")
	
	for(i in .db$tables$names){
		if(.db$tables$count[i] != 0){
			.db$main$nbList$line$nbList[[i]] <- list()
			.db$main$nbList$line$nbList[[i]]$group <- ggroup(hor = F,
			    cont = .db$main$nbList$line$notebook,
			    label = i,
			    expand = T
			)
			tablename <- paste("Table Name: ",
			    i,
				sep = ""
			)
			tablerows <- paste(
			    "Number of Records: ",
			    .db$tables$count[i],
			    sep = ""
			)
			.db$main$nbList$line$nbList[[i]]$infogroup <- ggroup(hor = T,
			    cont = .db$main$nbList$line$nbList[[i]]$group,
			    expand = F
			)
			.db$main$nbList$line$nbList[[i]]$tablename <- glabel(
			    tablename,
			    cont = .db$main$nbList$line$nbList[[i]]$infogroup
			)
			addSpring(.db$main$nbList$line$nbList[[i]]$infogroup)
			
			.db$main$nbList$line$nbList[[i]]$tablerows <- glabel(
			    tablerows,
			    cont = .db$main$nbList$line$nbList[[i]]$infogroup
			)

			.dbhistory <<- paste(.dbhistory,"\n#dbGetQuery(con, \"SELECT * FROM ",i," LIMIT 5\")", sep = "")
			
		}else{
			warn.message <- paste("Failed to create Line Viewer for:", i,"\nReason: Empty Table\n")
			warning(warn.message, immediate. = TRUE, call. = FALSE)
		}
	}
	.dbhistory <<- paste(.dbhistory, "\n", sep = "")
	
	.db$tables$usable <- names(.db$main$nbList$line$nbList)
	.db$main$nbList$line$activated <- rep(FALSE, length(.db$tables$usable))
	names(.db$main$nbList$line$activated) <- .db$tables$usable

	#Add handler to the line viewer
	lineHandler <- function(h, ...){
		j <- h$pageno
		i <- names(.db$main$nbList$line$nbList)[j]

		# If this tab isn't activated yet then create it.
		if(.db$main$nbList$line$activated[i] == FALSE){
			lineViewerQuery <- paste("SELECT * FROM",i,"LIMIT 5")
			lineTable <- dbGetQuery(con, lineViewerQuery)
			.db$main$nbList$line$nbList[[i]]$table <<- gtable(lineTable, label = i)
				#container = .db$main$nbList$line$nbList[[i]]$group,
			        #label = i,
			        #expand = T
		        #)
			#Yay this finally worked
			add(.db$main$nbList$line$nbList[[i]]$group,
			    .db$main$nbList$line$nbList[[i]]$table,
			     expand = T
			)
			.db$main$nbList$line$activated[i] <<- TRUE
		}
	}
	addHandlerChanged(.db$main$nbList$line$notebook, handler = lineHandler)	

	#Change to the first tab and make a note that it gets activated.
	if(.db$tables$totalNumber == 1){
		h <- list(pageno = 1)
		lineHandler(h)
	}	
	svalue(.db$main$nbList$line$notebook) <- 1
	.db$main$nbList$line$activated[1] <- T


	###############################################
	#### Make Variable Browser Tab
	###############################################

	.db$main$nbList$variable$notebook <- gnotebook(
	    cont = .db$main$nbList$variable$group,
	    expand = T
	)
	.db$main$nbList$variable$nbList <- list()

	n <- .db$tables$usable
	.db$main$nbList$variable$activated <- rep(FALSE, length(n))
	names(.db$main$nbList$variable$activated) <- n
	
	.dbhistory <<- paste(.dbhistory, "\n\n# Preview of variable types for each table", sep = "")
	for(i in n){
		.db$main$nbList$variable$nbList[[i]] <- list()
		.db$main$nbList$variable$nbList[[i]]$group <- ggroup(hor = F,
		    cont = .db$main$nbList$variable$notebook,
		    label = i,
		    expand = T
		)
		tablename <- paste("Table Name: ",
		    i,
			sep = ""
		)
		tablerows <- paste(
			"Number of Records: ",
			.db$tables$count[i],
			sep = ""
		)
		.db$main$nbList$variable$nbList[[i]]$infogroup <- ggroup(hor = T,
			cont = .db$main$nbList$variable$nbList[[i]]$group,
			expand = F
		)
		.db$main$nbList$variable$nbList[[i]]$tablename <- glabel(
			tablename,
			cont = .db$main$nbList$variable$nbList[[i]]$infogroup
		)
		addSpring(.db$main$nbList$variable$nbList[[i]]$infogroup)
			
		.db$main$nbList$variable$nbList[[i]]$tablerows <- glabel(
			tablerows,
			cont = .db$main$nbList$variable$nbList[[i]]$infogroup
		)		
		
		q <- paste("SELECT COLUMN_NAME, DATA_TYPE ",
		    "FROM INFORMATION_SCHEMA.COLUMNS ",
		    "WHERE table_name = '",i,"'",
		    sep = ""
		)
		
		.dbhistory <<- paste(.dbhistory, "\n#dbGetQuery(con, \"", q, "\")", sep = "")
	}

	variableHandler <- function(h, ...){
		j <- h$pageno
		i <- .db$tables$usable[j]
		if(.db$main$nbList$variable$activated[i] == FALSE){
		
			# This is just the code from the line viewer duplicated
			#Must be a better way but this is easier for now.
			
		
			#Works fine but let's try something else
			#q <- paste("DESCRIBE",i)
			
			# q <- paste("SELECT COLUMN_NAME, DATA_TYPE, IS_NULLABLE, COLUMN_DEFAULT ",
			    # "FROM INFORMATION_SCHEMA.COLUMNS ",
			    # "WHERE table_name = '",i,"'",
				# sep = ""
			# )
			q <- paste("SELECT COLUMN_NAME, DATA_TYPE ",
			    "FROM INFORMATION_SCHEMA.COLUMNS ",
			    "WHERE table_name = '",i,"'",
				sep = ""
			)			
			
			description <- dbGetQuery(con, q)
			R_Data_Type = sqlToR(description[,"DATA_TYPE"])
			
			#q <- paste("SELECT",insertFunction(description$COLUMN_NAME, "min"),"FROM",i)
			#mins <- t(dbGetQuery(.con,q))
			
			#q <- paste("SELECT",insertFunction(description$COLUMN_NAME, "max"),"FROM",i)
			#maxs <- t(dbGetQuery(.con,q))

			dat <- data.frame(description,
			    R_Data_Type = R_Data_Type)#,
				# MIN = as.character(mins),
			    # MAX = as.character(maxs)
			# )
			
			.db$main$nbList$variable$nbList[[i]]$table <<- gtable(dat,
			    label = i,
			    expand = T
			)

			add(.db$main$nbList$variable$nbList[[i]]$group,
			    .db$main$nbList$variable$nbList[[i]]$table,
			    expand = T
			)
			.db$main$nbList$variable$activated[i] <<- TRUE
		}
	}
	addHandlerChanged(.db$main$nbList$variable$notebook, handler = variableHandler)
	.dbhistory <<- paste(.dbhistory, "\n", sep = "")
	
	if(.db$tables$totalNumber == 1){
		h <- list(pageno = 1)
		variableHandler(h)
	}	
	svalue(.db$main$nbList$variable$notebook) <- 1
	.db$main$nbList$variable$activated[1] <- T
	
	#~ #################################################
	#~ #### Make SubSample Tab
	#~ #################################################

	
	.db$noLimit <- ""
	.db$initialValue <- ""
	#.db$compareOptions <-  c("=","<>",">","<",">=","<=","BETWEEN","LIKE","IN")
	.db$compareOptions <- data.frame(
		options = c("=","<>",">","<",">=","<=","BETWEEN","LIKE","IN"),
		icon = rep("blank", 9),
		tooltip = c("is equal to", "is not equal to", "is greater than", "is less than", "is greater than or equal to",
		                 "is less than or equal to", "is between the following values", "is like", "is in the following range")
		)
	.db$optionComboBoxWidth <- 80

	#Create Tables frame and radio buttons for usable tables
	if(.db$tables$totalNumber > 1){
		.db$main$nbList$subsample$tableRadio <- gradio(.db$tables$usable, hor = F)
		.db$main$nbList$subsample$tableRadioFrame <- gframe("Tables",
		hor = F,
		expand = F,
		cont = .db$main$nbList$subsample$group
		)
		add(.db$main$nbList$subsample$tableRadioFrame, .db$main$nbList$subsample$tableRadio)
	}
	.db$main$nbList$subsample$right <- list()

	#Overall group on the right
	.db$main$nbList$subsample$right$group <- ggroup(hor = F, 
	    expand = F, 
	    cont = .db$main$nbList$subsample$group
	)
			
	#Group for everything except the last line
	.db$main$nbList$subsample$right$groupTop <- ggroup(hor = F,
	    expand = F,
	    cont = .db$main$nbList$subsample$right$group
	)
			
	#Group for the last line
	.db$main$nbList$subsample$right$groupBottom <- ggroup(hor = T,
	    expand = F,
	    cont = .db$main$nbList$subsample$right$group
	)
			
	.db$main$nbList$subsample$right$rows <- list()
	
	# Set up Top Row with extra options
	.db$main$nbList$subsample$right$top <- list()
	.db$main$nbList$subsample$right$top$group <- ggroup(hor = T,
	    expand = F, 
	    cont = .db$main$nbList$subsample$right$groupTop
	)

	# Set up bottom row
	.db$main$nbList$subsample$right$top$addButton <- gbutton("Add Row",
	    cont = .db$main$nbList$subsample$right$groupBottom
	)
	.db$main$nbList$subsample$right$top$removeButton <- gbutton("Remove Row",
	    cont = .db$main$nbList$subsample$right$groupBottom
	)
	.db$main$nbList$subsample$right$top$resetButton <- gbutton("Reset to defaults",
	    cont = .db$main$nbList$subsample$right$groupBottom
	)
	addSpring(.db$main$nbList$subsample$right$groupBottom)

	.db$main$nbList$subsample$right$top$limitGroup <- ggroup(hor = T,
	    cont = .db$main$nbList$subsample$right$groupBottom
	)	
	.db$main$nbList$subsample$right$top$limitLabel <- glabel("Limit:",
	    cont = .db$main$nbList$subsample$right$top$limitGroup
	)
	.db$main$nbList$subsample$right$top$limitEdit <- gedit(.db$noLimit,
	    cont = .db$main$nbList$subsample$right$top$limitGroup
	)
	.db$main$nbList$subsample$right$top$helpButton <- gbutton("Help",
	    cont = .db$main$nbList$subsample$right$groupBottom
	)

	# Set up first real row
	.db$main$nbList$subsample$right$rows[[1]] <- list()
	.db$main$nbList$subsample$right$rows[[1]]$group <- ggroup(
	    hor = T,
	    expand = F,
	    cont = .db$main$nbList$subsample$right$groupTop
	)
	
	.db$main$nbList$subsample$right$top$executeButton <- gbutton("Execute Query",
	    cont = .db$main$nbList$subsample$right$rows[[1]]$group
	)
	addSpring(.db$main$nbList$subsample$right$rows[[1]]$group)		
	
	.db$main$nbList$subsample$right$rows[[1]]$variableComboBox <- gcombobox(
	    .db$fields[[.db$tables$usable[1] ]],
	    cont = .db$main$nbList$subsample$right$rows[[1]]$group
	)
	.db$main$nbList$subsample$right$rows[[1]]$optionComboBox <- gcombobox(.db$compareOptions,
	    cont = .db$main$nbList$subsample$right$rows[[1]]$group,
	    width = .db$optionComboBoxWidth
	)
	.db$main$nbList$subsample$right$rows[[1]]$valueEdit <- gedit(
	    .db$initialValue,
	    cont = .db$main$nbList$subsample$right$rows[[1]]$group,
		expand = T
	)

	#Add handler for the 'Add Row Button'
	addRowHandler <- function(h, ...){
		currentTable <- .db$tables$names[1]
		if(.db$tables$totalNumber > 1){
			currentTable <- svalue(.db$main$nbList$subsample$tableRadio)
		}
		i <- length(.db$main$nbList$subsample$right$rows) + 1
		.db$main$nbList$subsample$right$rows[[i]] <<- list()
		
		#Create group
		.db$main$nbList$subsample$right$rows[[i]]$group <<- ggroup(
		    hor = T,
		    expand = T,
		)

		#Create radio buttons
		.db$main$nbList$subsample$right$rows[[i]]$andor <<- gradio(
		    c("AND","OR"),
		    hor = T,
		    cont = .db$main$nbList$subsample$right$rows[[i]]$group
		)
		addSpring(.db$main$nbList$subsample$right$rows[[i]]$group)
		
		#Create combo box for variables
		.db$main$nbList$subsample$right$rows[[i]]$variableComboBox <<- gcombobox(
		    .db$fields[[currentTable ]],
		    cont = .db$main$nbList$subsample$right$rows[[i]]$group
		)
		
		#Create combo box for options
		#width of 80 is required or else the Windows version sucksssss
		.db$main$nbList$subsample$right$rows[[i]]$optionComboBox <<- gcombobox(
		    .db$compareOptions,
		    cont = .db$main$nbList$subsample$right$rows[[i]]$group,
		    width = .db$optionComboBoxWidth
		)
		
		#Create editable line for values
		.db$main$nbList$subsample$right$rows[[i]]$valueEdit <<- gedit(
		    .db$initialValue,
		    cont = .db$main$nbList$subsample$right$rows[[i]]$group,
			expand = T
		)
		add(.db$main$nbList$subsample$right$groupTop,
		    .db$main$nbList$subsample$right$rows[[i]]$group
		)
	} # End of addRowHandler
	addHandlerClicked(.db$main$nbList$subsample$right$top$addButton, handler = addRowHandler)
	
	removeRowHandler <- function(h, ...){
		i <- length(.db$main$nbList$subsample$right$rows)
		if(i <=1){
			warning("Not allowed to remove ALL of the records", immediate. = TRUE, call. = FALSE)
		}else{
			delete(.db$main$nbList$subsample$right$groupTop,
			       .db$main$nbList$subsample$right$rows[[i]]$group
			)
			.db$main$nbList$subsample$right$rows <<- .db$main$nbList$subsample$right$rows[-i]
		}
	}
	addHandlerClicked(.db$main$nbList$subsample$right$top$removeButton, handler = removeRowHandler)

	resetHandler <- function(h, ...){
		i <- length(.db$main$nbList$subsample$right$rows)
		if(i > 1){
			for(j in 2:i){
				removeRowHandler(h,...)
			}		
		}
		svalue(.db$main$nbList$subsample$right$rows[[1]]$valueEdit) <- ""
		svalue(.db$main$nbList$subsample$right$rows[[1]]$optionComboBox, index = T) <- 1
		svalue(.db$main$nbList$subsample$right$rows[[1]]$variableComboBox, index = T) <- 1
		svalue(.db$main$nbList$subsample$right$top$limitEdit) <- .db$noLimit
	}
	addHandlerClicked(.db$main$nbList$subsample$right$top$resetButton, handler = resetHandler)

	helpDialog <- function(message, wrap = 150){
		window <- gwindow("Help!")
		group <- ggroup(cont = window,hor = F)
		top.group <- ggroup(hor = F, expand = T, cont = group)
		gimage("home",dirname = "stock",size = "dialog", container = top.group)	
		gtext(#message,
		    strwrap(message, wrap),
		    cont = top.group,
		    expand = T, editable = F
		)
		bottom.group <- ggroup(hor = T, expand = F, cont = group)
		gbutton("ok",
		    expand = T,
		    handler = function(h,...){dispose(window)},
		    cont = bottom.group
		)
	}
	
	
	helpHandler <- function(h, ...){

		subsamplehelp <-"
'=' : Checks for equality\n
Ex:  'Year = 1986' specifies only the data from 1986
\n---\n
'<>' : Checks for inequality\n
Ex:  'Year <> 1986' specifies any information except data from 1986
\n---\n
'>' : Checks for greater than\n
Ex:  'Year > 1986' specifies data from 1987 to the present (or future)
\n---\n
'<' : Checks for less than\n
Ex:  'Year < 1986' species data from 1985 or earlier
\n---\n
'>=' : Checks for greater than or equal to\n
Ex:  'Year >= 1986' specifies data from 1986 to the present (or future)
\n---\n
'<=' : Checks for less than or equal to\n
Ex:  'Year <= 1986' species data from 1986 or earlier
\n---\n
'BETWEEN' : Checks for data between two values\n
Ex: 'Year BETWEEN 2000 and 2002' specifies data in the years 2000, 2001, and 2002\n
Ex: Which is equivalent to 'Year >= 2000 AND Year <= 2002'
\n---\n
'LIKE' : Used to match with wildcards, '%' matchs any # of characters (even 0), '_' matches exactly 1\n
Ex: 'Year LIKE 200_' will match 2000, 2001, ..., 2009\n
Ex: 'Year LIKE 20%' will match 20, 200, 201, ..., 209, 2000, 2001, ...,
2010, ..., 2098, 2099, ..., 20000,...
\n---\n
'IN : Used to match with a specific set of values\n
Ex: 'Year IN (2001, 2003, 2007)' will match only years 2001, 2003, 2007
"

		helpDialog(subsamplehelp)
	}
	addHandlerClicked(
	    .db$main$nbList$subsample$right$top$helpButton,
	    handler = helpHandler
	)
	if(.db$tables$totalNumber > 1){
		tableRadioHandler <- function(h, ...){
			newTable <- svalue(h$obj)
			newVariables <- .db$fields[[newTable]]
			#For however many rows there are change the variable combo boxes
			for(i in 1:length(.db$main$nbList$subsample$right$rows)){
				.db$main$nbList$subsample$right$rows[[i]]$variableComboBox[] <<- newVariables
				svalue(.db$main$nbList$subsample$right$rows[[i]]$variableComboBox, index = T) <<- 1
			}
		}
		addHandlerChanged(.db$main$nbList$subsample$tableRadio, tableRadioHandler)
	}
	
	#Add the execute handler later so that we can send to the query line

	################################################
	#### Make Query Tab
	################################################
	
	# Query Line
	.db$main$nbList$query$queryLine <- list()
	.db$main$nbList$query$queryLine$group <- ggroup(
	    hor = T,
		cont = .db$main$nbList$query$group
	)
	.db$main$nbList$query$queryLine$leftGroup <- ggroup(
	    hor = T,
	    cont = .db$main$nbList$query$queryLine$group
	)
	.db$main$nbList$query$queryLine$button <- gbutton(
	    "Run following query:",
	    cont = .db$main$nbList$query$queryLine$leftGroup,
	)

	 starttext <- paste("SELECT * FROM",.db$tables$usable[1],"LIMIT 10")
	.db$main$nbList$query$queryLine$line <- gedit(
	    starttext,
	    cont = .db$main$nbList$query$queryLine$group,
	    expand = T
	)
	
	# Tab name Line
	.db$main$nbList$query$tabLine$group <- ggroup(hor = T,
		cont = .db$main$nbList$query$group
	)
	.db$main$nbList$query$tabLine$leftGroup <- ggroup(
	    hor = T,
	    cont = .db$main$nbList$query$tabLine$group
	)
	.db$main$nbList$query$tabLine$label <- glabel(
	    "Tab Name: ",
		#cont = .db$main$nbList$query$tabLine$group,
		cont = .db$main$nbList$query$tabLine$leftGroup
	)

	.db$main$nbList$query$tabLine$line <- gedit("Tab",
	    cont = .db$main$nbList$query$tabLine$group,
	    expand = T
	)	
	
	# Save in R Line
	.db$main$nbList$query$saveLine$group <- ggroup(hor = T,
		cont = .db$main$nbList$query$group
	)
	.db$main$nbList$query$saveLine$leftGroup <- ggroup(hor = T,
		cont = .db$main$nbList$query$saveLine$group
	)
	.db$main$nbList$query$saveLine$checkbox <- gcheckbox(
	    "Save in R?  ",
		cont = .db$main$nbList$query$saveLine$leftGroup,
	)
	.db$main$nbList$query$saveLine$label <- glabel("Save as:")
	.db$main$nbList$query$saveLine$line <- gedit("", expand = T)


	saveCheckBoxHandler <- function(h, ...){
		save <- svalue(h$obj)
		if(save){
			add(.db$main$nbList$query$saveLine$group, .db$main$nbList$query$saveLine$label)
			add(.db$main$nbList$query$saveLine$group, .db$main$nbList$query$saveLine$line, expand = T)
		}else{
			delete(.db$main$nbList$query$saveLine$group, .db$main$nbList$query$saveLine$label)
			delete(.db$main$nbList$query$saveLine$group, .db$main$nbList$query$saveLine$line)
		}
	}
	addHandlerChanged(.db$main$nbList$query$saveLine$checkbox, saveCheckBoxHandler)
	
	
	#make notebook
	.db$main$nbList$query$notebook <- gnotebook(cont = .db$main$nbList$query$group, expand = T)
	notebookHandler <- function(h, ...){print(h)}
	addHandlerRightclick(.db$main$nbList$query$notebook, notebookHandler)
	
	queryHandler <- function(h, ...){
		query <- svalue(.db$main$nbList$query$queryLine$line)
		tabname <- svalue(.db$main$nbList$query$tabLine$line)
		output <- dbGetQuery(con, query)
		
		
		numberofrows <- dim(output)[1]
		
		#Add column for ids
		output <- cbind(1:numberofrows, output)
		names(output)[1] <- "ID"
		gtable(as.data.frame(output),
		    cont = .db$main$nbList$query$notebook,
			label = tabname
		)
		
		if(svalue(.db$main$nbList$query$saveLine$checkbox)){
			saveName <- svalue(.db$main$nbList$query$saveLine$line)
			tmpSave <- gsub(" ", "", saveName)
			if(tmpSave != saveName){
				print("Extra Spaces were removed.  Results saved as")
				print(tmpSave)
			}else if(tmpSave == ""){
				print("No name given: Saving to .tmp")
				print("Please rename to something else")
				tmpSave <- ".tmp"
			}
			temp <- paste(tmpSave,"<<- output")
			.dbhistory <<- paste(.dbhistory, "\n", tmpSave, " <- dbGetQuery(con, \"",query, "\")", sep = "")
			eval(parse(text=temp))
		}else{
			tmpSave <- "temp"
			.dbhistory <<- paste(.dbhistory, "\n", tmpSave, " <- dbGetQuery(con, \"",query, "\")", sep = "")	
		}
	}
	addHandlerClicked(.db$main$nbList$query$queryLine$button, queryHandler)
	

	
	# corresponds to the execute query in the subsample tab
	executeHandler <- function(h, ...){
		table <- .db$tables$names[1]
		if(.db$tables$totalNumber > 1){
			table <- svalue(.db$main$nbList$subsample$tableRadio)
		}
		pageno <- which(.db$tables$usable == table)
		variableHandler(list(pageno = pageno))
		tablevariables <- .db$main$nbList$variable$nbList[[pageno]]$table[]
		
		n <- length(.db$main$nbList$subsample$right$rows)
		out <- paste("SELECT * FROM", table, "WHERE")
		line <- ""
		for(i in 1:n){		
			#Need to reset added and andor each time through loop
			added <- ""
			andor <- ""
			value <- svalue(.db$main$nbList$subsample$right$rows[[i]]$valueEdit)
			variable <- svalue(.db$main$nbList$subsample$right$rows[[i]]$variableComboBox)
			typ <- tablevariables[tablevariables$COLUMN_NAME == variable, "R_Data_Type"]
			
			# Kludge solution to adding/removing quotes that they
			# either need or don't need
			quoteifcharacter <- "'"
			if(typ == "numeric"){
				# if it's numeric we don't want to put a quote
				quoteifcharacter <- ""
			}
			if(substr(value,1,1) %in% c("\"", "'")){
				# if they put their own quotes remove them
				# and have it add them manually later
				value <- substr(value,2,nchar(value)-1)
			}
			
			#Create the string to add to what we have so far
			added <- paste(
			    "( ",
				variable,
				" ",
				svalue(.db$main$nbList$subsample$right$rows[[i]]$optionComboBox),
				" ",
				quoteifcharacter,
				value,
				quoteifcharacter,
				" )",
				sep = ""
			)
			
			# We don't want to add an and/or on the first time
			if(i != 1){
				andor <- svalue(.db$main$nbList$subsample$right$rows[[i]]$andor)
				line <- paste("(", line, andor, added, ")")
			}else{
				line <- added
			}
		} # End of for

		# Correctly set the limit
		limittext <- svalue(.db$main$nbList$subsample$right$top$limitEdit)
		if(limittext == "" | limittext == .db$noLimit){
			limittext <- ""
		}else{
			limittext <- paste("LIMIT",limittext)
		}
		output <- paste(out,line,limittext)

		# Send output to the Query Tab
		svalue(.db$main$nbList$query$queryLine$line) <- output
		svalue(.db$main$notebook) <- which(names(.db$main$notebook) == "Query")
	} # End of executeHandler
	addHandlerClicked(.db$main$nbList$subsample$right$top$executeButton, executeHandler)

	#size(.db$main$nbList$query$tabLine$leftGroup) <- size(.db$main$nbList$query$queryLine$leftGroup)
}

	dlgConnect(toolkit)
}

#TODO: Put this inside the other functions.
#~ confirmDialog <- function(message){
	#~ window <- gwindow("Confirm")
	#~ group <- ggroup(cont = window,hor = F)
	#~ gimage("home",dirname = "stock",size = "dialog", container = group)	
	#~ rgrp <- ggroup(hor = F, cont = group, expand = T)
	#~ inner.group <- ggroup(hor = F, cont = rgrp)
	#~ tmp <- ggroup(cont = inner.group, hor = F, expand = T)
	#~ glabel(message, cont = tmp, expand = T)
	
	#~ button.group <- ggroup(cont = inner.group, hor = T)

	#~ gbutton("ok",expand = T, handler = function(h,...){dispose(window)}, cont = button.group)
	#~ return()
#~ }

#DatabaseConnect()
