.packageName <- "prettyR"
createIndxFile <- 
function(HTMLfile, navfile, listfile, title = "R listing") {
	# Create a connection
	basecon <- file(HTMLfile, "w")
	
	# Create the file with framesets
	cat("<html>\n <head>\n", file = basecon)
	cat("  <title>", title, "</title>\n </head>\n", sep = "", file = basecon)
	cat(" <frameset cols=200,* border=1>\n", file = basecon)
	cat("  <frame src=\"", navfile, sep = "", file = basecon)
	cat("\" name=\"nav\" scrolling=yes>\n", file = basecon)
	cat("  <frame src=\"", listfile, sep = "", file = basecon)
	cat("\" name=\"list\" scrolling=yes>\n", file = basecon)
	cat(" </frameset>\n</html>\n", file = basecon)
	close(basecon)
}

StartNav <-
function(navcon, title = "R listing", bgcolor = "#dddddd") {
	# Beginning of the nav file
	cat("<html>\n <body bgcolor=\"", bgcolor, "\">\n", sep = "", file = navcon)
	cat("<center><h2>", title, "</h2></center>\n", sep = "", file = navcon)
}

AddNav <- 
function(navcon, Rcommand, listname) {
	# Add an entry in the nav file
	if (nzchar(Rcommand) > 20)
		Rcommand <- paste(substring(Rcommand, 1, 18), "...", sep = "")
	nametag <- make.names(paste("ni", Rcommand, sep = ""))
	cat(" <a href=\"", listname, "#", nametag, "\" target=\"list\">\n",
		sep = "", file = navcon)
	cat("  ", Rcommand, "</a>\n<br>\n", sep = "", file = navcon)
	return(nametag)
}

HTMLgraph <-
function(listfile, graphfile = NULL, type = "png", ...) {
	# Create a graphic device associated with graphfile
	# and add the corresponding entry in the listfile
	
	# Make sure to close all previously opened device(s)
	graphics.off()
	
	# Get the path from listfile (we put all figures in the same directory)
	HTMLdir <- dirname(listfile)
	
	# If graphfile == NULL or graphfile == "", compute fig[#] automatically
	if (is.null(graphfile) || graphfile == "") {
		# Got a name for the graphfile that is not used yet
		i <- 0
		repeat({
			i <- i + 1;
			graphfile <- file.path(HTMLdir, paste("fig", i, ".", type, sep = ""));
			if (!file.exists(graphfile)) break
		})
	} else { # A name is provided for this graphfile
		# Make sure it is located in the same directory as HTMLfile
		graphfile <- file.path(HTMLdir, basename(graphfile))
		# Try to guess the format of the file from its extension
		graphtype <- tolower(sub("^.*[.](.+)$", "\\1", graphfile))
		if (graphtype != graphfile) type <- graphtype
	}
		
	# Create the graph device	
	res <- switch(type,
		png = png(graphfile, ...),
		bmp = bmp(graphfile, ...),
		jpeg = jpeg(graphfile, ...),
		stop("'type' must be one of: 'png', 'bmp' or 'jpeg'"))
	
	# Create the corresponding entry in the listfile
	cat("<img src=\"", basename(graphfile), "\">\n", sep = "")
} 

R2html <- 
function(Rfile, HTMLfile, echo = TRUE, split = FALSE, browse = TRUE,
	title = "R listing", bgcolor = "#dddddd", ...) {
	
	# Rfile must exist!
	if (!file.exists(Rfile))
		stop("You must specify an existing 'Rfile'!")
	
	# Make sure that HTMLfile has .htm[l] extension
	if (length(grep("\\.[hH][tT][mM][lL]?$", HTMLfile)) == 0)
		stop("'HTMLfile' must be a filename with .htm[l] extension!")
	
	# Make sure one has the full path to HTMLfile
	if (basename(HTMLfile) == HTMLfile)
		HTMLfile <- file.path(getwd(), HTMLfile)
	
	# Create names and connections for nav and list files
	baseHTMLfile <- basename(HTMLfile)	# Only the file name witrhout path
	HTMLdir <- dirname(HTMLfile)
	basenavfile <- sub("\\.[hH][tT][mM][lL]?", "_nav.html", baseHTMLfile)
	navfile <- file.path(HTMLdir, basenavfile)
	baselistfile <- sub("\\.[hH][tT][mM][lL]?", "_list.html", baseHTMLfile)
	listfile <- file.path(HTMLdir, baselistfile)
	
	# Create the file defining frames
	createIndxFile(HTMLfile, basenavfile, baselistfile, title)
	
	# Create the nav file
	navcon <- file(navfile, "w")
	StartNav(navcon, title = title, bgcolor = bgcolor)
	
	# Create the list file and redirect output and messages to it
	listcon <- file(listfile, "w")
	StartList(listcon, title = title, bgcolor = bgcolor)
	sink(listcon, append = TRUE, type = "output", split =split)
	sink(listcon, append = TRUE, type = "message", split =split)
	
	# Instructions to execute when we exit the function
	on.exit({
		sink(NULL, type = "output");
		sink(NULL, type = "message");
		close(navcon);
		close(listcon)
	})
	
	# Read the content of the R script file
	Rcon <- file(Rfile, "r")
	Rcommands <- readLines(Rcon)
	close(Rcon)
	
	# Make sure there is no opened graphic devices, or close them
 	graphics.off()
	
	# Process commands one by one
	cmd <- ""
	for(i in 1:length(Rcommands)) {
		# Possibly paste successives lines of a multiline command
		cmd <- paste(cmd, Rcommands[i], sep = "")
		cmdcon <- textConnection(cmd)
		commexp <- try(parse(cmdcon), silent = TRUE)
		close(cmdcon)
		
		# Do we got an error?
		if (inherits(commexp, "try-error")) {
			# Is it an incomplete code line?
			if (length(grep("\n2:", commexp)) == 0) {
				# This is a syntax error
		
				# Add a tag in the navfile
				nametag <- AddNav(navcon, cmd, baselistfile)
				cat("<a name=\"", nametag, "\"></a>\n", sep = "", file = listcon)
			
				# Do we echo the command?
				if(echo) {
					cat("<pre><font color =\"Navy\">\n> ", file = listcon)
					cat(cmd, "\n", file = listcon)
					cat("</font></pre>\n", file = listcon)
				}
				
				# Print the error message
				cat("<pre><font color =\"Red\">\n", file = listcon)
				cat(commexp, file = listcon)
				cat("</font></pre>\n", file = listcon)
			
				cmd <- ""
			}
		} else if (is.expression(commexp)) {	# This code appears syntactically correct
			# Add a tag in the navfile
			nametag <- AddNav(navcon, cmd, baselistfile)
			cat("<a name=\"", nametag, "\"></a>\n", sep = "", file = listcon)
			
			# Look if there is a forbidden instruction in the code
			forbidden <- c("connection", "fifo", "file", "sink")
			dont <- FALSE
			for (i in 1:length(forbidden))
				if (length(grep(forbidden[i], cmd)) > 0) {
					dont <- TRUE
					cmd <- paste("#", cmd, sep = "")
					break
				}
			if(echo) {
				cat("<pre><font color =\"Navy\">\n> ", file = listcon)
				cat(cmd, "\n", file = listcon)
				cat("</font></pre>\n", file = listcon)
			}
			
			# Do we need to create a new HTML graph?
			# We are looking for: # --FIG[:figurename]--
			if (length(grep("#\\s*--FIG(:[^-]+)?--", cmd)) > 0) {	# Create a new HTMLgraph
				# Get the name off the figure to create
				graphfile <- sub("^.*#\\s*--FIG:?([^-]+)?--.*$", "\\1", cmd)
				HTMLgraph(listfile, graphfile)	# Create a new graph
			}
			
			# Execute the command
			if (!dont) {
				cat("<pre>\n", file = listcon)
				## PhG: I suppose that you want to evaluate this in the global environment?
				eval(parse(text = cmd), envir = .GlobalEnv)
				cat("</pre>\n", file = listcon)
			}
			cmd <- ""
		}
 	}
 	
 	# Close all currently opened graph devices
 	graphics.off()
	 
	# Finalize files
	EndHTML(navcon)
	EndHTML(listcon)
	
	# Possibly
	if (browse)
		browseURL(paste("file://", HTMLfile, sep = ""))
}
