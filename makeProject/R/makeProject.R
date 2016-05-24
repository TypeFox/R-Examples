makeProject <-
function (name = "myProject",  path = getwd(), force=FALSE, author="Your Name", email="yourfault@somewhere.net" ){
  
	# Function to create directories safely
    safe.dir.create <- function(path) {
        dirTest <- function(x) !is.na(isdir <- file.info(x)$isdir) & 
            isdir
        if (!dirTest(path) && !dir.create(path)) 
            stop(gettextf("cannot create directory '%s'", path), 
                domain = NA)
    }
    
    # Sanity check on the project name
    if (!is.character(name)){ 
        stop("Project name must be a character vector")
    }
    
    # Sanity check on path
    if (!is.character(path)){ 
        stop("Path must be a character vector")
    }
    
    # Create the directories
    dir <- file.path(path, name)
    if (file.exists(dir) && !force) {
        stop(gettextf("directory '%s' already exists", dir), domain = NA)
    }
    
    message("Creating Directories ...")
    safe.dir.create(dir)
    safe.dir.create(code_dir <- file.path(dir, "code"))
    safe.dir.create(docs_dir <- file.path(dir, "data"))
    
    # Create a description file
    description <- file(file.path(dir, "DESCRIPTION"), "wt")
    cat("Project: ", name, "\n", "Title: What the project does (short line)\n", 
        "Version: 1.0\n", "Date: ", format(Sys.time(), format = "%Y-%m-%d"), 
        "\n", "Author: ", author, "\n", "Maintainer: Who to complain to <", email, ">\n", 
        "Description: More about what it does (maybe more than one line)\n",
        file = description, sep = "")    
    close(description)
    
    message("Creating Code Files ...")
     # Create a main file
    main <- file(file.path(dir, "main.R"), "wt")
    cat("# Project: ", name, "\n", "# Author: ", author, "\n", "# Maintainer: Who to complain to <", email, ">\n\n", 
    	"# This is the main file for the project\n",
    	"# It should do very little except call the other files\n",
    	"\n",
    	"### Set the working directory\n",
    	"setwd(", '"', dir, '")', "\n",
    	"\n\n",
    	"### Set any global variables here\n",
    	"####################\n",
    	"\n\n\n",
    	"####################\n",
    	"\n\n",
    	"### Run the code\n",
		'source("code/load.R")', "\n",
		'source("code/clean.R")', "\n",
		'source("code/func.R")', "\n",
		'source("code/do.R")', "\n",
        file = main, sep = "")    
    close(main)

     # Create a load file
    loadfile <- file(file.path(dir, "code/load.R"), "wt")
    cat("# Project: ", name, "\n", "# Author: ", author, "\n", "# Maintainer: Who to complain to <", email, ">\n\n", 
    	"# This file loads all the libraries and data files needed \n",
    	"# Don't do any cleanup here\n",
    	"\n",
    	"### Load any needed libraries\n",
    	"#load(LibraryName)\n",
    	"\n\n",
    	"### Load in any data files\n",
    	'#read.csv("data/FileName" as.is=T)', "\n",
        file = loadfile, sep = "")    
    close(loadfile)    

     # Create a clean file
    cleanFile <- file(file.path(dir, "code/clean.R"), "wt")
    cat("# Project: ", name, "\n", "# Author: ", author, "\n", "# Maintainer: Who to complain to <", email, ">\n\n", 
    	"# All the potentially messy data cleanup\n",
        file = cleanFile, sep = "")    
    close(cleanFile) 
    
     # Create a func file
    funcFile <- file(file.path(dir, "code/func.R"), "wt")
    cat("# Project: ", name, "\n", "# Author: ", author, "\n", "# Maintainer: Who to complain to <", email, ">\n\n", 
    	"# Functions for the project\n",
    	"\n",
    	"myFunc <- function(){\n\n",
    	"}\n",
        file = funcFile, sep = "")    
    close(funcFile)  
    
    # Create a do file
    doFile <- file(file.path(dir, "code/do.R"), "wt")
    cat("# The actual work\n",
        file = doFile, sep = "")    
    close(doFile) 
    
    message("Complete ...")

}
