#############################################################################################
# This file contains functions that are used to parse zipped csv files from Reuters.         
# After sourcing these functions (this script),                                             
# 1st run "configureTRTH" which will create an environment to hold parameter values         
# 2nd run "download_reut" to download the big zipped csv files to your archive directory    
# 3rd run "splitCSV" that will unzip the big zipped Reuters files that were downloaded.     
#   Then it will split the file such that there will be a file for each day for each symbol 
#   and it will put those split files in your csv directory                                 
# 4th run "FEreut2xts" (requires foreach) which will read the csv files into R, 
#   do a little data scrubbing, and save the xts data into your tick directory.  
#   Then it will convert the data to 1 second frequency data and save that into your sec 
#   directory.  Also, if tick.image and/or sec.image are TRUE, it will create plots of the 
#   data and store them.             
#                                                                                           
#############################################################################################
#                  Reuters Backfill Configuration Parameters                                #
#############################################################################################
## Arguments and their defaults
# config.file (character name of config file) is optional.  
# If provided, config.file will be sourced; 
# i.e you can use it instead of specifying all the parameters in the dots. 
# any arguments provided in dots will override arguments of same name from config.file
##
#path.output = '~/TRTH/'                # base directory for output
#tick_dir = '~/TRTH/tick'               # directory in which to store tick data
#archive_dir = '~/TRTH/archive'         # directory in which to store downloaded .gz files
#csv_dir = '~/TRTH/csv'                 # directory in which to store zipped csv files
#sec_dir = '~/TRTH/sec'                 # directory in which to store second data
#sec.image = TRUE                       # save a chart of the second data?
#tick.image = TRUE                      # save a chart of the tick data?
#default_type = 'guaranteed_spread'     # passed to instrument.auto if type cannot be inferred from RIC 
#default_currency = 'USD'               # passed to instrument.auto if type cannot be inferred from RIC
#digits.sec = 6                         # for options(digits.secs=digits.secs)
#width = 200                            # for options(width=width)
#use.instrument = FALSE                 # If TRUE, non-defined instruments will be defined, and 
                                        # negative prices will be removed from non-synthetics
#instrument_file = [searches path.output for filename containing 'instruments'] #name of instrument envir RData file
#job.name = ""                          # the Reuters TRTH job name (by default all files not on disk)
#no.cores = 4                           # number of cores for foreach
#overwrite = FALSE                      # will not redownload, or overwrite files unless this is TRUE
#username = stop("")                    #TRTH user name, usually your email address
#password = stop("")                    #TRTH password
#
#email_to <- 'someuser@somehost.com'    #NOT IN USE
#email_from <- 'someuser@somehost.com'  #NOT IN USE
#############################################################################################
#               Below is how you would typically use these functions                        #
#############################################################################################
##configureTRTH('~/TRTH/TRTH_config_file.R')
## OR
#configureTRTH(
#    path.output = '~/TRTH/',
#    width = 200,
#    digits.secs = 6,
#    use.instrument = FALSE,
#    instrument_file = '~/TRTH/instruments.RData',
#    username = 'email@domain.com',
#    password = 'password',
#    default_type = 'guaranteed_spread',
#    default_currency = "USD",
#    job.name = "ReutersJobName",
#    overwrite = FALSE,
#    tick.image = TRUE,
#    sec.image = TRUE,
#    no.cores = 20
#)
#
#download_reut(.TRTH)                   # Download big zipped CSV
#system.time(splitCSV(.TRTH))           # Split into daily CSVs
#system.time(Out <- FEreut2xts(.TRTH))  # Convert to xts data: tick and second
#############################################################################################

#TODO: if user changes the value of .TRTH$path.output, do we want to change the value of .TRTH$csv_dir, xts_dir, etc?

CleanUpArchive <- function(archive_dir) {
    # If a job is killed, or it fails, there will probably be files in the
    # archive directory that should be removed.  This will delete files
    # that do not end with .csv.gz and are not confirmation files.
    # This will be called from `configureTRTH` and on.exit in `splitCSV`
    if (substr(archive_dir, nchar(archive_dir), nchar(archive_dir)) != "/") {
        archive_dir <- paste(archive_dir, "/", sep="")
    }    
    archive.files <- list.files(archive_dir)
    to.remove <- archive.files[!grepl("\\.csv\\.gz", archive.files)]
    to.remove <- to.remove[!grepl("confirmation", to.remove)]
    if (length(to.remove) > 0) warning(paste("Cleaning up archive_dir: removing", 
                                       cat(paste(to.remove, collapse="\n"), '\n')))
    paste(archive_dir, to.remove, sep="")
    unlink(to.remove, force=TRUE)
}

## Some convenience functions
addslash <- function(x) {
    if (substr(x, nchar(x), nchar(x)) != '/') paste(x, "/", sep="")
    else x   
}
makeDir <- function(x) { #if directory does not exist, create it
     dir.create(x, showWarnings=FALSE, recursive=TRUE, mode="0775") #why not mode="0664" ???
}


configureTRTH <- function(config.file, path.output='~/TRTH/', ...) {
    ## Create environment to hold variables that more than one function needs to access    
    if (!exists('.TRTH', .GlobalEnv)) {
        .TRTH <- new.env(parent=.GlobalEnv)
    } else .TRTH <- get('.TRTH', pos=.GlobalEnv)
    dargs <- list(...)

    ## Load required packages
    #require(qmao)
    require(FinancialInstrument)
    require(doMC)
    #require(sendmailR) # for email on failure

    ## Source the config_file -- this will be overwritten by any arguments in dots
    if (!missing(config.file)) source(config.file)

    # There are some variables that we need that should be in the config file.
    # Anything passed in through dots will override arguments of the same name that were in config_file
    # Some things (subdirectory names) we will create if they aren't in dots or config_file
    pickDirArg <- function(x) {
        if (!is.null(dargs[[x]])) return(dargs[[x]]) #passed through dots
        if (!is.null(.TRTH[[x]])) return(.TRTH[[x]])        
        if (exists(x)) return(addslash(get(x)))
        addslash(paste(path.output, sub("_dir", "", x), sep=""))
    }

    #if (!is.null(dargs$path.output)) 
    .TRTH$path.output <- path.output <- addslash(path.output)
    .TRTH$archive_dir <- pickDirArg("archive_dir")
    .TRTH$csv_dir <- pickDirArg("csv_dir")
    .TRTH$tick_dir <- pickDirArg("tick_dir")
    .TRTH$sec_dir <- pickDirArg("sec_dir")

    # Make sure the directories we need exist.
    makeDir(.TRTH$path.output)
    makeDir(paste(.TRTH$path.output, "tmp", sep=""))
    makeDir(.TRTH$archive_dir)
    makeDir(.TRTH$csv_dir)
    makeDir(.TRTH$tick_dir)
    makeDir(.TRTH$sec_dir)

    pickArg <- function(x, default=NULL) {
        # if argument "x" was passed through dots, use that
        # otherwise, if it was in config_file, use that
        # if it's neither in dots, nor in config_file, use default
        if (!is.null(dargs[[x]])) return(dargs[[x]]) #passed through dots
        if (!is.null(.TRTH[[x]])) return(.TRTH[[x]])
        if (exists(x)) return(get(x))
        default
    }

    ## Set some options/preferences
    .TRTH$width <- pickArg('width', 200)
    options(width=.TRTH$width)
    .TRTH$digits.sec <- pickArg('digits.sec', 6)
    options(digits.secs=.TRTH$digits.secs)

    .TRTH$username <- pickArg('username', stop("Please provide your username"))
    .TRTH$password <- pickArg('password', stop("Please provide your password"))
    .TRTH$job.name <- pickArg('job.name', "")
    .TRTH$default_type <- pickArg('default_type', 'guaranteed_spread')
    .TRTH$default_currency <- pickArg('default_currency', 'USD')
    .TRTH$overwrite <- pickArg('overwrite', FALSE)
    .TRTH$tick.image <- pickArg('tick.image', TRUE)
    .TRTH$sec.image <- pickArg('sec.image', TRUE)
    .TRTH$no.cores <- pickArg('no.cores', 4)

    .TRTH$use.instrument <- pickArg('use.instrument', FALSE) 
    # if `use.instrument` is TRUE, then the code will remove negative prices from
    # instruments that are not of type='synthetic'.  An instrument_file that contains an .instrument is required
    # and will be loaded. If the RIC cannot be found in the .instrument environment, it will be auto defined
    # using instrument.auto.
    # Using use.instrument=TRUE can take a while.

    if (isTRUE(.TRTH$use.instrument)) {
        instr.file.bak <- tail(list.files(path.output)[grep("instruments", list.files(path.output))], 1)
        .TRTH$instrument_file <- pickArg('instrument_file', instr.file.bak)
        if (length(.TRTH$instrument_file) == 0 || is.na(.TRTH$instrument_file) || !file.exists(.TRTH$instrument_file))
            stop("Please specify a valid filepath for instrument_file or move a file with 'instruments' in its name to 'path.output'")
    }

    registerDoMC(.TRTH$no.cores)
    # registerDoSEQ()

    # create a text file that contains the job.name so that if a job is running,
    # someone other than the user that started the job can find out which job 
    # is running.
    system(paste('echo "', Sys.time(), ' configureTRTH job.name: ', .TRTH$job.name, '" > ', 
                 paste(.TRTH$path.output, "current.job.txt", sep=""), 
                 sep=""))

    .TRTH$tmp <- path.expand(paste(addslash(.TRTH$path.output), "tmp", sep=""))
    makeDir(.TRTH$tmp)

    assign('.TRTH', .TRTH, pos=.GlobalEnv)

    if (Sys.getenv("TMPDIR") == "") {
        # if the TMPDIR environment variable was not set before R was started, 
        # we need to set it.  That requries restarting R.
        wd <- getwd()
        setwd(.TRTH$path.output)
        assign(".First", function() {
            #Maybe this could just recursively call configureTRTH..
            require(FinancialInstrument)
            require(doMC)
            registerDoMC(.TRTH$no.cores)
            file.remove(".RData") # already been loaded
            rm(".Last", pos=.GlobalEnv) #otherwise, won't be able to quit R without it restarting
            setwd(wd)
            .TRTH
        }, pos=.GlobalEnv)
        assign(".Last", function() {
            system(paste('TMPDIR=', .TRTH$tmp, 
                         ' R --no-site-file --no-init-file --quiet', sep=""))
        }, pos=.GlobalEnv)
        save.image() # so we can load it back when R restarts
        #TODO: make copy of .RData if it already exists instead of clobbering.
        q("no")
    }
    # If the previous if-block executed, the rest of this function will be 
    # ignored. If more code needs to be run after the restart, it should be 
    # added to the .First function above.
    .TRTH
}


download_reut <- function(.TRTH) {
    if (missing(.TRTH)) {
        .TRTH <- try(get('.TRTH', pos=.GlobalEnv))
        if (inherits(.TRTH, 'try-error')) stop("Run configureTRTH function first")
    }

    # edit text file so others can see what job we're working on
    system(paste('echo "', Sys.time(), ' download_reut job.name: ', .TRTH$job.name, '" > ', 
                 paste(.TRTH$path.output, "current.job.txt", sep=""), 
                 sep=""))

    Sys.umask("0002")

    Archive.output <- list.files(.TRTH$archive_dir)
    Archive.output <- Archive.output[grep("\\.gz",Archive.output)]
    omit <- c(grep("confirmation",Archive.output),grep("report",Archive.output))
    if (length(omit) > 0) Archive.output <- Archive.output[-omit]

    listflag=FALSE
    while(!listflag)#try to download file list 
    {
        clear <- warnings() #currency loads from oanda alway generate warnings, clear them out
        Reuters <- system(paste("curl ftp://tickhistory-ftp.thomsonreuters.com:15500/results/ -u ",
                                .TRTH$username,":",.TRTH$password," --ftp-ssl -k -l",sep=""),intern=TRUE)
        cat("\n")
        w <- ''
        w <- warnings()[!warnings() %in% clear]
        if(!as.logical(length(Reuters)) || isTRUE(grep('curl',names(w))))
        {
            tmpmsg<-paste("curl returned error code", names(w),'\n',
                        'while attempting to download file list','\n',
                        'script will wait and retry in 30 min')
            #sendmail(email_to,email_from,"error downloading Reuters file list",msg=tmpmsg)
            Sys.sleep(180)
        } else listflag=TRUE
        
    }
    # now we're past the while loop, so we have a file list
    Reuters.report <- Reuters[grep("report",Reuters)]

    Reuters.output <-  Reuters[-c(grep("report",Reuters),grep("confirmation",Reuters))]
    Reuters.output <-  Reuters.output[grep(.TRTH$job.name, Reuters.output)]

    .TRTH$files.gz <- Reuters.output[!(Reuters.output %in% Archive.output)]
    #files.gz <- paste(username, "-", job.name, ".csv.gz", sep="")
    
    if (length(.TRTH$files.gz) == 0) .TRTH$files.gz <- Reuters.output
    if (length(.TRTH$files.gz) == 0) stop('Cannot find .gz file containing "job.name" Maybe it has already been purged?')
    assign(".TRTH", .TRTH, pos=.GlobalEnv)
    
    for(i in 1:length(.TRTH$files.gz))
    {	
        filename.gz <- .TRTH$files.gz[i]
        filename.csv <- substr(filename.gz,1,(nchar(filename.gz)-3))
	
        alias <- unlist(strsplit(filename.gz,"-"))[3]
        alias <- unlist(strsplit(alias,".csv.gz"))

	    ## Download New Datasets
        print(paste("Downloading ",filename.gz,sep=""))
        fileflag=FALSE
        Reuters2 <- 0
        while(!fileflag) #try to download individual files
        {
            if (!file.exists(paste(.TRTH$archive_dir, filename.gz, sep="")) || .TRTH$overwrite) {
                Reuters2 <- system(paste("curl -m 10800 --max-filesize 1610612736 ftp://tickhistory-ftp.thomsonreuters.com:15500/results/", 
                                    filename.gz, " -u ", .TRTH$username, ":", .TRTH$password, " --ssl -k > ", .TRTH$archive_dir, filename.gz, sep=""))
            } #else cat(paste(filename.gz, 'already exists, and overwrite==FALSE; not re-downloading.'), "\n")
            if(Reuters2 != 0)
            {
                w2 <- ''
                w2 <- warnings()
                tmpmsg <- paste("curl returned error code", Reuters2,"\n",
                                w2,'\n','while attempting to download',filename.gz,'\n',
                                'will wait and retry in 10 min')
                #sendmail(email_to,email_from,paste("error downloading Reuters file",filename.gz),msg=tmpmsg)
                Sys.sleep(600)
            } else fileflag=TRUE
        }
        
	    ## Download Report s
        if (!file.exists(paste(.TRTH$archive_dir, Reuters.report[grep(alias,Reuters.report)], sep="")) || .TRTH$overwrite) { 
    	    system(paste("curl ftp://tickhistory-ftp.thomsonreuters.com:15500/results/",
                        Reuters.report[grep(alias,Reuters.report)], " -u ", .TRTH$username, ":", .TRTH$password,
                        " --ftp-ssl -k > ", .TRTH$archive_dir, Reuters.report[grep(alias,Reuters.report)], sep=""))
	        #system(paste("gzip -d -f ",archive_dir,"Report/",Reuters.report[grep(alias,Reuters.report)],sep=""))
	        cat("\n")
        } #else cat(paste(Reuters.report[grep(alias,Reuters.report)], 
          #      "already exists, and overwrite==FALSE; not re-downloading.\n"))
    }

    #save(files.gz, file=paste(.TRTH$archive_dir, 'files.gz.tmp.rda', sep=""))
    #files.gz
#    detach(.TRTH)
    assign(".TRTH", .TRTH, pos=.GlobalEnv)
    .TRTH
}


get_files.gz <- function(archive_dir, job.name){
    # Don't _really_ need this function now that .TRTH envir is being passed around
    # but might as well use it since it's already written
    if (!file.exists(archive_dir)) stop("archive_dir does not exist")

    Archive.output <- list.files(archive_dir)
    Archive.output <- Archive.output[grep("\\.gz",Archive.output)]
    omit <- c(grep("confirmation",Archive.output),grep("report",Archive.output))
    if (length(omit) > 0) Archive.output <- Archive.output[-omit]
    Reuters.output <-  Archive.output[grep(job.name, Archive.output)]
    #if (length(Reuters.output) == 0) Reuters.output <- Archive.output
    Reuters.output
}

splitCSV <- function(.TRTH) {
    if (missing(.TRTH) && !exists(".TRTH")) stop("Run configureTRTH function first")
    if (Sys.getenv("TMPDIR") == "") {
        stop(paste("TMPDIR environment variable must be set",  
                   "either manually or by running configureTRTH"))
    }
    # edit text file so others can see what job we're working on
    system(paste('echo "', Sys.time(), ' splitCSV job.name: ', .TRTH$job.name, '" > ', 
                 paste(.TRTH$path.output, "current.job.txt", sep=""), 
                 sep=""))

    # make a temp dir to use for splitting so that (fingers crossed)
    # more than one instance can be run at a time in separate R sessions.
    dir.create(.TRTH$tmp_archive_dir <- addslash(tempdir()), showWarnings=FALSE, mode='0775')
    on.exit(unlink(.TRTH$tmp_archive_dir, recursive=TRUE, force=TRUE))
        
    if (substr(.TRTH$path.output, nchar(.TRTH$path.output), nchar(.TRTH$path.output)) != "/") {
        .TRTH$path.output <- paste(.TRTH$path.output, "/", sep="")
    }

    # get files.gz if it is NULL, or if the job.name was changed
    if (is.null(.TRTH$files.gz) || identical(integer(0), grep(.TRTH$job.name, .TRTH$files.gz))) 
        .TRTH$files.gz <- get_files.gz(.TRTH$archive_dir, .TRTH$job.name)

    if (is.null(.TRTH$instrument_file)) { #Don't need this anymore
        tmp <- list.files(paste(.TRTH$path.output))
        instrument_file <- paste(.TRTH$path.output, tail(tmp[grep("instruments", tmp)], 1), sep="")
        if (!file.exists(instrument_file)) {
            stop("Could not find instrument_file; please specify")
        } else .TRTH$instrument_file <- instrument_file
    }

    if (isTRUE(.TRTH$use.instrument)) loadInstruments(.TRTH$instrument_file)
    registerDoMC(.TRTH$no.cores)

    ## unzip to tempdir and split (new unzip method does not require rezip; keeps original gz file)
    setwd(.TRTH$archive_dir)

    foreach(i = 1:length(.TRTH$files.gz)) %dopar% 
    { # unzip in parallel
        filename.gz <- .TRTH$files.gz[i]
        filename.csv <- substr(filename.gz,1,(nchar(filename.gz)-3))
	    #unzip the file
	    print(paste("unzipping ",filename.gz, sep=""))
        #system(paste("gzip -d -f ",archive_dir,filename.gz,sep=""))
        system(paste("gunzip -f < ", .TRTH$archive_dir, filename.gz, " > ", .TRTH$tmp_archive_dir, filename.csv, sep=""))
        gc()
        print(paste(filename.gz, " unzipped.", sep=""))
    }
    ignored.csvs <- NULL #this will hold the names of CSVs that already have a header    
    setwd(.TRTH$tmp_archive_dir) #this directory contains the big CSVs that were unzipped
    .TRTH$big.csv <- list.files(.TRTH$tmp_archive_dir)
    for (i in 1:length(.TRTH$big.csv)) 
    {
        filename.csv <- .TRTH$big.csv[i]
        # Use awk to split the big CSV into daily CSVs.  Each CSV will have a single
        # row which we will then overwrite with the column headers.  Then we'll
        # use awk again to put the data into the split files
        #
        # First, make empty files (er, 1-row files that will be overwritten with header)
        # awk string says to make a file and put this row in it if the RIC or date are different than the previous row's RIC/date
        print(paste('Making headers from', filename.csv))
        system(paste('awk -v f2="" -F "," '," '",'{f1 = $1"."$2".csv";if(f1 != f2) { print >> f1; close(f2); f2=f1; } }',"' ",filename.csv, sep=""))
        
        tmpfiles <- list.files(.TRTH$tmp_archive_dir)
        files.header <- tmpfiles[grep("RIC",tmpfiles)]

        big.files <- tmpfiles[grep("@", tmpfiles)] #Big zipped CSVs from Reuters have e-mail address in name
        #big.files <- tmpfiles[grep(job.name, tmpfiles)]
        # csv files will be everthing that is not in "ignore" below
        # all these things we're ignoring should actually be in path.output, not here
        ignore <- c(big.files, files.header, 'NA', "Report", 
                    tmpfiles[grep("Tick2Sec|TRTH_config_file", tmpfiles)],
                    tmpfiles[grep("\\.rda", tmpfiles)], 
                    tmpfiles[grep("\\.RData", tmpfiles)],
                    tmpfiles[grep("missing_instruments", tmpfiles)],
                    ignored.csvs
                    )
        tmp.files.csv <- tmpfiles[!tmpfiles %in% ignore]
        # Extract single row and overwritting file with it 
        #(shouldn't be necessary anymore because it shouldn't have more than 1 row)
        system(paste('tail -1 "', files.header, '" > header.csv', sep="")) # extract one line
        # head -1 "RIC.Date[G].csv" > header.csv
        system(paste('mv header.csv "', files.header, '"', sep="")) # replace files.header with only one line
        # mv header.csv "RIC.Date[G].csv"

        for (fl in tmp.files.csv) { # make files with header that awk will later populate
            system(paste('cp "', files.header, '" ', paste(.TRTH$tmp_archive_dir, fl, sep=""), sep=""))
            #cp "#RIC.Date[G].csv" /home/garrett/TRTH/archive/GEM1-U1.01-APR-2008.csv
        }
        # after we've put a header in a file, we need to ignore that file the
        # next time through the loop so that we don't overwrite data.
        ignored.csvs <- c(ignored.csvs, tmp.files.csv)

        # If all of the files that awk just created already exist in csv_dir and overwrite==FALSE, then
        # there is no need to split this csv because we're not going to move any to csv_dir anyway.
        # a file in csv_dir might be "2012.02.10.AAPL.O.csv", (or *.csv.gz) but we have "AAPL.O.10-FEB-2012.csv"
        tmp <- gsub("\\.csv", "", tmp.files.csv)
        new.names <- do.call(c, lapply(strsplit(tmp, "\\."), function(x) {
            day <- gsub("-", ".", as.Date(x[length(x)], format='%d-%b-%Y'))
            fl <- make.names(paste(x[-length(x)], collapse="."))
            paste(day, "/", day, ".", fl, sep="")
        }))
        if (!all(file.exists(paste(paste(.TRTH$csv_dir, new.names, sep=""), ".csv", sep=""))) || 
            !all(file.exists(paste(paste(.TRTH$csv_dir, new.names, sep=""), ".csv.gz", sep=""))) || 
            isTRUE(.TRTH$overwrite)) {
            ## Split the Files 
            print(paste("Splitting ",filename.csv,sep=""))
            # The following awk will put data in our CSV files which currently only have column headers;
            #  Improved awk w/ file close to deal with 'too many open files', thanks to Josh Ulrich
            system(paste('awk -v f2="" -F "," '," '",'{f1 = $1"."$2".csv";print >> f1; if(f1 != f2) { close(f2); f2=f1; } }',"' ",filename.csv, sep=""))
            ## command line awkstring would look like this:
            # awk -v f2="" -F ","  '{f1 = $1"."$2".csv"; print >> f1; if(f1 != f2) { close(f2); f2=f1; } }' sourcefile.csv
            ## NOTE: if you get errors on 'too many open files' from awk, you'll need to adjust ulimit/nolimit
            print(paste('Done splitting ', filename.csv, sep=""))
        } else print('All CSVs created by awk already exist. Not re-splitting')
        # remove header file
        invisible(file.remove(paste(.TRTH$tmp_archive_dir, files.header, sep="")))
        # remove unzipped csv
        invisible(file.remove(paste(.TRTH$tmp_archive_dir, filename.csv, sep="")))
	    ## Zip the File
        # print(paste("zipping ",filename.csv,sep=""))
        # system(paste("gzip -f ",archive_dir,filename.csv,sep=""))
    }
    .TRTH$files.csv <- unique(c(tmpfiles[!tmpfiles %in% ignore], ignored.csvs))

    # Move split CSVs into csv_dir
    files.xts <- NULL
#    foreach (k = icount(length(.TRTH$files.csv))) %dopar%
    for (k in 1:length(.TRTH$files.csv))
    {
        #print(k)
        name.csv <- .TRTH$files.csv[k]                        # "ASBC.O.08-JAN-2011.csv"
        #name <- unlist(strsplit(name.csv,".",fixed=TRUE))[1]
        spl.name <- unlist(strsplit(name.csv, "\\."))   # "ASBC" "O" "08-JAN-2011" "csv" 
        last2 <- (length(spl.name) - 1):length(spl.name)# 3 4
        name <- paste(spl.name[-last2], collapse=".")   # "ASBC.O"
        #RIC.date <- try(as.Date(unlist(strsplit(name.csv,".",fixed=TRUE))[2], format="%d-%b-%Y"))
        RIC.date <- try(as.Date(spl.name[last2[1]], format="%d-%b-%Y"))
        date.format <- gsub("-",".",RIC.date)

        ## Handle leading digits and VIX and Cash
        name.new <- if(substr(name,1,1)==1){
            make.names(substr(name,2,nchar(name)))
        } else make.names(name)

        ## Create directory if it does not exist
        dir.create(paste(.TRTH$csv_dir, date.format, "/", sep=""), showWarnings=FALSE, recursive=TRUE, mode='0775') #mode='0664'

        ## Move files to appropriate place
        if (isTRUE(.TRTH$overwrite)) {
            system(paste("mv -fv ", name.csv, " ", .TRTH$csv_dir, date.format, "/", date.format, ".", name.new, ".csv", sep=""))
        } else if (!file.exists(paste(.TRTH$csv_dir, date.format, "/", date.format, ".", name.new, ".csv.gz", sep=""))) {
            system(paste("mv -nv ", name.csv, " ", .TRTH$csv_dir, date.format, "/", date.format, ".", name.new, ".csv", sep=""))
        } else print(paste(date.format, ".", name.new, ".csv.gz not overwritten.", sep=""))
        #print(paste(date.format, name.new, "moved", sep=" "))
        files.xts <- rbind(files.xts,as.data.frame(cbind(name.new,date.format),stringsAsFactors=FALSE))
    }
    files.xts$type <- rep(NA, NROW(files.xts))

    .TRTH$files.xts <- files.xts
    assign('.TRTH', .TRTH, pos=.GlobalEnv)
    setwd(.TRTH$archive_dir)
 
    if (isTRUE(.TRTH$use.instrument)) {
        missing_i <- NULL
        instr_s <- unique(files.xts[,'name.new'])
        alldefined <- unique(c(ls_instruments(), ls_instruments_by("X.RIC", NULL, in.slot='identifiers')))
        #FIXME: we really need a list of X.RICs, not a list of instrument_names that have X.RICs
        print(paste('Creating files.xts.  No more than', 
                    length(instr_s[!instr_s %in% alldefined]), 
                    'missing instruments will have to be defined'))
        missing_list <- list() # list to hold auto-defined missing instruments
        for(i in 1:length(instr_s)){
            instr <- getInstrument(instr_s[i], silent=TRUE)
            iauto <- NULL
            if(is.instrument(instr)){ 
                files.xts[files.xts$name.new ==instr_s[i],]$type <- paste(instr$type, collapse=";")
            } else {
                #NOTE: If we skip all this define-on-the-fly stuff, it would be much faster.
                pid <- try(parse_id(instr_s[i]))
                tmpid <- if(!inherits(pid, 'try-error') 
                            && !"" %in% c(pid$root, pid$suffix)) {
                    paste(pid$root, pid$suffix, sep="_")
                } else if (!inherits(pid, 'try-error')) {
                    pid$root
                } else instr_s[i]
                iauto <- instrument.auto(tmpid, currency=.TRTH$default_currency, 
                                        default_type=.TRTH$default_type, identifiers=list(X.RIC=instr_s[i]), assign_i=FALSE)
                if (!is.instrument(iauto)) {
                    warning(paste("Could NOT create ", .TRTH$default_type, " from ", 
                                instr_s[i], ". Creating _unknown_ instrument instead.", sep=""))
                    iauto <- try(suppressWarnings(instrument.auto(instr_s[i], currency=.TRTH$default_currency,
                                                default_type="unknown", assign_i=FALSE)))
                }
                missing_list[[iauto$primary_id]] <- iauto
                #assign(iauto$primary_id, iauto, pos=missing_i_envir) 
                files.xts[files.xts$name.new==instr_s[i],]$type <- paste(iauto$type, collapse=";")
                missing_i <- c(missing_i, instr_s[i])
            }
        }
      
        if (length(missing_list) > 0) {
            # Remove everything from .instrument, put back the auto-defined missing instruments and save them
            print("saving RData file with auto-defined missing instruments.")
            try(rm_instruments(), silent=TRUE)
            lapply(missing_list, function(x) {
                assign(x$primary_id, x, pos=FinancialInstrument:::.instrument)
            })
            saveInstruments(paste("missing_instr",  format(Sys.time(), "%Y.%m.%d_%H%M%S"), sep='_'), .TRTH$path.output)
            # now that we've saved only the newly defined instruments, we can load back our other instruments
            loadInstruments(.TRTH$instrument_file)
            if (!is.null(iauto)) {
                .TRTH$missing_i <- missing_i <- data.frame(symbol=missing_i,type=iauto$type[1]) #legacy
                write.csv(missing_i,file=paste(.TRTH$path.output,'missing_instruments.CSV',sep='')) 
            }
        }
    }
    .TRTH$files.xts <- files.xts
    gc()
    assign('.TRTH', .TRTH, pos=.GlobalEnv)
    .TRTH
}


FEreut2xts <- function(.TRTH) {
    if (missing(.TRTH) && !exists(".TRTH")) stop("Run configureTRTH function first")
    #attach(.TRTH)
    # Make sure csv_dir exists since it is where we read the data from
    if (!file.exists(.TRTH$csv_dir)) stop("There is no directory", paste(.TRTH$csv_dir))
    if (is.null(.TRTH$files.xts)) stop("Cannot find 'files.xts' -- Run splitCSV first")

    # edit text file so others can see what job we're working on
    system(paste('echo "', Sys.time(), ' FEreut2xts job.name: ', .TRTH$job.name, '" > ', 
                 paste(.TRTH$path.output, "current.job.txt", sep=""), 
                 sep=""))

    files.xts <- .TRTH$files.xts

    oldTZ <- Sys.getenv("TZ")
    Sys.setenv(TZ='GMT')

    write.tick <- TRUE #if the tickdata file already exists and overwrite==FALSE, this will be set to FALSE
    write.sec <- TRUE #if the secdata file already exists and overwrite==FALSE, this will be set to FALSE


    nc <- nchar(.TRTH$path.output) # make sure path.output ends with a forward slash
    if(substr(.TRTH$path.output, nc, nc) != "/") .TRTH$path.output <- paste(.TRTH$path.output, "/", sep="") 

    # Function that we'll use to save charts of the data
    makeImages <- function(Data, dir, RIC, date) {
        stopifnot(file.exists(paste(dir, RIC, sep="")))
        ## Bid
        dir.create(paste(dir, RIC, "/Bid.Image/", sep=""), showWarnings=FALSE, mode='0775')
        png(filename=paste(dir,RIC,"/Bid.Image/",date,".",RIC,".png",sep=""),width=1500,height=1000)
        try(chartSeries(to.minutes(Data$Bid.Price,1),type="bar"),silent=TRUE)
        dev.off()
        ## Ask
        dir.create(paste(dir,RIC,"/Ask.Image/",sep=""), showWarnings=FALSE, mode='0775')
        png(paste(dir,RIC,"/Ask.Image/",date,".",RIC,".png",sep=""),width=1500,height=1000)
        try(chartSeries(to.minutes(Data$Ask.Price,1),type="bar"),silent=TRUE)
        dev.off()
        ## Price
        Data.1 <- Data[!is.na(Data$Price),]
        if(dim(Data.1)[1]>50)
        {
            dir.create(paste(dir,RIC,"/Price.Image/",sep=""), showWarnings=FALSE, mode='0775')
            png(paste(dir,RIC,"/Price.Image/",date,".",RIC,".png",sep=""),width=1500,height=1000)
            try(chartSeries(to.minutes(na.omit(Data$Price),1),type="bar"),silent=TRUE)
            dev.off()
        }
    }

    Out <- foreach(ii=icount(NROW(files.xts)), .inorder=FALSE, .errorhandling='pass') %dopar% {
        RIC=files.xts[ii, 1]
        date=files.xts[ii, 2] 
        type=files.xts[ii, 3]        

        file.name.xts <- paste(.TRTH$tick_dir, RIC, "/", date, ".", RIC, ".RData", sep="")
        file.name.sec <- paste(.TRTH$sec_dir, RIC, "/", date, ".", RIC, ".RData", sep="")	
        if(!isTRUE(.TRTH$overwrite)) {
            if (file.exists(file.name.xts)){
	            cat(paste(file.name.xts, "already exists, not overwriting\n"))
                write.tick <- FALSE
                .TRTH$tick.image <- FALSE
            }
            if (file.exists(file.name.sec)) {
            	cat(paste(file.name.sec, "already exists, not overwriting\n"))
                write.sec <- FALSE
                .TRTH$sec.image <- FALSE
            }
        }

        print(paste(date, RIC, paste(c("tick", "sec")[c(write.tick, write.sec)], collapse=" "), ii, "of", NROW(files.xts), sep=" "))

        # if xts and sec data already exist for this product/Date, and overwrite == FALSE, 
        # there is nothing to be done -- return NULL
        if (!any(c(write.tick, write.sec))) return(NULL) 
        #TODO: unzip to a tempdir
        CSV.name <- paste(.TRTH$csv_dir, date, '/', date, '.', RIC, '.csv', sep="")
        if (!file.exists(CSV.name) && file.exists(paste(CSV.name, ".gz", sep=""))) {
            #only zipped file on disk. We'll have to unzip.
            system(paste("gzip -d -f ", CSV.name, ".gz", sep=""))
        }
        Data <- try(read.csv(CSV.name, stringsAsFactors=FALSE, header=TRUE), silent=TRUE)
        if (inherits(Data, 'try-error')) {
            Data <- read.csv(CSV.name, stringsAsFactors=FALSE, header=FALSE)
            colnames(Data) <- make.names(Data[1, ])
            Data <- Data[-1, ]
        }
        # Now that we've read the CSV, zip it and delete original to conserve disk space
        # print(paste("zipping ", CSV.name, sep=""))
        system(paste("gzip -f ", CSV.name, sep=""))

        OTC.remove <- grep("IRGCOND",Data$Qualifiers)
        #OTC.remove <- c(OTC.remove,grep("High[USER]",Data$Qualifiers,fixed=TRUE))
        #OTC.remove <- c(OTC.remove,grep("Low[USER]",Data$Qualifiers,fixed=TRUE))
        OTC.remove <- c(OTC.remove, grep("[USER]", Data$Qualifiers, fixed=TRUE))

        if(substr(RIC,1,(nchar(RIC)-2))=="ICF"){OTC.remove <- NULL}
        if(substr(RIC,1,(nchar(RIC)-2))=="DOL"){OTC.remove <- NULL}
        if(dim(Data)[1]<=25){return(NULL)}

        ## Remove block trades
        if(length(OTC.remove)){Data <- Data[-OTC.remove, ]}

        index.new <- as.POSIXct(paste(Data$Date.G.,Data$Time.G,sep=" "),format="%d-%b-%Y%H:%M:%OS",tz="GMT")

        ## Force Everything to numerics <-- should not be necessary, but I'll leave it as is
        Data <- Data[,c("Price","Volume","Bid.Price","Bid.Size","Ask.Price","Ask.Size")]
        Data$Price <- as.numeric(Data$Price)
        Data$Volume <- as.numeric(Data$Volume)
        Data$Bid.Price <- as.numeric(Data$Bid.Price)
        Data$Bid.Size <- as.numeric(Data$Bid.Size)
        Data$Ask.Price <- as.numeric(Data$Ask.Price)
        Data$Ask.Size <- as.numeric(Data$Ask.Size)

        Data <- xts(Data,order.by=index.new,tz="GMT")

        if (isTRUE(.TRTH$use.instrument)) {
            ## Turn bids/offers that are less than zero into NA for outrights
            type <- try(unlist(strsplit(type, ";")))
            if (inherits(type, 'try-error')) {
                warning('type is incorrect. Using "synthetic"')
                type <- 'synthetic'
            }        
            if(!any(c("unknown", "synthetic") %in% type))
            { #outrights
                Data$Bid.Price[Data$Bid.Price < 0, ] <- NA
	            Data$Ask.Price[Data$Ask.Price < 0, ] <- NA
	            Data$Price[Data$Price < 0, ] <- NA
            } 
        }
    
        ## If Bid.Price and Bid.Size are zero set both to NA
        zero.replace <- which(Data$Bid.Price == 0 & Data$Bid.Size == 0)
        if (length(zero.replace) != 0) {
            Data$Bid.Price[zero.replace] <- NA
            Data$Bid.Size[zero.replace] <- NA
        }
        ## Do same thing with Ask Price/Size
        zero.replace <- which(Data$Ask.Price == 0 & Data$Ask.Size == 0)
        if (length(zero.replace) != 0) {
            Data$Ask.Price[zero.replace] <- NA
            Data$Ask.Size[zero.replace] <- NA
        }

        ## Carry last bid/offer forward
        Data$Bid.Price <- na.locf(Data$Bid.Price)
        Data$Bid.Size <- na.locf(Data$Bid.Size)
        Data$Ask.Price <- na.locf(Data$Ask.Price)
        Data$Ask.Size <- na.locf(Data$Ask.Size)

        ## Remove Trades with Volume of zero
        Volume.remove <- which(Data$Volume == 0)
        if(length(Volume.remove) != 0) Data <- Data[-Volume.remove]

        ## Remove Bids with Size of zero
        Bid.remove <- which(Data$Bid.Size == 0)
        if(length(Bid.remove) != 0) Data <- Data[-Bid.remove]

        ## Remove Asks with Size of zero
        Ask.remove <- which(Data$Ask.Size == 0)
        if(length(Ask.remove) != 0) Data <- Data[-Ask.remove]

        if(dim(Data)[1]<=25){return(NULL)}

        ## Remove Price w/ Volume of NA and
        ## Volume w/ Price of NA	
        na.remove <- c(which(!is.na(Data$Price) & is.na(Data$Volume)),
        which(is.na(Data$Price) & !is.na(Data$Volume)))
        if (length(na.remove)!=0) {	Data <- Data[-na.remove] }

        ## not enough rows
        if(dim(Data)[1]<=10){return(NULL)}

        ## Remove leading NAs on Bid/Ask
        bid.remove <- which(is.na(Data$Bid.Price))
        ask.remove <- which(is.na(Data$Ask.Price))
        union.remove <- c(bid.remove,ask.remove)
        if(length(union.remove)>0){Data <- Data[-union.remove]}

        ## not enough rows
        if(dim(Data)[1]<=25){return(NULL)}

        if(write.tick) {
            dir.create(paste(.TRTH$tick_dir, RIC, sep=""), showWarnings=FALSE, mode='0775')
            assign(RIC, Data)  # Rename Data to RIC
            save(list=RIC, file=file.name.xts)
        }

        datarange <- range(index(Data),na.rm = TRUE)
        datarange.dif <- difftime(datarange[2],datarange[1],units="secs")
        if(isTRUE(.TRTH$tick.image) && datarange.dif>3600) makeImages(Data, .TRTH$tick_dir, RIC, date)

        # Convert to 1 second data and save
        if (write.sec) {
            dir.create(paste(.TRTH$sec_dir, RIC, "/", sep=""), showWarnings=FALSE, mode='0775')
            secData <- to_secBATV(Data)
            if (length(secData) == 0) return(NULL)
            assign(RIC, secData)
            save(list=RIC, file=file.name.sec)
        }
        if (isTRUE(.TRTH$sec.image) && datarange.dif > 3600) makeImages(Data, .TRTH$sec_dir, RIC, date)
        gc()
    } # End foreach loop
    #rm(list = 'RIC')
    gc()
    Sys.setenv(TZ=oldTZ)
    save(.TRTH, file=paste(.TRTH$path.output, 'config.env.RData', sep=""))
    assign('.TRTH', .TRTH, pos=.GlobalEnv)
    Out
}  ## End fn reut2xts 

#
## now clean up
#files.rm <- list.files(archive_dir)
#files.rm <- files.rm[-grep(".csv.gz",files.rm)]
#files.rm <- files.rm[grep(".csv",files.rm)]
#file.remove(files.rm)
#file.remove('files.xts.tmp.rda')
#
#rm(missing_i)
##rm(Out)

###############################################################################
# Copyright (c) 2009-2011
# Peter Carl,  Brian G. Peterson, Lance Levenson, Joshua Ulrich, Garrett See
#
# This code is distributed under the terms of the GNU Public License (GPL)
# for full details see the file COPYING
#
# $Id: TRTH_BackFill.R 944 2012-02-23 17:21:32Z gsee $
#
###############################################################################

