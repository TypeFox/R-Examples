#' Read nonmem table files produced.
#' 
#' The function reads in nonmem table files.  The files can be created from the 
#' \code{$EST} line or from the \code{$SIM} line in a NONMEM model file.
#' 
#' Currently the function searches the \code{$TABLE} for multiple header lines, 
#' and uses that to identify multiple simulations. The function expects at least
#' one header line  (\code{NOHEADER} option is not allowed in NONMEM table 
#' files).
#' 
#' @param nm_table The NONMEM table file to read. A text string. If \code{dplyr}
#'   and \code{readr} are available and \code{method="default"} or
#'   \code{method="readr_*"} then \code{nm_table} can either a path to a file, a
#'   connection, or literal data (either a single string or a raw vector). Files
#'   ending in .gz, .bz2, .xz, or .zip will be automatically uncompressed. Files
#'   starting with http://, https://, ftp://, or ftps:// will be automatically
#'   downloaded. Remote gz files can also be automatically downloaded &
#'   decompressed.
#' @param only_obs Should the non-observation lines in the data set be removed? 
#'   Currently filtered uisng the expected \code{MDV} column. \code{TRUE} or 
#'   \code{FALSE}.
#' @param method Can be one of \code{default}, \code{readr_1},  \code{readr_2}, 
#'   \code{readr_3}, \code{slow}. \code{readr_1} should be fastest.  All but 
#'   \code{slow} require \code{\link[dplyr]{dplyr}} and 
#'   \code{readr} (version >= 0.2.2).
#' @param quiet Should there be verbose messages about what the function is 
#'   doing?
#' @param sim_num Should the function add a column to the returned data frame 
#'   that identifies the simulation number (if present)?
#' @param sim_name The name of the resulting column in the returned data frame if \code{sim_num} is true.
#'   
#' @return Returns a data frame of the simulated table with an added column for 
#'   the simulation number. The data frame is given class \code{c("tbl_df", 
#'   "tbl", "data.frame")} for easy use with \code{\link[dplyr]{dplyr}}.
#' @export
#'   

read_nm_table <- function (nm_table, only_obs=FALSE, method="default",
                           quiet=TRUE,sim_num=FALSE,sim_name="NSIM"){
  
  # \code{\link[dplyr]{dplyr}} and \code{\link[readr]{readr}} are available and
  # \code{\link[readr]{readr}} (version >= 0.2.2).
  "MDV" <- "packageVersion" <- "write.table" <- "count.fields" <- "read.table" <- NULL
  rm(list=c("MDV","packageVersion","write.table","count.fields","read.table"))
  
  if(method=="default"){
    if (requireNamespace("readr", quietly = TRUE) && 
          (packageVersion("readr") >= "0.2.2") && 
          requireNamespace("dplyr", quietly = TRUE)) {
      method <- "readr_1" 
    } else {
      method <- "slow" 
    }
  }
  
  read_nm_tab_readr_1 <- function(nm_table,sim_num){
    "MDV" <- "packageVersion" <- "write.table" <- "count.fields" <- "read.table" <- NULL
    rm(list=c("MDV","packageVersion","write.table","count.fields","read.table"))
    
    #tab_dat <- read_table(nm_table, skip = 1) 
    #tab_dat <- tab_dat %>% mutate_each(funs(suppressWarnings(as.numeric(.))))
    
    ## get header names
    header_line <- readr::read_lines(nm_table,n_max=2)[2]
    comma_sep <- FALSE
    if(length(grep(",",header_line))!=0) comma_sep <- TRUE
    header_line <- sub("^\\s+","",header_line) 
    header_names <- strsplit(header_line,"\\s+,*\\s*")[[1]]
    
    if(!comma_sep){
      tab_dat <- readr::read_table(nm_table, col_names = header_names, 
                                   col_types=paste0(rep("d",length(header_names)),collapse = ""),
                                   skip = 2) 
    } else {
      tab_dat <- readr::read_csv(nm_table, col_names = header_names, 
                                 col_types=paste0(rep("d",length(header_names)),collapse = ""),
                                 skip = 2) 
    }
    
    # Handle multiple simulations
    if(any(is.na(tab_dat[1]))){ 
      if(sim_num){
        ## create simulation number
        args <- lazyeval::interp(~ cumsum(is.na(var))+1, var = as.name(names(tab_dat)[1]))
        tab_dat <- dplyr::mutate_(tab_dat,NSIM=args)
      }
      
      ## filter out NA columns
      args <- lazyeval::interp(~ !is.na(var), var = as.name(names(tab_dat)[1]))
      tab_dat <- dplyr::filter_(tab_dat,args)
    }
    return(tab_dat)
  }
  
  
  read_nm_tab_readr_2 <- function(nm_table){
    "MDV" <- "packageVersion" <- "write.table" <- "count.fields" <- "read.table" <- NULL
    rm(list=c("MDV","packageVersion","write.table","count.fields","read.table"))
    
    ## get header names
    header_line <- readr::read_lines(nm_table,n_max=2)[2]
    comma_sep <- FALSE
    if(length(grep(",",header_line))!=0) comma_sep <- TRUE
    header_line <- sub("^\\s+","",header_line) 
    header_names <- strsplit(header_line,"\\s+,*\\s*")[[1]]
    
    
    ## Check if we have unequal number of fields in the file
    ## used for multiple simulations
    tmp   <- readr::read_lines(nm_table)
    inds  <- grep("TABLE",tmp)
    if (length(inds)!=1){
      inds  <- inds[c(2:length(inds))]
      inds2 <- inds+1
      tempfile<- paste(nm_table,".xptmp",sep="")
      write.table(tmp[-c(inds,inds2)],file=tempfile,
                  row.names=FALSE,quote=FALSE,col.names = FALSE)
      #assign(paste("n.",filename,sep=""),read.table(tempfile,skip=2,header=T,sep=sep.char))
      if(!comma_sep){
        tab_dat <- readr::read_table(tempfile, col_names = header_names, 
                                     col_types=paste0(rep("d",length(header_names)),collapse = ""),
                                     skip = 2) 
      } else {
        tab_dat <- readr::read_csv(tempfile, col_names = header_names, 
                                   col_types=paste0(rep("d",length(header_names)),collapse = ""),
                                   skip = 2) 
      }
      unlink(tempfile)
    } else {
      if(!comma_sep){
        tab_dat <- readr::read_table(nm_table, col_names = header_names, 
                                     col_types=paste0(rep("d",length(header_names)),collapse = ""),
                                     skip = 2) 
      } else {
        tab_dat <- readr::read_csv(nm_table, col_names = header_names, 
                                   col_types=paste0(rep("d",length(header_names)),collapse = ""),
                                   skip = 2) 
      }
    }
    
    return(tab_dat)
  }
  
  read_nm_tab_readr_3 <- function(nm_table,sim_num){
    "MDV" <- "packageVersion" <- "write.table" <- "count.fields" <- "read.table" <- NULL
    rm(list=c("MDV","packageVersion","write.table","count.fields","read.table"))
    
    ## get header names
    header_line <- readr::read_lines(nm_table,n_max=2)[2]
    comma_sep <- FALSE
    if(length(grep(",",header_line))!=0) comma_sep <- TRUE
    header_line <- sub("^\\s+","",header_line) 
    header_names <- strsplit(header_line,"\\s+,*\\s*")[[1]]
    
    # read in all lines of file
    tmp   <- readr::read_lines(nm_table)
    #tmp_table <- nm_table
    tmp_table <- paste(tmp,collapse="\n")
    fun_name <- "readr::read_table"
    if(comma_sep) fun_name <- "readr::read_csv"
    skip=2
    
    ## Check for multiple table lines
    inds  <- grep("TABLE",tmp)
    if (length(inds)!=1){
      inds2 <- inds+1 # additional header lines
      if(sim_num){
        NSIM <- rep(1,length(tmp))
        NSIM[inds] <- NA
        NSIM <- cumsum(is.na(NSIM))
        tmp <- paste(tmp,NSIM) 
        header_names <- c(header_names,"NSIM")
      }
      tmp_table <- paste(tmp[-c(inds,inds2)],collapse="\n")
      skip=0
    }
    tab_dat <- do.call(eval(parse(text=paste0(fun_name))),args = list(tmp_table,col_names = header_names, 
                                                                      col_types=paste0(rep("d",length(header_names)),collapse = ""),
                                                                      skip = skip))
    
    return(tab_dat)
  }
  
  
  
  
  read_nm_tab_slow <- function (filename, quiet){
    "MDV" <- "packageVersion" <- "write.table" <- "count.fields" <- "read.table" <- NULL
    rm(list=c("MDV","packageVersion","write.table","count.fields","read.table"))
    
    ## Check which type of separator we have in our tables
    header.line = scan(file=filename,nlines=1,skip=1,what="character",sep="\n",quiet=T)
    sep.char = ""
    if(length(grep(",",header.line))!=0) sep.char = ","
    
    ## Check if we have unequal number of fields in the file
    ## used for multiple simulations
    fields.per.line <- count.fields(filename)
    fields.in.first.line <- fields.per.line[1]
    fields.in.rest <- fields.per.line[-1]
    if((length(unique(fields.in.rest))!=1) ||
         (all(fields.in.first.line==fields.in.rest))){ 
      if(!quiet) {
        cat(paste("Found different number of fields in ",filename,".\n",sep=""))
        cat("This may be due to multiple TABLE and header rows \n")
        cat("caused by running multiple simulations in NONMEM (NSIM > 1).\n")
        cat("Will try to remove these rows. It may take a while...\n")
      }
      tmp   <- readLines(filename, n = -1)
      inds  <- grep("TABLE",tmp)
      if (length(inds)!=1){
        inds  <- inds[c(2:length(inds))]
        inds2 <- inds+1
        tempfile<- paste(filename,".xptmp",sep="")
        write.table(tmp[-c(inds,inds2)],file=tempfile,
                    row.names=FALSE,quote=FALSE)
        #assign(paste("n.",filename,sep=""),read.table(tempfile,skip=2,header=T,sep=sep.char))
        tab_dat <- read.table(tempfile,skip=2,header=T,sep=sep.char)
        tab_dat$NSUB <- rep(1:(nrow(tab_dat)/(inds[1]-3)), each=(inds[1]-3))
        unlink(tempfile)
      } else {
        #assign(paste("n.",filename,sep=""),read.table(filename,skip=1,header=T,sep=sep.char))
        tab_dat <- read.table(filename,skip=1,header=T,sep=sep.char)
      }
    } else {
      #assign(paste("n.",filename,sep=""),read.table(filename,skip=1,header=T,sep=sep.char))
      tab_dat <- read.table(filename,skip=1,header=T,sep=sep.char)
    }
    return(tab_dat)
  }
  
  tab_dat <- switch(method,
                    readr_1 = read_nm_tab_readr_1(nm_table, sim_num=sim_num),
                    readr_2 = read_nm_tab_readr_2(nm_table),
                    readr_3 = read_nm_tab_readr_3(nm_table,sim_num=sim_num),
                    slow = read_nm_tab_slow(nm_table,quiet=quiet)
  )
  
  ## remove non-observation rows
  if(only_obs){
    if(any("MDV"==names(tab_dat))){
      tab_dat <- dplyr::filter(tab_dat,MDV==0)   
    } else {
      warning('\nMDV data item not listed in header, 
              Could not remove dose events!\n')
    }
  }
  
  if(sim_num) names(tab_dat)[match("NSIM",names(tab_dat))] <- sim_name
  
  tab_dat <- data.frame(tab_dat)
  tab_dat <- dplyr::as_data_frame(tab_dat)
  
  return(tab_dat)
}