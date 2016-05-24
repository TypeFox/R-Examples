# Project: RcppOctave
# 
# Author: Renaud Gaujoux
# Created: Jul 21, 2015
###############################################################################

slashes <- function(x) gsub("\\", "/", x, fixed = TRUE)

Sys.path <- function(normalize = TRUE){
    path <- strsplit(Sys.getenv('PATH'), .Platform$path.sep)[[1L]]
    if( normalize ) path <- normalizePath(path, mustWork = FALSE)
    path
}

ac_msg_checking <- function(..., appendLF = FALSE){
    message("Checking ", ..., "... ", appendLF = appendLF)
}

ac_msg_result <- function(...) message(...)

ac_msg_notice <- function(...) message("configure: ", ...)

#' Retrieve a Program Configuration Variable
ac_prog_var <- function(prog, varname, msg = '', intern = FALSE, slashes = FALSE){
	if( is.character(intern) ) intern <- nzchar(intern)
	if( !nzchar(msg) ) msg <- paste(basename(prog), varname)
	message(sprintf("Checking %s... ", msg), appendLF = FALSE)
    cmd <- sprintf('"%s" --print %s', prog, varname)
    value <- shell(cmd, intern = TRUE)
    if( slashes ) value <- slashes(value)
	if( !intern ){
        message(value)
        cat(value)
    }
    invisible(value)
}

#' Retrieve a Program Configuration Path
ac_prog_varpath <- function(...){
    ac_prog_var(..., slashes = TRUE)
}

ac_shell <- function(cmd, msg, intern = FALSE, slashes = TRUE){
  if( is.character(intern) ) intern <- nzchar(intern)
  message(sprintf("Checking %s... ", msg), appendLF = FALSE)
  value <- shell(cmd, intern = TRUE, ignore.stderr = TRUE)
  if( slashes ) value <- slashes(value)
  if( !intern ){
    message(value)
    cat(value)
  }
  invisible(value)
}

#' Finds a Program Location 
ac_path_prog <- function(prog, notfound = '', path = "", msg = "", mode = c('first', 'all', 'deep'), strip.flags = FALSE, intern = FALSE){
    
    strip.flags <- !missing(strip.flags) || !strip.flags
    if( strip.flags ) prog <- strsplit(prog, ' -', fixed = TRUE)[[1L]][1L]
        
    if( !nzchar(msg) ) msg <- prog
    path <- path[nzchar(path)]
    if( !length(path) ) path <- Sys.path()
    
    if( length(path) == 1L ) msg <- sprintf("%s [in %s]", msg, path)
    message(sprintf("Checking %s... ", msg), appendLF = FALSE)
    if( missing(mode) || !nzchar(mode) ) mode <- 'first'
    mode <- match.arg(mode)
    
    res <- sapply(path, function(p){
        if( mode == 'deep' ) f <- list.files(p, recursive = TRUE, pattern = prog, full.names = TRUE)
        else f <- Sys.which2(prog, p)
       
       # remove duplicates
       if( length(f) > 1L ){
            f <- f[order(nchar(basename(f)), decreasing = TRUE)]
            f <- f[!duplicated(tools::md5sum(f))]
        }
        f
    })

    # remove empty elements
    res <- unlist(res[lengths(res) > 0], use.names = FALSE)
    # slashes
    res <- slashes(res[nzchar(res)])
    
    if( mode == 'first' ) res <- head(res, 1L)
    
    res <- unname(res)
    if( !length(res) ){
        message('no')
        if( !intern ) cat(notfound)
    }else{
        if( length(res) == 1L ) message(res)
        else{
            message(sprintf('multiple [%i]', length(res)))
            lapply(paste(' *', res), message)
        }
        if( !intern ) cat(res, sep = "\n")
    }
    invisible(res)
    
}

strm <- function(x) message(paste0(capture.output(str(x)), collapse = "\n"))

Sys.which2 <- function(names, path = ""){
    
    if( !nzchar(path) ) path <- Sys.path()
    # backup and restore original PATH
    PATH <- Sys.getenv('PATH')
    on.exit( Sys.setenv(PATH = PATH) )
    
    res <- sapply(path, function(p){
                Sys.setenv(PATH=p)
                normalizePath(Sys.which(names), mustWork = FALSE)
        }
    )
    # slashes
    res <- slashes(res[nzchar(res)])
    unname(res)
}

#' Finds Rtools
ac_path_rtools <- function(gcc, intern = FALSE){
    
    message("Checking Rtools ... ", appendLF = FALSE)
    
    # look for file VERSION.txt that contains the string "rtools"
    parts <- strsplit(gcc, '/')[[1L]]
    path <- NULL
    for( i in 1:length(parts) ){
        p <- do.call(file.path, as.list(head(parts, i)))
        if( file.exists(version_file <- file.path(p, 'VERSION.txt')) && 
                length(grep('rtools', readLines(version_file), ignore.case = TRUE)) ){
            path <- p
            break
        }
    }
    
    if( !length(path) ){ # not found
        path <- slashes(normalizePath(file.path(dirname(gcc), "../..")))
      
        # if( !intern ) cat("")
        message("no")
        ac_msg_notice("using guessed Rtools path ", appendLF = FALSE)
    }
    # else
    { # found
        if( !intern ) cat(path)
        message(path)
        invisible(path)
    }
}

#' Finds a gcc version in Rtools/ that is compatible with given gcc version 
ac_cc_compatible <- function(target, rtools, absolute = TRUE, nolookup = TRUE){
    
    r_gcc <- NULL
    if( file_test("-f", rtools) ){
      r_gcc <- rtools
      rtools <- ac_path_rtools(rtools, intern = TRUE)
    }
    if( !length(rtools)){
      message("none")
      cat('')
      return(invisible())
    }
    
    if( !nolookup ) path <- ac_path_prog('gcc[-0-9.]*.exe$', path = rtools, msg = "all Rtools compiler(s)", mode = 'deep', intern = TRUE)
    else path <- r_gcc
      
    .cc_spec <- function(cc){
        v <- shell(sprintf("\"%s\" -dumpversion", cc), intern = TRUE)
        m <- shell(sprintf("\"%s\" -dumpmachine", cc), intern = TRUE)
        res <- list(version = package_version(v), machine = m, full = paste0(m,'-', v))
#        print(res)
    }
    
    ac_msg_checking("Octave compiler version")
    oct_spec <- .cc_spec(target)
    ac_msg_result(oct_spec$full)
    ac_msg_checking("compatible Rtools compiler")
    
    # give priority to current compiler if any
    if( !is.null(r_gcc) ){
      m <- match(normalizePath(path), normalizePath(r_gcc))
      path <- path[order(m)]
    }
    
    for( p in path ){
        spec <- .cc_spec(p)

        if( spec$full == oct_spec$full || 
              ( normalizePath(p) == normalizePath(r_gcc) && spec$version >= oct_spec$version) ||
              ( spec$machine == oct_spec$machine && spec$version >= oct_spec$version)
            ){
           bin <- dirname(p)
           message(sprintf("%s [%s]", basename(p), bin))
           
           # check compiler aliases exist therein
            # gcc
           message("Checking gcc alias... ", appendLF = FALSE)
           if( !file.exists(alias <- file.path(bin, 'gcc.exe')) || 
                   tools::md5sum(alias) !=  tools::md5sum(p) ){
               break
           }
           message(alias)
           
            # g++
           cxx <- ac_path_prog('g++.exe', path = bin, msg = "g++", intern = TRUE)
           if( !nzchar(cxx) ) break
           
           p <- gsub(rtools, '', bin, fixed = TRUE)
           cat(p)
           
           if( !absolute ) p <- gsub(rtools, "", p, fixed = TRUE)
           
           return(invisible(p))
        }
    }
    # not found
    message('none')
    cat('')
    invisible()
}

ac_cc_compatible_octave <- function(dir, cc, nolookup = FALSE){
    
  if( is.character(nolookup) ) nolookup <- nzchar(nolookup)
  
  r_gcc <- ac_path_prog(cc, strip.flags = TRUE, msg = 'current Rtools compiler', intern = TRUE)
  
  # check octave-config is in the directory
  octave_config <- ac_path_prog('octave-config', path = dir, intern = TRUE)
  if( !length(octave_config) ){
      return(invisible())
  }

  # get Octave version
  ac_msg_checking("Octave version")
  ver <- shell(sprintf("%s -p VERSION", octave_config), intern = TRUE)
  ac_msg_result(ver)
  lookup_dir <- dir
  if( compareVersion(ver, "3.6.4") <= 0 ){
    lookup_dir <- c(lookup_dir, file.path(dir, "../mingw/bin"))
  }
  
  octave_gcc <- ac_path_prog('gcc.exe', path = lookup_dir, msg = 'Octave compiler', intern = TRUE)
  ac_cc_compatible(octave_gcc, r_gcc, nolookup = nolookup)
}

