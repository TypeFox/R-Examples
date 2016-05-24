######################
## WINDOWS ROUTINES ##
######################

## Taken from pbatR
## getTimeStamp(...)                                                #
## Gets a unique time-stamp string for the output!                  #
getTimeStamp <- function() {
  zpad <- function(n, pad=2) {
    if( nchar(n)<pad )
      return( paste( rep("0",pad-nchar(n)), n, sep="" ) );
    return(n);
  }

  d <- as.POSIXlt( Sys.time() );
  return( paste( 1900+d$year, zpad(d$mon), zpad(d$mday), zpad(d$hour), zpad(d$min), zpad(floor(d$sec)), sep="" ) );  ## R2.3 change... seconds decide to have bloody decimal points... why can't we just be consistent between releases???? WHY??????
}

## Also taken from pbatR!
## OS functionality
#isWindows <- function()
#  return( Sys.info()["sysname"]=="Windows" );
## A more robust version of the above - I'm really not sure what to expect from vista
isWindows <- function()
  return(tolower(Sys.info()["sysname"]) == "windows" || tolower(R.Version()$os) == "mingw32")

## Returns the name of the file that will be touched upon completion
asyncRunCommands <- function(cmds){
  stamp <- paste("RBATCH", getTimeStamp(), runif(1)*10000000, sep="")
  rfile <- paste(stamp, ".R", sep="")
  touchfile <- paste(stamp, ".touch", sep="")

  f <- file(rfile, "w")
  for(cmd in cmds){
    ## handle the " character
    cmd <- paste(unlist(strsplit(cmd, "\"")), collapse = "\\\"")
    cat("system(\"", cmd, "\")\n", sep="", file=f)
    cat("cat(\"", cmd, " completed.\\n\")\n", sep="", file=f)
  }
  cat("f2 <- file(\"", touchfile, "\", \"w\")\ncat('touch',file=f2)\nclose(f2)\n", sep="", file=f)
  close(f)
  system(paste(Rwin(), "--vanilla <", rfile), wait=FALSE)

  return(touchfile)
}

parallelRunCommands <- function(cmds, n, sleep=10, wait=TRUE){ ## Sleep every second? every minute? every 10 seconds?
  cmdsList <- msplit(cmds, n)
  touchfiles <- rep("", n)
  for(i in 1:n)
    touchfiles[i] <- asyncRunCommands(cmdsList[[i]])

  #print("cmds\n")
  #print(cmds)
  #print("touchfiles\n")
  #print(touchfiles)

  if(wait){
    done <- rep(FALSE, n)
    while(!all(done)){
      for(i in 1:length(done))
        done <- file.exists(touchfiles[i])
      Sys.sleep(sleep)
    }

    ## And remove the files...
    #print(touchfiles)
    file.remove(touchfiles)
    file.remove(paste(substring(touchfiles, 1, nchar(touchfiles) - nchar(".touch")), ".R", sep=""))
  }

  return(invisible())
}

Rwin <- function(){
  r <- paste(Sys.getenv("R_HOME"), "/bin/R", sep="")
  if(!file.exists(r))
    r <- paste(r, ".exe", sep="")
  if(!file.exists(r))
    stop("Cannot find R executable/binary.")
  return(r)
}

#######################
rbatch.default <- function(){
  res <- parseCommandArgs(FALSE) #DF()
  #print(res)
  if(is.null(res$RBATCH)){
    return(rbatch.local())
  }else{
    res$RBATCH <- tolower(res$RBATCH)
    if(res$RBATCH == "local")
      return(rbatch.local())
    if(res$RBATCH == "lsf")
      return(rbatch.lsf())
    else if(res$RBATCH == "mosix")
      return(rbatch.mosix())
  }
  stop(paste("rbatch.default error, you supplied argument RBATCH='", res$RBATCH, "', but it should be one of 'local', 'lsf', or 'mosix' (case insensitive).", sep=""))
}

rbatch.parseCommandArgs <- function(BATCH, BATCHPOST, QUOTE, ARGQUOTE, RUN, MULTIPLIER, LOCAL=0){
  RES <- parseCommandArgs(FALSE) #DF()
  #cat("rbatch.parseCommandArgs, RES=\n")
  #print(RES)
  #print(RES$RUN)

  if(!is.null(RES$BATCH))
    BATCH <- RES$BATCH
  if(!is.null(RES$BATCHPOST))
    BATCHPOST <- RES$BATCHPOST
  if(!is.null(RES$QUOTE))
    QUOTE <- RES$QUOTE
  if(!is.null(RES$ARGQUOTE))
    ARGQUOTE <- RES$ARGQUOTE
  if(!is.null(RES$RUN)){
    #print("isn't null")
    #print(RES$RUN)
    RUN <- RES$RUN
    #print(RUN)
  }
  if(!is.null(RES$MULTIPLIER))
    MULTIPLIER <- RES$MULTIPLIER
  if(!is.null(RES$LOCAL))
    LOCAL <- RES$LOCAL

  #cat("rbatch.parseCommandArgs, RES$RUN=", RES$RUN, "\n")
  #cat("rbatch.parseCommandArgs, RUN=", RUN, "\n")

  return(c(BATCH=BATCH, BATCHPOST=BATCHPOST, QUOTE=QUOTE, ARGQUOTE=ARGQUOTE, RUN=RUN, MULTIPLIER=MULTIPLIER, LOCAL=LOCAL))
}

rbatch.local <- function(BATCH="ALLCORES", BATCHPOST="", QUOTE="", ARGQUOTE='"', RUN=1, MULTIPLIER=1){
  ncores <- 1
  if(is.character(BATCH)){
    if(BATCH == "ALLCORES"){
      if(isWindows()){
        ncores <- as.numeric(Sys.getenv("NUMBER_OF_PROCESSORS"))
      }else{
        require(parallel)
        ncores <- parallel:::detectCores()
      }
    }
  }else
    ncores <- as.integer(BATCH)
  BATCH <- ""
  return(rbatch.parseCommandArgs(BATCH=BATCH, BATCHPOST=BATCHPOST, QUOTE=QUOTE, ARGQUOTE=ARGQUOTE, RUN=RUN, MULTIPLIER=MULTIPLIER, LOCAL=ncores))
}

rbatch.lsf <- function(BATCH="bsub -q normal", BATCHPOST="", QUOTE='"', ARGQUOTE='""', RUN=1, MULTIPLIER=1){
  return(rbatch.parseCommandArgs(BATCH=BATCH, BATCHPOST=BATCHPOST, QUOTE=QUOTE, ARGQUOTE=ARGQUOTE, RUN=RUN, MULTIPLIER=MULTIPLIER))
}

rbatch.mosix <- function(BATCH="nohup mosrun -e -b -q", BATCHPOST=" &", QUOTE="", ARGQUOTE='"', RUN=1, MULTIPLIER=1){
  return(rbatch.parseCommandArgs(BATCH=BATCH, BATCHPOST=BATCHPOST, QUOTE=QUOTE, ARGQUOTE=ARGQUOTE, RUN=RUN, MULTIPLIER=MULTIPLIER))
}


####################################################################
# SET UP THE GLOBAL VARIABLES INTERFACE                            #
####################################################################
rbatch._env <- new.env();
rbatch._env.set <- function( x, value )
  assign( x, value, envir=rbatch._env );
rbatch._env.get <- function( x, mode="any" )
  get( x, envir=rbatch._env, mode=mode, inherits=FALSE );
rbatch._env.set("rbatch.local._queue", c())
rbatch._env.set("rbatch.local._numcores", 1)
## INTERNAL
rbatch.local.pushback <- function(cmdstr, cores){
  #print("rbatch.local")
  #print(cmdstr)
  rbatch._env.set("rbatch.local._queue", c(rbatch._env.get("rbatch.local._queue"), cmdstr))
  rbatch._env.set("rbatch.local._numcores", cores)
  return(invisible())
}

## EXPORTED, forks everything off!
## -- defaults to all cores, but sometimes this is detected wrong...
rbatch.local.run <- function(ncores=NA){
  cmdstrs <- rbatch._env.get("rbatch.local._queue")
  if(length(cmdstrs) < 1) {
    cat("rbatch.local.run: no commands have been batched.\n")
    return()
  }

  if(is.na(ncores))
    ncores <- as.numeric(rbatch._env.get("rbatch.local._numcores"))
  cat("Local, ncores=", ncores, ".\n", sep="")
  if(isWindows()){
    ## Fix up the cmd strs to get the path to R...
    rbin <- paste(Sys.getenv("R_HOME"), "/bin/", sep="")
    rbin <- paste(unlist(strsplit(rbin, "\\", fixed=TRUE)), collapse = "/")
    for(i in 1:length(cmdstrs)){
      if(substring(cmdstrs[i], 1, 2) == "R "){
        cmdstrs[i] <- paste(rbin, cmdstrs[i], sep="")
      }else if(substring(cmdstrs[i], 1, 3) == " R "){
        cmdstrs[i] <- paste(rbin, substring(cmdstrs[i], 2), sep="")
      }
    }
    #print(cmdstrs) ## Good to here!
    parallelRunCommands(cmdstrs, ncores)
  }else{
    ## NOW USE parallel, and go ahead and batch them all off!!!
    #print(cmdstrs) ## Debug
    require(parallel)
    parallel::mclapply(cmdstrs, function(i){system(i); cat(i,"completed.\n");}, mc.cores=ncores)
    #cat("RAN")
  }

  ## Clear out the old queue!!!
  rbatch._env.set("rbatch.local._queue", c())
  return(invisible())
}

rbatch <- function(rfile, seed, ..., rbatch.control=rbatch.default()){
  call <- match.call(expand.dots=TRUE)
  #syscall <- sys.call(sys.parent())

  ##print( batch )
  #print(rbatch.control)
  #stop()

  ## turn an obj into an R expression
  ## currently, only works with vectors, etc.
  pasted <- function( obj, argquote ) {
    #print( eval(obj) )
    #obj <- eval(obj, parent.frame()) ## done externally now..

    if( length(obj) == 1 )
      return( as.character(obj) )

    ret <- paste( argquote, 'c(', sep="" )
    for( i in 1:length(obj) ) {
      ###ret <- paste(ret,obj[i],sep="") ## 8/17/2010 -- vectors of strings having issues
      if(is.character(obj[i])){ ## need to handle character strings a little specially...
        ret = paste(ret, "'", obj[i], "'", sep="")
      }else{ ## same as before
        ret <- paste(ret,obj[i],sep="") ## 8/17/2010 -- vectors of strings having issues
      }
      if( i != length(obj) )
        ret <- paste(ret,",",sep="")
    }
    ret <- paste(ret,')',argquote,sep="")
    return( ret )
  }

  ## first member of call is 'bsub'
  ## second member is 'rfile', with the name of the R file to run (but no need to use)
  ## 'seed' might as well be treated as any other argument, so there is really no need to mess with anything else here...

  for( m in 1:rbatch.control["MULTIPLIER"] ) {
    argstr <- ""
    for( i in 3:length(call) ) {
      callname <- names(call)[i]
      #if(callname!="BATCH" && callname!="QUOTE" && callname!="ARGQUOTE" && callname!="RUN" && callname!="MULTIPLIER" && callname!="BATCHPOST"){
      if(callname!="BATCH" && callname!="QUOTE" && callname!="ARGQUOTE" && callname!="RUN" && callname!="MULTIPLIER" && callname!="BATCHPOST" && callname!="rbatch.control"){
        #wh <- which(callname == names(syscall))
        #if(length(wh) > 0){
        #  ## here's the funky addition
        #  argstr <- paste(argstr, names(call)[i], pasted(syscall[[wh]], argquote=rbatch.control["ARGQUOTE"]))
        #}else{
        #  argstr <- paste(argstr, names(call)[i], pasted(call[[i]], argquote=rbatch.control["ARGQUOTE"]))
        #}
        #print(eval(call[[i]], parent.frame()))
        evalCalli <- eval(call[[i]], parent.frame())
        if(!is.null(evalCalli))
          argstr <- paste(argstr, names(call)[i], pasted(evalCalli, argquote=rbatch.control["ARGQUOTE"]))
      }
    }
    ##cat( argstr, "\n" )

    cmdstr <- paste( rbatch.control["BATCH"], " ", rbatch.control["QUOTE"], "R --vanilla --args ", argstr, " < ", rfile, " > ", rfile, "out", seed, rbatch.control["QUOTE"], rbatch.control["BATCHPOST"], sep="" )
    cat( cmdstr, "\n" )

    #cat("LOCAL", rbatch.control["LOCAL"], "\n") ## DEBUG ONLY
    #cat("RUN", rbatch.control["RUN"], "\n")
    if(rbatch.control["RUN"] == 1){
      #print("inside the local...")
      if(rbatch.control["LOCAL"] == 0){
        system(cmdstr)
      }else{
        rbatch.local.pushback(cmdstr, rbatch.control["LOCAL"])
      }
    }

    seed <- seed + 1
    call[which(names(call) == "seed")] <- seed
  }

  return( seed )
}


## DEBUG
#source("parseCommandArgs.R")
#print(rbatch.default())

#rbatch("test.R", seed=10, test="hello", foo="bar")
#rbatch("test.R", seed=11, a="bye")
#rbatch.local.run()
