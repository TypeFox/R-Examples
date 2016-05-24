# Copyright (c) 2012 by Barb Banbury, University of Tennessee, 
# Update on Foundation API 2013 by Kurt Michels, University of Arizona
# Update to Agave API 2014-2015 by Kurt Michels, University of Arizona
# Removal of Foundation API 2015 by Kurt Michels, University of Arizona
# -- A note on removal.  I only removed any mention of the Foundation API
#    on the help pages.  If one does api="foundation" in the Validate 
#    function, it will work.  But the Foundation API is depracated.  I am
#    keeping the structure of two APIs supported because some day another
#    API will be created, then the developer only needs to replace all
#    of the Foundation urls.
#
# rPlant directly interacts with iplant's command-line API for the 
# Discovery Environment (DE)

# -- AUTHENTICATION FUNCTIONS -- #

# utils::globalVariables(c("rplant.env"))
# myfun <- function(x) print(x)
# environment(myfun) <- as.environment("rplant.env")

#####################
#####################
#### Create_Keys ####
#####################
#####################
rplant.env <- new.env()

MHmakeRandomString <- function(n=1, lenght=12)
{
    randomString <- c(1:n)                  # initialize vector
    for (i in 1:n)
    {
        randomString[i] <- paste(sample(c(0:9, letters, LETTERS),
                                 lenght, replace=TRUE),
                                 collapse="")
    }
    return(randomString)
}

Create_Keys <- function(user, pwd, print.curl=FALSE) {
  # Calls the Agave API, if the user already has a key and secret for rPlant, 
  #   then it fetches them, o/w it creates the keys, and subscribes to the 
  #   correct API's.  This key and secret are required for Validation, all
  #   of this is in the background.
  #
  # Args:
  #   user: Valid iPlant username
  #   pwd: Valid iPlant password, this combo's with the iPlant username
  #
  # Returns:
  #   If invalid credentials an error is shown, o/w it returns the rPlant
  #     key and secret associated with the username and password.
  web <- "https://agave.iplantc.org/clients/v2"
  curl.call <- getCurlHandle(userpwd        = paste(user, pwd, sep=":"), 
                             httpauth       = 1L, 
                             ssl.verifypeer = FALSE)
  client <- MHmakeRandomString()
  res <- tryCatch(expr  = fromJSON(postForm(web, 
                                            clientName  = client,
                                            tier        = "Unlimited", 
                                            description = "", 
                                            callbackUrl = "", 
                                            style       = "POST", 
                                            curl        = curl.call)), 
                  error = function(err) {
                            return(paste(err))
                          }
                  )

  if (print.curl){
    curl.string <- paste0("curl -sku ", user, " -X POST -d clientName=rPlant -d tier=Unlimited -d description='' -d callbackUrl='' ", web)
    print(curl.string)
  }

  Error(res)
  return(list(res$result$consumerKey, res$result$consumerSecret))
}

#####################
#####################
##### Validate ######
#####################
#####################

Validate <- function(user, pwd, api="agave", print.curl=FALSE) {
  # Calls the Agave API.  Used to call both Foundation API or the Agave
  #   API, but Foundation is deprecated.  This function simply validates
  #   a users credentials for the API.
  #
  # Args:
  #   user: Valid iPlant username
  #   pwd: Valid iPlant password, this combo's with the iPlant username
  #   api: "agave"
  #   print.curl: Prints the associated curl statment
  #
  # Returns:
  #   An error if not valid credentials, o/w nothing

# The Foundation API has been depracated, so it is no longer working.  I am
#   only commenting out the Foundation API because some day a new API will
#   probably be created and everything moved.  When that happens, this can
#   be used again.  The next section details all the different parts.

#  rplant.env <<- new.env()

  api <- match.arg(api, c("agave", "foundation"))

  if (api == "foundation"){

    web_BASE <- "https://foundation.iplantcollaborative.org/"
    web <- paste(web_BASE, "auth-v1/", sep="")
    curl.string <- paste("curl -sku '", user, "' ", web, sep="")
    curl.call <- getCurlHandle(userpwd        = paste(user, pwd, sep=":"), 
                               httpauth       = 1L, 
                               ssl.verifypeer = FALSE)
    res <- tryCatch(expr  = fromJSON(getURL(web, curl = curl.call)), 
                    error = function(err) {
                              return(paste(err))
                            }
                    ) 
    Error(res)
    if (res$status == "success"){
#      assign(x     = "rplant.env", 
#             value = new.env(hash = TRUE), 
#             envir = .LocalEnv)
      assign(x     = "api", 
             value = "f", 
             envir = rplant.env)
      assign(x     = "webio",  
             value = paste(web_BASE, "io-v1/io/", user, sep=""), 
             envir = rplant.env)
      assign(x     = "webio1",  
             value = paste(web_BASE, "io-v1/io", sep=""), 
             envir=rplant.env)
     assign(x     = "webcheck",  
             value = paste(web_BASE, "io-v1/io/list/", user, sep=""), 
             envir=rplant.env)
      assign(x     = "weblist",  
             value = paste(web_BASE, "io-v1/io/list", sep=""), 
             envir=rplant.env)
      assign(x     = "webshare",  
             value = paste(web_BASE, "io-v1/io/share/", user, sep=""), 
             envir=rplant.env)
      assign(x     = "webtransform",  
             value = paste(web_BASE, "io-v1/data/transforms/", user, sep=""),  
             envir = rplant.env)
      assign(x     = "webappslist",  
             value = paste(web_BASE, "apps-v1/apps/list", sep=""),  
             envir = rplant.env)
      assign(x     = "webappsname",  
             value = paste(web_BASE, "apps-v1/apps/name", sep=""),  
             envir = rplant.env)
      assign(x     = "webjob",  
             value = paste(web_BASE, "apps-v1/job", sep=""),  
             envir = rplant.env)
      assign(x     = "webjoblist",  
             value = paste(web_BASE, "apps-v1/jobs/list", sep=""),  
             envir=rplant.env)
      assign(x     = "webprofiles",  
             value = paste(web_BASE, "profile-v1/profile/search/username/", user, sep=""),  
             envir = rplant.env)
      assign(x     = "first",  
             value = paste("curl -sku '", user, "'", sep=""), envir=rplant.env)
      assign(x     = "user",  
             value = user,  
             envir = rplant.env)
      assign(x     = "pwd",  
             value = pwd,  
             envir = rplant.env) 
      assign(x     = "curl.call",  
             value = getCurlHandle(userpwd        = paste(user, pwd, sep=":"), 
                                   httpauth       = 1L, 
                                   ssl.verifypeer = FALSE),  
             envir = rplant.env)
    } else {
      return(res$message)
    }
  } else {

    keys <- Create_Keys(user,pwd)
    # Since is the first function using RCurl, I will say in detail 
    #   how it works.

    # First get a base url which will be the url to be called.  The url
    #   can be very complicated, and in rplant.env all of the url's
    #   that are called are stored in there.  It is just a matter of
    #   calling the correct one.
    web_BASE <- "https://agave.iplantc.org/"

    # A couple of the curl statements (PUT, POST) include the following
    #   options.  A GET statement does not have the options.  
    content <- c()
    content[1] <- "grant_type=client_credentials"
    content[2] <- "scope=PRODUCTION"
    content[3] <- paste("username=", user, sep="")
    content[4] <- paste("password=", pwd, sep="")

    # These options are separated by '&'
    string <- paste(content, collapse = "&")

    # Take that string and convert it essentially to integers.  If you look
    #   up ASCII characters they have an associated integer value, this is
    #   what the values are converted to
    val <- charToRaw(string)

    web <- paste(web_BASE, "token", sep="")

    curl.string <- paste("curl -sku '", keys[[1]], ":", keys[[2]], 
                         "' -X POST -d '", string, "' ", web, sep="")

    # For RCurl a curl call must be made.  Essentially this is the validation
    #   part of the curl statement.  For this Validate function on Agave
    #   it is taking the key and secret from the Agave store and using that
    #   for validation.  This type of validation is for this part only,
    #   for most validation on Agave an access token is used.
    curl.call <- getCurlHandle(userpwd        = paste(keys[[1]], keys[[2]],
                                                      sep=":"), 
                               httpauth       = 1L, 
                               ssl.verifypeer = FALSE)
    expire <- as.POSIXlt(format(Sys.time(),"%Y-%m-%d %k:%M:%OS"))
    expire$hour=expire$hour+2
    
    # This is the RCurl statement, getURLContent() is used, it should be noted
    #   that this way is not unique.  Now the RCurl statement of course starts
    #   with the url, then the curl call, then it reads in the value vector.
    #   remember this vector was originally the options in the curl statement
    #   VERY importantly separated by a '&'.  Then the customrequest = "POST"
    #   because this is a POST curl statement.
    res <- tryCatch(expr  = fromJSON(getURLContent(web, 
                                                   curl          = curl.call, 
                                                   infilesize    = length(val), 
                                                   readfunction  = val, 
                                                   upload        = TRUE, 
                                                   customrequest = "POST")),
                    error = function(err) {
                              return(paste(err))
                            }
                    )
  
    if (print.curl){
      print(curl.string)
    }

    Error(res)
    # As I said the RCurl statements are not unique, here is another way to do
    #   that exact statement, using postFORM(), which the RCurl creator, Dr.
    #   Lang said he liked to use.
    #
    #   res <- tryCatch(expr  = fromJSON(postForm(web, 
    #                                             grant_type = "client_credentials",
    #                                             scope      = "PRODUCTION", 
    #                                             username   = user, 
    #                                             password   = pwd, 
    #                                             style      = "POST", 
    #                                             curl       = curl.call)), 
    #                   error = function(err) {
    #                             return(paste(err))
    #                           }
    #                   )

    if (length(res) == 4){
#      assign(x     = "rplant.env",   
#             value = new.env(hash = TRUE),    
#             envir = .LocalEnv)
      assign(x     = "api",   
             value = "a",    
             envir = rplant.env)
      assign(x     = "consumer_key",   
             value = keys[[1]],    
             envir=rplant.env)
      assign(x     = "consumer_secret",   
             value = keys[[2]],    
             envir = rplant.env)
      assign(x     = "webio",   
             value = paste(web_BASE, "files/v2/media/", user, sep=""),    
             envir = rplant.env)
      assign(x     = "webio1",   
             value = paste(web_BASE, "files/v2/media/", sep=""),    
             envir = rplant.env)
      assign(x     = "webcheck",   
             value = paste(web_BASE, "files/v2/listings/", user, sep=""),    
             envir=rplant.env)
      assign(x     = "weblist",   
             value = paste(web_BASE, "files/v2/listings", sep=""),    
             envir = rplant.env)
      assign(x     = "webshare",   
             value = paste(web_BASE, "files/v2/pems/", user, sep=""),    
             envir=rplant.env)
      assign(x     = "webtransform",   
             value = paste(web_BASE, "transforms/v2/", sep=""),    
             envir = rplant.env)
      assign(x     = "webappslist",   
             value = paste(web_BASE, "apps/v2", sep=""),    
             envir=rplant.env)
      assign(x     = "webappsname",   
             value = paste(web_BASE, "apps/v2", sep=""),    
             envir = rplant.env)
      assign(x     = "webjob",   
             value = paste(web_BASE, "jobs/v2", sep=""),    
             envir = rplant.env)
      assign(x     = "webjoblist",   
             value = paste(web_BASE, "jobs/v2", sep=""),    
             envir = rplant.env)
      assign(x     = "webprofiles",   
             value = paste(web_BASE, "profiles/v2/search/username/", 
                           user, sep=""),
             envir=rplant.env)
      assign(x     = "webauth",   
             value = paste(web_BASE, "token", sep=""),    
             envir = rplant.env)
      assign(x     = "first",   
             value = paste("curl -sk -H 'Authorization: Bearer ", 
                           res$access_token, "'", sep=""),    
             envir = rplant.env)
      assign(x     = "user",   
             value = user,    
             envir=rplant.env)
      assign(x     = "pwd",   
             value = pwd,    
             envir = rplant.env) 
      assign(x     = "expire",   
             value = expire,    
             envir = rplant.env) 
      assign(x     = "access_token",   
             value = res$access_token,    
             envir = rplant.env)
      assign(x     = "refresh_token",   
             value = res$refresh_token,    
             envir = rplant.env)
      assign(x     = "curl.call",   
             value = getCurlHandle(httpheader      = c(paste("Authorization: Bearer ", 
                                                             get(x     = "access_token", 
                                                                 envir = rplant.env),
                                                             sep="")), 
                                   httpauth       = 1L, 
                                   ssl.verifypeer = FALSE),   
             envir=rplant.env)
    } else {
      sub <- substring(res$status,1,5)
      if (length(sub) == 0){
        return(stop("API Error, please retry", call. = FALSE))
      } else if (sub == "error"){
        return(stop(res$message, call. = FALSE))
      } else {
        return(stop(res$status, call. = FALSE))
      }
    }
  }
}

#####################
#####################
#### RenewToken #####
#####################
#####################

RenewToken <- function(print.curl=FALSE) {
  # Calls the Agave API.  It simply renews the tokens that were already
  #   acquired.
  #
  # Args:
  #   print.curl: Prints the associated curl statement
  #
  # Returns:
  #   An error if not valid credentials, o/w nothing
  if (rplant.env$api == "a"){
    content <- c()
    content[1] <- "grant_type=refresh_token"
    content[2] <- "scope=PRODUCTION"
    content[3] <- paste("refresh_token=", rplant.env$refresh_token, sep="")
    string <- paste(content, collapse = "&")
    val <- charToRaw(string)

    curl.call <- getCurlHandle(userpwd        = paste(get(x     = "consumer_key", 
                                                          envir = rplant.env), 
                                                      get(x     = "consumer_secret", 
                                                          envir = rplant.env), 
                                                      sep=":"), 
                               httpauth       = 1L, 
                               ssl.verifypeer = FALSE)

    res <- tryCatch(expr  = fromJSON(getURLContent(rplant.env$webauth, 
                                                   curl          = curl.call, 
                                                   infilesize    = length(val), 
                                                   readfunction  = val, 
                                                   upload        = TRUE, 
                                                   customrequest = "POST")), 
                    error = function(err) {
                              return(paste(err))
                           }
                    )

    if (print.curl){
      curl.string <- paste("curl -sku ", rplant.env$consumer_key, ":", 
                           rplant.env$consumer_secret, " -X POST -d '", string,
                           "' ", rplant.env$webauth, sep="")
      print(curl.string)
    }

    Error(res)

    if (length(res) == 4){
      assign(x     = "access_token",    
             value = res$access_token,     
             envir = rplant.env)
      assign(x     = "refresh_token",    
             value = res$refresh_token,     
             envir = rplant.env) 
      assign(x     = "first",    
             value = paste("curl -sk -H 'Authorization: Bearer ", res$access_token, "'", sep=""),
             envir = rplant.env)
      assign(x     = "curl.call",     
             value = getCurlHandle(httpheader     = c(paste("Authorization: Bearer ", 
                                                            get("access_token", envir = rplant.env),
                                                            sep="")), 
                                   httpauth       = 1L, 
                                   ssl.verifypeer = FALSE),     
             envir = rplant.env)
    }
  }
}

#####################
#####################
####### Check #######
#####################
#####################

Check <- function(name, path="", suppress.Warnings=FALSE, 
                  shared.username=NULL, check=FALSE){
  # This takes a file name (or directory name) and it simply checks if that
  #   file (or directory) exist.  If not an error is returned. 
  #
  # Args:
  #   name: name of file (or directory) to be checked
  #   path: path to where file (or directory) is
  #   suppress.Warnings: Either TRUE or FALSE, if TRUE check will be skipped
  #   shared.username: A string of the username.  If there then the file
  #     (or directory) is in that users directory.
  #   check: Either TRUE or FALSE, if TRUE then check that object exists in
  #     directory, if FALSE check that object DOES NOT exist in directory.
  #
  # Returns:
  #   Returns an error if something about the path or name is incorrect.
  #     o/w returns nothing if file (or directory) does exist
  Time()
  Renew()
  if (suppress.Warnings == FALSE){
    if (is.null(shared.username)){# Not a shared user
      dir.exist <- fromJSON(getURL(url  = paste(rplant.env$webcheck, path, sep="/"), 
                                   curl = rplant.env$curl.call)) 
      if (length(dir.exist$result) != 0){# Path does exist, check object
        if (path==""){
          obj.exist <- fromJSON(getURL(url  = paste(rplant.env$webcheck, name, sep="/"), 
                                       curl = rplant.env$curl.call))
        } else {
          obj.exist <- fromJSON(getURL(url  = paste(rplant.env$webcheck, path, name, sep="/"),
                                       curl = rplant.env$curl.call))
        }
      } else {
        if (dir.exist$status == "error"){# Path does not exist, show appropriate error
          if ((dir.exist$message == "File does not exist") || 
              (dir.exist$message == "File/folder does not exist")){
            return(stop(paste("path '", path, "' not proper directory", sep=""),
                        call. = FALSE))
          } else {
            return(stop("improper username/password combination", call. = FALSE))
          }
        } else {# If no error, then no directory
          return(stop(paste("path '", path, "' not proper directory", sep=""), 
                      call. = FALSE))
        }
      }
    } else {# Shared username, get proper path and simply check it
      if (path == ""){
        web <- paste(rplant.env$weblist, shared.username, name, sep="/")
      } else {
        web <- paste(rplant.env$weblist, shared.username, path, name, sep="/")
      }
      obj.exist <- fromJSON(getURL(web, curl=rplant.env$curl.call))
    }
    # Check whether object exists or not
    if (check){# If check=TRUE and object IS in directory return error
      if (length(obj.exist$result) != 0){
        return(stop(paste("object '", name, "' already exists in '", path,
                          "' directory", sep=""), call. = FALSE))
      }
    } else {# If check=FALSE and object IS NOT in directory return error
      if (length(obj.exist$result) == 0){
        return(stop(paste("object '", name, "' doesn't exist in '", path, 
                          "' directory", sep=""), call. = FALSE))
      }
    }
  }
}

#####################
#####################
####### Wait ########
#####################
#####################

Wait <- function(job.id, minWaitsec, maxWaitsec, print=FALSE){
  # This function simply waits for the job to finish before proceeding.  It is
  #   used when result files from a job must be retrieved in order to do the
  #   next job.  It simply calls the API and checks the job status.  Once status
  #   is finished then it proceeds.
  #
  # Args:
  #   job.id: job id of job to be checked
  #   minWaitsec: The min wait time in seconds.
  #   maxWaitsec: The max wait time in seconds.  The job polls, and this is
  #     the maximum time you want to wait between polls.
  #   print: Prints job number and current status.
  #
  # Returns:
  #   Returns nothing unless printing
  currentStatus= ''
  currentWait = minWaitsec
  if (rplant.env$api == 'f'){
    # For the Foundation API the job isn't done until 'ARCHIVING_FINISHED'
    while (( currentStatus != 'FAILED' ) && (currentStatus != 'ARCHIVING_FINISHED')) {
      # cache the status from previous inquiry
      oldStatus = currentStatus
      currentStatus = CheckJobStatus( job.id )

      if (currentStatus == oldStatus) {  # Status hasn't changed from last time we asked
        currentWait = currentWait * 1.10 #   so wait 10% longer to poll in the future

        if (currentWait > maxWaitsec) {
          currentWait = maxWaitsec       #   but don't wait too long
        }
      } else {
        currentWait = minWaitsec # status changed so reset wait counter to min value
      }
      # Sit idle for proscribed time. If you are using an event-based programming 
      #   model, you could just schedule the next check currentWait sec in the future 
      Sys.sleep(currentWait) 
    }
  } else {
    # For the Agave API the job isn't done until 'FINISHED'
    while (( currentStatus != 'FAILED' ) && (currentStatus != 'FINISHED')) {
      oldStatus = currentStatus
      currentStatus = CheckJobStatus( job.id )

      if (currentStatus == oldStatus) {
        currentWait = currentWait * 1.10
        if (currentWait > maxWaitsec) {
          currentWait = maxWaitsec
        }
      } else {
        currentWait = minWaitsec
      }
      Sys.sleep(currentWait) 
    }
  }

  if (print == TRUE) {
    message(paste("Job number: '", job.id, "' has status: ", currentStatus, sep=""))
  }
}

#####################
#####################
####### Misc. #######
#####################
#####################

#############
### Renew ###
#############

Renew <- function(){
  # This is called before every call.  It simply refreshes the curl call
  #
  # Returns:
  #   Nothing
  if (rplant.env$api == "a") {
    # The type of curl call for RCurl in the Agave API uses the Bearer access
    #   token.  Because of this the curl call is slightly different from the
    #   Foundation API
    assign(x     = "curl.call", 
           value = getCurlHandle(httpheader     = c(paste("Authorization: Bearer ", 
                                                          get("access_token", envir=rplant.env), 
                                                          sep="")), 
                                 httpauth       = 1L, 
                                 ssl.verifypeer = FALSE), 
           envir = rplant.env)
  } else {
    # The curl call for RCurl on the Foundation API simply uses username and
    #   password.
    assign(x     = "curl.call",
           value = getCurlHandle(userpwd        = paste(get("user", envir = rplant.env), 
                                                        get("pwd", envir = rplant.env), 
                                                        sep=":"), 
                                 httpauth       = 1L, 
                                 ssl.verifypeer = FALSE), 
           envir = rplant.env)
  }
}

#############
### Time ####
#############

Time <- function(){
  # For the Agave API the access token expires after 2 hours.  This function
  #   is called before every curl call.  The time is only kept track of in
  #   the R workspace, and if the token expires then it is renewed.
  #
  # Returns:
  #   Nothing
  if (rplant.env$api != "f"){
    compare <- as.POSIXlt(format(Sys.time(),"%Y-%m-%d %k:%M:%OS"))
    if (compare > rplant.env$expire){ # If it does expire
      expire <- as.POSIXlt(format(Sys.time(),"%Y-%m-%d %k:%M:%OS"))
      expire$hour=expire$hour+2 # insert a new expire time
      assign("expire", expire, envir=rplant.env)
      RenewToken() # Renew the token
    }
  }
}

#############
## TestApp ##
#############

TestApp <- function(APP){
  # This application takes the application name, and returns a short description
  #   of the application
  #
  # Args:
  #   APP: name of application
  #
  # Returns:
  #   Short description of application
  if (rplant.env$api == "f"){
    first_string <- "res$result[[1]]"
      if (substring(APP,nchar(APP)-1,nchar(APP)-1) == "u"){
      priv.APP <- substring(APP,1,nchar(APP)-2)
    } else if (substring(APP,nchar(APP)-2,nchar(APP)-2) == "u"){
      priv.APP <- substring(APP,1,nchar(APP)-3)
    } else {
      priv.APP <- APP
    }
  } else {
    first_string <- "res$result"
    priv.APP <- APP
  }

  Renew()
  res <- tryCatch(expr  = fromJSON(getForm(uri          = paste(rplant.env$webappsname,
                                                                priv.APP, sep="/"),
                                           .checkparams = FALSE, 
                                           curl         = rplant.env$curl.call)), 
                  error = function(err) {
                            return(paste(err))
                          }
                  )

  if (length(res) == 1){
    return(list(NULL))
  } else {
    shortd <- eval(parse(text=paste(first_string, "$shortDescription", sep="")))
    shortn <- nchar(shortd)
    longd <- eval(parse(text=paste(first_string, "$longDescription", sep="")))
    if (is.null(longd)) {longn = 0} else {longn <- nchar(longd)}
    if (longn >= shortn) {
      description <- longd
    } else {
      description <- shortd
    }
    return(c(eval(parse(text=paste(first_string, "$id", sep=""))), description))
  }
}

#############
### Error ###
#############

Error <- function(ERR){
  # This is the error checking component of the package.  It takes in an 
  #   object.  If the object is a string then it is most likely an error, 
  #   if it's not a string then the status of the of the object needs to be 
  #   checked to make sure there were no errors.  If an error did occur then
  #   an appropriate error is returned.
  #
  # Args:
  #   ERR: object (could be anything)
  #
  # Returns:
  #   Nothing if there is no error, o/w it returns an appropriate error
  if (length(ERR) == 1){
    sub1 <- substring(ERR,8,8)
    if (sub1 == "B"){
      return(stop("Bad Request", call. = FALSE))
    } else if (sub1 == "U"){
      return(stop("Invalid username/password combination", call. = FALSE))
    } else if ((sub1 == "F") || (sub1 == "N")){
      return(stop("file or directory or job id does not exist", call. = FALSE))
    } else {
      len <- nchar(ERR)
      return(stop(substring(ERR,8,len-3), call. = FALSE))
    }
  } else {
    for (i in 1:length(ERR)){
      if (names(ERR)[i] == "status"){
        if (ERR$status == "error"){
          return(stop(ERR$message, call. = FALSE))
        }
        break;
      }
    }
  }
}

#############
## appINFO ##
#############

appINFO <- function(application, dep=FALSE, input=FALSE){
  # This is the error checking component of the package.  It takes in an 
  #   object.  If the object is a string then it is most likely an error, 
  #   if it's not a string then the status of the of the object needs to be 
  #   checked to make sure there were no errors.  If an error did occur then
  #   an appropriate error is returned.
  #
  # Args:
  #   application: application name (string)
  #   dep: Either TRUE or FALSE indicating to check if application is
  #     depracated.  If application is deprecated then an error is sent.
  #   input: Either TRUE or FALSE, if TRUE then include application info
  #
  # Returns:
  #   Returns different information depending on the inputs.  All information
  #     is about the application.
  Time()
  Renew()
  # For the Foundation API and Agave API, naming schemes are a litle bit different.
  #   This needs to be accounted for on every function.
  if (rplant.env$api == "f"){
    tmp_string <- "tmp$result[[len]]"
    tmp_str <- "$public"
    if (substring(application,nchar(application)-1,nchar(application)-1) == "u"){
      priv.APP <- substring(application,1,nchar(application)-2)
    } else if (substring(application,nchar(application)-2,nchar(application)-2) == "u"){
      priv.APP <- substring(application,1,nchar(application)-3)
    } else {
      priv.APP <- application
    }
  } else {
    tmp_string <- "tmp$result"
    tmp_str <- "$isPublic"
    priv.APP <- application
  }
  # This part depends on the application name ending in "u1" etc.   If the 
  #   naming scheme for the public applications still does this, which is
  #   true for current (2014) Agave and Foundation API.  Then this simply
  #   takes that part off.  So the name 'Muscleu2' becomes 'Muscle'
  if (substring(application,nchar(application)-1,nchar(application)-1) == "u"){
    version <- as.numeric(substring(application, nchar(application), nchar(application)))
    text <- "Public App"
  } else if (substring(application, nchar(application)-2, nchar(application)-2) == "u"){
    version <- as.numeric(paste(substring(application, nchar(application)-1, nchar(application)-1),
                                substring(application, nchar(application), nchar(application)),
                                sep=""))
    text <- "Public App"
  } else {
    text <- "Private App"
  }

  tmp <- tryCatch(expr  = fromJSON(getForm(uri          = paste(rplant.env$webappsname,
                                                                priv.APP, sep="/"),
                                           .checkparams = FALSE, 
                                           curl         = rplant.env$curl.call)), 
                  error = function(err) {
                            return(paste(err))
                          }
                  )
  Error(tmp)

  len <- length(tmp$result)
  if (eval(parse(text=paste(tmp_string, tmp_str, sep=""))) == FALSE) {
    text <- "Private App"
  } else if (length(tmp) == 0) {
    return(stop("No information on application: not valid", call. = FALSE))
  }
  # This depends on the naming scheme, for 'Muscleu2' the version number is 'u2'
  if (text == "Public App"){
    APP <- eval(parse(text=paste(tmp_string, "$id", sep="")))
    if (substring(APP,nchar(APP)-1,nchar(APP)-1) == "u"){
      priv.APP <- substring(APP,1,nchar(APP)-2)
      version.APP <- as.numeric(substring(APP,nchar(APP),nchar(APP)))
    } else if (substring(APP,nchar(APP)-2,nchar(APP)-2) == "u"){
      priv.APP <- substring(APP,1,nchar(APP)-3)
      version.APP <- as.numeric(paste(substring(APP, nchar(APP)-1, nchar(APP)-1),
                                      substring(APP,nchar(APP),nchar(APP)),
                                      sep=""))
    }
    # When the application is looked up under the name 'Muscle', it finds the
    #   newest version.  So 'Muscleu2' is compared to the most recent version
    #   which could be 'u3', if this is so the application is deprecated and
    #   you are told so.
    if (version.APP > version){
      v.text <- paste("Deprecated, the newest version is:", APP)
    } else {
      v.text <- "Newest Version"
    }
    # When submitting a job use dep=TRUE, that way if application is deprecated
    #   the job will not be submitted because an error is returned.
    if (dep){
      if (substring(v.text, 1, 1) == "D"){
        return(stop(paste("Application deprecated, should be:", APP), call. = FALSE))
      }
    }
  }
  # This finds if the application can be parallelized, and it is returned
  set <- eval(parse(text=paste(tmp_string, "$parallelism", sep="")))
  # Don't return verbose output
  if (!input){
    if (text == "Private App"){
      return(list("Private App", priv.APP, tmp, set))
    } else {
      return(list("Public App", priv.APP, v.text, APP, tmp, set))
    }
  } else { # Return verbose output
    app.info<-c()
    for (input in sequence(length(eval(parse(text=paste(tmp_string, "$inputs", sep="")))))) {
      app.info <- rbind(app.info, eval(parse(text=paste(tmp_string, "$inputs[[input]]$id", 
                                                        sep=""))))
    }
    if (text == "Private App"){
      return(list("Private App", priv.APP, tmp, app.info, set))
    } else {
      return(list("Public App", priv.APP, v.text, APP, tmp, app.info, set))
    }
  }
}

# -- END -- #




# -- MAIN FUNCTIONS -- #

# These functions are duplicate for both the File and Directory functions.
#   Therefore they are written as a single, then wrappers are used for
#   the File and Dir functions.

#####################
#####################
###### Rename #######
#####################
#####################

Rename <- function(name, new.name, path="", print.curl=FALSE, suppress.Warnings=FALSE) {
  # This function simply takes the object 'name' and renames it to 'new.name'
  #
  # Args:
  #   name: Current name of object
  #   new.name: New name of object
  #   path: Path to where object is
  #   print.curl: Prints the associated curl statement
  #   suppress.Warnings: Don't do any error checking (faster)
  #
  # Returns:
  #   Returns nothing unless an error
  content <- c()
  if (rplant.env$api == "f") {
    content[1] <- "action=rename"
    content[2] <- paste("newName=", new.name, sep="")
  } else {
    if (path == ""){
      content[1] <- paste("path=", rplant.env$user, new.name, sep="/")
    } else {
      content[1] <- paste("path=", rplant.env$user, path, new.name, sep="/")
    }
    content[2] <- "action=move"
  }

  if (path == ""){
    web <- paste(rplant.env$webio, name, sep="/")
  } else {
    web <- paste(rplant.env$webio, path, name, sep="/")
  }

  if (print.curl){
    curl.string <- paste(rplant.env$first, " -X PUT -d '", 
                         paste(content, collapse = "&"), "' ", web, sep="")
    print(curl.string)
  }

  val <- charToRaw(paste(content, collapse = "&"))
  Renew()
  res <- tryCatch(expr  = fromJSON(httpPUT(url     = web, 
                                           content = val, 
                                           curl    = rplant.env$curl.call)),
                  error = function(err) {
                            return(paste(err))
                          }
                  )
  if (!suppress.Warnings){Error(res)}
}

#####################
#####################
####### Copy ########
#####################
#####################

Copy <- function(name, org.path="", end.path="", print.curl=FALSE, suppress.Warnings=FALSE) {
  # This function simply moves the object 'name' from 'org.path'
  #   to 'end.path'
  #
  # Args:
  #   name: Name of object
  #   org.path: Original or current path where object is
  #   end.path: Path to where object will be moved
  #   print.curl: Prints the associated curl statement
  #   suppress.Warnings: Don't do any error checking (faster)
  #
  # Returns:
  #   Returns nothing unless an error
  if (rplant.env$api == "f") {
    path <- "newPath="
  } else {
    path <- "path="
  }

  content <- c()
  if (end.path == ""){
    content[1] <- paste(path, rplant.env$user, name, sep="/")
  } else {
    content[1] <- paste(path, rplant.env$user, end.path, name, sep="/")
  }
  content[2] <- "action=copy"

  if (org.path == ""){
    web <- paste(rplant.env$webio, name, sep="/")
  } else {
    web <- paste(rplant.env$webio, org.path, name, sep="/")
  }

  if (print.curl){
    curl.string <- paste(rplant.env$first, " -X PUT -d '", 
                         paste(content, collapse = "&"), "' ", 
                         web, sep="")
    print(curl.string)
  }

  val <- charToRaw(paste(content, collapse = "&"))
  Renew()
  res <- tryCatch(expr  = fromJSON(httpPUT(url     = web, 
                                           content = val, 
                                           curl    = rplant.env$curl.call)),
                  error = function(err) {
                            return(paste(err))
                          }
                  )
  if (!suppress.Warnings){Error(res)}
}

#####################
#####################
####### Move ########
#####################
#####################

Move <- function(name, org.path="", end.path="", print.curl=FALSE, suppress.Warnings=FALSE) {
  # This function simply moves the object 'name' from 'org.path'
  #   to 'end.path'
  #
  # Args:
  #   name: Name of object
  #   org.path: Original or current path where object is
  #   end.path: Path to where object will be moved
  #   print.curl: Prints the associated curl statement
  #   suppress.Warnings: Don't do any error checking (faster)
  #
  # Returns:
  #   Returns nothing unless an error
  if (rplant.env$api == "f") {
    path <- "newPath="
  } else {
    path <- "path="
  }

  content <- c()
  if (end.path == ""){
    content[1] <- paste(path, rplant.env$user, name, sep="/")
  } else {
    content[1] <- paste(path, rplant.env$user, end.path, name, sep="/")
  }
  content[2] <- "action=move"

  if (org.path == ""){
    web <- paste(rplant.env$webio, name, sep="/")
  } else {
    web <- paste(rplant.env$webio, org.path, name, sep="/")
  }

  if (print.curl){
    curl.string <- paste(rplant.env$first, " -X PUT -d '", 
                         paste(content, collapse = "&"), "' ", 
                         rplant.env$webio, sep="")
    print(curl.string)
  }

  val <- charToRaw(paste(content, collapse = "&"))
  Renew()
  res <- tryCatch(expr  = fromJSON(httpPUT(url     = web, 
                                           content = val, 
                                           curl    = rplant.env$curl.call)),
                  error = function(err) {
                            return(paste(err))
                          }
                  )
  if (!suppress.Warnings){Error(res)}
}

#####################
#####################
###### Delete #######
#####################
#####################

Delete <- function(name, path="", print.curl=FALSE, suppress.Warnings=FALSE) {
  # This function removes the object 'name' from 'path'
  #
  # Args:
  #   name: Name of object
  #   path: Path to current object
  #   print.curl: Prints the associated curl statement
  #   suppress.Warnings: Don't do any error checking (faster)
  #
  # Returns:
  #   Returns nothing unless an error
  if (path == "") {
    web <- paste(rplant.env$webio, name, sep="/")
  } else {
    web <- paste(rplant.env$webio, path, name, sep="/")
  }

  if (print.curl) {
    curl.string <- paste(rplant.env$first, " -X DELETE ", web, sep="")
    print(curl.string)
  }
  Renew()
  res <- tryCatch(expr  = fromJSON(httpDELETE(web, curl = rplant.env$curl.call)), 
                  error = function(err) {
                            return(paste(err))
                          }
                  )
  if (!suppress.Warnings){Error(res)}
}

#####################
#####################
####### Share #######
#####################
#####################

Share <- function(name, path="", shared.username, read=TRUE, execute=TRUE, 
                  write=TRUE, print.curl=FALSE, suppress.Warnings=FALSE, D=FALSE) {
  # This function shares the object 'name' with shared.username.  Also one can
  #   decide which permissions to give to the shared user.
  #
  # Args:
  #   name: Name of object
  #   path: Current path where object is
  #   shared.username: String, valid iPlant username with whom the object
  #     is being shared.
  #   read: Gives read permissions to object
  #   execute: Gives execute permissions to object
  #   write: Gives write permissions to object
  #   D: Either TRUE or FALSE, if TRUE then object is a directory o/w
  #     the object is a file.
  #   print.curl: Prints the associated curl statement
  #   suppress.Warnings: Don't do any error checking (faster)
  #
  # Returns:
  #   Returns nothing unless an error
  content <- c()

  if (rplant.env$api == "f") {
    if (read == TRUE) {content <- append(content, "can_read=true")} 
    if (execute == TRUE) {content <- append(content, "can_execute=true")}
    if (write == TRUE) {content <- append(content, "can_write=true")}
    if ((read == FALSE) && (execute == FALSE) && (write == FALSE)) {
      return(stop("Must select some permissions", call. = FALSE))
    }
  } else {
    if ((read == TRUE) && (execute == TRUE) && (write == TRUE)) {
      content[1] <- "permission=all"
    } else if ((read == TRUE) && (execute == TRUE) && (write == FALSE)) {
      content[1] <- "permission=read_execute"
    } else if ((read == TRUE) && (execute == FALSE) && (write == TRUE)) {
      content[1] <- "permission=read_write"
    } else if ((read == TRUE) && (execute == FALSE) && (write == FALSE)) {
      content[1] <- "permission=read"
    } else if ((read == FALSE) && (execute == TRUE) && (write == TRUE)) {
      content[1] <- "permission=write_execute"
    } else if ((read == FALSE) && (execute == TRUE) && (write == FALSE)) {
      content[1] <- "permission=execute"
    } else if ((read == FALSE) && (execute == FALSE) && (write == TRUE)) {
      content[1] <- "permission=write"
    } else {
      return(stop("Must select some permissions", call. = FALSE))
    }
  }

  content <- append(content,paste("username=", shared.username, sep=""))

  if(D) { # Directory, so add recursive=true so all contents in the directory have same perms
    content <- append(content,"recursive=true")
  }

  if (path == "") {
    web <- paste(rplant.env$webshare, name, sep="/")
  } else {
    web <- paste(rplant.env$webshare, path, name, sep="/")
  }

  if (print.curl){
    curl.string <- paste(rplant.env$first," -X POST -d '", 
                         paste(content, collapse = "&"), "' ", web, sep="")
    print(curl.string)
  }

  val <- charToRaw(paste(content, collapse = "&"))
  Renew()
  res <- tryCatch(expr  = fromJSON(getURLContent(web, 
                                                 curl          = rplant.env$curl.call, 
                                                 infilesize    = length(val), 
                                                 readfunction  = val, 
                                                 upload        = TRUE, 
                                                 customrequest = "POST")), 
                  error = function(err) {
                            return(paste(err))
                          }
                  )
  if (!suppress.Warnings){Error(res)}
}

#####################
#####################
#### Permissions ####
#####################
#####################

Pems <- function(name, path="", print.curl=FALSE, suppress.Warnings=FALSE) {
  # This function looks at the permissions on an object.  It will return
  #   the object name, users with whom the object is shared and their
  #   permissions.
  #
  # Args:
  #   name: Name of object
  #   path: Current path where object is
  #   print.curl: Prints the associated curl statement
  #   suppress.Warnings: Don't do any error checking (faster)
  #
  # Returns:
  #   Returns the object name, users with whom the object is shared and
  #     their permissions.
  if (rplant.env$api == "f"){
    tmp_string <- "tmp$result$permissions"
  } else {
    tmp_string <- "tmp$result"
  }
    
  if (path == ""){
    web <- paste(rplant.env$webshare, name, sep="/")
  } else {
    web <- paste(rplant.env$webshare, path, name, sep="/")
  }

  if (print.curl){
    curl.string <- paste(rplant.env$first, web)
    print(curl.string)
  }

  Renew()
  tmp <- tryCatch(expr  = fromJSON(getURL(web, curl=rplant.env$curl.call)), 
                  error = function(err) {
                            return(paste(err))
                          }
                  )
  if (!suppress.Warnings){Error(tmp)}

  if (rplant.env$api == "a"){
    used <- c()
    len <- length(eval(parse(text=tmp_string))) - 1
    first <- 2 # We don't want to return the user themself, so start at 2
    total <- first + len - 1
  } else {
    used <- c("you", "admin_proxy", "ipcservices", "rodsBoot", "QuickShare",
              "ibp-proxy", "ipcservices", "ipc_admin", "admin2", 
              "proxy-de-tools", "de-irods", "rodsadmin")
    len <- 0
    total <- length(eval(parse(text=tmp_string))) - 1
    first <- 2 # The other users start at position 11
    for (i in first:total) { # Check permissions
      if (!eval(parse(text=paste(tmp_string, "[[", i, "]]$username", sep=""))) %in% used) {
        len = len + 1
      }
    }

  }
  
  if (len == 0){# If the object is not shared, still return something
    res <- matrix(, len + 1, 3) 
  } else { 
    res <- matrix(, len, 3) 
  }
  colnames(res) <- c("Name", "Username", "Permissions")
  res[1, 1] <- name
  if (len == 0){# If the object is not shared, return "None"
    res[1, 2] <- "None"
    res[1, 3] <- "None"
  } else {
    cnt = 1
    for (i in first:total) { # Check permissions
      if (!eval(parse(text=paste(tmp_string, "[[", i, "]]$username", sep=""))) %in% used) {
        if (cnt != 1){res[cnt,1] <- ""}
        res[cnt, 2] <- eval(parse(text=paste(tmp_string, "[[", i, "]]$username", 
                                             sep="")))
        if (eval(parse(text=paste(tmp_string, "[[", i, "]]$permission$read", 
                                  sep=""))) == TRUE) {
          R <- TRUE
        } else {
          R <- FALSE
        }
        if (eval(parse(text=paste(tmp_string, "[[", i, "]]$permission$write", 
                                  sep=""))) == TRUE) {
          W <- TRUE        
        } else {
          W <- FALSE
        }
        if (eval(parse(text=paste(tmp_string, "[[", i, "]]$permission$execute", 
                                  sep=""))) == TRUE) {
          E <- TRUE        
        } else {
          E <- FALSE
        }

        if ((R == TRUE) && (E == TRUE) && (W == TRUE)) {
          str <- "All"
        } else if ((R == TRUE) && (E == TRUE) && (W == FALSE)) {
          str <- "R/E"
        } else if ((R == TRUE) && (E == FALSE) && (W == TRUE)) {
          str <- "R/W"
        } else if ((R == TRUE) && (E == FALSE) && (W == FALSE)) {
          str <- "R"
        } else if ((R == FALSE) && (E == TRUE) && (W == TRUE)) {
          str <- "W/E"
        } else if ((R == FALSE) && (E == TRUE) && (W == FALSE)) {
          str <- "E"
        } else if ((R == FALSE) && (E == FALSE) && (W == TRUE)) {
          str <- "W"
        }

        res[cnt, 3] <- str
        cnt = cnt + 1
      }
    }
  }
  return(res)
}

# -- END -- #




# -- FILE FUNCTIONS -- #

#####################
#####################
#### UploadFile #####
#####################
#####################

UploadFile <- function(local.file.name, local.file.path="", filetype=NULL,
                       print.curl=FALSE, suppress.Warnings=FALSE) {
  # This function takes a file in a users local directory and it
  #   uploads the file onto iPlant's servers.
  #
  # Args:
  #   local.file.name: Name of file
  #   local.file.path: Current path where object is on the local directory
  #   filetype: Not required, but we can assign the file a file type
  #   print.curl: Prints the associated curl statement
  #   suppress.Warnings: Don't do any error checking (faster)
  #
  # Returns:
  #   Returns nothing unless there is an error.  Errors include, file already
  #     on iPlant servers.
  if (local.file.path == ""){
    file.path = paste(getwd(), local.file.name, sep="/")
  } else {
    file.path = paste(local.file.path, local.file.name, sep="/")
  }

  if (rplant.env$api == "f") {
    options <- list(userpwd        = paste(rplant.env$user, rplant.env$pwd, sep=":"), 
                    ssl.verifypeer = FALSE, 
                    httpauth       = AUTH_BASIC, 
                    useragent      = "R", 
                    followlocation = TRUE)
  } else {
    options <- list(httpheader=c(paste("Authorization: Bearer ", rplant.env$access_token, sep="")), 
                    ssl.verifypeer = FALSE, 
                    httpauth       = AUTH_BASIC, 
                    useragent      = "R", 
                    followlocation = TRUE)
  }
  # Check that file is not in iPlant directory
  Check(local.file.name, suppress.Warnings=suppress.Warnings, check=TRUE)

  if (!is.null(filetype)){
    res <- tryCatch(expr  = fromJSON(postForm(rplant.env$webio, 
                                              style        = "httppost", 
                                              fileToUpload = fileUpload(file.path), 
                                              fileType     = filetype, 
                                              .opts        = options)), 
                    error = function(err) {
                              return(paste(err))
                            }
                    )
    if (!suppress.Warnings){Error(res)}
    curl.string <- paste(rplant.env$first, " -F 'fileToUpload=@", file.path, 
                         "' -F 'fileType=", filetype, "' ", rplant.env$webio, 
                         sep="")
  } else {
    res <- tryCatch(expr  = fromJSON(postForm(rplant.env$webio, 
                                              style        = "httppost", 
                                              fileToUpload = fileUpload(file.path),
                                              .opts        = options)), 
                    error = function(err) {
                              return(paste(err))
                            }
                    )
    if (print.curl==TRUE){
      curl.string <- paste(rplant.env$first," -F 'fileToUpload=@", file.path, 
                           "' ", rplant.env$webio, sep="")
      print(curl.string)
    }

    if (!suppress.Warnings){Error(res)}

  }
}

#####################
#####################
##### ShareFile #####
#####################
#####################

ShareFile <- function(file.name, file.path="", shared.username, read=TRUE, 
                      execute=TRUE, write=TRUE, print.curl=FALSE, 
                      suppress.Warnings=FALSE) {

  # This function shares the 'file.name' with shared.username.  Also one can
  #   decide which permissions to give to the shared user.
  #
  # Args:
  #   file.name: Name of file
  #   file.path: Current path where file is
  #   shared.username: String, valid iPlant username with whom the object
  #     is being shared.
  #   read: Gives read permissions to file
  #   execute: Gives execute permissions to file
  #   write: Gives write permissions to file
  #   print.curl: Prints the associated curl statement
  #   suppress.Warnings: Don't do any error checking (faster)
  #
  # Returns:
  #   Returns nothing unless an error, file does not exist

  # Check 'file.name'
  Check(file.name, file.path, suppress.Warnings)

  Share(file.name, file.path, shared.username, read, execute, write, print.curl)
}

#####################
#####################
## PermissionsFile ##
#####################
#####################

PermissionsFile <- function(file.name, file.path="", print.curl=FALSE, 
                            suppress.Warnings=FALSE) {

  # This function looks at the permissions on a file.name.  It will return
  #   the files name, users with whom the file is shared and their
  #   permissions.
  #
  # Args:
  #   file.name: Name of file
  #   file.path: Current path where fiile is
  #   print.curl: Prints the associated curl statement
  #   suppress.Warnings: Don't do any error checking (faster)
  #
  # Returns:
  #   Returns the file name, users with whom the object is shared and
  #     their permissions o/w an error if file does not exist
    Check(file.name, file.path, suppress.Warnings)
    
    Pems(file.name, file.path, print.curl, suppress.Warnings)
}

#####################
#####################
#### RenameFile #####
#####################
#####################

RenameFile <- function(file.name, new.file.name, file.path="",
                       print.curl=FALSE, suppress.Warnings=FALSE) {
  # This function simply takes 'file.name' and renames it to 'new.file.name'
  #
  # Args:
  #   file.name: Current name of file
  #   new.file.name: New name of file
  #   file.path: Path to where file is
  #   print.curl: Prints the associated curl statement
  #   suppress.Warnings: Don't do any error checking (faster)
  #
  # Returns:
  #   Returns nothing unless an error if file does not exist, or if
  #     'new.file.name' does exist in 'file.path'.

  # Check file.name
  Check(file.name, file.path, suppress.Warnings)

  # Check new.file.name
  Check(new.file.name, file.path, suppress.Warnings, check=TRUE)

  Rename(file.name, new.file.name, file.path, print.curl) 
}

#####################
#####################
##### CopyFile ######
#####################
#####################

CopyFile <- function(file.name, file.path="", end.path="", 
                     print.curl=FALSE, suppress.Warnings=FALSE) {
  # This function simply moves the 'file.name' from 'file.path'
  #   to 'end.path'
  #
  # Args:
  #   file.name: Name of file
  #   file.path: Original or current path where file is
  #   end.path: Path to where file will be moved
  #   print.curl: Prints the associated curl statement
  #   suppress.Warnings: Don't do any error checking (faster)
  #
  # Returns:
  #   Returns nothing unless an error.  Error if 'file.name'
  #     does not exist or if 'file.name' in 'end.path' already
  #     does exist

  # Check 'file.name'
  Check(file.name, file.path, suppress.Warnings)

  # Check 'file.name' in 'end.path'
  Check(file.name, end.path, suppress.Warnings, check=TRUE)

  Copy(file.name, file.path, end.path, print.curl)
}

#####################
#####################
##### MoveFile ######
#####################
#####################

MoveFile <- function(file.name, file.path="", end.path="", 
                     print.curl=FALSE, suppress.Warnings=FALSE) {
  # This function simply moves the 'file.name' from 'file.path'
  #   to 'end.path'
  #
  # Args:
  #   file.name: Name of file
  #   file.path: Original or current path where file is
  #   end.path: Path to where file will be moved
  #   print.curl: Prints the associated curl statement
  #   suppress.Warnings: Don't do any error checking (faster)
  #
  # Returns:
  #   Returns nothing unless an error.  Error if 'file.name'
  #     does not exist or if 'file.name' in 'end.path' already
  #     does exist

  # Check 'file.name'
  Check(file.name, file.path, suppress.Warnings)

  # Check 'file.name' in 'end.path'
  Check(file.name, end.path, suppress.Warnings, check=TRUE)

  Move(file.name, file.path, end.path, print.curl)
}

#####################
#####################
#### DeleteFile #####
#####################
#####################

DeleteFile <- function(file.name, file.path="", print.curl=FALSE, 
                       suppress.Warnings=FALSE) {
  # This function removes the 'file.name' from 'file.path'
  #
  # Args:
  #   name: Name of file
  #   path: Path to current file
  #   print.curl: Prints the associated curl statement
  #   suppress.Warnings: Don't do any error checking (faster)
  #
  # Returns:
  #   Returns nothing unless an error if 'file.name' does not exist

  # Check 'file.name'
  Check(file.name, file.path, suppress.Warnings)

  Delete(file.name, file.path, print.curl)
}

#####################
#####################
#### SupportFile ####
#####################
#####################

SupportFile <- function(print.curl=FALSE, suppress.Warnings=FALSE) {  
  # This function lists all supported file types on either the
  #   Foundation API or Agave API.
  #
  # Args:
  #   print.curl: Prints the associated curl statement
  #   suppress.Warnings: Don't do any error checking (faster)
  #
  # Returns:
  #   Returns all supported file types.
  Time()
  Renew()
  res <- tryCatch(expr  = fromJSON(getForm(uri          = rplant.env$webtransform, 
                                           .checkparams = FALSE,
                                           curl         = rplant.env$curl.call)), 
                  error = function(err) {
                            return(paste(err))
                          }
                  )

  if (print.curl) {
    curl.string <- paste(rplant.env$first, "-X GET", rplant.env$webtransform)
    print(curl.string)
  }

  if (!suppress.Warnings){Error(res)}

  file.types <- c()
  for(i in 1:length(res$result)) {
    file.types <- c(file.types, res$result[[i]]$name)
  }
  return(file.types)
}

# -- END -- #




# -- DIRECTORY FUNCTIONS -- #

#####################
#####################
###### ListDir ######
#####################
#####################

ListDir <- function(dir.name="", dir.path="", print.curl=FALSE, 
                    shared.username=NULL, suppress.Warnings=FALSE,
                    show.hidden=FALSE) {
  # This function lists all files in the 'dir.name' contained in 'dir.path'
  #   A user can also list files shared with them from the shared user.
  #
  # Args:
  #   dir.name: Name of directory
  #   dir.path: Path to current directory
  #   shared.username: String, valid iPlant username with whom the object
  #     is being shared.
  #   print.curl: Prints the associated curl statement
  #   suppress.Warnings: Don't do any error checking (faster)
  #
  # Returns:
  #   Returns nothing unless an error if 'dir.name' does not exist
  if (is.null(shared.username)){
    web <- paste(rplant.env$weblist, rplant.env$user, sep="/")
  } else {
    web <- paste(rplant.env$weblist, shared.username, sep="/")
  }

  if (dir.path == ""){
    web <- paste(web, dir.name, sep="/")
  } else {
    web <- paste(web, dir.path, dir.name, sep="/")
  }
  # Check 'dir.name'
  Check(dir.name, dir.path, suppress.Warnings, shared.username) 

  if (print.curl){
    curl.string <- paste(rplant.env$first, " ", web, sep="")
    print(curl.string)
  }
  Renew()
  tmp <- tryCatch(expr  = fromJSON(getURL(web, curl=rplant.env$curl.call)), 
                  error = function(err) {return(paste(err))})
  if (!suppress.Warnings){Error(tmp)}
  nms <- NULL  #create names vector to parse hiddens
  type <- NULL
  for(i in sequence(length(tmp$result))){
    nms <- c(nms, tmp$result[[i]]$name)
    type <- c(type, tmp$result[[i]]$type)
  }
    
  toIgnore <- grep("^\\.+$", nms) # always ignore home directory
  whichHiddens <- grep("^[.]\\D+", nms)
  if(show.hidden)
    toIgnore <- union(toIgnore, whichHiddens)

  nms <- nms[-toIgnore]
  type <- type[-toIgnore]

  # This portion is necessary to weed out the artifact folders.  It probably
  #   won't be implemented because it is slow.
  # newnms <- NULL
  # newtype <- NULL
  # for (i in 1:length(nms)) {
  #   path <- paste(dir.path, dir.name, nms[i], sep="/")
  #   dir.exist <- fromJSON(getURL(url  = paste(rplant.env$webcheck, path, sep=""), curl = rplant.env$curl.call))
  #   if (dir.exist$status != 'error') {
  #     newnms <- append(newnms, nms[i])
  #     newtype <- append(newtype, type[i])
  #   }
  # }

  res <- matrix(nrow=length(nms), ncol=2)
  colnames(res) <- c("name", "type")

  for (i in sequence(dim(res)[1])) {
    res[i, 1] <- nms[i]
    res[i, 2] <- type[i]
  }
  return(res)
}

#####################
#####################
##### ShareDir ######
#####################
#####################

ShareDir <- function(dir.name, dir.path="", shared.username, read=TRUE, 
                     execute=TRUE, write=TRUE, print.curl=FALSE, suppress.Warnings=FALSE) {
  # This function shares the 'dir.name' with shared.username.  Also one can
  #   decide which permissions to give to the shared user.  All contents
  #   within the directory are shared with the shared.username.
  #
  # Args:
  #   dir.name: Name of directory
  #   dir.path: Current path where directory is
  #   shared.username: String, valid iPlant username with whom the object
  #     is being shared.
  #   read: Gives read permissions to directory
  #   execute: Gives execute permissions to directory
  #   write: Gives write permissions to directory
  #   print.curl: Prints the associated curl statement
  #   suppress.Warnings: Don't do any error checking (faster)
  #
  # Returns:
  #   Returns nothing unless an error, directory does not exist

  # Check 'dir.name'
  Check(dir.name, dir.path, suppress.Warnings)

  Share(dir.name, dir.path, shared.username, read, execute, write, print.curl, TRUE)
}

#####################
#####################
## PermissionsDir ###
#####################
#####################

PermissionsDir <- function(dir.name, dir.path="", print.curl=FALSE, suppress.Warnings=FALSE) {
  # This function looks at the permissions on 'dir.name'.  It will return
  #   the directories name, users with whom the directory is shared and their
  #   permissions.
  #
  # Args:
  #   dir.name: Name of directory
  #   dir.path: Current path where directory is
  #   print.curl: Prints the associated curl statement
  #   suppress.Warnings: Don't do any error checking (faster)
  #
  # Returns:
  #   Returns the directory name, users with whom the directory is shared and
  #     their permissions; o/w an error if directory does not exist

  # Check 'dir.name'
  Check(dir.name, dir.path, suppress.Warnings)
    
  Pems(dir.name, dir.path, print.curl, suppress.Warnings)
}

#####################
#####################
##### RenameDir #####
#####################
#####################

RenameDir <- function(dir.name, new.dir.name, dir.path="", print.curl=FALSE, suppress.Warnings=FALSE) {
  # This function simply takes 'dir.name' and renames it to 'new.dir.name'
  #
  # Args:
  #   dir.name: Current name of directory
  #   new.dir.name: New name of directory
  #   dir.path: Path to where directory is
  #   print.curl: Prints the associated curl statement
  #   suppress.Warnings: Don't do any error checking (faster)
  #
  # Returns:
  #   Returns nothing unless an error if directory does not exist, or if
  #     'new.dir.name' does exist in 'dir.path'.

  # Check 'dir.name' in 'dir.path'
  Check(dir.name, dir.path, suppress.Warnings)

  # Check 'new.dir.name' in 'dir.path'
  Check(new.dir.name, dir.path, suppress.Warnings, check=TRUE)

  Rename(dir.name, new.dir.name, dir.path, print.curl) 
}

#####################
#####################
######CopyDir ######
#####################
#####################

CopyDir <- function(dir.name, dir.path="", end.path="", print.curl=FALSE, suppress.Warnings=FALSE) {
  # This function simply moves the 'dir.name' from 'dir.path'
  #   to 'end.path'
  #
  # Args:
  #   dir.name: Name of directory
  #   dir.path: Original or current path where directory is
  #   end.path: Path to where directory will be moved
  #   print.curl: Prints the associated curl statement
  #   suppress.Warnings: Don't do any error checking (faster)
  #
  # Returns:
  #   Returns nothing unless an error.  Error if 'dir.name'
  #     does not exist or if 'dir.name' in 'end.path' already
  #     does exist

  # Check 'dir.name'
  Check(dir.name, dir.path, suppress.Warnings)

  # Check 'dir.name' in 'end.path'
  Check(dir.name, end.path, suppress.Warnings, check=TRUE)

  Copy(dir.name, dir.path, end.path, print.curl)
}

#####################
#####################
###### MoveDir ######
#####################
#####################

MoveDir <- function(dir.name, dir.path="", end.path="", print.curl=FALSE, suppress.Warnings=FALSE) {
  # This function simply moves the 'dir.name' from 'dir.path'
  #   to 'end.path'
  #
  # Args:
  #   dir.name: Name of directory
  #   dir.path: Original or current path where directory is
  #   end.path: Path to where directory will be moved
  #   print.curl: Prints the associated curl statement
  #   suppress.Warnings: Don't do any error checking (faster)
  #
  # Returns:
  #   Returns nothing unless an error.  Error if 'dir.name'
  #     does not exist or if 'dir.name' in 'end.path' already
  #     does exist

  # Check 'dir.name'
  Check(dir.name, dir.path, suppress.Warnings)

  # Check 'dir.name' in 'end.path'
  Check(dir.name, end.path, suppress.Warnings, check=TRUE)

  Move(dir.name, dir.path, end.path, print.curl)
}

#####################
#####################
##### DeleteDir #####
#####################
#####################

DeleteDir <- function(dir.name, dir.path="", print.curl=FALSE, suppress.Warnings=FALSE) {
  # This function removes the 'dir.name' from 'dir.path'
  #
  # Args:
  #   dir.name: Name of directory
  #   dir.path: Path to current directory
  #   print.curl: Prints the associated curl statement
  #   suppress.Warnings: Don't do any error checking (faster)
  #
  # Returns:
  #   Returns nothing unless an error if 'dir.name' does not exist

  # Check 'dir.name'
  Check(dir.name, dir.path, suppress.Warnings)

  Delete(dir.name, dir.path, print.curl)
}

#####################
#####################
###### MakeDir ######
#####################
#####################

MakeDir <- function(dir.name, dir.path="", print.curl=FALSE, 
                    suppress.Warnings=FALSE) {
  # This function simply makes the directory 'dir.name' in 'dir.path'
  #
  # Args:
  #   dir.name: Name of directory
  #   dir.path: Current path where directory is
  #   print.curl: Prints the associated curl statement
  #   suppress.Warnings: Don't do any error checking (faster)
  #
  # Returns:
  #   Returns nothing unless an error.  Error if 'dir.name'
  #     already exists in 'dir.path'
  content <- c()
  if (rplant.env$api == "f") {
    content[1] <- paste("dirName=", dir.name, sep="")
    if (dir.path==""){
      web <- rplant.env$webio
    } else {
      web <- paste(rplant.env$webio, dir.path, sep="/")
    }
  } else {
    web <- rplant.env$webio
    if (dir.path==""){
      content[1] <- paste("path=", dir.name, sep="") 
    } else {
      content[1] <- paste("path=", dir.path, "/", dir.name, sep="") 
    }
  }

  Check(dir.name, dir.path, suppress.Warnings, check=TRUE)

  content[2] <- "action=mkdir"

  if (print.curl) {
    curl.string <- paste(rplant.env$first, " -d '", 
                         paste(content, collapse = "&"), "' ", 
                         web, sep="")
    print(curl.string)
  }

  val <- charToRaw(paste(content, collapse = "&"))
  Renew()
  res <- tryCatch(expr  = fromJSON(httpPUT(url     = web, 
                                           content = val, 
                                           curl    = rplant.env$curl.call)),
                  error = function(err) {
                            return(paste(err))
                          }
                  )
  if (!suppress.Warnings){Error(res)}
}

# -- END -- #




# -- APPLICATION FUNCTIONS -- #

#####################
#####################
##### ListApps ######
#####################
#####################

ListApps<- function (description=FALSE, print.curl=FALSE, suppress.Warnings=FALSE) 
{

  # This function simply lists all of the public applications available to a
  #   user.
  #
  # Args:
  #   description: Either TRUE or FALSE, if TRUE then a short description of
  #     the application is included.
  #   print.curl: Prints the associated curl statement
  #   suppress.Warnings: Don't do any error checking (faster)
  #
  # Returns:
  #   Returns a list of applications, the list contains no duplicates and
  #     only the most current version of an app. o/w errors

  Time()
  Renew()
  tmp <- tryCatch(expr  = fromJSON(getForm(rplant.env$webappslist, 
                                           .checkparams = FALSE, 
                                           curl         = rplant.env$curl.call)), 
                  error = function(err) {
                            return(paste(err))
                          }
                  )

  if (print.curl) {
    curl.string <- paste(rplant.env$first, "-X GET", rplant.env$webappslist)
    print(curl.string)
  }

  if (!suppress.Warnings){Error(tmp)}

  Apps <- list()
  for (j in 1:length(tmp$result)){
    ans <- TestApp(tmp$result[[j]]$id) # returns App id and description
    # The loop is to make sure there are no duplicates
    if ((j != 1) & (!is.null(ans[[1]]))){ # If NOT the first, and not NULL
      for (k in 1:length(Apps)){ # Go through entire App list
        if (ans[[1]] == Apps[[k]][1]){ # if equality, then return NULL
          ans <- list(NULL, NULL)
          break
        }
      }
    }
    if (!is.null(ans[[1]])){ # Now no duplicate so add to list
      Apps <- append(Apps,list(c(ans)))
    } else {
      Apps <- append(Apps,list(c(tmp$result[[j]]$id, "Private Application")))
    }
  }
  if (description == TRUE){ # If description
    res <- matrix(, length(Apps))
    colnames(res) <- "Application" # Below, list both
    for (i in 1:length(Apps)) res[i, 1] <- paste(Apps[[i]], collapse=" - ")
  } else {
    res <- matrix(, length(Apps))
    colnames(res) <- "Application"
    for (i in 1:length(Apps)) res[i, 1] <- Apps[[i]][1] # Just list the ids
  }
  return(sort(res))
}

#####################
#####################
#### GetAppInfo #####
#####################
#####################

GetAppInfo <- function(application, return.json=FALSE, print.curl=FALSE) {
  # This takes the application name and returns basic information about the
  #   app.  For Foundation API and Agave API the applications are inherently
  #   private to the user only, or are public, where all iPlant users can
  #   use it.  An application is private until the user makes it public.
  #
  # Args:
  #   application: A string, application name
  #   return.json: Return json of app, contains all information
  #   print.curl: Prints the associated curl statement
  #
  # Returns:
  #   Returns a description about the application, whether the application
  #     is public, newest version, then vital information about the app,
  #     including input, output, etc. o/w errors
  if (rplant.env$api == "f"){
    tmp_string <- "tmp$result[[len]]"
  } else {
    tmp_string <- "tmp$result"
  }
  # Get the basic application ifo
  result <- appINFO(application)
  text <- result[[1]]
  priv.APP <- result[[2]]

  if (text == "Public App"){
    v.text <- result[[3]]
    APP <- result[[4]]
    tmp <- result[[5]]
  } else {
    tmp <- result[[3]]
  }

  if (print.curl) {
    curl.string <- paste(rplant.env$first, " -X GET ", rplant.env$webappsname,
                         "/", priv.APP, sep="")
    print(curl.string)
  }
  
  if (return.json) {
    return(tmp)
  } else {
    # Go through the json and tease out all of the info we need.
    app.info<-c()
    len <- length(tmp$result)
    for (input in sequence(length(eval(parse(text=paste(tmp_string, "$inputs",
                                                        sep="")))))) {
      app.info <- rbind(app.info, 
                        c("input", 
                          eval(parse(text=paste(tmp_string, 
                                                "$inputs[[input]]$id",
                                                sep=""))), 
                          eval(parse(text=paste(tmp_string, 
                                                "$inputs[[input]]$semantics$fileTypes[1]",
                                                sep=""))), 
                          eval(parse(text=paste(tmp_string, 
                                                "$inputs[[input]]$details$label", 
                                                sep="")))))
    }

    for (parameter in sequence(length(eval(parse(text=paste(tmp_string, "$parameters", 
                                                            sep="")))))) {
      app.info <- rbind(app.info, 
                        c("parameters", 
                          eval(parse(text=paste(tmp_string, 
                                                "$parameters[[parameter]]$id", 
                                                sep=""))), 
                          eval(parse(text=paste(tmp_string, 
                                                "$parameters[[parameter]]$value$type", 
                                                sep=""))), 
                          eval(parse(text=paste(tmp_string, 
                                                "$parameters[[parameter]]$details$label",
                                                sep="")))))


#      if (eval(parse(text=paste(tmp_string,"$parameters[[parameter]]$value$type", sep=""))) == "enumeration") {
#          for (i in c(1:length(eval(parse(text=paste(tmp_string,"$parameters[[parameter]]$value$enum_values", sep="")))))) {
#app.info <- rbind(app.info, c("", paste("enum-choice",i), names(eval(parse(text=paste(tmp_string,"$parameters[[parameter]]$value$enum_values[[i]]", sep="")))), eval(parse(text=paste(tmp_string,"$parameters[[parameter]]$value$enum_values[[i]]", sep="")))))
#}
#      }
      if (eval(parse(text=paste(tmp_string,"$parameters[[parameter]]$value$type", sep=""))) == "enumeration") {
          for (i in c(1:length(eval(parse(text=paste(tmp_string,"$parameters[[parameter]]$value$enum_values", sep="")))))) {
app.info <- rbind(app.info, c("", paste("enum-choice",i), names(eval(parse(text=paste(tmp_string,"$parameters[[parameter]]$value$enum_values[[i]]", sep="")))), eval(parse(text=paste(tmp_string,"$parameters[[parameter]]$value$enum_values[[i]][[1]]", sep="")))))
}
      }
  }
    
    for (output in sequence(length(eval(parse(text=paste(tmp_string, "$outputs",
                                                         sep="")))))) {


        
      app.info <- rbind(app.info, 
                        c("output", 
                          eval(parse(text=paste(tmp_string, 
                                                "$outputs[[output]]$id", 
                                                sep=""))), 
                          eval(parse(text=paste(tmp_string, 
                                                "$outputs[[output]]$semantics$fileTypes[1]", 
                                                sep=""))), 
                          eval(parse(text=paste(tmp_string, 
                                                "$outputs[[output]]$details$label", 
                                                sep=""))))) 
    }


    shortd <- eval(parse(text=paste(tmp_string, "$shortDescription", sep="")))
    shortn <- nchar(shortd)
    longd <- eval(parse(text=paste(tmp_string, "$longDescription", sep="")))
    if (is.null(longd)) {longn = 0} else {longn <- nchar(longd)}
    if (longn >= shortn) {
      description <- longd
    } else {
      description <- shortd
    }
    colnames(app.info)<-c("kind", "id", "fileType/value", "details")
    if (text == "Private App"){
      return(list(Description=description, Application=c(application, text),
                  Information=app.info))
    } else {
      return(list(Description=description, Application=c(application, text, v.text),
                  Information=app.info))
    }
  }
}


# -- END -- #
 

# -- JOB FUNCTIONS -- #

#####################
#####################
##### SubmitJob #####
#####################
#####################

SubmitJob <- function(application, file.path="", file.list=NULL, input.list, 
                      args.list=NULL, job.name, nprocs=1, private.APP=FALSE, 
                      suppress.Warnings=FALSE, shared.username=NULL,
                      print.curl=FALSE) {
  
  if (private.APP) {
    suppress.Warnings=TRUE
  }
  
  # This takes the application name and returns basic information about the
  #   app.
  #
  # Args:
  #   application: A string, application name
  #   file.path: Path to where ALL input files are located
  #   file.list: List of input files, can be many input
  #   input.list: List corresponding to file list, is type of input
  #     see help(SubmitJob) for details.  Use GetAppList to find input list
  #   args.list: List of arguments for the specific application.  This list
  #     has a very specific format that is included in the help(SubmitJob) file
  #   job.name: Job name adds a time stamp to make them unique
  #   nprocs: Number of processors allocated to job.  This number depends
  #     on if application is parallelizable.
  #   private.APP: Either TRUE or FALSE, if TRUE the application is private
  #     to the user, o/w the app is public
  #   email: Either TRUE or FALSE, if TRUE the user is sent an email when
  #     jov is finished.
  #   shared.username: String, valid iPlant username with whom the object
  #     is being shared.
  #   print.curl: Prints the associated curl statement
  #   suppress.Warnings: Don't do any error checking (faster)
  #
  # Returns:
  #   Returns the job id (number).  o/w an error

  # Job Name is automatically time stamped
  job.name <- paste(unlist(strsplit(paste(job.name, "_", format(Sys.time(), 
                    "%Y-%m-%d_%k-%M-%OS3"), sep=""), " ")), collapse="")

  if (nchar(job.name) > 128) {
    total = nchar(job.name) - 127
    job.name = substring(job.name, total, nchar(job.name))
  }

  Time()
  Renew()
  content <- c()
  # Create the options
  if (rplant.env$api == "f") {
    tmp_string <- "tmp$result[[len]]"
    eml_string <- "callbackUrl="
    content[1] <- paste("jobName=", job.name, sep="")
    content[2] <- paste("softwareName=", application, sep="")
    content[3] <- "requestedTime=24:00:00"
  } else {
    tmp_string <- "tmp$result"
    eml_string <- "callbackURL="
    content[1] <- paste("name=", job.name, sep="")
    content[2] <- paste("appId=", application, sep="")
    content[3] <- "maxRunTime=24:00:00"
  }
  # Check that all of the files exist
  for (i in 1:length(file.list)){
    Check(file.list[[i]], file.path, suppress.Warnings, shared.username)
  }

  if (suppress.Warnings == FALSE){
    # Output includes info about inputs of app
    result <- appINFO(application, FALSE, TRUE)
   
    if (result[[1]] == "Public App"){
      input <- result[[6]] # Input information, compare to input.list
      set <- result[[7]] # Parallelization information of app
    } else {
      input <- result[[4]] # Input information, compare to input.list
      set <- result[[5]] # Parallelization information of app
    }
    # Compare the input.list to actual inputs of the application.
    #   If they don't match throw an error.
    test.input <- rep(0, length(input.list))

#   for (j in 1:length(input.list)){
#     if (!input.list[[j]] % in % input){
#       test.input[j] <- 1
#       break;
#     }
#   }

    for (j in 1:length(input.list)){
      cnt = 0
      for (i in 1:length(input)){
        if (input.list[[j]] == input[i]){
          cnt = cnt + 1
        }
      }
      if (cnt != 1){
        test.input[j] <- 1
        break;
      }
    }

    # Throw an error if one of the inputs in input.list is incorrect
    if (sum(test.input) > 0 ){
      return(stop(paste("At least one of the inputs in 'input.list' is incorrect,",
                        "check GetAppInfo function for proper inputs", sep = ""), 
                  call. = FALSE))
    }
  }
  # If suppressing Warnings then the look up about the application didn't 
  #   happen.  Need to look up the application to get parallelization
  #   information.
  if (suppress.Warnings==TRUE){
    result <- appINFO(application)
    if (result[[1]] == "Public App"){
      set <- result[[6]]
    } else {
      set <- result[[4]]
    }
  }
  # If the application is a private app, then stop, because we can't run
  #   a job on someone elses private app
  if (private.APP==FALSE){
    if (result[[1]] == "Private App"){
      return(stop(paste("Private application, not valid for SubmitJob.  If it", 
                        "is your own private application use private.APP=TRUE"),
                  call. = FALSE))
    }
  }

  # If the application can be run in parallel, then increase the number of
  #   processors to be used.  A user cannot use more than 512 processors
  if (set == "PARALLEL"){
    if (nprocs < 2){
      nprocs = 12
    } else if (nprocs > 512){
      nprocs = 512
    }
  } else { # If the application is not parallel nprocs is set to one
    if (nprocs != 1){
      nprocs = 1
    }
  }

  if (rplant.env$api == "f") {
    content[4] <- paste("processorCount=", nprocs, sep="")
  } else {
    content[4] <- paste("nodeCount=", nprocs, sep="")
  }

  # Automatically makes analyses directory; will not overwrite if already present
  MakeDir("analyses", suppress.Warnings=TRUE)

  # Set archivePath, where job output will be returned
  content[5] <- "archive=1"
  content[6] <- paste("archivePath=", rplant.env$user, "analyses", 
                      job.name, sep="/"); x <- 6; # x tells the length of options

  ### I had the email working, but it currently does not, if someone in the 
  ###   future gets it running again, that would be swell.  Currently I will just
  ###   comment it out.
  # If email is TRUE
#  if (email==TRUE){
#    Renew()
#    res <- tryCatch(expr  = fromJSON(getURLContent(url  = rplant.env$webprofiles, 
#                                                   curl = rplant.env$curl.call)), 
#                    error = function(err) {
#                              return(paste(err))
#                            }
#                    )
#    Error(res)
#    content[7] <- paste(eml_string, res$result[[1]]$email, sep=""); x <- 7;
#  }

  # For the loop below n needs to be initialized
  if(!is.null(file.list)){n <- length(file.list)} else {n <- 0}

  # For all of the files in file.list the options needed to line up with the
  #   input.list, and also go to the correct directory.
  if (n > 0){
    for (i in c(1:n)){
      if (file.path=="") {
        if (is.null(shared.username)){
          content[x+i] <- paste(input.list[[i]],"=/", rplant.env$user, 
                                "/", file.list[[i]], sep="")
        } else {
          content[x+i] <- paste(input.list[[i]],"=/", shared.username, 
                                "/", file.list[[i]], sep="")
        }
      } else {
        if (is.null(shared.username)){
          content[x+i] <- paste(input.list[[i]],"=/", rplant.env$user, "/",
                                file.path, "/", file.list[[i]], sep="")
        } else {
          content[x+i] <- paste(input.list[[i]],"=/", shared.username, "/",
                                file.path, "/", file.list[[i]], sep="")
        }
      }
    }
  }

  # For the loop below m needs to be initialized
  if(!is.null(args.list)){m <- length(args.list)} else {m <- 0}

  # For the specific format of args.list, if there are arguments
  #   add them to the options.
  if (m > 0){
    for (i in c(1:m)){
      content[x+n+i] <- paste(args.list[[i]][1],"=", 
                               args.list[[i]][2], sep="")
    }
  }

  if (print.curl) {
    curl.string <- paste(rplant.env$first," -X POST -d '", 
                         paste(content, collapse = "&"), "' ", 
                         rplant.env$webjob, sep="")
    print(curl.string)
  }

  val <- charToRaw(paste(content, collapse = "&"))
  Renew()
  res <- tryCatch(expr  = fromJSON(getURLContent(rplant.env$webjob, 
                                                 curl          = rplant.env$curl.call,
                                                 infilesize    = length(val), 
                                                 readfunction  = val, 
                                                 upload        = TRUE, 
                                                 customrequest = "POST")),
                  error = function(err) {
                            return(paste(err))
                          }
                  )
  if (!suppress.Warnings){Error(res)}
  cat("Job submitted. \n")
  cat(paste("You can check your job using CheckJobStatus(", 
            res$result$id, ")", sep=""), "\n")
  # return(res$result$id)
  output <- vector("list", 2)
  names(output) <- c("id", "name")
  output$id <- res$result$id
  output$name <- job.name
  return(output)
}

#####################
#####################
## CheckJobStatus ###
#####################
#####################

CheckJobStatus <- function(job.id, history=FALSE, print.curl=FALSE) {
  # This function checks the job status of the job with that job number
  #
  # Args:
  #   job.id: Job number of job to be checked
  #   history:  Either TRUE or FALSE, only for Agave API, if TRUE
  #     then will show entire history of job.
  #   print.curl: Prints the associated curl statement
  #
  # Returns:
  #   Returns the status of the job:
  #     PENDING            
  #     STAGING_INPUTS     
  #     CLEANING_UP      
  #     ARCHIVING         
  #     STAGING_JOB        
  #     FINISHED          
  #     KILLED            
  #     FAILED             
  #     STOPPED            
  #     RUNNING            
  #     PAUSED             
  #     QUEUED             
  #     SUBMITTING         
  #     STAGED             
  #     PROCESSING_INPUTS  
  #     ARCHIVING_FINISHED 
  #     ARCHIVING_FAILED  
 
  Time()
  Renew()

  web <- paste(rplant.env$webjob, job.id, sep="/")

  if (!(((rplant.env$api == "f") && (history == TRUE)) || (history == FALSE))){
     web <- paste(web, "history", sep="")
  }

  res <- tryCatch(expr  = fromJSON(getForm(web, 
                                           .checkparams = FALSE, 
                                           curl         = rplant.env$curl.call)), 
                  error = function(err) {
                            return(paste(err))
                          }
                  )

  if (print.curl) {
    curl.string <- paste(rplant.env$first, web)
    print(curl.string)
  }

  Error(res)

  if (!(((rplant.env$api == "f") && (history == TRUE)) || (history == FALSE))){
    return(res$result)
  } else { 
    return(res$result$status)
  }
}

#####################
#####################
###### KillJob ######
#####################
#####################

KillJob <- function(job.id, print.curl=FALSE) {
  # This function stops the job with the job number
  #
  # Args:
  #   job.id: Job number of job to be checked
  #   print.curl: Prints the associated curl statement
  #
  # Returns:
  #   Returns nothing unless an Error
  Time()
  Renew()

  web <- paste(rplant.env$webjob, job.id, sep="/")

  content <- c()
  content[1] <- "action=stop"

  val <- charToRaw(paste(content, collapse = "&"))

  res <- tryCatch(expr  = fromJSON(getURLContent(web, 
                                                 curl          = rplant.env$curl.call,
                                                 infilesize    = length(val), 
                                                 readfunction  = val, 
                                                 upload        = TRUE, 
                                                 customrequest = "POST")),
                  error = function(err) {
                            return(paste(err))
                          }
                  )

  if (print.curl) {
    curl.string <- paste(rplant.env$first, " -X POST -d '", 
                         paste(content, collapse = "&"), "' ", 
                         web, sep="")
    print(curl.string)
  }

  Error(res)

  DeleteJob(job.id)
}

#####################
#####################
##### DeleteOne #####
#####################
#####################

DeleteOne <- function(job.id, print.curl=FALSE) {
  # This function deletes the job with the job number, and it deletes
  #   the folder the result files are in
  #
  # Args:
  #   job.id: Job number of job to be checked
  #   print.curl: Prints the associated curl statement
  #
  # Returns:
  #   Returns nothing unless an Error
  Time()
  Renew()

  web <- paste(rplant.env$webjob, job.id, sep="/")

  # Check that job id exists
  JS <- tryCatch(expr  = fromJSON(getForm(web, 
                                          .checkparams = FALSE, 
                                          curl         = rplant.env$curl.call)), 
                 error = function(err) {
                           return(paste(err))
                         }
                 )
  Error(JS)

  # If the job is finished or stopped then it can be deleted, it it's running
  #   an appropriate error is returned
  if ((JS$result$status == "FINISHED") || (JS$result$status == "STOPPED") || 
      (JS$result$status == "ARCHIVING_FINISHED") || (JS$result$status == "FAILED")){
    # DeleteJob deletes the directory the result files are in, this finds that folder
    if (JS$result$archive == TRUE){
      dir.name <- unlist(strsplit(JS$result$archivePath, "/"))[
                         length(unlist(strsplit(JS$result$archivePath, "/")))]

      dir.path <- substr(JS$result$archivePath, nchar(rplant.env$user) + 3, 
                         nchar(JS$result$archivePath)-nchar(dir.name)-1)

    # Delete the directory
      tmp <- tryCatch(expr  = fromJSON(httpDELETE(paste(rplant.env$webio, dir.path, 
                                                        dir.name, sep="/"), 
                                                  curl = rplant.env$curl.call)),
                      error = function(err) {
                                return(paste(err))
                              }
                      )
    }
    # Delete the job
    tmp <- tryCatch(expr  = fromJSON(httpDELETE(web, curl = rplant.env$curl.call)),
                    error = function(err) {
                              return(paste(err))
                            }
                    )
  
    if (print.curl) {
      curl.string <- paste(rplant.env$first, "-X DELETE", web)
      print(curl.string)
    }

    Error(tmp)

  } else {
    return(stop(paste("Error: Could not delete, job status:", 
                      JS$result$status), call. = FALSE))
  }
}
#####################
#####################
##### DeleteALL #####
#####################
#####################

DeleteALL <- function() {
  # Deletes all of the jobs in the job history, and their associated
  #   directories
  #
  # Args:
  #   Nothing
  #
  # Returns:
  #   Returns nothing unless an Error
  Time()
  Renew()

  # Get the entire job list
  res <- tryCatch(expr  = fromJSON(getForm(rplant.env$webjoblist, 
                                           .checkparams = FALSE, 
                                           curl         = rplant.env$curl.call)), 
                  error = function(err) {
                            return(paste(err))
                          }
                  )
  Error(res)

  if (length(res$result) == 0) {
    message("No jobs in job history")
  } else { # Go through each job and delete it
    for (i in 1:length(res$result)){
      JS <- tryCatch(expr = DeleteOne(res$result[[i]]$id), 
                 error = function(err) {
                           return(paste(err))
                         }
                 )
    }
  }
}



#####################
#####################
##### DeleteJob #####
#####################
#####################

DeleteJob <- function(job.id, print.curl=FALSE, ALL=FALSE) {
  # Deletes the job with the job number, also have the option to delete all
  #   jobs
  #
  # Args:
  #   job.id: Job number of job to be deleted
  #   print.curl: Prints the associated curl statement
  #   ALL: Delete all jobs in the job history
  #
  # Returns:
  #   Returns nothing unless an Error
  if (ALL==TRUE){
    DeleteALL()
    if (print.curl) {
      message("No curl statement to print")
    }
  } else {
    DeleteOne(job.id, print.curl)
  }
}

#####################
#####################
#### RetrieveOne ####
#####################
#####################


RetrieveOne <- function(file, archive.path, file.path, print.curl) {  
  # This function takes the file in the archive.path from the iPlant servers
  # to the file.path on the local computer.
  #
  # Args:
  #   file: String of the file name
  #   archive.path: Path to the file on the iPlant server side
  #   file.path: Path to where file will be downloaded on local computer
  #   print.curl: Prints the associated curl statement
  #
  # Returns:
  #   Returns nothing unless an Error
  Time()
  Renew()

  web <- paste(rplant.env$webio1, archive.path, "/", file, sep="")

  curlPerform(url       = web, 
              curl      = rplant.env$curl.call, 
              writedata = CFILE(file.path(file.path,file), 
              mode      = "wrb")@ref)
  gc()

  if (print.curl) {
    curl.string <- paste(rplant.env$first, "-X GET", web)
    print(curl.string)
  }
}

#####################
#####################
#### RetrieveJob ####
#####################
#####################

RetrieveJob <- function(job.id, file.vec=NULL, print.curl=FALSE, verbose=FALSE) {  
  # From the job number, retrieve the files in the file.vec and download
  # to the current directory on the local computer.  A folder the name
  # of the job is created, and all files put inside.
  #
  # Args:
  #   job.id: Job number of job whose files will be retreived
  #   file.vec: Vector containing file names, if NULL all files will be 
  #     downloaded
  #   print.curl: Prints the associated curl statement
  #   verbose: Either TRUE or FALSE, if TRUE print a statment listing
  #     which files have been downloaded to which folder
  #
  # Returns:
  #   Returns nothing unless verbose=TRUE or an Error
  Time()
  Renew()

  web <- paste(rplant.env$webjob, job.id, sep="/")

  JS <- tryCatch(expr  = fromJSON(getForm(web, 
                                          .checkparams=FALSE, 
                                           curl=rplant.env$curl.call)), 
                 error = function(err) {
                           return(paste(err))
                         }
                 )
  Error(JS)

  if ((JS$res$status == "ARCHIVING_FINISHED") || (JS$res$status == "FINISHED")) {
    # Create a local folder in current R working directory, the folder's name
    #   is the job name given to the job previously.
    dir.path <- file.path(getwd(), JS$result[[2]])
    if(!file.exists(dir.path)){
      if (.Platform$OS.type=="windows") {
        invisible(shell(paste("mkdir ", JS$result[[2]], sep="")))
      } else {
        dir.create(dir.path)
      }
    }

    # If file.vec is NULL, get all result file name for job number
    if(is.null(file.vec)){
      file.vec <- ListJobOutput(job.id, print.total=FALSE)
    }  

    # Get all result file name for job number
    fileList <- ListJobOutput(job.id, print.total=FALSE)

    # Go through each file in file.vec
    for (file in 1:length(file.vec)) {
      # if file exists in output then download
      if (file.vec[file] %in% fileList) {

        RetrieveOne(file.vec[file], JS$result$archivePath, dir.path, print.curl)

        if (verbose==TRUE) { # If TRUE, verbose output
          message(paste("Downloaded", file.vec[file], "to", dir.path))
        }
      } else { # If file not there, return an error
        return(stop(paste("`",file.vec[file], "' is not found within `", 
                          job.id,"'", sep=""), call. = FALSE))
      }
    }
  } else { # If job is not finished
    return(stop(paste("Job is", JS$res$status), call. = FALSE))
  }
}

#####################
#####################
### ListJobOutput ###
#####################
#####################

ListJobOutput <- function(job.id, print.curl=FALSE, print.total=TRUE) {
  # List the names of the result files with the job number
  #
  # Args:
  #   job.id: Job number of job whose files will be retreived
  #   print.curl: Prints the associated curl statement
  #   print.total: Either TRUE or FALSE, if TRUE prints number of files
  #
  # Returns:
  #   Returns the file list, else an error
  Time()
  Renew()

  # Check status of the job, and that it exists
  JS <- tryCatch(expr  = fromJSON(getForm(paste(rplant.env$webjob, job.id, sep="/"), 
                                          .checkparams = FALSE, 
                                          curl         = rplant.env$curl.call)), 
                 error = function(err) {
                           return(paste(err))
                         }
                 )
  Error(JS)

  file.vec <- c()
  if ((JS$res$status == "FINISHED") || (JS$res$status == "ARCHIVING_FINISHED")) {
    # Knowing it does, and knowing the archivePath, where the result files
    #   are located, get the list of files in the folder
    web <- paste(rplant.env$weblist, JS$result$archivePath, sep="")
    res <- fromJSON(getURLContent(web, curl=rplant.env$curl.call))

    if (print.curl) {
      curl.string <- paste(rplant.env$first, "-X GET", web)
      print(curl.string)
    }

    len <- length(res$result)
    if (len == 0){
      return(paste("There are ", len, " output files for job '", 
                   job.id,"'", sep=""))
    }

    if (print.total == TRUE) {
      message(paste("There are ", len-1, " output files for job '", 
                    job.id,"'", sep=""))
    }
    # Add each name into 'file.vec'
    for (i in 2:length(res$result)) {
      file.vec <- append(file.vec, res$result[[i]]$name)
    }
    return(file.vec)
  } else {
    return(stop(paste("Job is", JS$res$status), call. = FALSE))
  }
}

#####################
#####################
### GetJobHistory ###
#####################
#####################

GetJobHistory <- function(return.json=FALSE, print.curl=FALSE) {
  # List all of the jobs in the job history
  #
  # Args:
  #   return.json: Returns a json containing all information
  #   print.curl: Prints the associated curl statement
  #
  # Returns:
  #   Returns a list of jobs, with id, name and current status.  If no
  #     job then return "No jobs in history"
  Time()
  Renew()
  if (rplant.env$api == "f") {
    tmp_string <- "res$result[[i]]$software"
  } else {
    tmp_string <- "res$result[[i]]$appId"
  }

  jobList <- c()

  res <- tryCatch(expr  = fromJSON(getForm(rplant.env$webjoblist, 
                                           .checkparams = FALSE, 
                                           curl         = rplant.env$curl.call)), 
                  error = function(err) {
                            return(paste(err))
                          }
                  )

  if (print.curl) {
    curl.string <- paste(rplant.env$first, "-X GET", rplant.env$webjoblist)
    print(curl.string)
  }

  Error(res)

  if (length(res$result) == 0){
    return("No jobs in history")
  } else {
    if (return.json) 
      return(res)

    for (i in 1: length(res$result)) {
      job <- c(res$result[[i]]$id, res$result[[i]]$name, 
               eval(parse(text=tmp_string)), res$result[[i]]$status) 
      jobList <- rbind(jobList, job)
      colnames(jobList) <- c("job.id", "job.name", "application", "status")
    } 
    return(jobList)
  }
}

# -- END -- #
