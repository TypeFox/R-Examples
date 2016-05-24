#
#    _____      _ ____  ____
#   / ___/_____(_) __ \/ __ )
#   \__ \/ ___/ / / / / __  |
#  ___/ / /__/ / /_/ / /_/ / 
# /____/\___/_/_____/_____/  
#
#
#
# BEGIN_COPYRIGHT
#
# This file is part of SciDB.
# Copyright (C) 2008-2014 SciDB, Inc.
#
# SciDB is free software: you can redistribute it and/or modify
# it under the terms of the AFFERO GNU General Public License as published by
# the Free Software Foundation.
#
# SciDB is distributed "AS-IS" AND WITHOUT ANY WARRANTY OF ANY KIND,
# INCLUDING ANY IMPLIED WARRANTY OF MERCHANTABILITY,
# NON-INFRINGEMENT, OR FITNESS FOR A PARTICULAR PURPOSE. See
# the AFFERO GNU General Public License for the complete license terms.
#
# You should have received a copy of the AFFERO GNU General Public License
# along with SciDB.  If not, see <http://www.gnu.org/licenses/agpl-3.0.html>
#
# END_COPYRIGHT
#

# This file contains general utility routines including most of the shim
# network interface.

# An environment to hold connection state
.scidbenv = new.env()

# Force evaluation of an expression that yields a scidb or scidbdf object,
# storing result to a SciDB array when eval=TRUE.
# name: (character) optional SciDB array name to store to
# gc: (logical) optional, when TRUE tie result to R garbage collector
# temp (logical, optional): when TRUE store as a SciDB temp array
scidbeval = function(expr, eval=TRUE, name, gc=TRUE, temp=FALSE)
{
  ans = eval(expr)
  if(!(inherits(ans,"scidb") || inherits(ans,"scidbdf"))) return(ans)
# If expr is a stored temp array, then re-use its name
  if(!is.null(ans@gc$temp) && ans@gc$temp && missing(name)) name=ans@name
  .scidbeval(ans@name, `eval`=eval, name=name, gc=gc, schema=ans@schema, temp=temp)
}

# Create a new scidb reference to an existing SciDB array.
# name (character): SciDB expression or array name
# gc (logical, optional): Remove backing SciDB array when R object is
#     garbage collected? Default is FALSE.
# data.frame (logical, optional): If true, return a data.frame-like object.
#   Otherwise an array.
scidb = function(name, gc, `data.frame`)
{
  if(missing(name)) stop("array or expression must be specified")
  if(missing(gc)) gc=FALSE
  if(is.scidb(name) || is.scidbdf(name))
  {
    query = name@name
    return(.scidbeval(name@name, eval=FALSE, gc=gc, `data.frame`=`data.frame`, depend=list(name)))
  }
  escape = gsub("'","\\\\'",name,perl=TRUE)
  query = sprintf("join(show('filter(%s,true)','afl'), explain_logical('filter(%s,true)','afl'))", escape,escape)
  query = iquery(query, `return`=TRUE)
  logical_plan = query$logical_plan
  schema = gsub("^.*<","<",query$schema, perl=TRUE)
  obj = scidb_from_schemastring(schema, name, `data.frame`)
  obj@logical_plan = logical_plan
  if(gc)
  {
    if(length(grep("\\(",name))==0)
    {
      obj@gc$name = name
    }
    obj@gc$remove = TRUE
    reg.finalizer(obj@gc, function(e)
        {
          if (e$remove && exists("name",envir=e))
            {
              tryCatch(scidbremove(e$name,warn=FALSE), error = function(e) invisible())
            }
        }, onexit = TRUE)
  } else obj@gc = new.env()
  obj
}

# TRUE if the SciDB object 'name' or character name of a SciDB array is
# a temporary SciDB array, FALSE otherwise.
is.temp = function(name)
{
  if(is.scidb(name) || is.scidbdf(name))
  {
    if(!is.null(name@gc$temp)) return(name@gc$temp)
    name = name@name
  }
  query = sprintf("filter(list('arrays'),name='%s')", name)
  ans = iquery(query,return=TRUE)
  if(nrow(ans)<1) return(FALSE)
  if(is.null(ans$temporary)) return(FALSE)
  ans$temporary == "true"
}

# An important internal convenience function that returns a scidb object.  If
# eval=TRUE, a new SciDB array is created the returned scidb object refers to
# that.  Otherwise, the returned scidb object represents a SciDB array promise.
#
# INPUT
# expr: (character) A SciDB expression or array name
# eval: (logical) If TRUE evaluate expression and assign to new SciDB array.
#                 If FALSE, infer output schema but don't evaluate.
# name: (optional character) If supplied, name for stored array when eval=TRUE
# gc: (optional logical) If TRUE, tie SciDB object to  garbage collector.
# depend: (optional list) An optional list of other scidb or scidbdf objects
#         that this expression depends on (preventing their garbage collection
#         if other references to them go away).
# data.frame: (optional, logical) If TRUE, return a data.frame object, false
#             return a scidb object. Default is missing, in which case an
#             automatic decision is made about the object return class.
# schema, temp: (optional) used to create SciDB temp arrays
#               (requires scidb >= 14.8)
#
# OUTPUT
# A `scidb` or `scidbdf` array object.
#
# NOTE
# Only AFL supported.
`.scidbeval` = function(expr,eval,name,gc=TRUE, depend, `data.frame`, schema, temp)
{
  ans = c()
  if(missing(depend)) depend=c()
  if(missing(schema)) schema=""
  if(missing(temp)) temp=FALSE
  if(!is.list(depend)) depend=list(depend)
# Address bug #45. Try to cheaply determine if expr refers to a named array
# or an AFL expression. If it's a named array, then eval must be set TRUE.
  if(!grepl("\\(", expr, perl=TRUE)) eval = TRUE
  if(`eval`)
  {
    if(missing(name) || is.null(name))
    {
      newarray = tmpnam()
      if(temp)
      {
# SciDB temporary array syntax varies with SciDB version
        TEMP = "'TEMP'"
        if(compare_versions(options("scidb.version")[[1]],14.12)) TEMP="true"
        query   = sprintf("create_array(%s, %s, %s)", newarray, schema, TEMP)
        iquery(query, `return`=FALSE)
      }
    }
    else newarray = name
    query = sprintf("store(%s,%s)",expr,newarray)
    scidbquery(query, interrupt=TRUE, stream=0L)
    ans = scidb(newarray,gc=gc,`data.frame`=`data.frame`)
    if(temp) ans@gc$temp = TRUE
# This is a fix for a SciDB issue that can unexpectedly change schema
# bounds. And another fix to allow unexpected dimname and attribute name
# changes. Arrgh.
    if(schema!="" && !compare_schema(ans, schema, ignore_attributes=TRUE, ignore_dimnames=TRUE))
    {
      ans = repart(ans, schema)
    }
  } else
  {
    ans = scidb(expr,gc=gc,`data.frame`=`data.frame`)
# Assign dependencies
    if(length(depend)>0)
    {
      assign("depend",depend,envir=ans@gc)
    }
  }
  ans
}


# store the connection information and obtain a unique ID
scidbconnect = function(host=options("scidb.default_shim_host")[[1]],
                        port=options("scidb.default_shim_port")[[1]],
                        username, password,
                        auth_type=c("pam","digest"), protocol=c("http","https"))
{
  scidbdisconnect()
  auth_type = match.arg(auth_type)
  protocol = match.arg(protocol)
  assign("host",host, envir=.scidbenv)
  assign("port",port, envir=.scidbenv)
  assign("protocol",protocol,envir=.scidbenv)
  if(missing(username)) username=c()
  if(missing(password)) password=c()
# Check for login using either pam or basic authentication
  if(!is.null(username))
  {
    if(auth_type=="pam")
    {
      auth = GET(resource="login",list(username=username, password=password),header=FALSE)
      if(nchar(auth)<1) stop("Authentication error")
      assign("auth",auth,envir=.scidbenv)
    } else # HTTP basic digest auth
    {
      assign("digest",paste(username,password,sep=":"),envir=.scidbenv)
    }
    assign("authtype",auth_type,envir=.scidbenv)
    assign("authenv",new.env(),envir=.scidbenv)
    reg.finalizer(.scidbenv$authenv, function(e) scidblogout(), onexit=TRUE)
  }

# Use the query ID from a query as a unique ID for automated
# array name generation.
  x = tryCatch(
        scidbquery(query="setopt('precision','16')",release=1,resp=TRUE,stream=0L),
        error=function(e) stop("Connection error"))
  if(is.null(.scidbenv$uid))
  {
    id = tryCatch(strsplit(x$response, split="\\r\\n")[[1]],
           error=function(e) stop("Connection error"))
    id = id[[length(id)]]
    assign("uid",id,envir=.scidbenv)
  }
# Try to load the dense_linear_algebra library
  tryCatch(
    scidbquery(query="load_library('dense_linear_algebra')",
               release=1,resp=FALSE, stream=0L),
    error=invisible)
# Try to load the example_udos library (>= SciDB 13.6)
  tryCatch(
    scidbquery(query="load_library('example_udos')",release=1,resp=FALSE,stream=0L),
    error=invisible)
# Try to load the superfunpack
  tryCatch(
    scidbquery(query="load_library('superfunpack')",release=1,resp=FALSE,stream=0L),
    error=invisible)
# Try to load the P4 library
  tryCatch(
    scidbquery(query="load_library('linear_algebra')",release=1,resp=FALSE,stream=0L),
    error=invisible)
# Save available operators
  assign("ops",iquery("list('operators')",return=TRUE),envir=.scidbenv)
# Update the scidb.version option
  v = scidbls(type="libraries")[1,]
  if("major" %in% names(v))
  {
    options(scidb.version=paste(v$major,v$minor,sep="."))
  }
# Update the gemm_bug option
  if(compare_versions(options("scidb.version")[[1]],14.3))
  {
    options(scidb.gemm_bug=FALSE) # Yay
  }
  invisible()
}

scidblogout = function()
{
  if(exists("auth",envir=.scidbenv) && exists("host",envir=.scidbenv))
  {
    GET(resource="logout",header=FALSE)
  }
}

scidbdisconnect = function()
{
# forces call to scidblogout via finalizer, see scidbconnect:
  if(exists("authenv",envir=.scidbenv)) rm("authenv",envir=.scidbenv)
  gc()
  rm(list=ls(envir=.scidbenv), envir=.scidbenv)
}

make.names_ = function(x)
{
  gsub("\\.","_",make.names(x, unique=TRUE),perl=TRUE)
}

# x is vector of existing values
# y is vector of new values
# returns a set the same size as y with non-conflicting value names
make.unique_ = function(x,y)
{
  z = make.names(gsub("_",".",c(x,y)),unique=TRUE)
  gsub("\\.","_",tail(z,length(y)))
}

# Make a name from a prefix and a unique SciDB identifier.
tmpnam = function(prefix="R_array")
{
  salt = basename(tempfile(pattern=prefix))
  if(!exists("uid",envir=.scidbenv)) stop("Not connected...try scidbconnect")
  paste(salt,get("uid",envir=.scidbenv),sep="")
}

# Return a shim session ID or error
# This will also return an authenticaion token string if one is available.
getSession = function()
{
  session = GET("/new_session",header=FALSE)
  if(length(session)<1) stop("SciDB http session error; are you connecting to a valid SciDB host?")
  session = gsub("\r","",session)
  session = gsub("\n","",session)
  session
}

# Returns the eponymous package option setting used to globally overried
# interrupt handling. TRUE means abide by interrupts, FALSE ignore.
scidb.interrupt = function()
{
  ans = options("scidb.interrupt")[[1]]
  if(is.null(ans)) ans=TRUE
  tryCatch(as.logical(ans), error=function(e) TRUE)
}

# Supply the base SciDB URI from the global host, port and auth
# parameters stored in the .scidbenv package environment.
# Every function that needs to talk to the shim interface should use
# this function to supply the URI.
# Arguments:
# resource (string): A URI identifying the requested service
# args (list): A list of named query parameters
URI = function(resource="", args=list())
{
  if(!exists("host",envir=.scidbenv)) stop("Not connected...try scidbconnect")
  if(exists("auth",envir=.scidbenv))
    args = c(args,list(auth=get("auth",envir=.scidbenv)))
  prot = paste(get("protocol",envir=.scidbenv),"//",sep=":")
  if("username" %in% names(args) || "auth" %in% names(args)) prot = "https://"
  ans  = paste(prot, get("host",envir=.scidbenv),":",get("port",envir=.scidbenv),sep="")
  ans = paste(ans, resource, sep="/")
  if(length(args)>0)
    ans = paste(ans,paste(paste(names(args),args,sep="="),collapse="&"),sep="?")
  ans
}

# Send an HTTP GET message to the shim service
# resource: URI service name
# args: named list of HTTP GET query parameters
# header: TRUE=return the HTTP response header, FALSE=don't return it
# async: Doesn't really do anything, should be removed XXX
# interrupt: Set to FALSE to disable SIGINT, otherwise handle gracefully.
# ...: additional curl options passed to getURI
#
# returns NULL if async=TRUE or the return value of RCurl getURI otherwise
GET = function(resource, args=list(), header=TRUE, async=FALSE, interrupt=FALSE, ...)
{
  if(!(substr(resource,1,1)=="/")) resource = paste("/",resource,sep="")
  uri = URI(resource, args)
  uri = oldURLencode(uri)
  uri = gsub("\\+","%2B",uri,perl=TRUE)
  on.exit(sigint(SIG_DFL)) # Reset signal handler on exit
  callback = curl_signal_trap
  interrupt = interrupt && scidb.interrupt()   # Check package override
  if(interrupt)
  {
    sigreset ()
    sigint(TRAP())          # Handle SIGINT
  } else
  {
    sigint(SIG_IGN)          # Ignore SIGINT
  }
#  digest = curlopts()
  digest = digest_auth("GET",uri)
  if(async)
  {
    getURI(url=uri,
      .opts=c(list(header=header,
                 'ssl.verifyhost'=as.integer(options("scidb.verifyhost")),
                 'ssl.verifypeer'=0,
                 noprogress=!interrupt, progressfunction=curl_signal_trap),
                 list(...), httpheader=digest),
      async=TRUE)
    return(NULL)
  }
# The gsub is here because basic digest authentication might return TWO
# concatenated responses: the first a 401 auth error as the digest
# authentication is negotiated and the 2nd potentially the desired response (or
# maybe another error). We strip everything up to the last response with gsub.
  gsub(".*\\r\\n\\r\\nHTTP","HTTP",getURI(url=uri, .opts=c(list(header=header,'ssl.verifyhost'=as.integer(options("scidb.verifyhost")),'ssl.verifypeer'=0, noprogress=!interrupt, progressfunction=curl_signal_trap),list(...),httpheader=digest)))
}

# Check if array exists
.scidbexists = function(name)
{
  Q = scidblist()
  return(name %in% Q)
}

# scidblist: A convenience wrapper for list().
# Input:
# type (character), one of the indicated list types
# verbose (boolean), include attribute and dimension data when type="arrays"
# n: maximum lines of output to return
# Output:
# A list.
scidblist = function(pattern,
type= c("arrays","operators","functions","types","aggregates","instances","queries","libraries"),
              verbose=FALSE, n=Inf)
{
  type = match.arg(type)
  Q = iquery(paste("list('",type,"')",sep=""), return=TRUE, n=n)

  if(dim(Q)[1]==0) return(NULL)
  z=Q[,-1,drop=FALSE]
  if(type=="arrays" && !verbose) {
    z=z[,1]
    if(!missing(pattern))
      z = grep(z,pattern=pattern,value=TRUE)
  }
  z
}
scidbls = function(...) scidblist(...)

# Basic low-level query. Returns query id. This is an internal function.
# query: a character query string
# afl: TRUE indicates use AFL, FALSE AQL
# async: Does not do anything right now. Maybe in the future.
# save: Save format query string or NULL. If async=TRUE, save is ignored.
# release: Set to zero preserve web session until manually calling release_session
# session: if you already have a SciDB http session, set this to it, otherwise NULL
# resp(logical): return http response
# interrupt: Set to TRUE to enable user interrupt
# stream: Set to 0L or 1L to control streaming, otherwise use package options
# Example values of save:
# save="dcsv"
# save="lcsv+"
# save="(double NULL, int32)"
#
# Returns the HTTP session in each case
scidbquery = function(query, afl=TRUE, async=FALSE, save=NULL,
                      release=1, session=NULL, resp=FALSE, interrupt=FALSE, stream)
{
  DEBUG = FALSE
  STREAM = 0L
  if(!is.null(options("scidb.debug")[[1]]) && TRUE==options("scidb.debug")[[1]]) DEBUG=TRUE
  if(missing(stream))
  {
    if(!is.null(options("scidb.stream")[[1]]) && TRUE==options("scidb.stream")[[1]]) STREAM=1L
  } else STREAM = as.integer(stream)
  sessionid=session
  if(is.null(session))
  {
# Obtain a session from shim
    sessionid = getSession()
  }
  if(is.null(save)) save=""
  if(DEBUG)
  {
    cat(query, "\n")
    t1=proc.time()
  }
  ans = tryCatch(
    {
      if(is.null(save))
        GET("/execute_query",list(id=sessionid,release=release,
             query=query,afl=as.integer(afl), stream=0L),
             interrupt=interrupt)
      else
        GET("/execute_query",list(id=sessionid,release=release,
            save=save,query=query,afl=as.integer(afl),stream=STREAM), 
            interrupt=interrupt)
    }, error=function(e)
    {
# User cancel? Note not interruptable.
      GET("/cancel", list(id=sessionid,async='1'), async=TRUE)
      GET("/release_session", list(id=sessionid))
      "HTTP/1.0 206 ERROR"
    }, interrupt=function(e)
    {
      GET("/cancel", list(id=sessionid,async='1'), async=TRUE)
      GET("/release_session", list(id=sessionid))
      "HTTP/1.0 206 ERROR"
    })
  w = ans
  err = 200
  if(nchar(w)>9)
    err = as.integer(strsplit(substr(w,1,20)," ")[[1]][[2]])
  if(err==206)
  {
    stop("HTTP request canceled\n")
  } else if(err>206)
  {
    stop(strsplit(w,"\r\n\r\n")[[1]][[2]])
  }
  if(DEBUG) cat("Query time",(proc.time()-t1)[3],"\n")
  if(resp) return(list(session=sessionid, response=ans))
  sessionid
}

# scidbremove: Convenience function to remove one or more scidb arrays
# Input:
# x (character): a vector or single character string listing array names
# error (function): error handler. Use stop or warn, for example.
# async (optional boolean): If TRUE use expermental shim async option for speed
# force (optional boolean): If TRUE really remove this array, even if scidb.safe_remove=TRUE
# recursive (optional boolean): If TRUE, recursively remove this array and its dependency graph
# Output:
# null
scidbremove = function(x, error=warning, async, force, warn=TRUE, recursive=FALSE) UseMethod("scidbremove")
scidbremove.default = function(x, error=warning, async, force, warn=TRUE, recursive=FALSE)
{
  if(is.null(x)) return(invisible())
  if(missing(async)) async=FALSE
  if(missing(force)) force=FALSE
  errfun=error

  safe = options("scidb.safe_remove")[[1]]
  if(is.null(safe)) safe = TRUE
  if(!safe) force=TRUE
  uid = get("uid",envir=.scidbenv)
  if(is.scidb(x) || is.scidbdf(x)) x = list(x)
  for(y in x)
  {
# Depth-first, a bad way to traverse this XXX improve
    if(recursive && (is.scidb(y) || is.scidbdf(y)) && !is.null(unlist(y@gc$depend)))
    {
      for(z in y@gc$depend) scidbremove(z,error,async,force,warn,recursive)
    }
    if(is.scidb(y) || is.scidbdf(y)) y = y@name
    if(grepl("\\(",y)) next  # Not a stored array
    query = sprintf("remove(%s)",y)
    if(grepl(sprintf("^R_array.*%s$",uid),y,perl=TRUE))
    {
      tryCatch( scidbquery(query,async=async, release=1, stream=0L),
                error=function(e) if(!recursive && warn)errfun(e))
    } else if(force)
    {
      tryCatch( scidbquery(query,async=async, release=1, stream=0L),
                error=function(e) if(!recursive && warn)errfun(e))
    } else if(warn)
    {
      warning("The array ",y," is protected from easy removal. Specify force=TRUE if you really want to remove it.")
    }
  }
  invisible()
}
scidbrm = function(x,error=warning,...) scidbremove(x,error,...)

# df2scidb: User function to send a data frame to SciDB
# Returns a scidbdf object
df2scidb = function(X,
                    name=tmpnam(),
                    dimlabel="row",
                    chunkSize,
                    rowOverlap=0L,
                    types=NULL,
                    nullable,
                    schema_only=FALSE,
                    gc, start)
{
  if(!is.data.frame(X)) stop("X must be a data frame")
  if(missing(start)) start=1
  start = as.numeric(start)
  if(missing(gc)) gc=FALSE
  if(missing(nullable)) nullable = TRUE
  if(length(nullable)==1) nullable = rep(nullable, ncol(X))
  if(length(nullable)!=ncol(X)) stop ("nullable must be either of length 1 or ncol(X)")
  if(!is.null(types) && length(types)!=ncol(X)) stop("types must match the number of columns of X")
  n1 = nullable
  nullable = rep("",ncol(X))
  if(any(n1)) nullable[n1] = "NULL"
  anames = make.names(names(X),unique=TRUE)
  anames = gsub("\\.","_",anames,perl=TRUE)
  if(length(anames)!=ncol(X)) anames=make.names(1:ncol(X))
  if(!all(anames==names(X))) warning("Attribute names have been changed")
# Check for attribute/dimension name conflict
  old_dimlabel = dimlabel
  dimlabel = tail(make.names(c(anames,dimlabel),unique=TRUE),n=1)
  dimlabel = gsub("\\.","_",dimlabel,perl=TRUE)
  if(dimlabel!=old_dimlabel) warning("Dimension name has been changed")
  if(missing(chunkSize)) {
    chunkSize = min(nrow(X),10000)
  }
  chunkSize = as.numeric(chunkSize)
  m = ceiling(nrow(X) / chunkSize)

# Default type is string
  typ = rep(paste("string",nullable),ncol(X))
  args = "<"
  if(!is.null(types)) {
    for(j in 1:ncol(X)) typ[j]=paste(types[j],nullable[j])
  } else {
    for(j in 1:ncol(X)) {
      if("numeric" %in% class(X[,j])) 
        typ[j] = paste("double",nullable[j])
      else if("integer" %in% class(X[,j])) 
        typ[j] = paste("int32",nullable[j])
      else if("logical" %in% class(X[,j])) 
        typ[j] = paste("bool",nullable[j])
      else if("character" %in% class(X[,j])) 
        typ[j] = paste("string",nullable[j])
      else if("factor" %in% class(X[,j])) 
      {
        if("scidb_factor" %in% class(X[,j]))
        {
          typ[j] = paste("int64 NULL")
          Xj = X[,j]
          levels(Xj) = levels_scidb(Xj)
          X[,j] = as.vector(Xj)
        }
        else
          typ[j] = paste("string",nullable[j])
      }
      else if("POSIXct" %in% class(X[,j])) 
      {
        warning("Converting R POSIXct to SciDB datetime as UTC time. Subsecond times rounded to seconds.")
#        X[,j] = as.integer(X[,j])
        X[,j] = format(X[,j],tz="UTC")
        typ[j] = paste("datetime",nullable[j])
      }
    }  
  }
  for(j in 1:ncol(X)) {
    args = paste(args,anames[j],":",typ[j],sep="")
    if(j<ncol(X)) args=paste(args,",",sep="")
    else args=paste(args,">")
  }

  SCHEMA = paste(args,"[",dimlabel,"=",sprintf("%.0f",start),":",sprintf("%.0f",(nrow(X)+start-1)),",",sprintf("%.0f",chunkSize),",", rowOverlap,"]",sep="")

  if(schema_only) return(SCHEMA)

# Obtain a session from the SciDB http service for the upload process
  session = getSession()
  on.exit(GET("/release_session",list(id=session)) ,add=TRUE)

# Create SciDB input string from the data frame
  scidbInput = .df2scidb(X,chunkSize,start)

# Post the input string to the SciDB http service
  uri = URI("upload_file",list(id=session))
  hdr = digest_auth("POST",uri)
  tmp = tryCatch(postForm(uri=uri, uploadedfile=fileUpload(contents=scidbInput,filename="scidb",contentType="application/octet-stream"),.opts=curlOptions(httpheader=hdr,'ssl.verifyhost'=as.integer(options("scidb.verifyhost")),'ssl.verifypeer'=0)), error=invisible)
  tmp = tmp[[1]]
  tmp = gsub("\r","",tmp)
  tmp = gsub("\n","",tmp)

# Load query
  query = sprintf("store(input(%s, '%s'),%s)",SCHEMA, tmp, name)
  scidbquery(query, async=FALSE, release=1, session=session, stream=0L)
  scidb(name,`data.frame`=TRUE,gc=gc)
}

# Sparse matrix to SciDB
.Matrix2scidb = function(X,name,rowChunkSize=1000,colChunkSize=1000,start=c(0,0),gc=TRUE,...)
{
  D = dim(X)
  N = Matrix::nnzero(X)
  rowOverlap=0L
  colOverlap=0L
  if(length(start)<1) stop ("Invalid starting coordinates")
  if(length(start)>2) start = start[1:2]
  if(length(start)<2) start = c(start, 0)
  start = as.integer(start)
  type = .scidbtypes[[typeof(X@x)]]
  if(is.null(type)) {
    stop(paste("Unupported data type. The package presently supports: ",
       paste(.scidbtypes,collapse=" "),".",sep=""))
  }
  if(type!="double") stop("Sorry, the package only supports double-precision sparse matrices right now.")
  schema = sprintf(
      "< val : %s null>  [i=%.0f:%.0f,%.0f,%.0f, j=%.0f:%.0f,%.0f,%.0f]", type, start[[1]],
      nrow(X)-1+start[[1]], min(nrow(X),rowChunkSize), rowOverlap, start[[2]], ncol(X)-1+start[[2]],
      min(ncol(X),colChunkSize), colOverlap)
  schema1d = sprintf("<i:int64 null, j:int64 null, val : %s null>[idx=0:*,100000,0]",type)
# Obtain a session from shim for the upload process
  session = getSession()
  if(length(session)<1) stop("SciDB http session error")
  on.exit(GET("/release_session",list(id=session)) ,add=TRUE)

# Compute the indices and assemble message to SciDB in the form
# double,double,double for indices i,j and data val.
  dp = diff(X@p)
  j  = rep(seq_along(dp),dp) - 1
# Upload the data
# XXX I couldn't get RCurl to work using the fileUpload(contents=x), with 'x'
# a raw vector. But we need RCurl to support SSL. As a work-around, we save
# the object. This copy sucks and must be fixed. XXX problem is in RCurl, see
# scidb-functions.R.
  fn = tempfile()
  bytes = writeBin(.Call("scidb_raw",as.vector(t(matrix(c(X@i + start[[1]],j + start[[2]], X@x),length(X@x)))),PACKAGE="scidb"),con=fn)
  url = URI("/upload_file",list(id=session))
  hdr = digest_auth("POST",url)
  ans = tryCatch(postForm(uri = url, uploadedfile = fileUpload(filename=fn), .opts=curlOptions(httpheader=hdr,'ssl.verifyhost'=as.integer(options("scidb.verifyhost")),'ssl.verifypeer'=0)), error=invisible)
  unlink(fn)
  ans = ans[[1]]
  ans = gsub("\r","",ans)
  ans = gsub("\n","",ans)

  query = sprintf("store(redimension(input(%s,'%s',0,'(double null,double null, double null)'),%s),%s)",schema1d, ans, schema, name)
  iquery(query)
  scidb(name,gc=gc)
}


iquery = function(query, `return`=FALSE,
                  afl=TRUE, iterative=FALSE,
                  n=10000, excludecol, binary=FALSE, ...)
{
  DEBUG = FALSE
  if(!is.null(options("scidb.debug")[[1]]) && TRUE==options("scidb.debug")[[1]]) DEBUG=TRUE
  if(!afl && `return`) stop("return=TRUE may only be used with AFL statements")
  if(iterative && !`return`) stop("Iterative result requires return=TRUE")
  if(is.scidb(query) || is.scidbdf(query))  query=query@name
  if(missing(excludecol)) excludecol=NA
  if(iterative)
  {
    sessionid = scidbquery(query,afl,async=FALSE,save="lcsv+",release=0,interrupt=TRUE)
    return(iqiter(con=sessionid,n=n,excludecol=excludecol,...))
  }
  qsplit = strsplit(query,";")[[1]]
  m = 1
  if(n==Inf) n = -1    # Indicate to shim that we want all the lines of output
  for(query in qsplit)
  {
# Only return the last query, mimicking command-line iquery.
    if(`return` && m==length(qsplit))
    {
      if(binary)
      {
        return(scidb_unpack_to_dataframe(query,...))
      }
      ans = tryCatch(
       {
        sessionid = scidbquery(query,afl,async=FALSE,save="lcsv+",release=0, interrupt=TRUE)
        dt1 = proc.time()
        result = tryCatch(
          {
            GET("/read_lines",list(id=sessionid,n=as.integer(n+1)),header=FALSE,interrupt=TRUE)
          },
          error=function(e)
          {
             GET("/cancel",list(id=sessionid))
             GET("/release_session",list(id=sessionid))
             stop(e)
          })
        GET("/release_session",list(id=sessionid))
        if(DEBUG) cat("Data transfer time",(proc.time()-dt1)[3],"\n")
        dt1 = proc.time()
# Handle escaped quotes
        result = gsub("\\\\'","''",result, perl=TRUE)
        result = gsub("\\\\\"","''",result, perl=TRUE)
# Map SciDB missing (aka null) to NA, but preserve DEFAULT null.
# This sucks, need to avoid this parsing and move on to binary xfer.
        result = gsub("DEFAULT null","@#@#@#kjlkjlkj@#@#@555namnsaqnmnqqqo",result,perl=TRUE)
        result = gsub("null","NA",result, perl=TRUE)
        result = gsub("@#@#@#kjlkjlkj@#@#@555namnsaqnmnqqqo","DEFAULT null",result,perl=TRUE)
        val = textConnection(result)
        ret = tryCatch({
                read.table(val,sep=",",stringsAsFactors=FALSE,header=TRUE,...)},
                error=function(e){ warning(e);c()})
        close(val)
        if(DEBUG) cat("R parsing time",(proc.time()-dt1)[3],"\n")
        ret
       }, error = function(e)
           {
             stop(e)
           })
    } else
    {
      ans = scidbquery(query, afl, async=FALSE, release=1, interrupt=TRUE, stream=0L)
    }
    m = m + 1
  }
  if(!(`return`)) return(invisible())
  ans
}

iqiter = function (con, n = 1, excludecol, ...)
{
  if(missing(excludecol)) excludecol=NA
  dostop = function(s=TRUE)
  {
    GET("/release_session",list(id=con))
    if(s) stop("StopIteration",call.=FALSE)
  }
  if (!is.numeric(n) || length(n) != 1 || n < 1)
    stop("n must be a numeric value >= 1")
  init = TRUE
  header = c()
  nextEl = function() {
    if (is.null(con)) dostop()
    if(init) {
      ans = tryCatch(
       {
        result = GET("/read_lines",list(id=con,n=n),header=FALSE)
        val = textConnection(result)
        ret=read.table(val,sep=",",stringsAsFactors=FALSE,header=TRUE,nrows=n,...)
        close(val)
        ret
       }, error = function(e) {dostop()},
          warning = function(w) {dostop()}
      )
      header <<- colnames(ans)
      init <<- FALSE
    } else {
      ans = tryCatch(
       {
        result = GET("/read_lines",list(id=con,n=n),header=FALSE)
        val = textConnection(result)
        ret = read.table(val,sep=",",stringsAsFactors=FALSE,header=FALSE,nrows=n,...)
        close(val)
        ret
       }, error = function(e) {dostop()},
          warning = function(w) {dostop()}
      )
      colnames(ans) = header
    }
    if(!is.na(excludecol) && excludecol<=ncol(ans))
    {
      rownames(ans) = ans[,excludecol]
      ans=ans[,-excludecol, drop=FALSE]
    }
    ans
  }
  it = list(nextElem = nextEl, gc=new.env())
  class(it) = c("abstractiter", "iter")
  it$gc$remove = TRUE
  reg.finalizer(it$gc, function(e) if(e$remove) 
                  tryCatch(dostop(FALSE),error=function(e) stop(e)),
                  onexit=TRUE)
  it
}


# Utility csv2scidb function
.df2scidb = function(X, chunksize, start=1)
{
  if(missing(chunksize)) chunksize = min(nrow(X),10000L)
#  scipen = options("scipen")
#  options(scipen=20)
#  buf = capture.output(
#         write.table(X, file=stdout(), sep=",",
#                     row.names=FALSE,col.names=FALSE,quote=FALSE)
#        )
#  options(scipen=scipen)
#  x = sapply(1:length(buf), function(j)
#    {
#      chunk = floor(j/chunksize)
#      if(j==length(buf))
#      {
#        if((j-1) %% chunksize == 0)
#          tmp = sprintf("{%.0f}[\n(%s)];",start + chunk*chunksize, buf[j])
#        else
#          tmp = sprintf("(%s)];",buf[j])
#      } else if((j-1) %% chunksize==0)
#      {
#        tmp = sprintf("{%.0f}[\n(%s),",start + chunk*chunksize, buf[j])
#      } else if((j) %% chunksize == 0)
#      {
#        tmp = sprintf("(%s)];",buf[j])
#      } else
#      {
#        tmp = sprintf("(%s),",buf[j])
#      }
#      tmp
#    }
#  )
#  x = paste(x,collapse="")
#  ans = paste(x,collapse="")

  ans = .Call("df2scidb",X,as.integer(chunksize), as.double(round(start)), "%.15f",PACKAGE="scidb")
  ans
}

# Internal utility function, make every attribute of an array nullable
make_nullable = function(x)
{
  cast(x,sprintf("%s%s",build_attr_schema(x,nullable=TRUE),build_dim_schema(x)))
}

noE = function(w) sapply(w,
  function(x)
  {
    if(is.character(x)) return(x)
    sprintf("%.0f",x)
  })

# If a scidb object has the "sparse" attribute set, return that. Otherwise
# interrogate the backing array to determine if it's sparse or not.  Set
# count=FALSE to skip the count and return TRUE if the attribute is not set.
is.sparse = function(x, count=TRUE)
{
  ans = attr(x,"sparse")
  if(is.null(ans) && !count) return(TRUE)
  if(is.null(ans))
  {
    return(count(x) < prod(dim(x)))
  }
  ans
}

# Check for scidb missing flag
is.nullable = function(x)
{
  any(scidb_nullable(x))
}

# Returns TRUE if version string x is greater than or equal to than version y
compare_versions = function(x,y)
{
  b = as.numeric(strsplit(as.character(x),"\\.")[[1]])
  a = as.numeric(strsplit(as.character(y),"\\.")[[1]])
  ans = b[1] > a[1]
  if(b[1] == a[1]) ans = b[2] >= a[2]
  ans
}

# Reset array coordinate system to zero-indexed origin
origin = function(x)
{
  N = paste(rep("null",2*length(dimensions(x))),collapse=",")
  GEMM.BUG = ifelse(is.logical(options("scidb.gemm_bug")[[1]]),options("scidb.gemm_bug")[[1]],FALSE)
  if(GEMM.BUG) query = sprintf("sg(subarray(%s,%s),1,-1)",x@name,N)
  else query = sprintf("subarray(%s,%s)",x@name,N)
  .scidbeval(query,`eval`=FALSE,depend=list(x),gc=TRUE)
}


chunk_map = function(`array`)
{
  if(!missing(`array`))
  {
    if(is.scidb(`array`) || is.scidbdf(`array`))
    {
      `array` = schema(`array`)
    }
    q1="redimension( apply(list('chunk map'), element_count, int64(nelem), iid, int64(inst), aid, int64(arrid)), <element_count:int64, iid:int64>[aid=0:*,1000,0,i=0:*,1000000000,0])"
    q2=sprintf("redimension( filter(apply( list('arrays', true), aid, int64(id), match, regex(name,'%s@*[0-9]*')), match=true), <name: string null> [aid = 0:*,1000,0])", `array`)
    query=sprintf("project( cross_join(%s as X, %s as Y, X.aid, Y.aid), element_count, iid)", q1, q2)
    return(iquery(query,`return`=TRUE)[,2:4])
  }


iquery("
project(
 cross_join(
  redimension(
   apply(filter(list('chunk map'), inst=instn and attid=0), iid, int64(inst), aid, int64(arrid)),
   <nchunks:uint64 null,
    min_ccnt:uint32 null,
    avg_ccnt:double null,
    max_ccnt:uint32 null,
    total_cnt: uint64 null>
   [iid = 0:*,1000,0, aid= 0:*,1000,0],
   count(*) as nchunks,
   min(nelem) as min_ccnt,
   avg(nelem) as avg_ccnt,
   max(nelem) as max_ccnt,
   sum(nelem) as total_cnt
  ) as A,
  redimension(
   apply( list('arrays', true), aid, int64(id)),
   <name: string null>
   [aid = 0:*,1000,0]
  ) as B,
  A.aid, B.aid
 ),
 name, nchunks, min_ccnt, avg_ccnt, max_ccnt, total_cnt
)",`return`=TRUE)
}

# Map scidbdf object column classes into R, adding an extra integer at the start for the index!
scidbdfcc = function(x)
{
  if(!is.null(options("scidb.test")[[1]]))
  {
    cat("Using old method for data.frame import")
    return(NA)
  }
  c("integer",as.vector(unlist(lapply(.scidbdftypes[scidb_types(x)],function(x) ifelse(is.null(x),NA,x)))))
}

.scidbstr = function(object)
{
  name = substr(object@name,1,35)
  if(nchar(object@name)>35) name = paste(name,"...",sep="")
  dn = "\nDimension:"
  if(length(dimensions(object))>1) dn = "\nDimensions:"
  cat("SciDB expression: ",name)
  cat("\nSciDB schema: ",schema(object))
  cat("\n\nAttributes:\n")
  cat(paste(capture.output(print(data.frame(attribute=object@attributes,type=scidb_types(object),nullable=scidb_nullable(object)))),collapse="\n"))
  cat(dn,"\n")
  bounds = scidb_coordinate_bounds(object)
  cat(paste(capture.output(print(data.frame(dimension=dimensions(object),start=bounds$start,end=bounds$end,chunk=scidb_coordinate_chunksize(object),row.names=NULL,stringsAsFactors=FALSE))),collapse="\n"))
  cat("\n")
}

# The RCurl package does not handle interruption well, sometimes segfaulting.
#
# Here is an example:
#
#    sigint(TRAP())
#    x=getBinaryURL(url,
#        .opts=list(noprogress=FALSE, progressfunction=curl_signal_trap))
#    sigint(SIG_DFL)
#
# Note that users might have to hold down CTRL+C a while because of the tiny
# time window it's available.
sigint = function(val)
{
  .Call("sig",as.integer(val),PACKAGE="scidb")
}
sigstate = function()
{
  .Call("state",PACKAGE="scidb")
}
sigreset = function()
{
  .Call("reset",PACKAGE="scidb")
}
curl_signal_trap = function(down,up)
{
# Unforunately, RStudio doesn't seem to let us set up custom signal handlers.
  if(TRAP() == SIG_IGN)
  {
    k = tryCatch(
    {
      sigint(SIG_DFL)    # Enable SIGINT
      Sys.sleep(0.02)    # Just plain ugly
      sigint(SIG_IGN)    # Disable
      0L
    }, interrupt = function(e)
    {
      cat("Canceling...\n",file=stderr())
      1L
    })
    return(k)
  }
# Fast custom signal handler
  k = sigstate()
  if(k>0)
  {
    cat("Canceling...\n",file=stderr())
    return(1L)
  }
  0L
}

# Walk the dependency graph, setting all upstreams array to persist
# Define a generic persist
persist = function(x, remove=FALSE, ...) UseMethod("persist")
persist.default = function(x, remove=FALSE, ...)
{
  DEBUG = FALSE
  if(!is.null(options("scidb.debug")[[1]]) && TRUE==options("scidb.debug")[[1]]) DEBUG=TRUE
  if(!(is.scidb(x) || is.scidbdf(x))) return(invisible())
  if(DEBUG) cat("Persisting ",x@name,"\n")
  x@gc$remove = remove
  if(is.null(unlist(x@gc$depend))) return()
  for(y in x@gc$depend)
  {
    if(DEBUG) cat("Persisting ",y@name,"\n")
    y@gc$remove = remove
    if(!is.null(y@gc$depend)) persist(y, remove=remove)
  }
}

# A special persist function for complicated model objects
persist.glm_scidb = function(x, remove=FALSE, ...)
{
  .traverse.glm_scidb(x, persist, remove)
}
# A special remove function for complicated model objects
scidbremove.glm_scidb = function(x, error=warning, async=FALSE, force=FALSE, warn=TRUE, recursive=TRUE)
{
  .traverse.glm_scidb(x, scidbremove, error, async, force, warn, recursive)
}

# Show the repository log (not in namespace)
show_commit_log = function()
{
  log = system.file("misc/log",package="scidb")
  if(file.exists(log))
  {
    file.show(log)
  }
}

# This is a workhorse SciDB data munger
# query is a query string, not a scidb object.
scidb_unpack_to_dataframe = function(query, ...)
{
  DEBUG = FALSE
  if(!is.null(options("scidb.debug")[[1]]) && TRUE==options("scidb.debug")[[1]]) DEBUG=TRUE
  buffer = 100000L
  row_names = NULL
  args = list(...)
  interrupt = scidb.interrupt()
  if(!is.null(args$buffer))
  {
    argsbuf = tryCatch(as.integer(args$buffer),warning=function(e) NA)
# Potential problem here if args$buffer is bogus or '*' or whatever.
# Deal with it. XXX FIX ME
    if(!is.na(argsbuf) && argsbuf < 100e6) buffer = as.integer(argsbuf)
  }
  if(!is.null(args$row.names)) row_names = args$row.names
  up = TRUE
  if(!is.null(args$unpack)) up = args$unpack
  if(up) x = scidb(sprintf("unpack(%s,row)",query), data.frame=TRUE)
  else x = scidb(sprintf("project(unpack(%s,row),%s)",query,scidb_attributes(x)), data.frame=TRUE)
  N = scidb_nullable(x)
  TYPES = scidb_types(x)
  ns = rep("",length(N))
  ns[N] = "null"
  format_string = paste(paste(TYPES,ns),collapse=",")
  format_string = sprintf("(%s)",format_string)
  sessionid = scidbquery(x@name, save=format_string, async=FALSE, release=0, interrupt=TRUE)
# Release the session on exit
  on.exit( GET("/release_session",list(id=sessionid)) ,add=TRUE)

  dt2 = proc.time()
  r = URI("/read_bytes",list(id=sessionid,n=0))
  sigint(TRAP())
  BUF = tryCatch(
        { 
          getBinaryURL(r, .opts=c(curlopts(),list('ssl.verifyhost'=as.integer(options("scidb.verifyhost")),'ssl.verifypeer'=0, noprogress=!interrupt, progressfunction=curl_signal_trap)))
        }, error=function(e)
        { 
          GET("/release_session",list(id=sessionid))
          sigint(SIG_DFL)
          stop("canceled")
        })
  sigint(SIG_DFL)
  if(DEBUG) cat("Data transfer time",(proc.time()-dt2)[3],"\n")
  dt1 = proc.time()
  len = length(BUF)
  p = 0
  ans = c()
  cnames = c(names(x),"lines","p")
  rnames = c()
  while(p<len)
  {
dt2 = proc.time()
    n = ncol(x)
    tmp   = .Call("scidb_parse", as.integer(buffer), TYPES, N, BUF, as.double(p))
    names(tmp) = cnames
    lines = tmp[[n+1]]
    p_old = p
    p     = tmp[[n+2]]
if(DEBUG) cat("  R buffer ",p,"/",len," bytes parsing time",(proc.time()-dt2)[3],"\n")
dt2 = proc.time()
    if(lines>0)
    {
      len_out = length(tmp[[1]])
      if(!is.null(row_names))
      {
        if(is.numeric(row_names) && length(row_names)==1)
        {
          rnames = c(rnames,tmp[[row_names]][1:lines])
          tmp = tmp[-row_names]
          n = n -1
        }
      }
      if(lines < len_out) tmp = lapply(tmp[1:n],function(x) x[1:lines])
# Let's adaptively re-estimate a buffer size
      avg_bytes_per_line = ceiling((p - p_old)/lines)
      buffer = ceiling(1.3*(len - p)/avg_bytes_per_line) # Engineering factor
# Assemble the data frame
      if(is.null(ans)) ans = data.frame(tmp[1:n],stringsAsFactors=FALSE)
      else ans = rbind(ans,data.frame(tmp[1:n],stringsAsFactors=FALSE))
    }
if(DEBUG) cat("  R rbind/df assembly time",(proc.time()-dt2)[3],"\n")
  }
dt2 = proc.time()
  if(!is.null(row_names))
  {
    if(is.numeric(row_names) && length(row_names)==1)
    {
      rownames(ans) = rnames
    } else
    {
      rownames(ans) = row_names
    }
  }
if(DEBUG) cat("  R rownames assignment",(proc.time()-dt2)[3],"\n")
  if(DEBUG) cat("Total R parsing time",(proc.time()-dt1)[3],"\n")
  ans
}

# Convenience function for digest authentication.
digest_auth = function(method, uri, realm="", nonce="123456")
{
  if(exists("authtype",envir=.scidbenv))
  {
   if(.scidbenv$authtype != "digest") return(c(Expect=""))
  }
  uri = gsub(".*/","/",uri)
  userpwd = .scidbenv$digest
  if(is.null(userpwd)) userpwd=":"
  up = strsplit(userpwd,":")[[1]]
  user = up[1]
  pwd  = up[2]
  if(is.na(pwd)) pwd=""
  ha1=digest(sprintf("%s:%s:%s",user,realm,pwd,algo="md5"),serialize=FALSE)
  ha2=digest(sprintf("%s:%s", method,  uri, algo="md5"),serialize=FALSE)
  cnonce="MDc1YmFhOWFkY2M0YWY2MDAwMDBlY2JhMDAwMmYxNTI="
  nc="00000001"
  qop="auth"
  response=digest(sprintf("%s:%s:%s:%s:%s:%s",ha1,nonce,nc,cnonce,qop,ha2),algo="md5",serialize=FALSE)
  sprintf('Authorization: Digest username="%s", realm=%s, nonce="%s", uri="%s", cnonce="%s", nc=%s, qop=%s, response="%s"', user, realm, nonce, uri, cnonce, nc, qop, response)
}

# Return extra curl options defined by the active connection
curlopts = function()
{
  ans = c()
  if(exists("digest",envir=.scidbenv))
    ans=c(ans,list(userpwd=get("digest",envir=.scidbenv)))
  ans
}

# Internal warning function
warnonce = (function() {
  state = list(
    unpack="The array contains multiple SciDB attributes, returning as an unpacked dataframe.",
    missing="An array might contain missing values (R NA/SciDB 'null'-able attribute).\n Missing values are not yet understood by SciDB multiplication operators\n and any missing values have been replaced with zero (see replaceNA).",
    nonum="Note: The R sparse Matrix package does not support certain value types like\ncharacter strings",
    toobig="Dimensions are too big! Returning data in unpacked data.frame form."
  )
  function(warn) {
    if(!is.null(state[warn][[1]])) {
      warning(state[warn], call.=FALSE)
      s <<- state
      s[warn] = c()
      state <<- s
    }
  }
}) () # oh boy


# Utility function that parses an aggregation expression
# makes sense of statements like
# max(a) as amax, avg(b) as bmean, ...
# Input
# x: scidb object
# expr: aggregation expression (as is required)
# Output
# New attribute schema
aparser = function(x, expr)
{
  map = list(
    ApproxDC    = function(x) "double null",
    avg         = function(x) "double null",
    count       = function(x) "uint64 null",
    first_value = function(x) sprintf("%s null",x),
    last_value  = function(x) sprintf("%s null",x),
    mad         = function(x) sprintf("%s null",x),
    max         = function(x) sprintf("%s null",x),
    min         = function(x) sprintf("%s null",x),
    median      = function(x) sprintf("%s null",x),
    prod        = function(x) sprintf("%s null",x),
    stdev       = function(x) "double null",
    sum         = function(x) sprintf("%s null",x),
    top_five    = function(x) sprintf("%s null",x),
    var         = function(x) "double null")

  p = strsplit(expr,",")[[1]]                      # comma separated
  p = lapply(p,function(v) strsplit(v,"as")[[1]])  # as is required
  new_attrs = lapply(p, function(v) v[[2]])        # output attribute names
  y = lapply(p,function(v)strsplit(v,"\\(")[[1]])  #
  fun = lapply(y, function(v) gsub(" ","",v[[1]])) # functions
  old_attrs = unlist(lapply(y,function(v)gsub("[\\) ]","",v[[2]])))
  i = 1:length(x@attributes)
  names(i)=x@attributes
  old_types = scidb_types(x)[i[old_attrs]]
  old_types[is.na(old_types)] = "int64"    # catch all
  z = rep("",length(new_attrs))
  for(j in 1:length(z))
  {
    z[j] = sprintf("%s:%s", new_attrs[j],map[unlist(fun)][[j]](old_types[j]))
  }
  sprintf("<%s>", paste(z, collapse=","))
}

# Some versions of RCurl seem to contain a broken URLencode function.
oldURLencode = function (URL, reserved = FALSE) 
{
    OK <- paste0("[^-ABCDEFGHIJKLMNOPQRSTUVWXYZ", "abcdefghijklmnopqrstuvwxyz0123456789$_.+!*'(),", 
        if (!reserved) 
            ";/?:@=&", "]")
    x <- strsplit(URL, "")[[1L]]
    z <- grep(OK, x)
    if (length(z)) {
        y <- sapply(x[z], function(x) paste0("%", as.character(charToRaw(x)), 
            collapse = ""))
        x[z] <- y
    }
    paste(x, collapse = "")
}

