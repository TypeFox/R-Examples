#' Open connection to a RethinkDB server
#'
#' Opens connection to a given RethinkDB server.
#'
#' @param host Host to connect to.
#' @param port Port to connect to.
#' @param authKey Authentication key. Not supported yet.
#' @param v Protocol version; \code{"V0_3"} and \code{"V0_4"} supported, the default should be used.
#' @return Object of a class \code{RethinkDB_connection}, which can be passed to \code{r()$...$run} and \code{r()$...$runAsync} functions.
#' @export
openConnection<-function(host='localhost',port=28015,authKey=NULL,v="V0_4"){
 JSON<-0x7e6970c7
 V0_3<-0x5f75e83e
 V0_4<-0x400c2d20

 if(!is.null(authKey))
  stop("Auth key not supported yet!");

 header<-NULL;
 if(identical(v,"V0_3")){
  header<-c(V0_3,0L,JSON);
 }
 if(identical(v,"V0_4")){
  header<-c(V0_4,0L,JSON);
 }
 if(is.null(header)) stop(sprintf("Unknown protocol version %s!",v));

 #Open connection; R will throw error in case of trouble here
 socketConnection(host,port,
  open='w+b',blocking=TRUE)->con;

 #Send handshake
 writeBin(as.integer(header),con,size=4,endian='little');

 while(TRUE){
  response<-readBin(con,character(),1);
  if(length(response)>0) break;
 }
 if(!identical(response,"SUCCESS"))
  stop(sprintf(
   "Could not establish connection with RethinkDB server.\nIt says: %s",
   response));

 ans<-new.env();
 ans$con<-con;
 ans$host<-host;
 ans$port<-port;
 ans$ver<-header[1];
 ans$handlers<-list();
 class(ans)<-"RethinkDB_connection";
 ans
}

#' Check if connection is opened
#'
#' Check whether a given connection is opened.
#' Closed connections cannot be used, and will throw errors on such an attempt; their associated callbacks and/or sync cursor are dead and won't fire/produce any more data.
#'
#' @param x Connection to to check.
#' @return \code{TRUE} if connection is opened and can be used for queries, \code{FALSE} otherwise.
#' @author Miron B. Kursa
#' @export
isOpened<-function(x){
 stopifnot(inherits(x,"RethinkDB_connection"));
 opened<-try(isOpen(x$con),silent=TRUE);
 if(inherits(opened,"try-error")){
  FALSE;
 }else{
  opened;
 }
}

#' Print RethinkDB connection
#'
#' Prints a RethinkDB connection details, including a number of pending callbacks.
#'
#' @param x Connection to print.
#' @param ... Ignored.
#' @method print RethinkDB_connection
#' @note Never blocks.
#' @export
print.RethinkDB_connection<-function(x,...){
 opened<-ifelse(isOpened(x),"Opened","Lost");
 cat(sprintf("\n %s connection to RethinkDB @ %s:%s.\n",
  opened,x$host,x$port));
 if(length(x$handlers)>0){
  cat(sprintf(" %s handlers await draining.\n",length(x$handlers)))
 }
 cat("\n");
}

#' Close RethinkDB connection
#'
#' Closes connection and stops all associated callbacks and/or sync cursor.
#'
#' @param con Connection to close.
#' @param ... Ignored.
#' @method close RethinkDB_connection
#' @export
close.RethinkDB_connection<-function(con,...){
 #Closing socket
 try(close(con$con),silent=TRUE);
 try(close(con$currentCursor),silent=TRUE);
 #Removing handlers
 con$handlers<-list();
 return(con);
}

START<-1;
CONTINUE<-2;
STOP<-3;
sendQuery<-function(x,type,query,token,options){
 stopifnot(inherits(x,"RethinkDB_connection"));

 if(missing(token))
  token<-as.integer(floor(stats::runif(1)*2^31));

 if(nchar(query)>0){
  if(missing(options)||(length(options)==0)){
   query<-sprintf("[%s,%s]",type,query);
  }else{
   query<-sprintf("[%s,%s,%s]",type,query,rjson::toJSON(options));
  }
 }else{
  query<-sprintf("[%s]",type);
  if(!missing(options)) stop("Can't have options without query.");
 }

 len<-as.integer(nchar(query)+1); #Also \0 at the end!
 writeBin(token,x$con,size=8,endian='little');
 writeBin(len,x$con,size=4,endian='little');
 writeBin(query,x$con,endian='little');

 return(token);
}

fetchResponseRaw<-function(x,syncToken){
 stopifnot(inherits(x,"RethinkDB_connection"));

 token<-readBin(x$con,integer(),size=8);
 if(length(token)==0){
  Sys.sleep(1);
  return(x);
 }
 len<-readBin(x$con,integer(),size=4);
 U<-readBin(x$con,raw(),size=1,n=len);
 stuff<-rjson::fromJSON(readBin(U,character()));

 SUCCESS_PARTIAL<-3;

 if(length(stuff$t)!=1) stop("Invalid response!");

 #Parse response type
 if(stuff$t>15){
  close(x);
  stop(sprintf("Error: %s",stuff$r));
 }

 if(!is.null(stuff$p)){
  message("Saved profile in the connection's lastProfile element.");
  x$lastProfile<-stuff$p;
 }

 token<-token;
 if(!is.null(x$handlers[[sprintf("t%s",token)]])){
  hadCallback<-TRUE;
  isMore<-stuff$t==SUCCESS_PARTIAL;
  cba<-TRUE;
  for(element in stuff$r){
   x$handlers[[sprintf("t%s",token)]](element)->cba;
   if(!identical(cba,TRUE)) break;
  }
  if(isMore){
   if(identical(cba,TRUE)){
    #Request more
    sendQuery(x,CONTINUE,'',token);
   }else{
    #Drop this query
    x$handlers[[sprintf("t%s",token)]]<-NULL;
    sendQuery(x,STOP,'',token);
   }
  }else{
   x$handlers[[sprintf("t%s",token)]]<-NULL;
  }
 }else{
  #No handler...
  if((!missing(syncToken))&&(token==syncToken)){
   #...maybe this is a syncToken?
   return(list(r=stuff$r,cc=stuff$t==SUCCESS_PARTIAL));
  }else{
   #...maybe this is a currentCursor syncToken?
   if((!is.null(x$currentCursor))&&(x$currentCursor$token==token)){
    #Add response to its cache
    x$currentCursor$contents<-c(x$currentCursor$contents,stuff$r);
    x$currentCursor$active<-(stuff$t==SUCCESS_PARTIAL);
   }
  }
  #if not, probably a dangling query stop ack; ignore.
 }
 return(invisible(NULL));
}

hookQuery<-function(x,query,cb,options){
 token<-sendQuery(x,START,query,options=options);
 x$handlers[[sprintf("t%s",token)]]<-cb;
 return(x);
}

syncQuery<-function(x,query,options){
 token<-sendQuery(x,START,query,options=options);
 while(TRUE){
  fetchResponseRaw(x,token)->ans;
  if(is.null(ans)) next; #That was an async response
  break;
 }
 if((length(ans$r)==1)&&(!ans$cc)){
  return(ans$r[[1]]);
 }else{
  if((!is.null(x$currentCursor))&&(x$currentCursor$active)){
   warning("Stopping existing, active sync cursor. Use async API to work with multiple streams!");
   close(x$currentCursor);
  }
  x$currentCursor<-makeCursor(x,ans$r,token,ans$cc);
  return(x$currentCursor);
 }
 return(ans);
}

#' Drain RethinkDB connection
#'
#' Drains a given RethinkDB connection, i.e. pull query responses and both call their associated callbacks (for async queries) and/or filling sync cursor local cache.
#' Draining ends when all async queries end; the function blocks for the entire time this is happening.
#'
#' The async query callback will only fire during \code{drainConnection} or (opportunistically) \code{\link{cursorNext}}; consequently this function must be run to guarantee that installed callbacks will have a chance to fire.
#'
#' @param x Connection to drain.
#' @export
drainConnection<-function(x){
 while(length(x$handlers)>0){
  if(isOpened(x)){
   fetchResponseRaw(x);
  }else{
   close(x); #Also clears x$handlers
  }
 }
 return(invisible(NULL));
}
