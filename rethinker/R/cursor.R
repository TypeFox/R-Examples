makeCursor<-function(con,contents,token,active){
 new.env()->ans;
 class(ans)<-"RethinkDB_cursor";
 ans$con<-con;
 ans$active<-active; # i.e. makes sense to pull more after contents are used
 ans$contents<-contents;
 ans$token<-token;
 ans
}

#' Pull next object from a cursor
#'
#' Pulls a datum from a given cursor, sending continuation queries when needed.
#'
#' @param cursor Cursor to pull from; a result of \code{r()$...$run(...)}.
#' @param inBatch If set to \code{TRUE}, enables batch mode, i.e., returning the whole local cache (this is usually NOT the whole data available under cursor) rather than a single result.
#' Values other than \code{TRUE} or \code{FALSE} are invalid.
#' @return In a default mode, a list representing the returned response JSON, or \code{NULL} if no data is available.
#' In a batch mode, list of such lists representing the whole cache (which may be empty, corresponding to default mode's \code{NULL}).
#' @note When this function empties local cache, it may ask RethinkDB for more data and hence block.
#' Use \code{\link{isCursorEmpty}} to decide if it makes sense to call \code{cursorNext}.
#' In case you don't need any more answers for the query, close cursor with \code{close} method.
#' @author Miron B. Kursa
#' @export
cursorNext<-function(cursor,inBatch=FALSE){
 stopifnot(inherits(cursor,"RethinkDB_cursor"));
 cursor$active<-(cursor$active)&&(!is.null(cursor$con))&&isOpened(cursor$con);
 if((length(cursor$contents)==0)&&(cursor$active)){
  #Empty cache? Query for more!
  sendQuery(cursor$con,CONTINUE,'',cursor$token);
  while(TRUE){
   fetchResponseRaw(cursor$con,cursor$token)->ans;
   if(is.null(ans)) next;
   break;
  }
  cursor$contents<-ans$r;
  cursor$active<-ans$cc;
 }
 #Use cache
 if(length(cursor$contents)>0){
  if(inBatch){
   ans<-cursor$contents;
   cursor$contents<-list();
  }else{
   ans<-cursor$contents[[1]];
   cursor$contents[-1]->cursor$contents;
  }
  return(ans);
 }

 #Nothing to return
 if(inBatch) return(list());
 return(NULL);
}

#' Convert cursor into a list
#'
#' Converts cursor into a list.
#' For convenience, when given anything other than cursor returns this object unchanged; this way can be used to wrap the result of \code{$run}, so that it is never a cursor.
#'
#' @param x RethinkDB cursor or any object.
#' @param maxResults Number of results sufficient to stop pulling from cursor.
#' @return A list of elements pulled from \code{x} if it is a cursor, \code{x} otherwise.
#' @note The lenght of a list may be larger than \code{maxResults} because RethinkDB transmits results in batches.
#' @author Miron B. Kursa
#' @export
cursorToList<-function(x,maxResults=10000){
 if(!inherits(x,"RethinkDB_cursor")) return(x);
 ans<-list();
 while((!isCursorEmpty(x))&&(length(ans)<maxResults)){
  ans<-c(ans,cursorNext(x,inBatch=TRUE));
 }
 return(ans);
}

#' Check if cursor is empty
#'
#' Check whether a given cursor is fully drained and will output no more datum.
#' The function never blocks; also verifies that the underlying connection is useful.
#'
#' @param cursor Cursor to check; a result of \code{r()$...$run(...)}.
#' @return \code{TRUE} if cursor has no more data to return.
#' @note It is possible that \code{\link{cursorNext}} will return \code{NULL} just after \code{isCursorEmpty} returns \code{FALSE}.
#' Changefeeds cursors (made with \code{r()$...$changes()$...}) will never become empty (provided that connection won't become broken).
#' @author Miron B. Kursa
#' @export
isCursorEmpty<-function(cursor){
 stopifnot(inherits(cursor,"RethinkDB_cursor"));
 cursor$active<-(cursor$active)&&(!is.null(cursor$con))&&isOpened(cursor$con);
 return((length(cursor$contents)==0)&&(!cursor$active));
}

#' Print cursor
#'
#' Prints a given cursor's status.
#'
#' @param x Cursor to print.
#' @param ... Ignored.
#' @method print RethinkDB_cursor
#' @note Never blocks; also checks whether the underlying connection is alive.
#' @export
print.RethinkDB_cursor<-function(x,...){
 empty<-isCursorEmpty(x);
 if(empty){
  cat("\n Empty RethinkDB cursor.\n")
 }else{
  cat(sprintf("\n Active RethinkDB cursor;\n %s response(s) cached",
   length(x$contents)));
  if(x$active){
   cat(", more to come.\n");
  }else{
   cat(", no more on the server.\n");
  }
 }
 cat("\n");
 return(invisible(x));
}

#' Close cursor
#'
#' Closes a given cursor and stops its associated query.
#' Should be called on the current cursor before a new sync query is invoked on the same connection.
#'
#' @param con Cursor to close.
#' @param ... Ignored.
#' @method close RethinkDB_cursor
#' @export
close.RethinkDB_cursor<-function(con,...){
 x<-con;
 x$active<-(x$active)&&(!is.null(x$con))&&isOpened(x$con);
 if(x$active)
  sendQuery(x$con,STOP,'',x$token);
 x$con<-NULL; #GC may use it to kill connection
 x$active<-FALSE;
 x$contents<-list();
 return(x);
}
