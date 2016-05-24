# $Id: rename.vars.R 625 2005-06-09 14:20:30Z nj7w $

rename.vars <- function(data,from='',to='',info=TRUE) {

   dsn <- deparse(substitute(data))
   dfn <- names(data)

   if ( length(from) != length(to)) {
     cat('--------- from and to not same length ---------\n')
     stop()
   }

   if (length(dfn) < length(to)) {
     cat('--------- too many new names ---------\n')
     stop()
   }

   chng <- match(from,dfn)

   frm.in <- from %in% dfn
   if (!all(frm.in) ) {
     cat('---------- some of the from names not found in',dsn,'\n')
     stop()
   }

   if (length(to) != length(unique(to))) {
     cat('---------- New names not unique\n')
     stop()
   }

   dfn.new <- dfn
   dfn.new[chng] <- to
   if (info) cat('\nChanging in',dsn)
   tmp <- rbind(from,to)
   dimnames(tmp)[[1]] <- c('From:','To:')
   dimnames(tmp)[[2]] <- rep('',length(from))
   if (info)
     {
       print(tmp,quote=FALSE)
       cat("\n")
     }
   names(data) <- dfn.new
   data
}


# GRW 2004-04-01
remove.vars <- function( data, names, info=TRUE)
  {
    for( i in names )
      {
        if(info)
          cat("Removing variable '", i, "'\n", sep="")
        data[[i]] <- NULL
      }
    data
  }
