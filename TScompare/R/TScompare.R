testequaltf <- function(s1, s2) {
    #r <- try(all( c(earliestEnd(s1) == earliestEnd(s2), 
    #                latestStart(s1) == latestStart(s2))), silent=TRUE)
    r <- try(all(c(tfend(s1) == tfend(s2), tfstart(s1) == tfstart(s2))), 
             silent=TRUE)
    if (inherits(r, "try-error")) r <- FALSE #includes failed as well as unequal
    # failure may be because of different time representations, 
    r
    }

TScompare <- function(ids, con1, con2, na.rm=FALSE, fuzz=1e-14) {
	#would be good to check if the cons are still both alive (not timed out)
	if(1 == NCOL(ids)) ids <- cbind(ids, ids)
	if(2 != NCOL(ids)) stop("ids must not have more than 2 columns.")
	rw <- rv <- rep(NA, NROW(ids))
	na1 <- na2 <- NULL
	for (i in 1:NROW(ids)){
	   s1 <- try(TSget(ids[i,1], con1), silent=TRUE)
	   if (inherits(s1, "try-error")) na1 <- c(na1, ids[i,1])
	   s2 <- try(TSget(ids[i,2], con2), silent=TRUE)
	   if (inherits(s2, "try-error")) na2 <- c(na2, ids[i,2])
	   if ((!inherits(s1, "try-error")) & (!inherits(s2, "try-error"))) {
	      if(na.rm) {
		   s1 <- trimNA(s1)
		   s2 <- trimNA(s2)
		   }
	      rw[i] <- testequaltf(s1, s2)
	   
	      ii <- ! is.na(s1)
	      #if window is FALSE values are FALSE automatically
	      if(rw[i]) rv[i] <- all( c(ii == !is.na(s2)) && 
		                 !any(abs(c(s1)[ii] - c(s2)[ii]) > fuzz))
	      else rv[i] <- FALSE
	      }
	   }
	r <- list(window=rw, value=rv, ids=ids, na1=na1, na2=na2,
	         con1=con1, con2=con2, na.rm=na.rm, fuzz=fuzz)
	class(r) <- "TScompare"
	r
	}

summary.TScompare  <- function(object, ...){
	x <- list(n=length(object$window),
		  na1=length(object$na1),
		  na2=length(object$na2),
		  na =sum(is.na(object$window)),
		  window=sum(object$window, na.rm=TRUE),
		  value=sum(object$value, na.rm=TRUE))
	class(x) <- "summary.TScompare"
	x
	}

print.summary.TScompare  <- function(x, digits=getOption("digits"), ...){
	cat(x$n - x$na1, " of ", x$n, "are available on con1.\n")
	cat(x$n - x$na2, " of ", x$n, "are available on con2.\n")
	cat(x$window, " of ", x$n - x$na, "remaining have the same window.\n")
	cat(x$value,  " of ", x$n - x$na, "remaining have the same window and values.\n")
	invisible(x)
	}

doubleCheck <- function(x, con1=x$con1, con2=x$con2, na.rm=FALSE, fuzz=1e-10) {
	# if the con has expired default will fail
	# go through a TScompare object with a relaxed fuzz and return a new
	# TScompare object with modified $value
	# other elements should not be changed
	if(x$fuzz > fuzz) 
	   stop("doubleCheck only works if fuzz is bigger than original fuzz.")
	ids <- x$ids
	rw  <- x$window
	rv  <- x$value

	for (i in 1:NROW(ids)){
	   if((!is.na(rv[i])) && (!rv[i]) && rw[i] ) {
	   s1 <- try(TSget(ids[i,1], con1), silent=TRUE)
	   s2 <- try(TSget(ids[i,2], con2), silent=TRUE)
	   if (inherits(s1, "try-error"))
	      warning( ids[i,1], "is now NA on con1 (but not change result).")
	   else if (inherits(s2, "try-error")) 
	      warning( ids[i,2], "is now NA on con2 (but not change result).")
	   else {
	      if(na.rm) {
		   s1 <- trimNA(s1)
		   s2 <- trimNA(s2)
		   }
	   
	      ii <- ! is.na(s1)
	      # only check 
	      if(rw[i] && !rv[i]) rv[i] <- all( c(ii == !is.na(s2)) && 
		                 !any(abs(c(s1)[ii] - c(s2)[ii]) > fuzz))
	      }
	   }}
	r <- list(window=rw, value=rv, ids=ids, na1=x$na1, na2=x$na2,
	       con1=con1, con2=con2, na.rm=na.rm, fuzz=fuzz)
	class(r) <- "TScompare"
	r
	}

tfDetails <- function(x, con1=x$con1, con2=x$con2, na.rm=FALSE) {
	# if the con has expired default will fail
	# go through a TScompare object and give details of series with 
	# different $window
	ids <- x$ids
	rw  <- x$window
	r   <- matrix("ok", NROW(ids), 2)
	for (i in 1:NROW(ids)){
	   if(is.na(rw[i]))  r[i,] <- "NA"
	   if(!rw[i]) {
	      s1 <- try(TSget(ids[i,1], con1), silent=TRUE)
	      s2 <- try(TSget(ids[i,2], con2), silent=TRUE)
	      if (inherits(s1, "try-error"))
	   	 warning( ids[i,1], "is now NA on con1 (not comparing windows).")
	      else if (inherits(s2, "try-error")) 
	   	 warning( ids[i,2], "is now NA on con2 (not comparing windows).")
	      else {
	   	 if(na.rm) {
	   	      s1 <- trimNA(s1)
	   	      s2 <- trimNA(s2)
	   	      }
	      
	   	 ii <- ! is.na(s1)
	   	 r[i,1] <- paste(start(s1), " to ", end(s1),collapse="")
	   	 r[i,2] <- paste(start(s2), " to ", end(s2),collapse="")
	   	 }
	      }}
	dimnames(r) <- list(ids[,1], c("con1", "con2"))
	r[!rw,]
	}

tfplot.TScompare  <- function(x, con1=x$con1, con2=x$con2, diff=FALSE, ...){
	# if the con has expired default will fail
	v <- x$value
	v[is.na(v)] <- FALSE
	ids <- x$ids[!(v & x$window),]
	if(0 == NROW(ids)) message("No differences to plot.")
	else for (i in 1:NROW(ids)){
	   if(diff) tfplot(TSget(ids[i,1], con1) - TSget(ids[i,2], con2),
	                   Title=ids[i,1])
	   else     tfplot(TSget(ids[i,1], con1),  TSget(ids[i,2], con2),
	                   Title=ids[i,1])
	   }
	invisible(x)
	}

AllIds <- function(con, vintage=getOption("TSvintage")){
  if (is.null(vintage)) Q <- "select distinct id from Meta;" 
  else {
    if( vintage == "current") vintage <- DBI::dbGetQuery(con,
         "select vintage from vintageAlias where alias = 'current';")$vintage
    Q <- paste(
       "select distinct id from Meta where vintage = '", vintage, "';", sep="")
    }
  DBI::dbGetQuery(con,Q)$id 
  }

AllPanels <- function(con){	
  if(!con@hasPanels) NULL
  else {
  if (is.null(vintage)) Q <- "select distinct panel from Meta;" 
  else {
    if( vintage == "current") vintage <- DBI::dbGetQuery(con,
         "select vintage from vintageAlias where alias = 'current';")$vintage
    Q <- paste(
       "select distinct panel from Meta where vintage = '", vintage, "';", sep="")
    }
  DBI::dbGetQuery(con, Q)$panel
  }
  }


AllVintages <- function(con){
  if(!con@hasVintages) NULL
  else DBI::dbGetQuery(con, "select distinct vintage from Meta;")$vintage
  }
