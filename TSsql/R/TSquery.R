
TSquery <- function (select, dateField, table, where = NULL, 
       frequency = "monthly", na.as=0, names=NULL, con = options()$connection){
    if(missing(con)&& is.null(con)) 
	  stop("con should be specified or set with options(connection=con). See ?TSquery.") 
    if (is.null(con)) stop("argument 'con' cannot be NULL")
    if (frequency == "monthly") 
        dates <- paste("EXTRACT(YEAR from ",  dateField,
	            "), EXTRACT(MONTH from ", dateField, ")")
    else if (frequency == "annual") 
        dates <- paste("EXTRACT(YEAR from ", dateField, ")")
    # something like this if dates have time? (or can zoo ignore?
    #else if (frequency == "daily") 
    #    dates <- paste("EXTRACT(YEAR  from ", dateField,
    #            "), EXTRACT(MONTH from ", dateField, 
    #	         "), EXTRACT(DAY   from ", dateField, ")")
    else if (frequency == "daily") 
        dates <- dateField
    else stop("frequency not supported.")
    q <- paste("SELECT ", dates, ", ", select, " FROM ", table)
    if (!is.null(where)) 
        q <- paste(q, " WHERE ", where)
    q <- paste(q, " GROUP BY ", dates, " ORDER BY ", dates, " ;")
    reso <- res <- dbGetQuery(con, q)
    if (any(dim(res) == 0)) 
        stop("empty query result.")
    res <- as.matrix(res)
    if (frequency == "monthly") {
        res <- res[!is.na(res[, 1]), , drop = FALSE]
        res <- res[!is.na(res[, 2]), , drop = FALSE]
        y <- res[1, 1]
        m <- res[1, 2]
        sampleT <- 1 + (12 * res[nrow(res), 1] + res[nrow(res), 
            2]) - (12 * res[1, 1] + res[1, 2])
        r <- matrix(numeric(3 * sampleT), sampleT, 3)
        j <- 1
        for (i in seq(sampleT)) {
            if ((res[j, 1] == y) & (res[j, 2] == m)) {
                r[i, ] <- res[j, ]
                j <- j + 1
            }
            m <- m + 1
            if (m == 13) {
                y <- y + 1
                m <- 1
            }
        }
    }
    else if (frequency == "annual") {
        res <- res[!is.na(res[, 1]), , drop = FALSE]
        y <- res[1, 1]
        sampleT <- 1 + res[nrow(res), 1] - y
        r <- matrix(numeric(2 * sampleT), sampleT, 2)
        j <- 1
        for (i in seq(sampleT)) {
            if (res[j, 1] == y) {
                r[i, ] <- res[j, ]
                j <- j + 1
            }
            y <- y + 1
        }
    }
    if (frequency == "monthly") 
        r <- ts(r[, -(1:2)], start = r[1, 1:2], frequency = 12)
    else if (frequency == "annual") 
        r <- ts(r[, -1], start = r[1, 1], frequency = 1)
    else if (frequency == "daily") {
        #expand with NA (or 0) in missing days
	st <- as.Date(res[1,1])
	en <- as.Date(res[NROW(res),1])
 	dex <- st + 0 : (en - st)
	ind <- dex %in% as.Date(res[,1]) 
	r   <- rep(na.as, length(ind))
	r[ind] <- as.numeric(res[,2])
	r <- zoo::zoo(r, dex)
	}
    if(!is.null(names)) seriesNames(r) <- names 
    r
    }

