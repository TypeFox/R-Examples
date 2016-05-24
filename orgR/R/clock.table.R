
map.clock.heading <- function(clock.ind, heading.vec){#
    d <- clock.ind - heading.vec#
    neg.ind <- which(d <= 0)[1]#
    return(neg.ind - 1 )#
}

##' Parse clock entry to ISO date 
##'
##' 
##' @title clock.table
##' @param clock.entries a standard clock entry from org-mode
##' @return POXICt object
##' @author Yi Tang
##' @export
##' @examples 
##' str <- c("CLOCK: [2014-11-26 Wed 09:36]--[2014-11-26 Wed 10:04] =>  0:28",
##'          "CLOCK: [2014-12-04 Thu 15:24]--[2014-12-04 Thu 16:25] =>  1:01")
##' ToISOdate(str)
ToISOdate <- function(clock.entries){
    s <- str_locate_all(clock.entries, "\\[")
    e <- str_locate_all(clock.entries, "\\]")
    tmp1 <- rep(NA, len = length(clock.entries))
    tmp2 <- rep(NA, len = length(clock.entries))
    for (i in seq_along(clock.entries)){
        s.i <- s[[i]]
        e.i <- e[[i]]
        tmp1[i] <- substr(clock.entries[i], s.i[1, 1] + 1 , e.i[1, 1] - 1)
        tmp2[i] <- substr(clock.entries[i], s.i[2, 1] + 1, e.i[2, 1] - 1)
    }
    res <- list(lubridate::ymd_hm(tmp1),
                lubridate::ymd_hm(tmp2))
    return(res)
}

##' Parse org file
##'
##' scan a org file and return the headlines and associated clock entries 
##' @title clock.table 
##' @export
##' @param org.file a org file 
##' @return a data.table 
##' @author Yi Tang
GetClockTable <- function(org.file = "~/tmp.org"){
    ## org.file <- "~/git/org/tmp.org"
    str <- readLines(org.file) 
    heading.lines <- grep("^\\*", str)
    clock.entry.lines <- grep("CLOCK: \\[", str)



    clock.entries <- str[clock.entry.lines]
    headings <- str_trim(str[heading.lines])

    ind <- sapply(clock.entry.lines, function(i) map.clock.heading(i, heading.lines))
    clock.table <- data.table(clock.entries,
                              headings = str[heading.lines[ind]],
                              head.ind = heading.lines[ind])
    clock.table$clock.closed <- grepl("--", clock.table$clock.entries)

ind <- clock.table$clock.closed == TRUE
    clock.table[ind, c("start", "end") := {
        ToISOdate(clock.entries)
    }]

    clock.table[, clock.entries := NULL]
    return(clock.table)

}
