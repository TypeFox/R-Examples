
##' ad descrition 
##'
##' lalla details
##' @rdname time.stampe
##' @title function 
##' @param str content of a .org file
##' @param ts.format format of time stamps used in the .org file. It is equivalent to \code{org-time-stamp-formats} in Emacs
##' @return a 
##' @export
##' @author yitang
GetTS <- function(str = "a", ts.format = c("<%Y-%m-%d %a>", "<%Y-%m-%d %a %H:%M>")){
    inactive.ts.format <- "\\[[0-9]{4}-[0-9]{2}-[0-9]{2} [[:alpha:]]{3} [0-9]{2}:[0-9]{2}\\]"
    
    ts.category <- c(close = "CLOSED: ", dealine = "DEADLINE: ", scheduel = "SCHEDULED: ")
    ts.entries <- lapply(ts.category, function(i){
        i.str <- str_extract(str, paste0(i , inactive.ts.format))
        i.ts <- gsub(i, "", i.str) 
        return(i.ts)
    })
    setDT(ts.entries)
    return(ts.entries)
}


##' org headlines 
##'
##' A function to parse org files, will return headlines and associated attributes, including tag, clock entries, shedules, deadlines, closed date, todo states,
##' @title Headlines 
##' @param org.file a file path point to a .org file 
##' @return a table of headlines and attributes 
##' @export
##' @author Yi Tang
GetHeadlines <- function(org.file = "~/tmp.org"){
    str <- readLines(org.file) 
    heading.lines <- grep("^\\*", str)#
    headings <- str_trim(str[heading.lines])
    
    levels <- str_count(headings, "\\*")
    ## remove stars from the heading
    headings <- gsub("^\\*{1, } ", "", headings)



    todo.keywords <- c("TODO", "NEXT", "DONE", "WAITING", "HOLD", "CANCELLED", "PHONE", "MEETING")
    ## if (is.null(todo.keywords)){
    ##     str <- strsplit(content, " ")
    ##     str <- sapply(str, "[", 2)
    ## }
    ## must provide todo lists

#### todo states 
    first.word.in.headings <- str_extract(headings, "[[:alpha:]]{1, }")
    todo.ind <- first.word.in.headings %in% todo.keywords
    todo.state <- NA
    todo.state[todo.ind] <- first.word.in.headings[todo.ind]

#### tags
    all.tags <- str_extract(headings, ":[[:alpha:]]{+}:")
    archive.tag <- ":ARCHIVE:" == all.tags 
    archive.tag[is.na(archive.tag)] <- FALSE


    the.line.after.heading <- str[heading.lines + 1]
    plan.ts <- GetTS(str = the.line.after.heading)

    head.info <- data.table(id = 1:length(headings),
                            heading = headings,
                            head.ind = heading.lines,
                            level = levels,
                            todo = todo.state,
                            tag = all.tags,
                            archive = archive.tag,
                            plan.ts)
    return(head.info)
}


##' Search for parent headlines 
##'
##' Given a headlines table and headline id, it will return the parent headlines.
##' @title search.parent 
##' @param head.info a head table from GetHeadlings
##' @param heading.id a unique id from head.info 
##' @return a data.table 
##' @export
##' @author Yi Tang
search.parent <- function(head.info, heading.id){
    level.1.ind <- which(head.info$level == 1)
    dist <- (level.1.ind - heading.id)^2
    ind <- which.min(dist)
    nearest.level.1 <- level.1.ind[ind] 
    head.info[nearest.level.1 : heading.id]
}


## search.parent(head.info, 40)


#### tree view of strcutre
##' Visualise org-mode headings 
##'
##' tree structure of org headlines
##' @title org-headings-tree
##' @param head.info a data.tabl returned by GetHeadlines()
##' @param output file to save the results, default setting is to print to scree 
##' @param plantuml TRUE/FALSE, for plantuml program?
##' @return a string that can be used in plantuml program
##' @export
##' @author Yi Tang
tree.headings <- function(head.info, output = "screen", plantuml = TRUE){
    tree.prep <- sapply(head.info$level, function(i) paste(rep("+", i), collapse = ""))
    tree <- paste(tree.prep, head.info$heading)
    if (plantuml) 
        tree <- c("@startuml", "salt", "{", "{T", tree, "}", "}", "@enduml")
    if (output == "screen")
        cat(tree, sep = "\n")
    else
        cat(tree, sep = "\n", file = output)
}

## tree.headings(head.info, output = "screen")
