# To confuse the code checker
Mem.Usage <- 
CPU.Time <- 
Session.Name <-
NULL


win_processes <- 
function(){
    "list machine processes"
    stopifnot(requireNamespace("lubridate"))
    output <- system2("tasklist", c("/v", "/fo csv"), stdout=T)
    file <- textConnection(output)
    on.exit(close(file))
    wp <- utils::read.csv( file, na.strings="N/A")
    names(wp) <- gsub("\\.$", "", names(wp))
    Session.Name
    transform( wp
             , Mem.Usage = as.numeric(gsub(" K", "", gsub(",", "", wp$Mem.Usage)))
             , CPU.Time  = lubridate::dhours(  as.integer(gsub("(\\d+):(\\d+):(\\d+)", "\\1", CPU.Time)))
                         + lubridate::dminutes(as.integer(gsub("(\\d+):(\\d+):(\\d+)", "\\2", CPU.Time)))
                         + lubridate::dseconds(as.integer(gsub("(\\d+):(\\d+):(\\d+)", "\\3", CPU.Time)))
             , Session.Name = tolower(Session.Name)
             )
}
win_users <-
function(){
    "List logged on users"
    stopifnot(requireNamespace("stringr"))
    output <- system("query user", intern=TRUE)
    output <- stringr::str_sub(output, 2)
    
    x <- stringr::str_split(output, "  +")
    
    y <- do.call(rbind, x[-1])
    colnames(y) <- x[[1]]
    z <- as.data.frame(y)
    z[["LOGON TIME"]] <- as.POSIXct(z[["LOGON TIME"]], format="%m/%d/%Y %I:%M %p")
    return(z)
}
if(FALSE){ #testing
    library(plyr)
    processes <- win_processes()
    head(processes)
    
    arrange(processes, ImageName)
    "NppToR.exe" %in% win_processes()$ImageName
    
    using(plyr)
    ddply( win_processes(), .(Session), summarize
         , Total.Mem = sum(Mem.Usage)
         )
}

win_load <- function(){
    "list the load on the CPU."
    system2("wmic", c("cpu", "get", "loadpercentage"))
}

win_kill <- 
function( ...       	#< Thrown Away, used to force user to specify full argument name.
        , image     	#< The executable name, such as 'Rgui.exe'
	, pid       	#< process ID, a number
	, force=TRUE	#< force close?
	, title=NULL	#< Filter by title
	){
    "Kill processes"
    if(!missing(image))
        for(i in image)
            system2("taskkill", c(if(force)"/F", "/IM", i))
    if(!missing(pid))
        for(p in pid)
            system2("taskkill", c(if(force)"/F", "/IM", p))    
    if(!missing(title))
        system2("taskkill", c(if(force)"/F", "/FI", shQuote(paste0("Windowtitle eq ", title))))
}

# @Suggests plyr
# @Suggests lubridate
whos_the_hog <- function(){
    "Lists the resource useage summarized by user."
    stopifnot(requireNamespace("plyr"), requireNamespace("lubridate"))
    wp <- win_processes()
    users <- win_users()
    names(users) <- c("Username", "Session.Name", "Session.ID", "State", "Idle.Time", "Logon.Time")
    joined <- merge( wp[setdiff(names(wp), "User.Name")]
                   , users, by="Session.Name")
    
    plyr::arrange(
	plyr::ddply( joined, plyr::as.quoted("Username"), plyr::summarize
                   , Mem.Usage = sum(Mem.Usage)
                   , N.Processes = NROW(Mem.Usage)
                   , CPU.Time = sum(CPU.Time)
                   ), Mem.Usage, CPU.Time)
}
