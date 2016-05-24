
#' LogBook DataBase Creation GUI
#' 
#' The \code{gui_lb_db_create} function implements the graphical user interface for the
#' creation of a logbook database.
#' This function, with an edited logbook dataset (see \code{\link{gui_lb_editraw}}),
#'  creates a new logbook database.
#'
#' @return This function does not return a value. 
#' After the execution a logbook database will be deployed.
#' 
#' @usage gui_lb_db_create()
#' 
#' @export gui_lb_db_create
#'
#'@seealso \code{\link{gui_lb_editraw}}

### LogBook DataBase manager gui


gui_lb_db_create <- function()
{
  lb_DB <- log_DB$new()
  
  lb_db_win <- gwindow("Create New logbook DataBase", visible = FALSE)
  one <- ggroup(horizontal = FALSE, container = lb_db_win)
  
  c_name <- gframe(text = "Enter a name for the new DataBase and press Enter",
                   horizontal = TRUE, container = one)
  addSpring(c_name)
  dbname <- gedit(initial.msg = "New logbook DB name...",
                  container = c_name, handler = function(h,...)
                  {
                    enabled(dats) <- TRUE
                  }) 
  addSpring(c_name)
  gimage(system.file("ico/go-down-3.png", package="vmsbase"), container = one)
  dats <- gbutton(text = "Select logbook data", container = one, handler = function(h,...)
  {
    svalue(dats) <- gfile(text = "Select logbook data",
                          type = "open",
                          filter = list("LB data" = list(patterns = c("*.logbook")),
                                        "All files" = list(patterns = c("*"))))
    enabled(go) <- TRUE
  })
  
  gimage(system.file("ico/go-down-3.png", package="vmsbase"), container = one)
  
  go <- gbutton("Create DB", container = one, handler = function(h,...)
  {
    enabled(lb_db_win) <- FALSE
    
    lb_DB$dir <- gfile(text = "Select Logbook DB destination",
                       type = "selectdir",
                       filter = list("LB data" = list(patterns = c("*.lb.sqlite"))))
    
    lb_DB$db <- paste(lb_DB$dir, "/", svalue(dbname), ".lb.sqlite", sep = "")  
    
    dbConnect(SQLite(), dbname = lb_DB$db)
    
    firstline <- readLines(svalue(dats), n = 1)
    if(nchar(firstline) > 45)
    {
      read.csv.sql(svalue(dats), sql = "CREATE TABLE elobo AS SELECT * FROM file", dbname = lb_DB$db, eol = "\n")
      sqldf("CREATE TABLE lb_cla(vessel INT, log_num INT, met_fo INT, met_des CHAR)", dbname = lb_DB$db)
      
    }else{
    read.csv.sql(svalue(dats), sql = "CREATE TABLE logbook AS SELECT * FROM file", dbname = lb_DB$db, eol = "\n")
    
    sqldf("CREATE TABLE elobo(vessUE INT, s_utc REAL, e_utc REAL, FAO_SPECIES INT)", dbname = lb_DB$db)
    sqldf("CREATE TABLE lb_cla(vessel INT, log_num INT, met_fo INT, met_des CHAR)", dbname = lb_DB$db)
    
    ##############
    
    distin <- sqldf("select distinct vessUE, s_utc, e_utc, gear, metier from logbook order by s_utc, e_utc", dbname = lb_DB$db)
    
    cat("\n\n   ---   Logbook Editing Started   ---\n", sep = "")
    
    numrow <- nrow(distin)
    species <- sqldf("select distinct specie from logbook order by specie", dbname = lb_DB$db)
    tocn <- which(species[,1] == "")
    ifelse(length(tocn) > 0, species <- species[-tocn,1], species <- species[,1])
    fao_lst <- read.table(system.file("extdata/FAO_CODES_SPECIES.csv", package="vmsbase"), sep = ",", dec = ".")
    tocl <- which(!(species %in% fao_lst[,1]))
    if(length(tocl) > 0){species <- species[-tocl]}
    spclst <- matrix(data = 0, nrow = 1, ncol = length(species))
    colnames(spclst) <- paste("FAO_", species, sep = "")
    proto_log <- data.frame("vessUE" = numeric(1),
                            "s_utc" = numeric(1),
                            "e_utc" = numeric(1),
                            "gear" = character(1),
                            "metier" = character(1))
    proto_log <- cbind(proto_log, spclst)
    hea_str <- paste("vessUE INT, ", paste(colnames(proto_log[2:3]), "REAL", sep = " ", collapse = ", "), ", ", paste(colnames(proto_log[4:5]), "CHAR", sep = " ", collapse = ", "), sep = "")
    cre_str <- paste(colnames(proto_log[6:ncol(proto_log)]), "INT", sep = " ", collapse = ", ")
    fin_str <- paste(hea_str, cre_str, sep = ", ")
    
    sqldf("drop table if exists elobo", dbname = lb_DB$db)
    query <- paste("CREATE TABLE elobo(", fin_str,")", sep = "")
    sqldf(query, dbname = lb_DB$db)
    
    for(n in 1:numrow)
    {
      cat("\nLogbook ", n, " of ", numrow, sep = "")
      if(distin[n,"s_utc"] == "NA" | distin[n,"e_utc"] == "NA" | distin[n,"s_utc"] > distin[n,"e_utc"])
      {
        cat(" - Skipped, invalid Date")
        next
      }
      vess_id <- distin[n,"vessUE"]
      utc_rs  <- distin[n,"s_utc"]
      utc_re  <- distin[n,"e_utc"]
      temp <- fn$sqldf("select specie, qty from logbook where vessUE = `vess_id` and s_utc >= `utc_rs` and e_utc <= `utc_re`", dbname = lb_DB$db)
      if(nrow(temp) == 0)
      {
        cat(" - Skipped", sep = "")
        next
      }
      tort <- which(temp[,1] == "")
      if(length(tort) > 0){temp <- temp[-tort,]}
      if(nrow(temp) == 0)
      {
        cat(" - Skipped!", sep = "")
        next
      }
      nclm <- (which(species %in% temp[,1]))+5
      if(max(as.numeric(table(temp[,1])))>1) temp <- as.data.frame(cbind(names(tapply(as.numeric(temp[,2]),temp[,1],sum)),as.numeric(tapply(as.numeric(temp[,2]),temp[,1],sum))))
      new_log <- data.frame("vessUE" = numeric(1), "s_utc" = numeric(1), "e_utc" = numeric(1), "gear" = character(1), "metier" = character(1))
      new_log <- cbind(new_log, spclst)
      new_log["vessUE"] <- vess_id
      new_log["e_utc"] <- utc_re
      new_log["s_utc"] <- utc_rs
      new_log["gear"] <- distin[n,"gear"]
      new_log["metier"] <- distin[n,"metier"]
      new_log[nclm] = new_log[nclm]+abs(as.numeric(temp[,2]))
      sqldf("insert into elobo select * from `new_log`", dbname = lb_DB$db)
    }
    
    cat("\n\n   ---   End Logbook Editing   ---\n\n", sep = "")
    }
    gconfirm(" LogBook DataBase \n \"Deploy & Editing\" \n Complete!",
             title = "Confirm",
             icon = "info",
             parent = lb_db_win,
             handler = function(h,...){dispose(lb_db_win)})
  })
  
  enabled(dats) <- FALSE
  enabled(go) <- FALSE
  visible(lb_db_win) <- TRUE
}