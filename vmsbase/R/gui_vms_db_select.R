
#' VMS DataBase Select GUI
#'  
#' The \code{gui_vms_db_sel} function implement the graphical user interface for the
#'  VMS DataBase Select routine.
#' 
#' This function, with a VMS DataBase (see \code{\link{gui_vms_db_stat}}),
#'  enables the user to perform queries on, and extract data from, the submitted VMS DataBase.
#'  
#' @param vms_db_name The path of a VMS DataBase
#' 
#' @return This function does not return a value. 
#' 
#' @usage gui_vms_db_sel(vms_db_name = "")
#' 
#' @export gui_vms_db_sel
#'
#'@seealso \code{\link{gui_vms_db_stat}}

gui_vms_db_sel <- function(vms_db_name = "")
{
  vms_DB <- vms_DB$new()
  que_vms_DB <- que_vms_DB$new()
  
  if(vms_db_name == "")
  {
    vms_db_name <- gfile(text = "Select VMS DataBase file",
                         type = "open",
                         filter = list("VMS DB file" = list(patterns = c("*.vms.sqlite"))))
  }
  
  vms_DB$db <- vms_db_name  
  
  min_date <- as.numeric(sqldf("select min(DATE) from ping", dbname = vms_DB$db)[,1])
  max_date <- as.numeric(sqldf("select max(DATE) from ping", dbname = vms_DB$db)[,1])
  
  
  
  min_date_l <- paste("First: ", days(min_date), "/", months(min_date), "/", years(min_date), sep = "")
  max_date_l <- paste("Last: ", days(max_date), "/", months(max_date), "/", years(max_date), sep = "")
  
  
  sel_vms_win <- gwindow("VMS Select Tool", height = 600, visible = FALSE)
  sel_big_g <- ggroup(container = sel_vms_win, horizontal = FALSE)
  sel_vms_top <- ggroup(horizontal = TRUE, container = sel_big_g)
  sel_vms_bot <- ggroup(horizontal = TRUE, container = sel_big_g, expand = TRUE)
  
  ### FROM
  sel_tab_f <- gframe(text = "Table", container = sel_vms_top)
  addSpring(sel_tab_f)
  from_cb <- gdroplist(c("Tracks", "Pings"), container = sel_tab_f,
                       handler = function(h,...)
                       {
                         enabled(sel_area_f) <- !enabled(sel_area_f)
                         enabled(sel_met_f) <- !enabled(sel_met_f)
                         enabled(sel_fish_f) <- !enabled(sel_fish_f)
                       })
  addSpring(sel_tab_f)
  
  ### DATE
  sel_date_f <- gframe(text = "Date", horizontal = FALSE, container = sel_vms_top)
  addSpring(sel_date_f)
  temp_bar <- gbutton(text = "View Date\nFrequencies", container = sel_date_f, handler = function(h,..){
    
    temp_win <- gwindow("VMS Date Frequencies Tool", width = 500, visible = FALSE)
    theplot <- ggraphics(container = temp_win, expand = TRUE)
    
    utc <- sqldf("select DATE from ping", dbname = vms_DB$db)
    utc <- as.numeric(utc[,1])    
    nbars <- 20
    min_utc <- floor(min(utc))
    max_utc <- ceiling(max(utc))
    bre_utc <- seq(min_utc,max_utc,length=nbars)
    bar_utc <- as.numeric(table(findInterval(utc,bre_utc)))
    col_utc <- plotrix::color.scale(1:length(bre_utc), extremes=c("olivedrab1","forestgreen"))
    bar_names <- chron(floor(bre_utc[-length(bre_utc)]))
    
    visible(temp_win) <- TRUE
    
    par(lwd=0.01, las=2, cex=0.7, mar=c(5,7,5,5))
    barplot(bar_utc, horiz = TRUE, axes=T, col=col_utc, names.arg=bar_names, xlab = "Num. of Pings")
    
    enabled(temp_bar) <- FALSE
    addHandlerDestroy(temp_win, handler = function(h,..){enabled(temp_bar) <- TRUE})
  })
  fro_to_g <- ggroup(horizontal = FALSE , container = sel_date_f)
  addSpring(fro_to_g)
  min_date_vms <- glabel(min_date_l, container = fro_to_g)
  sel_fro_f <- gframe(text = "From", container = fro_to_g)
  fro_d <- gdroplist(levels(days(min_date)), container = sel_fro_f)
  svalue(fro_d) <- days(min_date)
  fro_m <- gdroplist(levels(months(min_date)), container = sel_fro_f)
  svalue(fro_m) <- months(min_date)
  fro_y <- gdroplist(levels(years(min_date)), container = sel_fro_f)
  svalue(fro_y) <- years(min_date)
  
  sel_to_f <- gframe(text = "To", container = fro_to_g )
  to_d <- gdroplist(levels(days(max_date)), container = sel_to_f)
  svalue(to_d) <- days(max_date)
  to_m <- gdroplist(levels(months(max_date)), container = sel_to_f)
  svalue(to_m) <- months(max_date)
  to_y <- gdroplist(levels(years(max_date)), container = sel_to_f)
  svalue(to_y) <- years(max_date)
  max_date_vms <- glabel(max_date_l, container = fro_to_g)
  addSpring(fro_to_g)
  #use_date_g <- ggroup(container = fro_to_g, horizontal = TRUE)
  #glabel("Specify Date?", container = use_date_g)
  #use_date_r <- gradio(c("Yes", "No"), container = use_date_g, horizontal = TRUE)
  addSpring(sel_date_f)
  
  
  ### AREA
  sel_area_f <- gframe(text = "Area", horizontal = FALSE, container = sel_vms_top)
  addSpring(sel_area_f)
  
  sel_gsa_f <- gframe(text = "", horizontal = FALSE, container = sel_area_f)
  #gsa
  if(sqldf("select count(*) from p_area", dbname = vms_DB$db) == 0)
  {
    g_are_0 <- ggroup(horizontal = TRUE, container = sel_gsa_f)
    g_are_l <- glabel("No Area data\nRun Assign Area")
    enabled(sel_area_f) <- FALSE
    
  }else{
    gsas <- sqldf("select distinct(AREA) from p_area order by area", dbname = vms_DB$db)[,1]
    n_gsa <- length(gsas)
    gsas_li <- vector("list", n_gsa)
    for(i in 1:n_gsa)
    {
      if(((i-1) %% 4) == 0 )
      {
        new_gr <- paste("g_are_", i, sep = "")
        assign(new_gr, ggroup(horizontal = TRUE, container = sel_gsa_f))
        addSpring(get(new_gr))
      }
      
      gsas_li[[i]] <- gcheckbox(gsas[i], container = get(new_gr))
      addSpring(get(new_gr))
      
    }
  }
  #in_area <- gdroplist(c( "All", paste("GSA - ", 1:27, sep = "")), container = sel_area_f)
  addSpring(sel_area_f)
  use_area_g <- ggroup(container = sel_area_f, horizontal = TRUE)
  glabel("Specify Area?", container = use_area_g)
  use_area_r <- gradio(c("Yes", "No"), container = use_area_g, horizontal = TRUE, 
                       handler = function(h,...)
                       {
                         enabled(sel_gsa_f) <- !enabled(sel_gsa_f)
                       })
  
  ### METIER
  
  sel_met_f <- gframe(text = "Metier", horizontal = FALSE, container = sel_vms_top)
  addSpring(sel_met_f)
  sel_met_g <- ggroup(horizontal = FALSE, container = sel_met_f)
  #gsa
  
  if(sqldf("select count(*) from vms_lb", dbname = vms_DB$db) == 0)
  {
    g_met_0 <- ggroup(horizontal = TRUE, container = sel_met_g)
    g_are_l <- glabel("No Metier data\nRun VMS-LB Match")
    enabled(sel_met_f) <- FALSE
    
  }else{
    met_fou <- sqldf("select distinct(met_des) from vms_lb order by met_des", dbname = vms_DB$db)[,1]
    n_met <- length(met_fou)
    mets_li <- vector("list", n_met)
    for(i in 1:n_met)
    {
      if(((i-1) %% 2) == 0 )
      {
        new_gr <- paste("g_met_", i, sep = "")
        assign(new_gr, ggroup(horizontal = TRUE, container = sel_met_g))
        addSpring(get(new_gr))
      }
      
      mets_li[[i]] <- gcheckbox(met_fou[i], container = get(new_gr))
      addSpring(get(new_gr))
      
    }
  }
  addSpring(sel_met_f)
  use_met_g <- ggroup(container = sel_met_f, horizontal = TRUE)
  glabel("Specify Metier?", container = use_met_g)
  use_area_r <- gradio(c("Yes", "No"), container = use_met_g, horizontal = TRUE, 
                       handler = function(h,...)
                       {
                         enabled(sel_met_g) <- !enabled(sel_met_g)
                       })
  
  
  ### FISH
  sel_fish_f <- gframe(text = "Fishing Point", horizontal = FALSE, container = sel_vms_top)
  addSpring(sel_fish_f)
  glabel("Fishing Only", container = sel_fish_f)
  fishin_r <- gradio(c("Yes", "No"), container = sel_fish_f, horizontal = TRUE)
  #   addSpring(sel_fish_f)
  
  
  ##############
  sel_go_g <- ggroup(horizontal = FALSE, expand = TRUE, container = sel_vms_top)
  sel_go_b <- gbutton(text = "\n\n\nGO\n\n\n", container = sel_go_g, handler = function(h,...)
  {
    from <- switch(svalue(from_cb), "Tracks" = "intrp", "Pings" = "ping")
    date_fro <- as.numeric(chron(dates. = dates(paste(which(levels(months(min_date)) == (svalue(fro_m))), svalue(fro_d), svalue(fro_y), sep = "/")), times. = times("00:00:00"), format = c(dates = "m/d/y", times = "h:m:s")))
    date_to <- as.numeric(chron(dates. = dates(paste(which(levels(months(max_date)) == (svalue(to_m))), svalue(to_d), svalue(to_y), sep = "/")), times. = times("00:00:00"), format = c(dates = "m/d/y", times = "h:m:s")))
    
    if(from == "ping")
    {

      que_vms_DB$que <- paste("select * from ", from, 
                              " where DATE > ", date_fro, 
                              " and DATE < ", date_to,
                              sep = "")
      
    }
    
    if(from == "intrp")
    {
      if(svalue(use_area_r) == "No")
      {
        
        if(svalue(fishin_r) == "Yes")
        {
          que_vms_DB$que <- paste("select * from ", from, 
                                  ", p_fish where DATE > ", date_fro, 
                                  " and DATE < ", date_to,
                                  " and FISH = 1 and intrp.rowid = i_id",
                                  sep = "")
        }else{
          que_vms_DB$que <- paste("select * from ", from, 
                                  " where DATE > ", date_fro, 
                                  " and DATE < ", date_to,
                                  sep = "")
        }
        
      }
      
      
      if(svalue(use_area_r) == "Yes")
      {
        date_fro <- as.numeric(chron(dates. = dates(paste(which(levels(months(min_date)) == (svalue(fro_m))), svalue(fro_d), svalue(fro_y), sep = "/")), times. = times("00:00:00"), format = c(dates = "m/d/y", times = "h:m:s")))
        date_to <- as.numeric(chron(dates. = dates(paste(which(levels(months(max_date)) == (svalue(to_m))), svalue(to_d), svalue(to_y), sep = "/")), times. = times("00:00:00"), format = c(dates = "m/d/y", times = "h:m:s")))
        
        num_are <- length(gsas_li)
        selec <- character(num_are)
        for(j in 1:num_are)
        {
          if(svalue(gsas_li[[j]]) == TRUE)
          {
            selec[j] <- gsas_li[[j]][]
          }
        }
        
        
        num_met <- length(mets_li)
        selem <- character(num_met)
        for(l in 1:num_met)
        {
          if(svalue(mets_li[[l]]) == TRUE)
          {
            selem[l] <- mets_li[[l]][]
          }
        }
        
        num_sel <- length(which(selec != ""))
        num_sem <- length(which(selem != ""))
        
        ############ SEL 0
        if(num_sel == 0) 
        {
          if(num_sem == 0)
          {
            if(svalue(fishin_r) == "Yes")
            {
              from <- paste(from, ", p_fish", sep = "")
              plu_fis <- " and FISH = 1 and intrp.rowid = i_id"
              que_vms_DB$que <- paste("select * from ", from, 
                                      " where DATE > ", date_fro, 
                                      " and DATE < ", date_to,
                                      plu_fis,
                                      sep = "")
            }else{
              que_vms_DB$que <- paste("select * from ", from, 
                                      " where DATE > ", date_fro, 
                                      " and DATE < ", date_to,
                                      sep = "")
            }
          }
          
          if(num_sem == 1)
          {
            mete <- paste(" and met_des = '",mets_li[[which(selem != "")]][],"'", sep = "")
            from <- paste(from, ", vms_lb", sep = "")
            plu_met <- " and intrp.I_NCEE = vms_lb.vessel and intrp.T_NUM = vms_lb.track"
            if(svalue(fishin_r) == "Yes")
            {
              from <- paste(from, ", p_fish", sep = "")
              plu_fis <- " and FISH = 1 and intrp.rowid = i_id"
              que_vms_DB$que <- paste("select * from ", from, 
                                      " where DATE > ", date_fro, 
                                      " and DATE < ", date_to,
                                      plu_fis,
                                      plu_met,
                                      mete,
                                      sep = "")
            }else{
              que_vms_DB$que <- paste("select * from ", from, 
                                      " where DATE > ", date_fro, 
                                      " and DATE < ", date_to,
                                      plu_met,
                                      mete,
                                      sep = "")
            }
          }
          
          if(num_sem > 1)
          {
            mete = " and( "
            for(k in 1:num_sem)
            {
              value <- which(selem != "")[k]
              if(k == num_sem){
                mete <- paste(mete, "met_des = '",
                              mets_li[[value]][],
                              "') ", sep = "")
              }else{
                mete <- paste(mete, "met_des = '",
                              mets_li[[value]][],
                              "' or ", sep = "")
              }
            }
            
            from <- paste(from, ", vms_lb", sep = "")
            plu_met <- " and intrp.I_NCEE = vms_lb.vessel and intrp.T_NUM = vms_lb.track"
            if(svalue(fishin_r) == "Yes")
            {
              from <- paste(from, ", p_fish", sep = "")
              plu_fis <- " and FISH = 1 and intrp.rowid = i_id"
              que_vms_DB$que <- paste("select * from ", from, 
                                      " where DATE > ", date_fro, 
                                      " and DATE < ", date_to,
                                      plu_fis,
                                      plu_met,
                                      mete,
                                      sep = "")
            }else{
              que_vms_DB$que <- paste("select * from ", from, 
                                      " where DATE > ", date_fro, 
                                      " and DATE < ", date_to,
                                      plu_met,
                                      mete,
                                      sep = "")
            }
            
          }
        }
        
        ############ SEL 1
        if(num_sel == 1) 
        {
          
          area <- paste(" and AREA = ",gsas_li[[which(selec != "")]][], sep = "")
          from <- paste(from, ", p_area", sep = "")
          plu_are <- " and intrp.T_NUM = p_area.T_NUM and I_NCEE = p_area.vess_id"
          
          if(num_sem == 0)
          {
            if(svalue(fishin_r) == "Yes")
            {
              from <- paste(from, ", p_fish")
              plu_fis <- " and FISH = 1 and intrp.rowid = i_id"
              que_vms_DB$que <- paste("select * from ", from, 
                                      " where DATE > ", date_fro, 
                                      " and DATE < ", date_to,
                                      plu_fis,
                                      plu_are,
                                      area,
                                      sep = "")
            }else{
              que_vms_DB$que <- paste("select * from ", from, 
                                      " where DATE > ", date_fro, 
                                      " and DATE < ", date_to,
                                      plu_are,
                                      area,
                                      sep = "")
            }
            
          }
          
          if(num_sem == 1)
          {
            mete <- paste(" and met_des = '",mets_li[[which(selem != "")]][], "'", sep = "")
            from <- paste(from, ", vms_lb", sep = "")
            plu_met <- " and intrp.I_NCEE = vms_lb.vessel and intrp.T_NUM = vms_lb.track"
            
            if(svalue(fishin_r) == "Yes")
            {
              from <- paste(from, ", p_fish", sep = "")
              plu_fis <- " and FISH = 1 and intrp.rowid = i_id"
              que_vms_DB$que <- paste("select * from ", from, 
                                      " where DATE > ", date_fro, 
                                      " and DATE < ", date_to,
                                      plu_fis,
                                      plu_met,
                                      mete,
                                      plu_fis,
                                      plu_are,
                                      area,
                                      sep = "")
            }else{
              que_vms_DB$que <- paste("select * from ", from, 
                                      " where DATE > ", date_fro, 
                                      " and DATE < ", date_to,
                                      plu_met,
                                      mete,
                                      plu_are,
                                      area,
                                      sep = "")
            }
          }
          
          if(num_sem > 1)
          {
            mete = " and( "
            for(k in 1:num_sem)
            {
              value <- which(selem != "")[k]
              if(k == num_sem){
                mete <- paste(mete, "met_des = '",
                              mets_li[[value]][],
                              "') ", sep = "")
              }else{
                mete <- paste(mete, "met_des = '",
                              mets_li[[value]][],
                              "' or ", sep = "")
              }
            }
            
            from <- paste(from, ", vms_lb", sep = "")
            plu_met <- " and intrp.I_NCEE = vms_lb.vessel and intrp.T_NUM = vms_lb.track"
            
            
            if(svalue(fishin_r) == "Yes")
            {
              from <- paste(from, ", p_fish", sep = "")
              plu_fis <- " and FISH = 1 and intrp.rowid = i_id"
              que_vms_DB$que <- paste("select * from ", from, 
                                      " where DATE > ", date_fro, 
                                      " and DATE < ", date_to,
                                      plu_fis,
                                      plu_met,
                                      mete,
                                      plu_fis,
                                      plu_are,
                                      area,
                                      sep = "")
            }else{
              que_vms_DB$que <- paste("select * from ", from, 
                                      " where DATE > ", date_fro, 
                                      " and DATE < ", date_to,
                                      plu_met,
                                      mete,
                                      plu_are,
                                      area,
                                      sep = "")
            }
            
            
          }
          
        }
        ############ SEL 2
        if(num_sel > 1)
        {
          area = " and( "
          for(k in 1:num_sel)
          {
            value <- which(selec != "")[k]
            if(k == num_sel){
              area <- paste(area, "AREA = ",
                            gsas_li[[value]][],
                            ") ", sep = "")
            }else{
              area <- paste(area, "AREA = ",
                            gsas_li[[value]][],
                            " or ", sep = "")
            }
          }
          from <- paste(from, ", p_area", sep = "")
          plu_are <- " and I_NCEE = p_area.vess_id and intrp.T_NUM = p_area.T_NUM"
          
          
          if(num_sem == 0)
          {
            if(svalue(fishin_r) == "Yes")
            {
              from <- paste(from, ", p_fish")
              plu_fis <- " and FISH = 1 and intrp.rowid = i_id"
              que_vms_DB$que <- paste("select * from ", from, 
                                      " where DATE > ", date_fro, 
                                      " and DATE < ", date_to,
                                      plu_fis,
                                      plu_are,
                                      area,
                                      sep = "")
            }else{
              que_vms_DB$que <- paste("select * from ", from, 
                                      " where DATE > ", date_fro, 
                                      " and DATE < ", date_to,
                                      plu_are,
                                      area,
                                      sep = "")
            }
            
          }
          
          if(num_sem == 1)
          {
            mete <- paste(" and met_des = '",mets_li[[which(selem != "")]][],"'", sep = "")
            from <- paste(from, ", vms_lb", sep = "")
            plu_met <- " and intrp.I_NCEE = vms_lb.vessel and intrp.T_NUM = vms_lb.track "
            
            if(svalue(fishin_r) == "Yes")
            {
              from <- paste(from, ", p_fish", sep = "")
              plu_fis <- " and FISH = 1 and intrp.rowid = i_id"
              que_vms_DB$que <- paste("select * from ", from, 
                                      " where DATE > ", date_fro, 
                                      " and DATE < ", date_to,
                                      plu_fis,
                                      plu_met,
                                      mete,
                                      plu_fis,
                                      plu_are,
                                      area,
                                      sep = "")
            }else{
              que_vms_DB$que <- paste("select * from ", from, 
                                      " where DATE > ", date_fro, 
                                      " and DATE < ", date_to,
                                      plu_met,
                                      mete,
                                      plu_are,
                                      area,
                                      sep = "")
            }
          }
          
          if(num_sem > 1)
          {
            mete = " and( "
            for(k in 1:num_sem)
            {
              value <- which(selem != "")[k]
              if(k == num_sem){
                mete <- paste(mete, "met_des = '",
                              mets_li[[value]][],
                              "') ", sep = "")
              }else{
                mete <- paste(mete, "met_des = '",
                              mets_li[[value]][],
                              "' or ", sep = "")
              }
            }
            
            from <- paste(from, ", vms_lb", sep = "")
            plu_met <- " and intrp.I_NCEE = vms_lb.vessel and intrp.T_NUM = vms_lb.track "
            
            
            if(svalue(fishin_r) == "Yes")
            {
              from <- paste(from, ", p_fish", sep = "")
              plu_fis <- " and FISH = 1 and intrp.rowid = i_id"
              que_vms_DB$que <- paste("select * from ", from, 
                                      " where DATE > ", date_fro, 
                                      " and DATE < ", date_to,
                                      plu_fis,
                                      plu_met,
                                      mete,
                                      plu_fis,
                                      plu_are,
                                      area,
                                      sep = "")
            }else{
              que_vms_DB$que <- paste("select * from ", from, 
                                      " where DATE > ", date_fro, 
                                      " and DATE < ", date_to,
                                      plu_met,
                                      mete,
                                      plu_are,
                                      area,
                                      sep = "")
            }
            
            
          }
          
        }
      }
    }
    
    #############
    
    #cat("\n", que_vms_DB$que, "\n", sep = "")
    
    que_res_df <- sqldf(que_vms_DB$que, dbname = vms_DB$db)
    que_res <- gtable(que_res_df, expand = T)
    add(query_nb, que_res, label = "Query")
    
    enabled(export_que) <- TRUE
    enabled(save_que) <- TRUE
  })
  
  
  ################
  export_que <- gbutton(text = "Export Query\n to CSV", container = sel_go_g, handler = function(h,...)
  {
    
    write.table(sqldf(que_vms_DB$que, dbname = vms_DB$db),
                file = paste(gfile(text = "Save VMS query",
                                   type = "save"),
                             ".csv", sep = ""),
                sep = ";", dec = ".")
    
  })
  enabled(export_que) <- FALSE
  
  save_que <- gbutton(text = "Save Query\nto new DB", container = sel_go_g, handler = function(h,...)
  {
    
    vms_db_win <- gwindow("Create New VMS DataBase", visible = TRUE)
    one <- ggroup(horizontal = FALSE, container = vms_db_win)
    glabel("To proceed with the creation of a VMS Database,
  first enter a name for the new DB and press enter...", container = one)
    c_name <- ggroup(horizontal = TRUE, container = one)
    addSpring(c_name)
    dbname <- gedit(initial.msg = "New VMS DB name...",
                    container = c_name, handler = function(h,...)
                    {
                      enabled(go) <- TRUE
                    }) 
    addSpring(c_name)
    gimage(system.file("ico/go-down-3.png", package="vmsbase"), container = one)
    glabel("Provide a destination for the VMS DB 
      clicking on the \'Create DB\' button.", container = one)
    go <- gbutton("Create DB", container = one, handler = function(h,...)
    {
      
      que_vms_DB$dir <- gfile(text = "Select VMS DB destination",
                              type = "selectdir",
                              filter = list("VMS DB data" = list(patterns = c("*.vms.sqlite"))))
      
      que_vms_DB$db <- paste(que_vms_DB$dir, "/", svalue(dbname), ".vms.sqlite", sep = "")  
      
      request <- sqldf(que_vms_DB$que, dbname = vms_DB$db)
      
      dbConnect(SQLite(), dbname = que_vms_DB$db)
      
      cat("\nCreating tables... ", sep = "")
      sqldf("CREATE TABLE query AS SELECT * FROM `request`", dbname = que_vms_DB$db)
      
      sqldf("CREATE TABLE ping(I_NCEE INT, LAT REAL, LON REAL, DATE REAL, SPE REAL, HEA INT)", dbname = que_vms_DB$db)
      sqldf("CREATE TABLE warn(p_id INT, W_DUPL INT, W_HARB INT, W_LAND INT, W_COHE INT)", dbname = que_vms_DB$db)
      sqldf("CREATE TABLE track(I_NCEE INT, LAT REAL, LON REAL, DATE REAL, SPE REAL, HEA REAL, W_HARB INT, T_NUM INT, P_ID INT)", dbname = que_vms_DB$db)
      sqldf("CREATE TABLE intrp(I_NCEE INT, LAT REAL, LON REAL, DATE REAL, SPE REAL, HEA REAL, W_HARB INT, T_NUM INT, P_ID INT, P_INT INT, T_ID INT)", dbname = que_vms_DB$db)
      sqldf("CREATE TABLE p_depth(i_id INT, vess_id INT, DEPTH REAL)", dbname = que_vms_DB$db)
      sqldf("CREATE TABLE p_area(vess_id INT, t_num INT, AREA INT)", dbname = que_vms_DB$db)
      sqldf("CREATE TABLE vms_lb(vessel INT, track INT, logbook INT, log_id INT, met_des CHAR)", dbname = que_vms_DB$db)
      sqldf("CREATE TABLE p_fish(i_id INT, F_SPE INT, F_DEP INT, F_DIS INT, FISH INT)", dbname = que_vms_DB$db)
      
      if(svalue(from_cb) == "Pings")
      {
        sqldf("INSERT INTO ping SELECT * FROM `request`", dbname = que_vms_DB$db)
      }
      
      if(svalue(from_cb) == "Tracks")
      {
        cat("\nLoading interpolated... ", sep = "")
        ntr_dat <- request[,1:11]
        sqldf("INSERT INTO intrp SELECT * FROM `ntr_dat`", dbname = que_vms_DB$db)
        
        cat("\nLoading tracks... ", sep = "")
        tracks <- request[!is.na(request[,"T_ID"]), "T_ID"]
        tra_dat <- data.frame()
        for(i in 1:length(tracks))
        {
          #cat(".", sep = "")
          tra_dat <- rbind(tra_dat, 
                           fn$sqldf("select * from track where ROWID = `tracks[i]`", 
                                    dbname = vms_DB$db))
        }
        sqldf("INSERT INTO track SELECT * FROM `tra_dat`", dbname = que_vms_DB$db)
        cat(" Complete!", sep = "")
        cat("\nLoading pings... ", sep = "")
        pings <- request[!is.na(request[,"P_ID"]), "P_ID"]
        pin_dat <- data.frame()
        war_dat <- data.frame()
        for(j in 1:length(pings))
        {
          #cat(".", sep = "")
          pin_dat <- rbind(pin_dat, 
                           fn$sqldf("select * from ping where ROWID = `pings[j]`", 
                                    dbname = vms_DB$db))
          war_dat <- rbind(war_dat, 
                           fn$sqldf("select * from warn where ROWID = `pings[j]`", 
                                    dbname = vms_DB$db))
        }
        cat(" Complete!", sep = "")
        
        sqldf("INSERT INTO ping SELECT * FROM `pin_dat`", dbname = que_vms_DB$db)
        sqldf("INSERT INTO warn SELECT * FROM `war_dat`", dbname = que_vms_DB$db)
        
        cat("\nLoading depth, area, match... ", sep = "")
        vssl <- unique(request[,"I_NCEE"])
        for(l in 1:length(vssl))
        {
          cat(".", sep = "")
          are_dat <- fn$sqldf("select * from p_area where vess_id = `vssl[l]`", 
                   dbname = vms_DB$db)
          if(nrow(are_dat) != 0){
          sqldf("INSERT INTO p_area SELECT * FROM `are_dat`", dbname = que_vms_DB$db)
          }
          cou_dat <- fn$sqldf("select * from vms_lb where vessel = `vssl[l]`", 
                              dbname = vms_DB$db)
          if(nrow(cou_dat) != 0){
          sqldf("INSERT INTO vms_lb SELECT * FROM `cou_dat`", dbname = que_vms_DB$db)
          }
          dep_dat <- fn$sqldf("select * from p_depth where vess_id = `vssl[l]`", 
                              dbname = vms_DB$db)
          if(nrow(dep_dat) != 0){
          sqldf("INSERT INTO p_depth SELECT * FROM `dep_dat`", dbname = que_vms_DB$db)
          }
        }
        cat(" Complete!", sep = "")
        cat("\nLoading fishing points... ", sep = "")
        fisi <- request[which(!is.na(request[,"i_id"])), "i_id"]
        for(m in 1:length(fisi))
        {
          cat(".", sep = "")
          fis_dat <- fn$sqldf("select * from p_fish where i_id = `fisi[m]`", dbname = vms_DB$db)
          if(nrow(fis_dat) != 0){
          sqldf("INSERT INTO p_fish SELECT * FROM `fis_dat`", dbname = que_vms_DB$db)
          }
        }
        cat(" Complete!", sep = "")
        cat("\n\n   ---   VMS DataBase Deploy Complete   ---\n\n", sep = "")
      }
      
      gconfirm("VMS DataBase deploy complete!",
               title = "Confirm",
               icon = "info",
               parent = vms_db_win,
               handler = function(h,...){dispose(vms_db_win)})
    })
    
    enabled(go) <- FALSE
    
  })
  enabled(save_que) <- FALSE
  
  query_nb <- gnotebook(tab.pos = 3, closebuttons = TRUE, container = sel_vms_bot, expand = T)
  
  visible(sel_vms_win) <- TRUE
  
}






#         if(num_sel == 1) 
#         {
#           
#           area <- paste("and AREA = ",gsas_li[[which(selec != "")]][], sep = "")
#           from <- paste(from, ", p_area", sep = "")
#           plu_are <- " and intrp.T_NUM = p_area.T_NUM "
#           if(svalue(fishin_r) == "Yes")
#           {
#             from <- paste(from, ", p_fish")
#             plu_fis <- " and I_NCEE = vess_id and FISH = 1 and intrp.rowid = i_id"
#             que_vms_DB$que <- paste("select * from ", from, 
#                                     " where DATE > ", date_fro, 
#                                     " and DATE < ", date_to,
#                                     plu_fis,
#                                     plu_are,
#                                     area,
#                                     sep = "")
#           }else{
#             que_vms_DB$que <- paste("select * from ", from, 
#                                     " where DATE > ", date_fro, 
#                                     " and DATE < ", date_to,
#                                     plu_are,
#                                     area,
#                                     sep = "")
#           }
#           
#           
#           
#           
#           
#         }
#         



#  
#         if(num_sel > 1)
#         {
#           area = " and( "
#           for(k in 1:num_sel)
#           {
#             value <- which(selec != "")[k]
#             if(k == num_sel){
#               area <- paste(area, "AREA = ",
#                             gsas_li[[value]][],
#                             ") ", sep = "")
#             }else{
#               area <- paste(area, "AREA = ",
#                             gsas_li[[value]][],
#                             " or ", sep = "")
#             }
#           }
#           
#           if(svalue(fishin_r) == "Yes")
#           {
#             que_vms_DB$que <- paste("select * from ", from, 
#                                     ", p_area, p_fish where I_NCEE = vess_id and intrp.rowid = i_id and FISH = 1 and ",
#                                     from, ".T_NUM = p_area.T_NUM ",
#                                     area, 
#                                     " and DATE > ", date_fro, 
#                                     " and DATE < ", date_to,
#                                     sep = "")
#           }else{
#             
#             que_vms_DB$que <- paste("select * from ", from, 
#                                     ", p_area where I_NCEE = vess_id and ",
#                                     from, ".T_NUM = p_area.T_NUM ",
#                                     area, 
#                                     " and DATE > ", date_fro, 
#                                     " and DATE < ", date_to,
#                                     sep = "")
#           }
#         }
#         
#       }
#     }
#     