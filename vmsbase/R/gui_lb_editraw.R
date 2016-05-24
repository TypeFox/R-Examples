
#' Raw LogBook Editing GUI
#' 
#' The \code{gui_lb_editraw} function implements the graphical user interface for the
#' editing of a raw logbook dataset
#' 
#' This function, with a raw logbook dataset, produces an edited version of it. 
#' The user must select, for each mandatory field, the column name where the data is 
#' stored.
#'
#' @return This function does not return a value.
#'  After the execution the user is asked where to put the edited file.
#' 
#' @usage gui_lb_editraw()
#' 
#' @export gui_lb_editraw
#'

gui_lb_editraw <- function ()
{
  #   if (workDir$path == "")
  #   {
  #     
  #     gmessage("\nPlease, select a working directory!", title = "Alert", icon = "warning")
  #     stop("Select a working directory first!")
  #     
  #   }
  
  rawfile <- log_File$new()
  
  # Initialization of the main window
  lb_editraw_win <- gwindow("Raw LogBook file editor", visible = FALSE)
  
  # Definition of the main group
  group1 <- ggroup(horizontal = F, spacing = 10, container = lb_editraw_win)
  #gimage(system.file("img/ico/logbook.png", package="vmstest"), container = group1)
  group <- gframe(horizontal = F, container = group1)
  
  # Load button and file path group
  g1 <- ggroup(container = group)
  addSpring(g1)
  gsep <- gframe(text = "Separator", container = g1)
  separator <- gradio(c(";",","), container = gsep)
  gdec <- gframe(text = "Decimal", container = g1)
  decimal <- gradio(c(".",","), container = gdec)
  addSpring(g1)
  g1a <- ggroup(horizontal = FALSE, container = g1)
  addSpring(g1a)
  loadbutt <- gbutton("Load Raw LogBook data", container=g1a,
                      handler = function(h,...)
                      {
                        rawfile$path <- gfile(type = "open")
                        
                        svalue(logfilew) <- rawfile$path
                        
                        firstline <- readLines(rawfile$path, n = 1)
                        
                        col_name <- gsub("\"", "", gsub(" ", ".", unlist(strsplit(firstline, svalue(separator)))))
                        
                        set_wdgt_vals(wdgt_w_val, col_name)
                        
                        turn_wdgt_on(wdgt_lst)
                        
                      })
  addSpring(g1a)
  addSpring(g1a)
  loadefl <- gbutton("Load Raw EFLALO data", container = g1a, handler = function(h,...)
  {
    
    rawfile$path <- gfile(type = "open")
    
    svalue(logfilew) <- rawfile$path
    
    enabled(lb_editraw_win) <- FALSE
    
    if (svalue(saveModeSel) == "New")
    {
      rawfile$data <- read.table(file = rawfile$path,
                                 header = TRUE,
                                 sep = svalue(separator),
                                 dec = svalue(decimal))
      
      write.table(saveRawEflalo(rawfile),
                  
                  file = paste(gfile(text = "Save edited LogBook file",
                                     type = "save",
                                     filter = list("logBook files" = list(patterns = c("*.logbook")))), ".logbook", sep =""),
                  sep = ",",
                  dec = ".",
                  quote = FALSE,
                  row.names = FALSE)
    }
    
    if (svalue(saveModeSel) == "Append")
    {
      rawfile$data <- read.table(file = rawfile$path,
                                 header = TRUE,
                                 sep = svalue(separator),
                                 dec = svalue(decimal))
      
      write.table(saveRawEflalo(rawfile),
                  
                  file = paste(gfile(text = "Save edited LogBook file",
                                     type = "save",
                                     filter = list("logBook files" = list(patterns = c("*.logbook")))), ".logbook", sep =""),
                  sep = ",",
                  dec = ".",
                  quote = FALSE,
                  append = T,
                  col.names = F,
                  row.names = F)
      
    }
    
    gconfirm("Raw LogBook data editing complete!",
             title = "Confirm",
             icon = "info",
             parent = lb_editraw_win,
             handler = function(h,...){dispose(lb_editraw_win)})
  })
  addSpring(g1a)
  addSpring(g1a)
  loadOpt <- gbutton("Load Option File", container=g1a,handler = function(h,...)
  {
    
    header <- read.csv(file = gfile(text = "Select option file to read", type = "open"),
                       header = F,
                       sep = "\n")
    
    header <- as.character(header[,1])
    
    set_wdgt_vals(wdgt_w_val, header)
    
    set_wdgt_sel(wdgt_lst, header)
    
    turn_wdgt_on(wdgt_lst)
    
    if(header[2] == "No")
    {
      enabled(sdate) <- F
      enabled(stime) <- F
    }
    
    if(header[5] == "No")
    {
      enabled(edate) <- F
      enabled(etime) <- F
    }
    
  })
  logfilew <- glabel("no file loaded", container=g1a)
  addSpring(g1a)
  addSpring(g1)
  
  
  # Vessel label and dropdown list menu
  g2 <- gframe(text = "Vessel Informations", container = group, spacing = 5, horizontal = F)
  g2u <- ggroup(horizontal = T, container = g2)
  addSpring(g2u)
  glabel("Vessel Id", container = g2u)
  vessUE <- gcombobox(c("Select column"), container = g2u)
  addSpring(g2u)
  g3_meti <- gframe(text = "Metier data", container = g2, spacing = 5, horizontal = T)
  addSpring(g3_meti)
  glabel("Available", container = g3_meti) 
  g3_met_r <- gradio(c("Yes","No"), horizontal = T, container = g3_meti,
                     handler = function(h,...){
                       ifelse(svalue(g3_met_r) == "No",enabled(the_meti) <- F, enabled(the_meti) <- T)
                     })  
  addSpring(g3_meti)
  the_meti <- gcombobox(c("Select column"), container = g3_meti)
  addSpring(g3_meti)
  
  g3_gear <- gframe(text = "Gear data", container = g2, spacing = 5, horizontal = T)
  addSpring(g3_gear)
  glabel("Available", container = g3_gear) 
  g3_gea_r <- gradio(c("Yes","No"), horizontal = T, container = g3_gear,
                     handler = function(h,...){
                       ifelse(svalue(g3_gea_r) == "No",enabled(the_gear) <- F, enabled(the_gear) <- T)
                     })  
  addSpring(g3_gear)
  the_gear <- gcombobox(c("Select column"), container = g3_gear)
  addSpring(g3_gear)
  
  g3 <- gframe(text = "Logbook Data", container = group, spacing = 5, horizontal = F)
  g3_frm <- ggroup(horizontal = TRUE, container = g3)
  addSpring(g3_frm)
  date_frm <- gradio(c("DD/MM/YYYY","MM/DD/YYYY"), horizontal = T, container = g3_frm)
  addSpring(g3_frm)
  g3_datim <- ggroup(horizontal = TRUE, container = g3)
  g3start <- gframe(text = "Departure Date & Time", container = g3_datim, spacing = 5, horizontal = T)
  addSpring(g3start)
  #   glabel("Available", container = g3start) 
  #   g3sq <- gradio(c("Yes","No"), horizontal = T, container = g3start,
  #                  handler = function(h,...){
  #                    ifelse(svalue(g3sq) == "No",
  # {
  #   enabled(sdate) <- F
  #   enabled(stime) <- F
  # },
  # {
  #   enabled(sdate) <- T
  #   enabled(stime) <- T
  # })
  # })    
  #   addSpring(g3start)
  g3s <-ggroup(horizontal = F, container = g3start)
  g3s1 <-ggroup(horizontal = T, container = g3s)
  glabel("Start date", container = g3s1)   
  sdate <- gcombobox(c("Select column"), container = g3s1)
  g3s2 <-ggroup(horizontal = T, container = g3s)
  glabel("Start time", container = g3s2)   
  stime <- gcombobox(c("Select column"), container = g3s2)
  addSpring(g3start)
  
  g3end <- gframe(text = "Arrival Date & Time", container = g3_datim, spacing = 5, horizontal = T)
  addSpring(g3end)
  #   glabel("Available", container = g3end) 
  #   g3eq <- gradio(c("Yes","No"), horizontal = T, container = g3end,
  #                  handler = function(h,...){
  #                    ifelse(svalue(g3eq) == "No",
  # {
  #   enabled(edate) <- F
  #   enabled(etime) <- F
  # },
  # {
  #   enabled(edate) <- T
  #   enabled(etime) <- T
  # })})  
  #   addSpring(g3end)
  g3e <-ggroup(horizontal = F, container = g3end)
  g3e1 <-ggroup(horizontal = T, container = g3e)
  glabel("End date", container = g3e1)   
  edate <- gcombobox(c("Select column"), container = g3e1)
  g3e2 <-ggroup(horizontal = T, container = g3e)
  glabel("End time", container = g3e2)   
  etime <- gcombobox(c("Select column"), container = g3e2)
  addSpring(g3end)
  
  g3cat <- gframe(text = "Catch Data", container = g3, spacing = 5, horizontal = T)
  addSpring(g3cat)
#   g3cat1 <-ggroup(horizontal = T, container = g3cat)
#   addSpring(g3cat1)
  glabel("Species", container = g3cat)   
  species <- gcombobox(c("Select column"), container = g3cat)
#   addSpring(g3cat1)
  addSpring(g3cat)
#   g3cat2 <-ggroup(horizontal = T, container = g3cat)
#   addSpring(g3cat2)
  glabel("Quantity", container = g3cat)
  qty <- gcombobox(c("Select column"), container = g3cat)
#   addSpring(g3cat2)
  addSpring(g3cat)
  
#   g3_extra <- gframe(text = "Metier & Gear", container = g3, spacing = 5, horizontal = F)
  
  
  
  # Saving, new, append, ok, cancel
  g8 <- gframe(container = group, spacing = 5, text = "Save & Exit")
  addSpace(g8, 10)
  
  saveModeSel <- gradio(c("New","Append"), horizontal = T, container = g8)
  addSpring(g8)
  saveOpt <- gcheckbox(text = "Save format options?", checked = F, container = g8)
  gbutton("ok", container = g8, handler = function(h,...)
  {
    enabled(lb_editraw_win) <- FALSE
    
    if (svalue(saveOpt) == T)
    { 
      write(get_wdgt_vals(wdgt_lst), 
            file = gfile(text = "Save LogBook template file", type = "save", initialfilename = "rawformat.rawlogbook"),
            sep = ",")
    }
    
    if (svalue(saveModeSel) == "New")
    {
      rawfile$data <- read.table(file = rawfile$path,
                                 header = TRUE,
                                 sep = svalue(separator),
                                 dec = svalue(decimal))
      
      write.table(saveRawLogBook(rawfile,
                                 get_wdgt_vals(wdgt_lst)),
                  
                  file = paste(gfile(text = "Save edited LogBook file",
                                     type = "save",
                                     filter = list("logBook files" = list(patterns = c("*.logbook")))), ".logbook", sep =""),
                  sep = ",",
                  dec = ".",
                  quote = FALSE,
                  row.names = FALSE)
    }
    
    if (svalue(saveModeSel) == "Append")
    {
      rawfile$data <- read.table(file = rawfile$path,
                                 header = TRUE,
                                 sep = svalue(separator),
                                 dec = svalue(decimal))
      
      write.table(saveRawLogBook(rawfile,
                                 get_wdgt_vals(wdgt_lst)),
                  file = gfile(text = "Append Edited LogBook Data to existing file:",
                               type = "save",
                               filter = list("logBook files" = list(patterns = c("*.logbook")))),
                  sep = ",",
                  dec = ".",
                  quote = FALSE,
                  append = T,
                  col.names = F,
                  row.names = F)
      
    }
    
    gconfirm("Raw LogBook data editing complete!",
             title = "Confirm",
             icon = "info",
             parent = lb_editraw_win,
             handler = function(h,...){dispose(lb_editraw_win)})
  })
  gbutton("cancel", handler = function(h,...) dispose(lb_editraw_win), container=g8)
  
  wdgt_lst <- c(vessUE,
                g3_met_r,
                the_meti,
                g3_gea_r,
                the_gear,
#                 g3sq,
                sdate, stime,
#                 g3eq,
                edate, etime,
                species, qty,
                date_frm)
  wdgt_w_val <- c(vessUE, sdate, stime, edate, etime, species, qty, the_meti, the_gear)
  
  turn_wdgt_off(wdgt_lst)
  
  visible(lb_editraw_win) <- TRUE
  
}

