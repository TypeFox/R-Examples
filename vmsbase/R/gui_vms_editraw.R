
#' Raw VMS Editing GUI
#' 
#' The \code{gui_vms_editraw} function implements the graphical user interface for the
#' editing of a raw vms dataset
#' 
#' This function, with a raw vms dataset, produces an edited version of it. 
#' The user must select, for each mandatory field, the column name where the data is 
#' stored.
#'
#' @return This function does not return a value.
#'  After the execution the user is asked where to put the edited file.
#' 
#' @usage gui_vms_editraw()
#' 
#' @export gui_vms_editraw
#'

gui_vms_editraw <- function()
{
  
  # Initialization of a "vmsFile" object
  rawfile <- vms_File$new()
  
  # Initialization of the main window
  vms_editraw_win <- gwindow("Raw VMS file editor", visible = FALSE)
  
  # Definition of the main group
  group <- ggroup(horizontal = FALSE, container = vms_editraw_win)
  #group <- gframe(horizontal = FALSE, container = group1)
  
  # Load button and file path group
  g1 <- gframe(text = "Loading Options", container = group)
  addSpring(g1)
  gsep <- gframe(text = "Separator", container = g1)
  separator <- gradio(c(";",",","tab","space"), container = gsep)
  addSpring(g1)
  gdec <- gframe(text = "Decimal", container = g1)
  decimal <- gradio(c(".",","), container = gdec)
  addSpring(g1)
  g1a <- gframe(horizontal = FALSE, container = g1)
  addSpring(g1a)
  loadbutt <- gbutton("Load Raw VMS data", container = g1a,
                      handler = function(h,...)
                      {
                        rawfile$path <- gfile(type = "open")
                        
                        firstline <- readLines(rawfile$path, n = 1)
                        
                        col_name <- unlist(strsplit(firstline, switch(svalue(separator), ";" = ";", "," = ",", "tab" = "\t", "space" = " ")))
                        
                        set_wdgt_vals(wdgt_w_val, col_name)
                        
                        turn_wdgt_on(wdgt_lst)
                      })
  addSpring(g1a)
  loadtac <- gbutton("Load Raw TACSAT data", container = g1a,
                     handler = function(h,...)
                     {
                       svalue(separator) <- ","
                       svalue(decimal) <- "."
                       
                       rawfile$path <- gfile(type = "open")
                       
                       if(!is.na(rawfile$path))
                       {
                         firstline <- readLines(rawfile$path, n = 1)
                         
                         col_name <- unlist(strsplit(firstline, switch(svalue(separator), ";" = ";", "," = ",", "tab" = "\t", "space" = " ")))
                         
                         set_wdgt_vals(wdgt_w_val, col_name)
                         
                         header <- read.csv(file = system.file("extdata/tacsat_format.rawvms", package = "vmsbase"),
                                            header = F,
                                            sep = "\n")
                         
                         header <- as.character(header[,1])
                         
                         set_wdgt_vals(wdgt_w_val, header)
                         
                         set_wdgt_sel(wdgt_lst, header)
                         
                         turn_wdgt_on(wdgt_lst)
                       }
                     })
  addSpring(g1a)
  loadOpt <- gbutton("Load Option File", container = g1a,
                     handler = function(h,...)
                     {
                       
                       header <- read.csv(file = gfile(text = "Select option file to read", type = "open"),
                                          header = F,
                                          sep = "\n")
                       
                       header <- as.character(header[,1])
                       
                       set_wdgt_vals(wdgt_w_val, header)
                       
                       set_wdgt_sel(wdgt_lst, header)
                       
                       turn_wdgt_on(wdgt_lst)
                     })
  addSpring(g1a)
  addSpring(g1)
  
  
  # Vessel label and dropdown list menu
  g2 <- gframe(text = "Vessel ID", container = group)
  addSpring(g2)
  vessId <- gcombobox(c("Select column"), container=g2)
  addSpring(g2)
  
  # Latitude label, mode selector, dropdown list menu
  g3up <- ggroup(horizontal = TRUE, container = group)
  addSpring(g3up)
  g3 <- gframe(text = "Latitude", container = g3up)
  latModeSel <- gradio(c("sex","dec"), container=g3)
  g3a <- ggroup(horizontal = FALSE, container = g3,)
  g3b <- gframe(text = "Deg + Min + Sec + Dir", horizontal = TRUE, container = g3a)
  latDeg <- gcombobox(c("Degrees"), container=g3b)
  latMin <- gcombobox(c("Minutes"), container=g3b)
  latSec <- gcombobox(c("Seconds"), container=g3b)
  latDir <- gcombobox(c("Direction"), container=g3b)
  g3c <- gframe(text = "Decimal", horizontal = TRUE, container = g3a)
  latDec <- gcombobox(c("Dec Lat"), container=g3c)
  addSpring(g3up)
  # Longitude label, mode selector, dropdown list menu
  g4 <- gframe(text = "Longitude", container = g3up)
  lonModeSel <- gradio(c("sex","dec"), container=g4)
  g4a <- ggroup(horizontal = FALSE, container = g4)
  g4b <- gframe(text = "Deg + Min + Sec + Dir", horizontal = TRUE, container = g4a)
  lonDeg <- gcombobox(c("Degrees"), container=g4b)
  lonMin <- gcombobox(c("Minutes"), container=g4b)
  lonSec <- gcombobox(c("Seconds"), container=g4b)
  lonDir <- gcombobox(c("Direction"), container=g4b)
  g4c <- gframe(text = "Decimal", horizontal = TRUE, container = g4a)
  lonDec <- gcombobox(c("Dec Lon"), container=g4c)
  addSpring(g3up)
  
  # Time label, mode, dropdown list
  g5 <- gframe(text = "Date & Time", container = group)
  addSpring(g5)
  timeModeSel <- gradio(c("UTC", "Date + Time", "Date + H M S"), container = g5)
  g5a <- ggroup(horizontal = F, container = g5)
  g5b <- gframe(text = "Time UTC", horizontal = T, container = g5a)
  timeUtc <- gcombobox(c("UTC Time"), container = g5b)
  
  g5c <- gframe(text = "Date + Time", horizontal = T, container = g5a)
  timeDate2a <- gcombobox(c("Date"), container = g5c)
  timeDate2b <- gcombobox(c("Time"), container = g5c)
  
  g5d <- gframe(text = "Date + H + M + S", horizontal = T, container = g5a)
  timeDate <- gcombobox(c("Date"), container = g5d)
  timeHour <- gcombobox(c("Hour"), container = g5d)
  timeMinute <- gcombobox(c("Minute"), container = g5d)
  timeSecond <- gcombobox(c("Second"), container = g5d)
    addSpring(g5)

  g5aa <- ggroup(horizontal = F, container = g5)
  addSpring(g5aa)
  g5ab <- gframe(text = "Date Format", horizontal = T, container = g5aa)
  date_frm <- gradio(c("DD/MM/YYYY","MM/DD/YYYY"), horizontal = F, container = g5ab)

  addSpring(g5)
  
  # Speed label, mode selector, dropdown list menu
  g6up <- ggroup(horizontal = TRUE, container = group)
  addSpring(g6up)
  g6 <- gframe(text = "Speed", container = g6up)
  speedModeSel <- gradio(c("Knots","Km/h"), container=g6)
  speedCol <- gcombobox(c("Speed"), container=g6)
  addSpring(g6up)
  # Heading label, mode selector, dropdown list menu
  g7 <- gframe(text = "Heading", container = g6up)
  headModeSel <- gradio(c("Rad","Deg"), container=g7)
  headCol <- gcombobox(c("Heading"), container=g7)
  addSpring(g6up)
  
  # Saving, new, append, ok, cancel
  g8 <- gframe(text = "Save & Exit", container = group)
  help_b <- gbutton("Help", container = g8,
                    handler = function(h,...)
                    {
                      gmessage(
                        "Loading Options: select the character to consider as Separator and Decimal characters in the raw vms file.
                        
Load raw data: button to select and load the raw vms file.
Load options: button to select and load an option file generated by the raw vms edit tool.
                        
Vessel ID: dropdown list to select the column in the raw vms file with the vessel id information.
                        
Latitude and Longitude: dropdown list to select the column/s in the raw vms file with the lat/lon information.
  There are two possible modalities of input, sexagesimal (degree, minutes, seconds and direction) or decimal.
                        
Date & time: dropdown list to select the column/s in the raw vms file with the date and time information.
  There are two possible modalities of input, UTC (number of days since 01/01/1970) or date & times (date as d/m/y + hour + minutes + seconds).
                        
Speed: dropdown list to select the column in the raw vms file with the speed information and must be specified if this is in kmh or in knots.
                        
Heading: dropdown list to select the column in the raw vms file with the Heading information and must be specified if this is in radiants or degrees.
                        
",
                        title = "Edit Raw VMS Help", icon = c("info"))
                      
                    })
  addSpring(g8)
  saveModeSel <- gradio(c("New","Append"), horizontal = T, container = g8)
  addSpring(g8)
  saveOpt <- gcheckbox(text = "Save format options?", checked = F, container = g8)
  gbutton("ok", container = g8, handler = function(h,...)
  {
    enabled(vms_editraw_win) <- FALSE
    
    if (svalue(saveOpt) == T)
    { 
      
      write(get_wdgt_vals(wdgt_lst), 
            file = gfile(text = "Save VMS template file", type = "save", initialfilename = "rawformat.rawvms"),
            sep = ",")
    }
    
    if (svalue(saveModeSel) == "New")
    {
      rawfile$data <- read.table(file = rawfile$path,
                                 header = TRUE,
                                 sep = svalue(separator),
                                 dec = svalue(decimal))
      rep_txt <- ""
      write.table(saveRawVms(rawfile, get_wdgt_vals(wdgt_lst)),
                  
                  file = paste(gfile(text = "Save edited VMS file",
                                     type = "save",
                                     filter = list("VMS files" = list(patterns = c("*.vms")))), ".vms", sep = ""),
                  sep = ";",
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
      rep_txt = ""
      write.table(saveRawVms(rawfile,
                             get_wdgt_vals(wdgt_lst)),
                  
                  file = gfile(text = "Append VMS file to...",
                               type = "save",
                               filter = list("VMS Files" = list(patterns = c("*.vms")))),
                  sep = ";",
                  dec = ".",
                  quote = FALSE,
                  append = T,
                  col.names = F,
                  row.names = F)
    }
    
    gconfirm(paste("Raw VMS data editing complete!\n\n", rep_txt, sep = ""),
             title = "Confirm",
             icon = "info",
             parent = vms_editraw_win,
             handler = function(h,...){dispose(vms_editraw_win)})
  })
  
  
  gbutton("cancel", handler = function(h,...){dispose(vms_editraw_win)}, container=g8)
  addSpring(g8)
  
  wdgt_lst <- c(vessId,
                latModeSel, latDeg, latMin, latSec, latDir, latDec,
                lonModeSel, lonDeg, lonMin, lonSec, lonDir, lonDec,
                timeModeSel, timeUtc, timeDate2a, timeDate2b, timeDate, timeHour, timeMinute, timeSecond,
                speedModeSel, speedCol,
                headModeSel, headCol, date_frm)
  
  wdgt_w_val <- c(vessId,
                  latDeg, latMin, latSec, latDir, latDec,
                  lonDeg, lonMin, lonSec, lonDir, lonDec,
                  timeUtc, timeDate2a, timeDate2b, timeDate, timeHour, timeMinute, timeSecond,
                  speedCol,
                  headCol)
  
  turn_wdgt_off(wdgt_lst)
  
  visible(vms_editraw_win) <- TRUE
  
}