################################################################################
# TODO LIST
# TODO: ...

################################################################################
# CHANGE LOG (last 20 changes)
# 29.08.2015: Added importFrom.
# 14.12.2014: Fixed "Error in read.table..." when 'Save' is pressed without data.
# 14.12.2014: Fixed "Error in basename(x) : a character vector argument expected"
# 14.12.2014: Changed single gender marker to general sex markers (multiple).
# 11.10.2014: Added 'focus', added 'parent' parameter.
# 28.07.2014: Changed some gui text.
# 28.06.2014: Added help button and moved save gui checkbox.
# 10.11.2013: Fixed check that short name is provided.
# 22.09.2013: Fixed duplicate check.
# 21.09.2013: Fixed correct gender marker when reading from file.
# 21.09.2013: Fixed no gender marker by putting 'NA' as option.
# 22.06.2013: First version.

#' @title Make Kit
#'
#' @description
#' Add new kits or edit the kit file.
#'
#' @details A graphical user interface for reading information from 'bins' and 
#' 'panels' file for the creation of additional kits. It is also possible to
#' edit the short and full name of existing kits or removing kits.
#' The gender marker of each kits is auto detected but can be changed manually.
#' #' NB! Short name must be unique.
#' @param env environment in wich to search for data frames.
#' @param savegui logical indicating if GUI settings should be saved in the environment.
#' [Not currently used]
#' @param debug logical indicating printing debug information.
#' @param parent widget to get focus when finished.
#' 
#' @return TRUE
#' 
#' @export
#' 
#' @importFrom utils write.table read.delim help head tail  
#' 
#' @seealso \code{\link{readBinsFile}}, \code{\link{readPanelsFile}}, \code{\link{combineBinsAndPanels}}

makeKit_gui <- function(env=parent.frame(), savegui=NULL, debug=FALSE, parent=NULL){
  
  if(debug){
    print(paste("IN:", match.call()[[1]]))
  }
  
  # Global variables.
  .noSexMarkerString <- "<none>"
  .f3g1 <- NULL
  .separator <- .Platform$file.sep # Platform dependent path separator.
  .packagePath <- path.package("strvalidator", quiet = FALSE)
  .subFolder <- "extdata"
  .fileName <- "kit.txt"
  .filePath <- paste(.packagePath, .subFolder, .fileName, sep=.separator)
  .newKitInfo <- NULL
  .kitInfo <- NULL
  .binsFiles <- NULL
  .panelsFiles <- NULL
  
  if(debug){
    print("File path")
    print(.filePath)
  }
  
  # Load existing kit file.  
  if(!is.na(.filePath)){
    if(file.exists(.filePath)){
      # Read kit info file.
      .kitInfo <- read.delim(file=.filePath, header = TRUE, 
                                 sep = "\t", quote = "\"",
                                 dec = ".", fill = TRUE,
                                 stringsAsFactors=FALSE)
    } else {
      gmessage(message="The kit file was not found",
               title="File not found",
               icon = "error",
               parent = w) 
    }
  }        
  
  
  # Main window.  
  w <- gwindow(title="Manage kits", visible=FALSE)
  
  # Runs when window is closed.
  addHandlerDestroy(w, handler = function (h, ...) {
    
    # Save GUI state.
    #.saveSettings()
    
    # Focus on parent window.
    if(!is.null(parent)){
      focus(parent)
    }
    
  })
  
  # Vertical main group.
  gv <- ggroup(horizontal=FALSE,
               spacing=8,
               use.scrollwindow=FALSE,
               container = w,
               expand=TRUE) 

  # Help button group.
  gh <- ggroup(container = gv, expand=FALSE, fill="both")
  
  # No saved options yet!
  # savegui_chk <- gcheckbox(text="Save GUI settings", checked=FALSE, container=gh)
  
  addSpring(gh)
  
  help_btn <- gbutton(text="Help", container=gh)
  
  addHandlerChanged(help_btn, handler = function(h, ...) {
    
    # Open help page for function.
    print(help("makeKit_gui", help_type="html"))
    
  })
  
  # FRAME 0 ###################################################################
  
  f0 <- gframe(text = "STR Kits",
               horizontal=FALSE,
               spacing = 5,
               container = gv) 
  
  f0_opt <- gradio(items=c("Edit kit file", "Add new kits"), 
                                  selected = 2, container=f0)
  
  addHandlerChanged(f0_opt, handler = function (h, ...) {
    
    val_obj <- svalue(f0_opt, index=TRUE)
    
    if(debug){
      print(val_obj)
    }
    
     # Clear.
     if(!is.null(.f3g1)){
       delete(obj=.f3g1, widget=f3)
     }

    if(val_obj == 1){
      
      # Enable 'edit' objects.
      enabled(f1) <- TRUE
      
      # Disable 'new' objects.
      enabled(f2) <- FALSE
      
      # Autoselect 'Overwrite'.
      svalue(f4_opt, index=TRUE) <- 2
      
      # Update path.
      svalue(f1g1_file_edt) <- .filePath
      
    } else {
      
      # Enable 'new' objects.
      enabled(f2) <- TRUE
      
      # Disable 'edit' objects.
      enabled(f1) <- FALSE
      
      # Autoselect 'Overwrite'.
      svalue(f4_opt, index=TRUE) <- 1
      
    }
    
  } )  
  
  # FRAME 1 ###################################################################
  
  f1 <- gframe(text = "Kit file",
               horizontal=FALSE,
               spacing = 5,
               container = gv) 

  # This is disabled by default.
  enabled(f1) <- FALSE
  
  f1g1 <- glayout(container = f1, spacing = 1)
  
  f1g1[1,1] <- f1g1_file_edt <- gedit(text=.filePath, container=f1g1, width=80)
  f1g1[1,2] <- f1g1_file_btn <- gbutton(text="Load", container=f1g1)
  
  addHandlerChanged(f1g1_file_btn, handler = function (h, ...) {
    
    val_obj <- svalue(f1g1_file_edt)
    
    if(debug){
      print("Kit file:")
      print(val_obj)
    }

    if(!is.na(val_obj)){
      if(file.exists(val_obj)){
        
        # Read kit info file.
        .newKitInfo <<- read.delim(file=val_obj, header = TRUE, 
                               sep = "\t", quote = "\"",
                               dec = ".", fill = TRUE,
                               stringsAsFactors=FALSE)
        
        # Update GUI.
        .update(kitInfo=.newKitInfo, newKit=FALSE)
        
      } else {
        
        gmessage(message="The kit file was not found",
                 title="File not found",
                 icon = "error",
                 parent = w) 
      }
    }        
    
  } )  
  
  # FRAME 2 ###################################################################
  
  f2 <- gframe(text = "New kits",
               horizontal=FALSE,
               spacing = 5,
               container = gv) 
  
  f2g1 <- glayout(container = f2, spacing = 1)

  # BINS ----------------------------------------------------------------------
  
  f2g1[1,1] <- f2g1_bins_btn <- gbutton(text="Select Bins File",
                                        container=f2g1)

  f2g1[1,2] <- f2g1_bins_lbl <- glabel(text="", container=f2g1)

  addHandlerChanged(f2g1_bins_btn, handler = function (h, ...) {
    
    .binsFiles <<- gfile(text = "Select Bins file...",
                        type = "open", 
                        filter = list("text files" = list(mime.types = c("text/plain")),
                                      "All files" = list(patterns = c("*"))),
                        multi=FALSE)
    if(!is.na(.binsFiles)){
      svalue(f2g1_bins_lbl) <- basename(.binsFiles)
    }
    
  } )  
  
  # PANELS --------------------------------------------------------------------
  
  f2g1[2,1] <- f2g1_panels_btn <- gbutton(text="Select Panels File",
                                        container=f2g1)
  
  f2g1[2,2] <- f2g1_panels_lbl <- glabel(text="", container=f2g1)
  
  addHandlerChanged(f2g1_panels_btn, handler = function (h, ...) {
    
    .panelsFiles <<- gfile(text = "Select Panels file...",
                        type = "open", 
                        filter = list("text files" = list(mime.types = c("text/plain")),
                                      "All files" = list(patterns = c("*"))),
                        multi=FALSE)
    
    if(!is.na(.panelsFiles)){
      svalue(f2g1_panels_lbl) <- basename(.panelsFiles)
    }
    
  } )  
  
  # READ AND COMBINE ----------------------------------------------------------
  
  f2g1[3,1] <- f2g1_read_btn <- gbutton(text="Combine",
                                          container=f2g1)
  
  f2g1[3,2] <- glabel(text="Read and combine Bins and Panels files",
                      container=f2g1)
  
  addHandlerChanged(f2g1_read_btn, handler = function (h, ...) {
    
    if(debug){
      print("Bins file:")
      print(.binsFiles)
      print("Panels file:")
      print(.panelsFiles)
    }
    
    if(file.exists(.binsFiles) & file.exists(.panelsFiles)){

      # Read and combine files.
      .newKitInfo <<- combineBinsAndPanels(bin=readBinsFile(.binsFiles),
                           panel=readPanelsFile(.panelsFiles))    
      
      # Update GUI.
      .update(kitInfo=.newKitInfo, newKit=TRUE)
      
    } else {

      gmessage(message="One or several files was not found",
               title="File not found",
               icon = "error",
               parent = w) 
      
    }
    
  } )  
  
  # FRAME 3 ###################################################################
  
  f3 <- ggroup(horizontal=FALSE,
               use.scrollwindow=TRUE,
               expand=TRUE,
               fill="both",
               spacing = 5,
               container = gv) 
  
  # Function for updating GUI with kits.  
  .update <- function(kitInfo, newKit){
    # kitInfo - data.frame with kit information.
    # newKit - logical, if TRUE autodetect sex marker.
    #                   if FALSE read from file.
    
    # Get panels.
    panel <- unique(kitInfo$Panel)
    shortName <- unique(kitInfo$Short.Name)
    fullName <- unique(kitInfo$Full.Name)
    
    # Clear.
    if(!is.null(.f3g1)){
      delete(obj=.f3g1, widget=f3)
    }
    
    # Add container.
    .f3g1 <<- glayout(container = f3, spacing = 10, expand=TRUE, fill="both")
    
    # Add titles.
    .f3g1[1,1] <<- glabel(text="Remove", container=.f3g1)
    .f3g1[1,2] <<- glabel(text="Panel", container=.f3g1)
    .f3g1[1,3] <<- glabel(text="Short.Name", container=.f3g1)
    .f3g1[1,4] <<- glabel(text="Full.Name", container=.f3g1)
    .f3g1[1,5] <<- glabel(text="List sex markers (separate by comma)",
                          container=.f3g1)
    
    # Loop over panels and add objects.
    for(p in seq(along=panel)){
      
      # Get all markers.
      markers <- unique(kitInfo$Marker[kitInfo$Panel == panel[p]])
      
      # sex marker.
      sexMarkers <- NULL
      if(newKit){
        # Try to autodetect sex marker.
        sexMarkers <- grep("AM|Y", markers, ignore.case=TRUE, value=TRUE)
        if(length(sexMarkers) == 0){
          sexMarkers <- .noSexMarkerString
        }
      } else {
        # Get current sex marker.
        sexMarkers <- unique(kitInfo$Marker[kitInfo$Panel == panel[p] & kitInfo$Sex.Marker])
        if(length(sexMarkers) == 0){
          # If no matching marker, set to no sex marker string.
          sexMarkers <- .noSexMarkerString
        }
      }
      # Collapse to string.
      sexMarkers <- paste(sexMarkers, collapse=",")
      
      if(debug){
        print("newKit:")
        print(newKit)
        print("markers:")
        print(markers)
        print("sexMarkers:")
        print(sexMarkers)
      }

      # Add widgets, and populate.
      .f3g1[p + 1, 1] <<- gcheckbox(text="", checked=FALSE, container=.f3g1)
      .f3g1[p + 1, 2] <<- glabel(text=panel[p], container=.f3g1)
      .f3g1[p + 1, 3] <<- gedit(text=shortName[p], width = 20, container=.f3g1)
      .f3g1[p + 1, 4] <<- gedit(text=fullName[p], width = 40, container=.f3g1)
      .f3g1[p + 1, 5] <<- gedit(text=sexMarkers, width = 30, container=.f3g1)
      
    }
  
  }
  
  # FRAME 4 ###################################################################
  
  f4 <- gframe(text="Save as",
               horizontal=FALSE,
               expand=FALSE,
               spacing = 10,
               container = gv) 
  
  # SAVE ----------------------------------------------------------------------
  
  f4_opt <- gradio(items=c("Append to kit file",
                           "Overwrite kit file",
                           "Save as data frame"), 
                   selected = 1, 
                   horizontal = FALSE,
                   container=f4)
  
  f4_name_edt <- gedit(width=25, container=f4)
  enabled(f4_name_edt) <- FALSE
  
  f4_save_btn <- gbutton(text="Save", border=TRUE, expand=FALSE, container=f4)
  
  
  addHandlerChanged(f4_opt, handler = function (h, ...) {
    
    val_obj <- svalue(f4_opt, index=TRUE)
    
    if(debug){
      print(val_obj)
    }
    
    if(val_obj == 1){
      
      # Disable.
      enabled(f4_name_edt) <- FALSE
      
      # Update path.
      svalue(f1g1_file_edt) <- .filePath
      
    } else if (val_obj == 2) {
      
      # Disable.
      enabled(f4_name_edt) <- FALSE

    } else if (val_obj == 3) {
      
      # Enable.
      enabled(f4_name_edt) <- TRUE
      
    }
    
  } )  
  
  addHandlerChanged(f4_save_btn, handler = function (h, ...) {
    
    # Get variables.
    val_name <- svalue(f4_name_edt)
    val_opt <- svalue(f4_opt, index=TRUE)
    
    # Check if kit info exist.
    if(!is.null(.newKitInfo)){

      # Initiate vectors.
      removeKit <- logical()
      shortName <- character()
      fullName <- character()
      sexMarkers <- 0
      
      # Get panels.
      panel <- unique(.newKitInfo$Panel)
      
      # Get data.
      for(p in seq(along=panel)){
        removeKit[p] <- svalue(.f3g1[1 + p, 1])
        shortName[p] <- svalue(.f3g1[1 + p, 3])
        fullName[p] <- svalue(.f3g1[1 + p, 4])
        
        # Read sex marker value and replace any spaces.
        sexTmp <- svalue(.f3g1[1 + p, 5])
        sexTmp <- gsub(" ", "", sexTmp, fixed=TRUE)
        sexMarkers[p] <- sexTmp
      }
      
      if(debug){
        print("Short.Name:")
        print(shortName)
        print("Full.Name:")
        print(fullName)
        print("sexMarkers:")
        print(sexMarkers)
      }
      
      # Check that short name is provided for all kits not removed.
      missing <- shortName[!removeKit] %in% ""
      
      # Check if short name is missing.
      if(!any(missing)){
        
        if(val_opt == 1){ # Add new.
          # Check if short name exist in kit file.
          exist <- shortName[!removeKit] %in% getKit()
        } else { # Overwrite or save as data frame.
          # Check if short name exist in new kits.
          exist <- duplicated(toupper(shortName[!removeKit]))
        }
        
        # Check if short name exist.
        if(!any(exist)){
          
          # Change button.
          svalue(f4_save_btn) <- "Saving..."
          enabled(f4_save_btn) <- FALSE
          
          for(p in seq(along=panel)){
            
            if(debug){
              print("Panel:")
              print(panel[p])
            }
            
            # Crete columns and add data to kits:
            # Set short name.
            if(is.null(.newKitInfo$Short.Name)){
              .newKitInfo$Short.Name <<- NA
            }
            .newKitInfo$Short.Name[.newKitInfo$Panel == panel[p]] <<- shortName[p]
            
            # Set full name.
            if(is.null(.newKitInfo$Full.Name)){
              .newKitInfo$Full.Name <<- NA
            }
            .newKitInfo$Full.Name[.newKitInfo$Panel == panel[p]] <<- fullName[p]
            
            # Set sex marker flag.
            if(is.null(.newKitInfo$Sex.Marker)){
              .newKitInfo$Sex.Marker <<- NA
            }
            currentSexMarkers <- unlist(strsplit(sexMarkers[p], split=",", fixed=TRUE))
            selPanel <- .newKitInfo$Panel == panel[p]
            selMarker <- .newKitInfo$Marker %in% currentSexMarkers
            selection <- selPanel & selMarker
            .newKitInfo$Sex.Marker[selPanel] <<- FALSE # Reset current panel.
            .newKitInfo$Sex.Marker[selection] <<- TRUE # Flag sex markers.
            # Check for misspelled markers.
            ok <- currentSexMarkers %in% .newKitInfo$Marker[selPanel]
            if(!all(ok)){
              warning(paste("Given sex marker:",
                            paste(currentSexMarkers[!ok], collapse=","),
                            "not found in", panel[p]))
            }
            
          }
          
          # Remove kits if any.
          if(any(removeKit)){
            
            if(debug){
              print("panel:")
              print(panel)
              print("removeKit:")
              print(removeKit)
            }
            
            removePanel <- panel[removeKit]
            
            if(debug){
              print("removePanel:")
              print(removePanel)
            }
            
            for(p in seq(along=removePanel)){
              .newKitInfo <<- .newKitInfo[.newKitInfo$Panel != removePanel[p], ]
            }
            
          }  
          
          if(debug){
            print("newKitInfo:")
            print(head(.newKitInfo))
            print(tail(.newKitInfo))
          }
          
          # If append.
          if(val_opt == 1){
            if(debug){
              print("Append:")
              print(head(.kitInfo))
              print(head(.newKitInfo))
            }
            #.newKitInfo <<- rbind(.kitInfo, .newKitInfo)
            write.table(x=.newKitInfo, file=.filePath,
                        row.names = FALSE, 
                        append = TRUE,
                        col.names = FALSE,
                        quote = FALSE, sep = "\t")
          } else if(val_opt == 2){
            if(debug){
              print("Save as kit file")
            }
            write.table(x=.newKitInfo, file=.filePath,
                        row.names = FALSE,
                        append = FALSE,
                        col.names = TRUE,
                        quote = FALSE, sep = "\t")
          } else if(val_opt == 3){
            if(debug){
              print("Save as data frame")
            }
            assign(val_name, .newKitInfo, envir=env)
          }
          
          if(debug){
            print(paste("EXIT:", match.call()[[1]]))
          }
          
          # Close GUI.
          dispose(w)
          
        } else {
          
          message <- paste("A kit with short name",
                           shortName[exist][1], 
                           "already exist!\n\n",
                           "Short name must be unique!")
          
          gmessage(message, title="Duplicate short name",
                   icon = "error",
                   parent = w) 
          
        }
        
      } else {
        
        message <- "A short name must be provided for all new kits"
        
        gmessage(message, title="Missing short name",
                 icon = "error",
                 parent = w) 
        
      }
      
    }
    
  } )  

  # Show GUI.
  visible(w) <- TRUE
  focus(w)
  
}
