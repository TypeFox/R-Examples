
###########################################
# CRUD operation for perspectives

##' Add / List / Remove perspective to / of / from emuDB
##' 
##' Add / List / Remove perspective to / of / from emuDB. The EMU-webApp subdivides different ways 
##' to look at an emuDB into so called perspectives. These perspectives, 
##' between which you can switch in the web application, contain 
##' information on what levels are displayed, which ssffTracks are drawn, 
##' and so on. For more information on the structural elements of an emuDB 
##' see \code{vignette{emuDB}}.
##' @param emuDBhandle emuDB handle as returned by \code{\link{load_emuDB}}
##' @param name name of perspective
##' @name AddListRemovePerspective
##' @keywords emuDB database DBconfig Emu 
##' @examples
##' \dontrun{
##' 
##' ##################################
##' # prerequisite: loaded ae emuDB 
##' # (see ?load_emuDB for more information)
##' 
##' # add perspective called "justTones" to the ae emuDB
##' add_perspective(emuDBhandle = ae,
##'                 name = "justTones") 
##'                 
##' # add levelCanvasOrder so only the "Tone" level is displayed
##' set_levelCanvasesOrder(emuDBhandle = ae, 
##'                        perspectiveName = "justTones", 
##'                        order = c("Tone"))
##' 
##' # list perspectives of ae emuDB
##' list_perspectives(emuDBhandle = ae)
##' 
##' # remove newly added perspective
##' remove_perspective(emuDBhandle = ae,
##'                    name = "justTones")
##'                    
##' }
##' 
NULL

##' @rdname AddListRemovePerspective
##' @export
add_perspective <- function(emuDBhandle, 
                            name){
  
  DBconfig = load_DBconfig(emuDBhandle)
  
  curPersp = list_perspectives(emuDBhandle)
  # check if level defined
  if(name %in% curPersp$name){
    stop("Perspective with name: '", name, "' already exists")
  }
  
  persp = list(name = name, 
               signalCanvases = list(order = c("OSCI", "SPEC"), 
                                     assign = NULL, contourLims = NULL),
               levelCanvases = list(order = NULL),
               twoDimCanvases = list(order = NULL))
  
  l = length(DBconfig$EMUwebAppConfig$perspectives)
  
  DBconfig$EMUwebAppConfig$perspectives[[l + 1]] = persp
  
  # store changes
  store_DBconfig(emuDBhandle, DBconfig)
  
}


##' @rdname AddListRemovePerspective
##' @export
list_perspectives <- function(emuDBhandle){
  
  DBconfig = load_DBconfig(emuDBhandle)
  df = data.frame(name = character(),
                  signalCanvasesOrder = character(),
                  levelCanvasesOrder = character(),
                  stringsAsFactors = F)
  
  for(p in DBconfig$EMUwebAppConfig$perspectives){
    df = rbind(df , data.frame(name = p$name,
                               signalCanvasesOrder = paste(p$signalCanvases$order, collapse = "; "),
                               levelCanvasesOrder = paste(p$levelCanvases$order, collapse = "; "),
                               stringsAsFactors = F))
  }
  
  return(df)
}


##' @rdname AddListRemovePerspective
##' @export
remove_perspective <- function(emuDBhandle, 
                               name){
  
  DBconfig = load_DBconfig(emuDBhandle)
  
  curPersp = list_perspectives(emuDBhandle)
  
  # check if perspective defined
  if(!name %in% curPersp$name){
    stop("No perspective with name: '", name, "' found!")
  }
  
  for(i in 1:length(DBconfig$EMUwebAppConfig$perspectives)){
    if(DBconfig$EMUwebAppConfig$perspectives[[i]]$name == name){
      DBconfig$EMUwebAppConfig$perspectives[[i]] = NULL
    }
  }
  # store changes
  store_DBconfig(emuDBhandle, DBconfig)
  
}

###########################################
# CRUD operation for signalCanvasesOrder


##' Set / Get signalCanvasesOrder of / to / from emuDB
##' 
##' Set / Get signalCanvasesOrder array that specifies which signals are 
##' displayed in the according perspective by the EMU-webApp. An entry in this character vector 
##' refers to either the name of an ssffTrackDefinition or a predefined string: \code{"OSCI"} which 
##' represents the oscillogram or \code{"SPEC"} which represents the 
##' spectrogram. For more information on the structural elements of an emuDB 
##' see \code{vignette{emuDB}}.
##' 
##' @param emuDBhandle emuDB handle as returned by \code{\link{load_emuDB}}
##' @param perspectiveName name of perspective
##' @param order character vector containig names of ssffTrackDefinitions or "OSCI" / "SPEC"
##' @name SetGetSignalCanvasesOrder
##' @keywords emuDB database DBconfig Emu
##' @examples 
##' \dontrun{
##' 
##' ##################################
##' # prerequisite: loaded ae emuDB 
##' # (see ?load_emuDB for more information)
##' 
##' # get signal canvas order of the "default"
##' # perspective of the ae emuDB
##' get_signalCanvasesOrder(emuDBhandle = ae, 
##'                         perspectiveName = "default")
##'                         
##' }
##' 
NULL

##' @rdname SetGetSignalCanvasesOrder
##' @export
set_signalCanvasesOrder <- function(emuDBhandle,
                                    perspectiveName,
                                    order){
  
  DBconfig = load_DBconfig(emuDBhandle)
  
  curTracks = c("OSCI", "SPEC", list_ssffTrackDefinitions(emuDBhandle)$name)
  
  #check if tracks given are defined
  for(t in order){
    if(!t %in% curTracks){
      stop("No ssffTrackDefinition present with name '", t, "'!")
    }
  }
  
  for(i in 1:length(DBconfig$EMUwebAppConfig$perspectives)){
    if(DBconfig$EMUwebAppConfig$perspectives[[i]]$name == perspectiveName){
      DBconfig$EMUwebAppConfig$perspectives[[i]]$signalCanvases$order = as.list(order)
      break
    }
  }
  
  # store changes
  store_DBconfig(emuDBhandle, DBconfig)
}


##' @rdname SetGetSignalCanvasesOrder
##' @export
get_signalCanvasesOrder <- function(emuDBhandle,
                                    perspectiveName){
  
  DBconfig = load_DBconfig(emuDBhandle)
  
  order = NA
  for(p in DBconfig$EMUwebAppConfig$perspectives){
    if(p$name == perspectiveName){
      order = unlist(p$signalCanvases$order)
    }
  }
  return(order)
}

###########################################
# CRUD operation for levelCanvasesOrder


##' Set / Get level canvases order of emuDB
##' 
##' Set / Get which levels of an emuDB to display as level canvases (in a 
##' given perspective of the EMU-webApp),
##' and in what order. Level canvases refer to levels of 
##' the type "SEGMENT" or "EVENT" that are displayed by the EMU-webApp. Levels 
##' of type "ITEM" can always be displayed using the hierarchy view of the
##' web application but can not be displayed as level canvases.
##' For more information on the structural elements of an emuDB 
##' see \code{vignette{emuDB}}.
##' 
##' @param emuDBhandle emuDB handle as returned by \code{\link{load_emuDB}}
##' @param perspectiveName name of perspective
##' @param order character vector containig names of levelDefinitions
##' @name SetGetlevelCanvasesOrder
##' @keywords emuDB database DBconfig Emu 
##' @examples 
##' \dontrun{
##' 
##' ##################################
##' # prerequisite: loaded ae emuDB 
##' # (see ?load_emuDB for more information)
##' 
##' # get level canvases order of ae emuDB
##' order = get_levelCanvasesOrder(emuDBhandle = ae,
##'                                perspectiveName = "default")
##' 
##' # reverse the level canvases order of ae emuDB
##' set_levelCanvasesOrder(emuDBhandle = ae
##'                        perspectiveName = "default",
##'                        order = rev(order))
##'                        
##' # get level canvases order of ae emuDB                       
##' get_levelCanvasesOrder(emuDBhandle = ae,
##'                        perspectiveName = "default")
##' }
##' 
NULL

##' @rdname SetGetlevelCanvasesOrder
##' @export
set_levelCanvasesOrder <- function(emuDBhandle,
                                   perspectiveName,
                                   order){
  
  DBconfig = load_DBconfig(emuDBhandle)
  
  curLevelNames = list_levelDefinitions(emuDBhandle)$name
  curLevelTypes = list_levelDefinitions(emuDBhandle)$type
  
  #check if levels given are defined and of correct type
  for(t in order){
    if(!t %in% curLevelNames){
      stop("No levelDefinition present with name '", t, "'!")
    }
    lt = curLevelTypes[curLevelNames == t]
    if(!lt %in% c("SEGMENT", "EVENT")){
      stop("levelDefinition with name '", t, "' is not of type 'SEGMENT' or 'EVENT'")
    }
  }
  
  for(i in 1:length(DBconfig$EMUwebAppConfig$perspectives)){
    if(DBconfig$EMUwebAppConfig$perspectives[[i]]$name == perspectiveName){
      DBconfig$EMUwebAppConfig$perspectives[[i]]$levelCanvases$order = as.list(order)
      break
    }
  }  
  # store changes
  store_DBconfig(emuDBhandle, DBconfig)
}


##' @rdname SetGetlevelCanvasesOrder
##' @export
get_levelCanvasesOrder <- function(emuDBhandle,
                                   perspectiveName){
  
  DBconfig = load_DBconfig(emuDBhandle)
  
  order = NA
  for(p in DBconfig$EMUwebAppConfig$perspectives){
    if(p$name == perspectiveName){
      order = unlist(p$levelCanvases$order)
    }
  }
  return(order)
}


# FOR DEVELOPMENT 
# library('testthat') 
# test_file("tests/testthat/test_aaa_initData.R")
# test_file('tests/testthat/test_emuR-database.DBconfig.EMUwebAppConfig.R')