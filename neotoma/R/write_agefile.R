#' Write age control file to disk formatted for either Bacon or Clam
#'
#' Passing in a download object the function outputs a Bacon or Clam formatted file to a
#' user defined destination for age modelling with existing age-depth modeling software.
#'
#' @importFrom utils write.csv write.table
#' @param download A single site returned by \code{get_download}.
#' @param chronology Default is \code{1}, the default chronology for the core.  If a core has more than one chronology the user can define a different set of chronological controls.
#' @param path The location of the 'Cores' folder & working directory for Bacon.  Do not include "Cores" in the path name.
#' @param corename The intended handle for the core, to be used in writing to file.
#' @param cal.prog The method intended to build the age model, either \code{'Bacon'} or \code{'Clam'}.
#' 
#' @author Simon J. Goring \email{simon.j.goring@@gmail.com}
#' 
#' @return This command returns a file in location \code{path/Cores} containing all the relevant information required to build either the default or prior chronology for a core.
#' 
#' @examples \dontrun{
#' #  The point of pulling chronology tables is to re-build or examine the chronological 
#' #  information that was used to build the age-depth model for the core.
#' 
#' }
#' @references
#' Neotoma Project Website: http://www.neotomadb.org
#' API Reference:  http://api.neotomadb.org/doc/resources/contacts
#' @keywords Neotoma Palaeoecology API
#' @export
write_agefile <- function(download, chronology = 1, path, 
                          corename, cal.prog = 'Bacon'){

  if (!file.exists(paste0(path, '/Cores'))){
    stop(paste0('Core directory must exist.  ',
                'There is no directory at ', path, '/Cores'))
  }

  if (!'download' %in% class(download) | !(c('chronologies') %in% names(download) | c('chronologies') %in% names(download[[1]]))){
    stop(paste0('write_agefile can only operate on valid download ',
                'objects with valid chronologies'))
  }
  
  if ('download' %in% class(download) & (c('chronologies') %in% names(download) | c('chronologies') %in% names(download[[1]]))){
    
    if('download' %in% class(download[[1]])){download <- download[[1]]}
    
    if(class(download$chronologies[[chronology]]) == 'list'){
      # This is new.  We can push chroncontrols into the download object:
      chron.controls <- download$chronologies[[chronology]]$chroncontrol
    } else{
      # Otherwise, we get them from the download object.
      chron.controls <- get_chroncontrol(download,
                                         verbose = FALSE)
    }
    
    if (nrow(chron.controls$chron.control) < 2){
      stop('Chronology must have more than a single date for proper analysis.')
    }
    
    uncal <- c('Radiocarbon', 'Radiocarbon, reservoir correction', 
               'Radiocarbon, average of two or more dates')
    
    if (!tolower(cal.prog) %in% c('bacon', 'clam')){
      stop('You must define either Bacon or Clam as your model output.')
    }
    
    if (cal.prog == 'Bacon'){
      # If you're outputting to the Bacon family of models:
      chron <- data.frame(labid = paste0(chron.controls$chron.control$control.type, 
                                         "_",
                                         chron.controls$chron.control$chron.control.id),
                          age   = chron.controls$chron.control$age,
                          error = abs(chron.controls$chron.control$age - 
                                        chron.controls$chron.control$age.young),
                          depth = chron.controls$chron.control$depth,
                          cc    = ifelse(chron.controls$chron.control$control.type %in% uncal,
                                         1, 0), stringsAsFactors=FALSE)
      
      chron$labid[regexpr(',', chron$labid)>0] <- gsub(',', replacement='_', chron$labid[regexpr(',', chron$labid)>0])
    }
    
    if (cal.prog == 'Clam'){
      #  If you're using Clam:
      chron <- data.frame(ID        = paste0(chron.controls$chron.control$control.type, 
                                             "_",
                                             chron.controls$chron.control$chron.control.id),
                          C14_age   = chron.controls$chron.control$age,
                          cal_BP    = chron.controls$chron.control$age,
                          error     = abs(chron.controls$chron.control$age - 
                                          chron.controls$chron.control$age.young),
                          offset    = NA,
                          depth     = chron.controls$chron.control$depth,
                          thickness = chron.controls$chron.control$thickness,
                          stringsAsFactors=FALSE)
      
      chron$cal_BP [ chron.controls$chron.control$control.type %in% uncal] <- NA
      chron$C14_age[!chron.controls$chron.control$control.type %in% uncal] <- NA
      
      chron$ID[regexpr(',', chron$labid)>0] <- gsub(',', replacement='_', chron$labid[regexpr(',', chron$labid)>0])
    }

    depths <- download$sample.meta$depth
    
    #  There are a couple checks we need:
    if(any(regexpr('Core top', chron[,1])>0)){
      # Core tops need to have error added:
      if(any(is.na(chron$error[which(regexpr('Core top', chron[,1])>0)]))){
        warning('Core tops have no error.  By default we are setting the error to 2 yrs.')
        chron$error[which(regexpr('Core top', chron[,1])>0) & is.na(chron$error)] <- 2
      }
    }

    if(any(is.na(chron[,1]))){
      chron <- chron[!is.na(chron[,1]),]
      
      warning('Some samples in this chronology have no dates.  These samples are being removed.')
    }
    
    #  Now lets make sure the file path is there:
    if (!corename %in% list.files(paste0(path, '/Cores'))){
      works <- dir.create(path = paste0(path, '/Cores/', corename))
      if (!works) {
        stop(paste0('Could not create the directory.  ',
                    'Check the path, corename and your permissions.'))
      }
    }
    
    utils::write.csv(chron, paste0(path, '/Cores/', corename, '/', corename, '.csv'),
              row.names = FALSE, quote = TRUE)
    utils::write.table(depths, paste0(path, '/Cores/', corename, '/', corename, '_depths.txt'),
              row.names = FALSE, quote = TRUE,col.names=FALSE)
  }
}
