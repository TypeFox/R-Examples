
#' Import a single RSML file into a List. Work only if the roots in the rsml file have properties associated with them.
#' @param rsml.path    The path to the .rsml file. Mandatory.
#' @keywords rsml
#' @import XML 
#' @export
#' @examples
#' path <- "http://rootsystemml.github.io/images/examples/arabidopsis-simple.rsml"
#' pl.list <- rsmlToList(path) # import the file as a list
#' write.csv(pl.list$processed, "rsml-table.csv") # save it as a csv table
#' 


rsmlToList <- function(rsml.path){
    
  rsml <- xmlToList(xmlParse(rsml.path))
  plants <- rsml$scene
  
  # Get the properties used in this specific rsml file
  properties.list <- rsml$metadata$"property-definitions"
  properties.names <- vector(length = length(properties.list))
  for(i in 1:length(properties.list)){properties.names[i] <- properties.list[i]$"property-definition"$label}
  
  # Create the empty data frame
  rsml.summary <- data.frame(file_id = rsml.path,
                             plant_id = "-1",
                             root_id = "-1",
                             parent_id = "-1",
                             root_order = -1)
  for(p in properties.names) rsml.summary[[p]] <- -1
  rsml.summary$root_id <- as.character(rsml.summary$root_id)
  rsml.summary$plant_id <- as.character(rsml.summary$plant_id)
  rsml.summary$parent_id <- as.character(rsml.summary$parent_id)
  
  # Create the same data frame for temp storage
  template <- rsml.summary[1,]
  
  # Remove the data from the dataframe to have a clean start
  rsml.summary <- rsml.summary[!1,]
  
  # Get the data from the rsml structure
  k <- 1
  for(r0 in plants){
    for(r1 in r0){
      if("properties" %in% names(r1)){
        # Send data
        temp <- template
        temp[["plant_id"]][1] <- k
        if(!is.null(r0$.attr[1])) temp[["plant_id"]][1] <- r0$.attr[1]
        temp[["root_id"]][1] <- r1$.attr[1]
        temp[["root_order"]][1] <- 0
        for(p in properties.names){
          if(!is.null(r1$properties[[p]])){ temp[[p]][1] <- r1$properties[[p]]}
        }
        rsml.summary <- rbind(rsml.summary, temp)
        
        if("root" %in% names(r1)){
          for(r2 in r1){  
            if("properties" %in% names(r2)){
              
              temp <- template
              temp[["plant_id"]][1] <- k
              if(!is.null(r0$.attr[1])) temp[["plant_id"]][1] <- r0$.attr[1]
              temp[["root_id"]][1] <- r2$.attr[1]
              temp[["parent_id"]][1] <- r1$.attr[1]
              temp[["root_order"]][1] <- 1
              for(p in properties.names){
                if(!is.null(r2$properties[[p]])){ temp[[p]][1] <- r2$properties[[p]]}
              }
              rsml.summary <- rbind(rsml.summary, temp)
              
              if("root" %in% names(r2)){
                for(r3 in r2){
                  if("properties" %in% names(r3)){
                    
                    temp <- template
                    temp[["plant_id"]][1] <- k
                    if(!is.null(r0$.attr[1])) temp[["plant_id"]][1] <- r0$.attr[1]
                    temp[["root_id"]][1] <- r3$.attr[1]
                    temp[["parent_id"]][1] <- r2$.attr[1]
                    temp[["root_order"]][1] <- 2
                    for(p in properties.names){
                      if(!is.null(r3$properties[[p]])){ temp[[p]][1] <- r3$properties[[p]] }
                    }
                    rsml.summary <- rbind(rsml.summary, temp)
                    
                  }
                }
              }
            }
          }
        }
      }
    }
    k <- k+1
  }
  for(p in properties.names) rsml.summary[[p]] <- as.numeric(rsml.summary[[p]])
  
  rsml.data <- list()
  rsml.data$processed <- rsml.summary
  
  name <- tail(strsplit(rsml.path, "/")[[1]], n=1)
  rsml.data$files.list <- name  
  
  # Return the extended dataframe 
  return(rsml.data)
}
