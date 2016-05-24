
#' Import a single RSML file into a Plant object
#' @param rsml.path    The path to the .rsml file. Mandatory.
#' @param threed  Does the RSML file contains a 3D root system?
#' @keywords rsml
#' @keywords threed is the plant in 3d
#' @import XML 
#' @export
#' @examples
#' # 2D example
#' path <- "http://rootsystemml.github.io/images/examples/arabidopsis-simple.rsml"
#' pl <- rsmlToPlant(path) # import the file as a 2D plant
#' plot(pl) # plot the plant
#' 
#' # 3D example
#' path <- "http://rootsystemml.github.io/images/examples/anagallis.rsml"
#' pl <- rsmlToPlant(path, threed=TRUE) # import the file as a 2D plant
#' plot(pl, threed=TRUE) # plot the plant

rsmlToPlant <- function(rsml.path, threed = FALSE){
  
  # Create the plant object
  pl <- plant()
  rsml <- xmlToList(xmlParse(rsml.path))
  
  # Get the scale in pixels / cm
  meta <- rsml$meta
  scale <- 1
  if(!is.null(meta$unit) & !is.null(meta$resolution)){
    if(grepl("inch", meta$unit)) scale <- as.numeric(meta$resolution) / 2.54 
    if(grepl("cm", meta$unit)) scale <- as.numeric(meta$resolution)
    if(grepl("mm", meta$unit)) scale <- as.numeric(meta$resolution) * 10
  }
       
  plants <- rsml$scene
  inc  = 1
  for(r0 in plants){ # Each plant
    for(r1 in r0){ # Each first order root
      if(class(r1) == "list"){
        # Create an empty root
        id1 <- tryCatch({r1$.attrs[["ID"]]},
          error=function(cond) {r1$.attrs[["id"]]},
          finally={}
        )   
        r <- root(id = id1)
        # Add the nodes
        ns <- r1$geometry$polyline
        cumulDist <- 0
        for(i in 1:length(ns)){
          # Get the distance of the node from the base of the root
          if(i == 1){ dist <- 0
          }else{
            dx <- as.numeric(ns[[i]][[1]]) - as.numeric(ns[[i-1]][[1]]) # x
            dy <- as.numeric(ns[[i]][[2]]) - as.numeric(ns[[i-1]][[2]]) # y
            if(threed){
              dz <- as.numeric(ns[[i]][[3]]) - as.numeric(ns[[i-1]][[3]]) # z
              dist = sqrt(dx^2 + dy^2 + dz^2) / scale
            }else{
              dist = sqrt(dx^2 + dy^2) / scale
            }
          }
          cumulDist = cumulDist + dist
          
          # Create the node
          if(threed) {
            r <- addNodeToRoot(r, node(as.numeric(ns[[i]][[1]]) / scale, as.numeric(ns[[i]][[2]]) / scale, as.numeric(ns[[i]][[3]]) / scale, bLength=cumulDist))
          } else { 
            r <- addNodeToRoot(r, node(as.numeric(ns[[i]][[1]]) / scale, as.numeric(ns[[i]][[2]]) / scale, bLength=cumulDist))          
          }
        }  
        # Add the functions
        if("functions" %in% names(r1)){
          for(f in r1$functions){
            # Diameter
            if(grepl("diam", f$.attr[1])){
              for(i in 1:(length(f)-1)){
                r$nodes[[i]]$diameter = as.numeric(f[[i]]) / scale
              }
            }
            # Orientation
            if(grepl("orien", f$.attr[1])){
              for(i in 1:(length(f)-1)){
                r$nodes[[i]]$orientation = f[[i]]
              }
            }          
          }
        }
        # If their is lateral roots
        if("root" %in% names(r1)){
          for(r2 in r1){  
            if("geometry" %in% names(r2)){
              # Create an empty root
              id1 <- tryCatch({r2$.attrs[["ID"]]},
                              error=function(cond) {r2$.attrs[["id"]]},
                              finally={}
              )   
              rr <- root(id = id1, parent=r$id)
              # Add the nodes
              ns <- r2$geometry$polyline
              cumulDist <- 0
              for(i in 1:length(ns)){
                # Get the distance of the node from the base of the root
                if(i == 1){ dist <- 0
                }else{
                  dx <- as.numeric(ns[[i]][[1]]) - as.numeric(ns[[i-1]][[1]]) # x
                  dy <- as.numeric(ns[[i]][[2]]) - as.numeric(ns[[i-1]][[2]]) # y
                  if(threed){
                    dz <- as.numeric(ns[[i]][[3]]) - as.numeric(ns[[i-1]][[3]]) # z
                    dist = sqrt(dx^2 + dy^2 + dz^2) / scale
                  }else{
                    dist = sqrt(dx^2 + dy^2) / scale
                  }
                }
                cumulDist = cumulDist + dist            
                # Create the node
                if(threed) {
                  rr <- addNodeToRoot(rr, node(as.numeric(ns[[i]][[1]]) / scale, as.numeric(ns[[i]][[2]]) / scale, as.numeric(ns[[i]][[3]]) / scale, bLength=cumulDist))
                } else { 
                  rr <- addNodeToRoot(rr, node(as.numeric(ns[[i]][[1]]) / scale, as.numeric(ns[[i]][[2]]) / scale, bLength=cumulDist))          
                }
              }
              if("functions" %in% names(r2)){
                for(f in r2$functions){
                  if(grepl("diam", f$.attr[1])){
                    for(i in 1:(length(f)-1)){
                      rr$nodes[[i]]$diameter = as.numeric(f[[i]]) / scale
                    }
                  }
                  # Get the orientation and the insertion angel
                  if(grepl("orien", f$.attr[1])){
                    for(i in 1:(length(f)-1)){
                      rr$nodes[[i]]$orientation = f[[i]]
                    }
                    rr$insertion_angle <- getInsertionAngle(r, rr)                              
                  }          
                }
              }
              
              # If their is third order roots
              if("root" %in% names(r2)){
                for(r3 in r2){  
                  if("geometry" %in% names(r3)){
                    # Create an empty root
                    id1 <- tryCatch({r3$.attrs[["ID"]]},
                                    error=function(cond) {r3$.attrs[["id"]]},
                                    finally={}
                    )   
                    rrr <- root(id =id1, parent=rr$id)
                    # Add the nodes
                    ns <- r3$geometry$polyline
                    cumulDist <- 0
                    for(i in 1:length(ns)){
                      # Get the distance of the node from the base of the root
                      if(i == 1){ dist <- 0
                      }else{
                        dx <- as.numeric(ns[[i]][[1]]) - as.numeric(ns[[i-1]][[1]]) # x
                        dy <- as.numeric(ns[[i]][[2]]) - as.numeric(ns[[i-1]][[2]]) # y
                        if(threed){
                          dz <- as.numeric(ns[[i]][[3]]) - as.numeric(ns[[i-1]][[3]]) # z
                          dist = sqrt(dx^2 + dy^2 + dz^2) / scale
                        }else{
                          dist = sqrt(dx^2 + dy^2) / scale
                        }
                      }
                      cumulDist = cumulDist + dist                    
                      # Create the node
                      if(threed) {
                        rrr <- addNodeToRoot(rrr, node(as.numeric(ns[[i]][[1]]) / scale, as.numeric(ns[[i]][[2]]) / scale, as.numeric(ns[[i]][[3]]) / scale, bLength=cumulDist))
                      } else { 
                        rrr <- addNodeToRoot(rrr, node(as.numeric(ns[[i]][[1]]) / scale, as.numeric(ns[[i]][[2]]) / scale, bLength=cumulDist))          
                      }                  
                    }
                    
                    # Get the diameter an orientation, if
                    if("functions" %in% names(r3)){
                      for(f in r3$functions){
                        if(grepl("diam", f$.attr[1])){
                          for(i in 1:(length(f)-1)){
                            rrr$nodes[[i]]$diameter = as.numeric(f[[i]]) / scale
                          }
                        }
                        # Get the orientation and the orientation ange
                        if(grepl("orien", f$.attr[1])){
                          for(i in 1:(length(f)-1)){
                            rrr$nodes[[i]]$orientation = f[[i]]
                          }
                          rrr$insertion_angle <- getInsertionAngle(rr, rrr)                                    
                        }          
                      }
                    }
                    # Get the insertion position
                    rrr$insertion <- getInsertionPosition(rr, rrr)
                    rr <- addChildToRoot(rr, rrr)
                  }
                }
              }
              # Get the insertion position
              rr$insertion <- getInsertionPosition(r, rr)
              r <- addChildToRoot(r, rr)
            }
          }
        }
        pl <- addRootToPlant(pl, r)
      }
    }
  }
  pl
}
