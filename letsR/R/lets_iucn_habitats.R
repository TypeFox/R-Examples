#' Download species' habitat information from the IUCN RedList online database
#' 
#' @author Bruno Vilela
#' 
#' @description Get species' habitat information from the IUCN RedList 
#' website(\url{http://www.iucnredlist.org/}) for one or more species.
#' 
#' @param input Character vector with one or more species names,
#' or an object of the \code{\link{PresenceAbsence}} class.
#' @param count Logical, if \code{TRUE} a counting window will open.
#' 
#' @return A data frame with species names in the first column and the 
#' habitats where it occurs in the remaining columns,
#' '1' if species is present in that habitat and '0' otherwise.
#'
#' @details Note that you must be connected to the internet to use this function. 
#' 
#' @import XML
#' 
#' @seealso \code{\link{lets.iucn}}
#' @seealso \code{\link{lets.iucn.his}} 
#' 
#' @examples \dontrun{
#' # Single species
#' lets.iucn.ha("Pongo pygmaeus")
#' 
#' # Multiple species
#' sp <- c("Musonycteris harrisoni", "Ailuropoda melanoleuca",
#'         "Cebus flavius")
#' lets.iucn.ha(sp)
#' }
#' 
#' @export

lets.iucn.ha <- function(input, count = FALSE) {
  #keep species name(s)
  
  sps <- .getnames(input) 
  
  
  #Habitat names (and the name "Species" that will be used in the matrix columns names)
  names <- c("Species", "Forest", "Savanna", "Shrubland", "Grassland", 
             "Wetlands", "Rocky areas", "Caves and Subterranean Habitats", 
             "Desert", "Marine Neritic", "Marine Oceanic", "Marine Deep Ocean Floor", 
             "Marine Intertidal", "Marine Coastal/Supratidal", "Artificial/Terrestrial", 
             "Artificial/Aquatic", "Introduced Vegetation", "Other", "Unknown")
  
  #create an empty matrix
  habitat <- matrix(0, nrow = length(input), ncol = length(names)) 
  
  #Adding the column names
  colnames(habitat) <- names
  
  n <- length(input)
  
  # With count window
  if (count) {
    
    # Do not set a new device in rstudio to avoid warnings()
    if (!"tools:rstudio" %in% search()) {
      dev.new(width = 2, height = 2, pointsize = 12)
      par(mar = c(0, 0, 0, 0))
    }
    
    for(i in 1:n){
      plot.new()
      text(0.5, 0.5, paste(paste("Total:", n, "\n",
                                 "Species to go: ",
                                 (n - i))))
      
      ncolumns <- ncol(habitat)
      habitat[i, 2:(ncolumns - 1)] <- .Habitat(input, i, 
                                               habitat,
                                               names)
      
      # If none of the habitat names has been found or if 
      # the species has not been found in IUCN archives, 
      # it will have 1 in the column Unknwon
      if (sum(habitat[i, ]) == 0) {
        habitat[i, ncolumns] <- 1
      }
    }
    dev.off()
  }
  
  if (!count) {
    for(i in 1:n){
      
      ncolumns <- ncol(habitat)
      habitat[i, 2:(ncolumns - 1)] <- .Habitat(input, i, 
                                               habitat,
                                               names)
      
      # If none of the habitat names has been found or if 
      # the species has not been found in IUCN archives, 
      # it will have 1 in the column Unknwon
      if (sum(habitat[i, ]) == 0) {
        habitat[i, ncolumns] <- 1
      }
    }
  }
  
  #Putting species' names in the first column
  habitat[, 1] <- sps
  
  #Return the resulting matrix
  return(as.data.frame(habitat))
}


# Auxiliar funciton inside the loop
.Habitat <- function(input, i, habitat, names) {
  # Taking the Website code from the internet
  
  c <- .getcode(input[i])
  httpclas <- "http://www.iucnredlist.org/details/classify/"
  h2 <- try(htmlParse(paste(httpclas, c, "/0", sep = "")),
            silent = TRUE)
  
  # Taking the specific parts that contain the habitat names
  b2 <- try(xpathSApply(h2, '//html', xmlValue), silent = TRUE)
  
  # Look for the habitat names inside the string 
  # (if the sting contains the name, it will be marked 1)
  Nnames <- length(names)
  habitatParcial <- numeric(Nnames - 2)
  
  for(t in 2:(Nnames - 1)) {
    if (sum(grep(names[t], b2)) > 0) {
      habitatParcial[(t - 1)] <- 1
    }
  }
  return(habitatParcial)
}


