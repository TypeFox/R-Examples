#' Download species' information from the IUCN RedList online database
#' 
#' @author Bruno Vilela
#' 
#' @description Get species' information from the IUCN
#' website(\url{http://www.iucnredlist.org/}) for one or more species.
#' 
#' 
#' @param input Character vector with one or more species names,
#' or an object of class \code{\link{PresenceAbsence}}.
#' @param count Logical, if \code{TRUE} a counting window will open.
#' 
#' @return Returns a data frame with the Species Name, Family, Conservation Status, 
#' Criteria used to estabilish the conservation status, Population Status, Year of 
#' Description (only for animals), and the Countries where it occurs. If species do
#' not have information (i.e. have not been evaluated), the result is: NE (Not
#' evaluated).
#' 
#' @details Note that you must be connected to the internet to use this function.
#' 
#' @import XML
#' 
#' @seealso \code{\link{lets.iucn.ha}}
#' @seealso \code{\link{lets.iucn.his}}
#' 
#' @examples \dontrun{
#' # Single species
#' lets.iucn("Pongo pygmaeus")
#' 
#' # Multiple species
#' sp <- c("Musonycteris harrisoni", "Ailuropoda melanoleuca",
#'         "Cebus flavius")
#' lets.iucn(sp)
#' }
#' 
#' @export

lets.iucn <- function(input, count = FALSE) {
  
  input <- .getnames(input)
  
  # Matrices to save the result
  ln <- length(input)
  matriz1 <- matrix(nrow = ln, ncol = 1)
  status <-  matriz1
  criterio <- matriz1
  populacao <- matriz1
  familia <- matriz1
  autor <- matriz1  
  pais <- matriz1
  
  # With count window
  if (count) {
    
    # Do not set a new device in rstudio to avoid warnings()
    if (!"tools:rstudio" %in% search()) {
      dev.new(width = 2, height = 2, pointsize = 12)
      par(mar = c(0, 0, 0, 0))
    }
    
    for (i in 1:ln){
      plot.new()
      text(0.5, 0.5, paste(paste("Total:", ln, "\n",
                                 "Species to go: ",
                                 (ln - i))))
      
      loopresu        <- .LoopIUCN(input[i])
      status[i, ]    <- loopresu$Status
      criterio[i, ]  <- loopresu$Criteria
      populacao[i, ] <- loopresu$Population
      familia[i, ]   <- loopresu$Family
      autor[i, ]     <- loopresu$Author
      pais[i, ]       <- loopresu$Country
    }
    dev.off()
  }
  
  if(!count){
    for (i in 1:ln){          
      
      loopresu        <- .LoopIUCN(input[i])
      status[i, ]    <- loopresu$Status
      criterio[i, ]  <- loopresu$Criteria
      populacao[i, ] <- loopresu$Population
      familia[i, ]   <- loopresu$Family
      autor[i, ]     <- loopresu$Author
      pais[i, ]       <- loopresu$Country
    }
  }
  
  # Final table
  resu <- cbind(input, familia, status, criterio, populacao, autor, pais)
  colnames(resu) <- c("Species", "Family", "Status", "Criteria",
                      "Population", "Description_Year", "Country")
  
  # Making a data frame and guarantee that it is numeric the Year
  resu <- as.data.frame(resu)
  resu[, 6] <- as.numeric(levels(resu[, 6]))[resu[, 6]]
  return(resu)
}





# Inside the loop function -----------------------------

.LoopIUCN <- function(input) {
  
  # Get species code
  speciescode <- .getcode(input)
  # Find the web
  h <- try(htmlParse(paste("http://api.iucnredlist.org/details/",
                           speciescode, "/0", sep = "")),
           silent=TRUE)
  # If not find, species set to NE
  if (class(h)[1] == "try-error") {
    statusi <- "NE"
    criterioi <- ""
    populacaoi <- "Unknown"
    familiai <- ""
    autori <- ""
    paisi <- ""
  } else {
    # Look for the informations
    statusi <- try(xpathSApply(h, '//div[@id="red_list_category_code"]', 
                               xmlValue),
                   silent = TRUE)
    criterioi <- try(xpathSApply(h, '//div[@id="red_list_criteria"]',
                                 xmlValue),
                     silent = TRUE)
    pop <- try(xpathSApply(h, '//div[@id="population_trend"]',
                           xmlValue),
               silent=TRUE)
    
    
    familiai <- try(xpathSApply(h, '//div[@id="family"]',
                                xmlValue),
                    silent = TRUE)
    
    autori <- try(xpathSApply(h, '//div[@id="species_authority"]',
                              xmlValue),
                  silent = TRUE)
    
    # Error control
    if (class(statusi)[1] == "try-error") {
      statusi <-  "NE"
    } else {
      statusi <- statusi
    }
    
    if (class(criterioi)[1] == "try-error") {
      criterioi <-  ""
    } else {
      criterioi <- criterioi
    }
    
    # Sometimes the pop is an empty list.
    if (class(pop)[1] == "try-error" | class(pop) == "list") {
      populacaoi <-  "Unknown"
    } else {
      populacaoi <- pop
    }
    
    if (class(familiai)[1] == "try-error") {
      familiai <-  ""
    } else {
      familiai <- familiai
    }
    
    if (class(autori)[1] == "try-error") {
      autori <-  ""
    } else {
      autori <- gsub("\\D", "", autori)
      autori <- as.numeric(substr(autori, 1, 4))
      if (is.na(autori)) {
        autori <- ""
      }
    }
    
    # Countries
    distr1 <- try(xpathSApply(h, '//ul[@class="countries"]', 
                              xmlValue),
                  silent = TRUE)
    # Error control
    if (class(distr1)[1] == "try-error") {
      paisi <- ""
    } else {
      # Remove countries with comma
      distr2 <- try(unlist(strsplit(distr1, "\n")), silent = TRUE)
      distr2 <- gsub(',', '', distr2)
      paisi <- paste(distr2, collapse = ", ")
    }
  }
  # Return the information
  return(list("Status" = statusi,
              "Criteria" = criterioi,
              "Population" = populacaoi,
              "Family" = familiai,
              "Author" = autori,
              "Country" = paisi))
}
