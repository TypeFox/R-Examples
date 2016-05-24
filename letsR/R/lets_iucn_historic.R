#' Download species' temporal trend in conservation status from the IUCN RedList online database
#' 
#' @author Bruno Vilela
#' 
#' @description Get species conservation status over time (i.e. from 1980 to the present date available)
#' from the IUCN website(\url{http://www.iucnredlist.org/}) for one or more species.
#' 
#' @param input character vector with one or more species names,
#' or an object of class \code{\link{PresenceAbsence}}.
#' @param count Logical, if \code{TRUE} a counting window will open.

#' @return A data frame with the species names in the first column rows and the years (1980 - present) in
#' the remaining columns, the code represents the species' conservation status (see the IUCN RedList 
#' website for details). If species do not have information (i.e. have not been evaluated), the result
#' is: NE (Not evaluated).
#' 
#' @return Codes and categories:
#' @return \strong{EX}: Extinct
#' @return \strong{EW}: Extinct in the Wild
#' @return \strong{VU}: Vulnerable
#' @return \strong{EN}: Endangered
#' @return \strong{CR}: Critically Endangered
#' @return \strong{LC}: Least Concern
#' @return \strong{NT}: Near Threatened
#' @return \strong{DD}: Data Deficient
#' @return \strong{CT}: Commercially Threatened
#' @return \strong{IN}: Indeterminate
#' @return \strong{IK}: Insufficiently Known
#' @return \strong{LR}: Lower Risk
#' @return \strong{RA}: Rare
#'
#'
#' @details Note that you must be connected to the internet to use this function. 
#'
#' @import XML
#' 
#' @seealso \code{\link{lets.iucn.ha}}
#' @seealso \code{\link{lets.iucn}}
#' 
#' @examples \dontrun{
#' # Single species
#' lets.iucn.his("Panthera onca")
#' 
#' # Multiple species
#' sp <- c("Rhincodon typus", "Ailuropoda melanoleuca")
#' lets.iucn.his(sp)
#' }
#' 
#' @export


lets.iucn.his <- function(input, count = FALSE) {  
  
  input <- .getnames(input)
  n <- length(input)
  
  # Automate date
  data <- date()
  anofinal <- substr(data, (nchar(data) - 3), nchar(data))
  anos <- 1980:anofinal
  
  #Empty matrix
  resus <- matrix(ncol = length(anos), nrow = n)
  colnames(resus) <- anos
  # With count window
  if (count) {
    
    # Do not set a new device in rstudio to avoid warnings()
    if (!"tools:rstudio" %in% search()) {
      dev.new(width = 2, height = 2, pointsize = 12)
      par(mar = c(0, 0, 0, 0))
    }
    
    
    for(i in 1:n) {
      plot.new()
      text(0.5, 0.5, paste(paste("Total:", n, "\n",
                                 "Species to go: ",
                                 (n - i))))
      resus[i, ] <- .Hist(input[i], anos)
    }   
    
    dev.off()
  }
  
  if (!count) {
    
    for(i in 1:n){              
      resu <- .Hist(input[i], anos)
      resus[i, ] <- .Hist(input[i], anos)
    }   
    
  }
  
  Species <- gsub(as.matrix(input), pattern = "-", 
                  replacement = " ")
  final <- data.frame(Species, resus)
  colnames(final) <- c("Species", anos)
  return(final)
}



#-------------------------------
# Automate the current year

.Hist <- function(input, anos) {
  
  lanos <- length(anos)
  matriz <- matrix("NE", ncol = lanos)
  colnames(matriz) <- anos
  
  
  c <- .getcode(input)
  if (is.null(c)) {
    return(matriz)
  } else {
    h2 <- htmlParse(paste("http://api.iucnredlist.org/details/", 
                          c, "/0", sep = ""))
    ano1 <- xpathSApply(h2, '//div[@id="modified_year"]', xmlValue)
    ameaca1 <- xpathSApply(h2, '//div[@id="red_list_category_code"]', xmlValue)
    
    h <- try(htmlParse(paste("http://www.iucnredlist.org/details/full/", 
                             c, "/0", sep = "")), 
             silent = TRUE)
    
    a <- try(xpathSApply(h, '//td[table]', xmlValue),
             silent = TRUE)
    a <- a[2]
    
    if (is.na(a)) {
      matriz[anos %in% ano1:anos[lanos]] <- ameaca1
      return(matriz)
    } else {
      
      a <- gsub("\n", "", a)
      a <- gsub("\t", "", a)
      b <- strsplit(a, "          ")[[1]]
      b <- strsplit(b, "      ")
      c <- do.call("rbind", b)
      c <- matrix(c, ncol = 1)
      dupc <- duplicated(c)
      if (any(dupc)) {
        c <- c[!dupc, , drop = FALSE]
      }
      ano <- substr(gsub("\\D", "", c), 1, 4)
      remano <- ano != ""
      ano <- ano[remano, , drop = FALSE]
      c <- c[remano, , drop = FALSE]
      if (nrow(ano) == 0) {
        return(matriz)
      } else {
        
        if (length(ano) >= 1) {
          d <- gsub("[0-9]", "", c)
          d <- gsub("[[:punct:]]", "", d)
          d2 <- gsub("\\W", "", d)
          EX <- grep("Extinct", d2)
          EW <- grep("ExtinctintheWild", d2)
          VU <- grep("Vulnerable", d2)
          EN <- grep("Endangered", d2)
          CR <- grep("CriticallyEndangered", d2)
          LC <- grep("LeastConcern", d2)
          NT <- grep("NearThreatened", d2)
          DD <- grep("DataDeficient", d2)
          CT <- grep("CommerciallyThreatened", d2)
          IN <- grep("Indeterminate", d2)
          IK <- grep("InsufficientlyKnown", d2)
          LR <- grep("LowerRisk", d2)
          RA <- grep("Rare", d2)
          RA2 <- grep("rare", d2)
          ameaca <- numeric(length(ano))
          
          ameaca[EX] <- "EX"
          ameaca[EW] <- "EW"
          ameaca[VU] <- "VU"
          ameaca[EN] <- "EN"
          ameaca[CR] <- "CR"
          ameaca[LC] <- "LC"
          ameaca[NT] <- "NT"
          ameaca[IK] <- "IK"
          ameaca[DD] <- "DD"
          ameaca[IN] <- "IN"
          ameaca[RA] <- "RA"
          ameaca[RA2] <- "RA"
          ameaca[CT] <- "CT"
          ameaca[LR] <- "LR"
          
          ameaca <- ameaca[!(duplicated(ano))]
          ano <- ano[!(duplicated(ano))]
          ano <- as.numeric(ano)
          ameaca[ameaca == "0"] <- d2[ameaca == "0"]
          
          for(i in 1:length(ano)) {      
            matriz[, anos %in% ano[i]] <- ameaca[i]  
          }
          
          ameaca <- c(ameaca, ameaca1)
          ano <- c(ano, ano1)
          ameaca <- ameaca[!(duplicated(ano))]
          ano <- ano[!(duplicated(ano))]
          ano <- ano[ano %in% anos]
          ameaca <- ameaca[ano %in% anos]
          pos <- which(anos %in% ano)
          pos2 <- sort(ano, index.return = TRUE)$ix
          ameaca <- ameaca[pos2]
          
          for(i in 1:(length(ameaca) - 1)) {
            subseq <- seq(from = (pos[i] + 1), (pos[i + 1] - 1))
            matriz[, subseq] <- ameaca[i]
          }
        }
        if (ano1 %in% anos) {
          pos3 <- which(anos %in% ano1)
          matriz[, pos3:ncol(matriz)] <- ameaca1
        }
        return(matriz)
      }
    }
  }
}
