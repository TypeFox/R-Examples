################################################################################
# TODO LIST
# TODO: Option to simulate 1 random profile in replicates..
# TODO: Add size?

################################################################################
# NOTES
# ...

################################################################################
# CHANGE LOG (10 last changes)
# 14.04.2016: Version 1.0.0 released.
# 13.08.2015: Updated details and description.
# 24.02.2014: First version.


#' @title DNA Profile Simulator
#'
#' @description Simulates DNA profiles.
#'
#' @details There are three ways to create DNA profiles:
#' 1) Simulate DNA profiles of the selected kit using allele frequencies from
#' the provided db.
#' 2) Simulate DNA profiles using all available markers and allele frequencies
#' from the provided db.
#' 3) Provide a data.frame with a fixed DNA profile.
#' \code{sim} random profiles or \code{sim} replicates are created.
#' The resulting data.frame can be used as input to \code{\link{simSample}}.
#' NB! Homozygous alleles must be specified two times e.g. 16, 16.
#' 
#' @param data data.frame with columns 'Marker' and 'Allele' (creates fixed profiles).
#' @param kit character string to specify the typing kit (if kit=NULL all markers in db will be used).
#' @param sim integer to specify the number of replicates or random profiles.
#' @param name character string giving the sample base name ('Sim' is appended).
#' @param db data.frame with the allele frequency database (used if data=NULL to create random profiles).
#' @param debug logical TRUE to indicate debug mode.
#' 
#' @return data.frame with columns 'Sample.Name', 'Marker', 'Allele', 'Sim'.
#' 
# @importFrom strvalidator getKit
#' @importFrom utils head tail str
#' 
#' @export
#' 
#' @seealso \code{\link{simSample}}
#' 
#' @examples
#' # Create a data frame with a DNA profile.
#' markers = rep(c("D3S1358","TH01","FGA"), each=2)
#' alleles = c(15,18,6,10,25,25)
#' df <- data.frame(Marker=markers, Allele=alleles)
#' 
#' # Simulate sample
#' res <- simProfile(data=df, sim=10, name="Test")
#' print(res)

simProfile <- function(data=NULL, kit=NULL, sim=1, name=NULL, db=NULL, debug=FALSE) {
  
  # Debug info.
  if(debug){
    print(paste(">>>>>> IN:", match.call()[[1]]))
    print("CALL:")
    print(match.call())
    print("###### PROVIDED ARGUMENTS")
    print("STRUCTURE data:")
    print(str(data))
    print("HEAD data:")
    print(head(data))
    print("TAIL data:")
    print(tail(data))
    print("sim")
    print(sim)
    print("name")
    print(name)
    print("STRUCTURE db:")
    print(str(db))
    print("HEAD db:")
    print(head(db))
    print("TAIL db:")
    print(tail(db))
  }
  
  # CHECK PARAMETERS ##########################################################
  
  if(is.null(data) & is.null(db)){
    stop(paste("Either 'data' or 'db' must be specified!"))
  }

  # Check db.
  if(!is.null(data)){
    
    if(!is.data.frame(data)){
      stop(paste("'data' must be of type data.frame."))
    }

    if(!"Marker" %in% names(data)){
      stop(paste("'data' must have a colum 'Marker'."))
    }
    
    if(!"Allele" %in% names(data)){
      stop(paste("'data' must have a colum 'Allele'."))
    }
    
  }  
  
  # Check db.
  if(!is.null(db)){
    
    if(!is.data.frame(db)){
      stop(paste("'db' must be of type data.frame."))
    }
    
    if(!"Database" %in% names(db)){
      stop(paste("'db' must have a colum 'Database'."))
    }
    
    if(!"N" %in% names(db)){
      stop(paste("'db' must have a colum 'N'."))
    }
    
    if(!"Allele" %in% names(db)){
      stop(paste("'db' must have a colum 'Allele'."))
    }
    
    if(is.factor(db$Allele) %in% names(db)){
      db$Allele <- as.character(db$Allele)
      message(paste("'Allele' converted to character (was factor)."))
    }
    
  }
  
  if(!is.logical(debug)){
    stop(paste("'debug' must be logical."))
  }
  
  # PREPARE ###################################################################

  message("SIMULATE PROFILE")
  
  # MULTIPLY SPECIFIED PROFILE ################################################
  
  # Prepare for manual profile.
  if(!is.null(data)){

    # Get number of rows.
    rows <- nrow(data)
    
    # Replicate profile for each simulation.
    data <- do.call("rbind", replicate(sim, data, simplify = FALSE))
    
    if("Sim" %in% names(data)){
      data$Sim <- NA
      message("The 'Sim' column was overwritten!")
    } else {
      data$Sim <- NA
      message("'Sim' column added.")
    }
    
    # Add a column indicating simulation.
    data$Sim <- rep(seq(1:sim), each=rows)
    
    if("Sample.Name" %in% names(data)){
      data$Sample.Name <- NA
      message("The 'Sample.Name' column was overwritten!")
    } else {
      data$Sample.Name <- NA
      message("'Sample.Name' column added.")
    }
    
    # Add a column with the sample name.
    data$Sample.Name <- paste(name, rep(seq(1:sim), each=rows))
    
  }
  
  # SIMULATE PROFILE PER SIMULATION ###########################################
  
  # Prepare for simulated profile.
  if(!is.null(db)){
   
    # Get first and last marker column.
    first <- grep("Allele", names(db), fixed=TRUE) + 1
    last <- length(names(db))
    
    # Extract alleles.
    alleles <- db$Allele
    
    # Extract marker names.
    marker <- names(db)[first:last]
    
    
    if(is.null(kit)){
      # Use all available markers in db.

      if(debug){
        print("database markers")
        print(marker)
      }
      
      # Remove columns with no frequency data (they are not in the kit).
      markersToRemove <- which(colSums(db[first:last], na.rm=TRUE)==0)
      if(length(markersToRemove) > 0){
        marker <- marker[-markersToRemove]
      }
      
      if(debug){
        print("Keep markers")
        print(marker)
      }
      
      
    } else {
      # Use markers in specified kit (if available in db).
      
      # Get markers in kit.
      kitMarker <- getKit(kit=kit, what="Marker")
      
      if(debug){
        print("kitMarker")
        print(kitMarker)
      }
      
      # Keep only markers which exist in the selected kit.
      markersToKeep <-  which(marker %in% kitMarker)
      if(length(markersToKeep) > 0){
        marker <- marker[markersToKeep]
      }

      if(debug){
        print("Keep markers")
        print(marker)
      }
      
      # Remove columns with no frequency data.
      markersToRemove <- which(colSums(db[first:last][markersToKeep], na.rm=TRUE)==0)
      if(length(markersToRemove) > 0){
        marker <- marker[-markersToRemove]
      }

      if(debug){
        print("After remove")
        print(marker)
      }
      
    }
    
#     if(debug){
#       # Start the clock!
#       ptm <- proc.time()
#     }
#      
#     profileName <- NULL
#     profileMarker <- NULL
#     profileAllele <- NULL
#     profileSim <- NULL
#     
#     # FAST, BUT MARKERS IN A SAMPLE IS NOT KEEPT TOGETHER.
#     # Loop over markers.
#     for(m in seq(along=marker)){
#       
#       # Extract marker column values.
#       loci <- db[, marker[m]]
#       
#       # Check if values.
#       if(!all(is.na(loci))){
#         
#         # Extract observed allele values for current marker.
#         allele <- alleles[!is.na(loci)]
#         
#         # Extract frequencies for current alleles and marker.
#         freq <- loci[!is.na(loci)]
#         
#         # Draw two random alleles.
#         # TODO: NB! THIS DOES NOT CHECK FOR (Y,Y) GENOTYPE!
#         genotype <- sample(x=allele, size=2*sim, replace = TRUE, prob = freq)
#         
#         profileName <- c(profileName, rep(paste(name, seq(1:sim), sep=""), each=2))
#         profileMarker <- c(profileMarker, rep(marker[m],2*sim))
#         profileAllele <- c(profileAllele, genotype)
#         profileName <- c(profileName, rep(seq(1:sim), each=2))
#         
#       }
#       
#     }
# 
#     if(debug){
#       # Stop the clock
#       proc.time() - ptm
#     }

    # Pre-allocate vectors.
    vectorSize <- 2 * length(marker) * sim
    profileName <- rep(NA, vectorSize)
    profileMarker <- rep(NA, vectorSize)
    profileAllele <- rep(NA, vectorSize)
    profileSim <- rep(NA, vectorSize)
    
    # SLOWER, BUT MARKERS IN A SAMPLE ARE KEEPT TOGETHER.
    # Loop over markers.
    i <- 1
    for(s in 1:sim){
  
      # Loop over markers.
      for(m in seq(along=marker)){
        
        # Extract marker column values.
        loci <- db[, marker[m]]
        
        # Check if values.
        if(!all(is.na(loci))){
          
          # Extract observed allele values for current marker.
          allele <- alleles[!is.na(loci)]
          
          # Extract frequencies for current alleles and marker.
          freq <- loci[!is.na(loci)]
          
          repeat{
            # Draw two random alleles.
            genotype <- sample(x=allele, size=2, replace = TRUE, prob = freq)

            # Check for (Y,Y) genotype.
            if(!all(genotype=="Y")){
              break
            }
          }
          
          # Save result in vectors.
          profileMarker[i] <- marker[m]
          profileAllele[i] <- genotype[1]
          profileName[i] <- paste(name, s, sep="")
          profileSim[i] <- s
          
          i <- i + 1
          
          profileMarker[i] <- marker[m]
          profileAllele[i] <- genotype[2]
          profileName[i] <- paste(name, s, sep="")
          profileSim[i] <- s
          
          i <- i + 1

        }
        
      }
      
    }

    # Save in data frame.
    data <- data.frame(Sim=profileSim, Sample.Name=profileName,
                       Marker=profileMarker, Allele=profileAllele,
                       stringsAsFactors=FALSE)

# TODO: Markers must be sorted after generated from db?

    # Sort samples and markers.
#     sample <- unique(res$Sample.Name)
# 
#     for(s in seq(along=sample)){
#       
#       
#       
#     }
# 
#     sortMarker(data=res, kit="SGMPlus", addMissingLevels = FALSE, debug = FALSE)
#   
  
  }
    
  # RETURN ####################################################################
  
  # Debug info.
  if(debug){
    print("RETURN")
    print("STRUCTURE:")
    print(str(data))
    print("HEAD:")
    print(head(data))
    print("TAIL:")
    print(tail(data))
    print(paste("<<<<<< EXIT:", match.call()[[1]]))
  }
  
  # Return result. 
  return(data)
  
}