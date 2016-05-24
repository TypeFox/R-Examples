#' @include getLineup.R
#' @include simRedCard.R
NULL

###############################################
# --------------------------------------------#
# Class formation                             #
# --------------------------------------------#
###############################################

# --------------------------------------------
# Function for validity check
# --------------------------------------------

# Validity check function for objects of the formation class
# 
# @param object \code{S4} object of the \code{formation} class
#
# @return Returns a \code{TRUE}, if the settings of the object are valid.
validateFormation <- function(object) {
  errors <- character()
  GK <- object@GK
  SW <- object@SW
  DF <- object@DF
  MF <- object@MF
  ST <- object@ST
  hard <- object@hardness
  homeAdv <- object@homeAdv
  
  if(!(any(is.na(SW)) || is.numeric(SW))) {
    stop("SW must be a numeric or NA.")
  }
  
  for (param in c("GK", "SW")) {
    obj <- eval(parse(text = param))
    if (!(length(obj) == 1)) {
      msg <- paste(param, " should have length 1. Has length ", length(obj),
                   ".", sep = "", collapse = ",")
      errors <- c(errors, msg)
    }
  }
  
  if (any(is.na(SW))) {
    for (param in c("GK", "DF", "MF", "ST", "hard", "homeAdv")) {
      obj <- eval(parse(text = param))
      if (!(all(round(obj) == obj))) {
        msg <- paste("All arguments of ", param, " should be integers.",
                     sep = "", collapse = ",")
        errors <- c(errors, msg)
      }
      
      if (!(all(obj >= 0))) {
        msg <- paste("All arguments of ", param, " should be greater than zero.",
                     sep = "", collapse = ",")
        errors <- c(errors, msg)
      }
      
    }  
    
    if (GK[1] > 13) {
      msg <- paste("First argument of GK should be smaller than 14. GK has strength ", 
                   GK[1], ".", sep = "", collapse = ",")
      errors <- c(errors, msg)
    }
    
  } else {
      for (param in c("GK", "SW", "DF", "MF", "ST", "hard", "homeAdv")) {
        obj <- eval(parse(text = param))
        if (!(all(round(obj) == obj))) {
          msg <- paste("All arguments of ", param, " should be integers.",
                       sep = "", collapse = ",")
          errors <- c(errors, msg)
        }
        
        if (!(all(obj >= 0))) {
          msg <- paste("All arguments of ", param, " should be greater than zero.",
                       sep = "", collapse = ",")
          errors <- c(errors, msg)
        }
      }  
        
      for (param in c("GK", "SW")) {
        obj <- eval(parse(text = param))
        if (obj[1] > 13) {
          msg <- paste("First argument of ", param, " should be smaller than 14.", 
                        param," has strength ", obj[1], ".", sep = "", collapse = ",")
          errors <- c(errors, msg)
        }
      }  
  }
    
  if (!(length(hard) == 5)) {
    msg <- paste("Hardness should have length 5. Has length ", length(hard),
                 ".", sep = "", collapse = ",")
    errors <- c(errors, msg)
  }
  
  if (sum(hard) > 10) {
    msg <- paste("Sum of hardness should be smaller than 11. Sum of hardness is ", 
                 sum(hard), ".", sep = "", collapse = ",")
    errors <- c(errors, msg)
  }

  if (!(length(homeAdv) == 5)) {
    msg <- paste("HomeAdv should have length 5. Has length ", length(homeAdv),
                 ".", sep = "", collapse = ",")
    errors <- c(errors, msg)
  }
  
  if (is.numeric(SW) && hard[2] %% 2 != 0 && length(hard) == 5) {
    msg <- "Hardness given to the sweeper must be divisible by two."
    errors <- c(errors, msg)
  }
  
  if (is.numeric(SW) && homeAdv[2] %% 2 != 0 && length(homeAdv) == 5) {
    msg <- "Home advantage given to the sweeper must be divisible by two."
    errors <- c(errors, msg)
  }
  
  if (hard[1] %% 2 != 0 && length(hard) == 5) {
    msg <- "Hardness given to the goalkeeper must be divisible by two."
    errors <- c(errors, msg)
  }
  
  if (homeAdv[1] %% 2 != 0 && length(homeAdv) == 5) {
    msg <- "Home advantage given to the goalkeeper must be divisible by two."
    errors <- c(errors, msg)
  }
  
  if (any(is.na(SW)) && hard[2] != 0 && length(hard) == 5) {
    msg <- "Hardness given to the sweeper, which is not set up."
    errors <- c(errors, msg)
  }
  
  if (any(is.na(SW)) && homeAdv[2] != 0 && length(homeAdv) == 5) {
    msg <- "Home advantage given to the sweeper, which is not set up."
    errors <- c(errors, msg)
  }
  
  if (homeAdv[1] != 0 && length(homeAdv) == 5) {
    warning("Home advantage given to the goalkeeper (not possible according to the rules).")
  }
  
  if (homeAdv[2] != 0 && length(homeAdv) == 5) {
    warning("Home advantage given to the sweeper (not possible according to the rules).")
  }
  
  if (length(GK) + length(SW) + length(DF) + length(MF) + length(ST) - is.na(SW[1]) > 11) {
    msg <- "Number of players should be 11 or less.  Number of players are too much."
    errors <- c(errors, msg)
  }
  
  if (sum(hard) > 1 && sum(hard) < 11) {
    for (param in c("DF", "MF", "ST")) {
      obj <- eval(parse(text = param))
      if (!(all(round(obj) <= 13))) {
        msg <- paste("All players in ", param, " should be weaker than 14.",
                     sep = "", collapse = ",")
        errors <- c(errors, msg)
      }
    }
  }
  
  if (length(errors) == 0) TRUE else errors
}

# --------------------------------------------
# Class definition for formation
# --------------------------------------------

# formation class
setClass("formation",
         slots = c(GK = "numeric", SW = "vector", DF = "numeric", 
                   MF = "numeric", ST = "numeric",
                   hardness = "numeric", homeAdv = "numeric"),
         validity = validateFormation)

# --------------------------------------------
# Constructor function for formation
# --------------------------------------------

#' Representing a formation
#' 
#' Represents a valid united formation.
#' 
#' @inheritParams overview
#' 
#' @return 
#' \code{S4} object of the class \code{formation}.
#'
#' @export
formation <- function(GK, SW, DF, MF, ST, hardness = c(0, 0, 0, 0, 0), 
                      homeAdv = c(0, 0, 0, 0, 0)) {
  new("formation", GK = GK, SW = SW, DF = DF, MF = MF, ST = ST, 
      hardness = hardness, homeAdv = homeAdv)
}


# --------------------------------------------
# Methods for formation
# --------------------------------------------

#' @rdname getLineup
setMethod("getLineup", 
          signature(obj = "formation"),
          function(obj) {
            GK <- obj@GK
            SW <- ifelse(is.na(obj@SW), 0, obj@SW)
            DF <- sum(obj@DF)
            MF <- sum(obj@MF)
            ST <- sum(obj@ST)
            hard <- obj@hardness
            homeAdv <- obj@homeAdv
            #3:1 rule before adding hardness and home advantage
            minStrength <- min(DF, MF, ST)[1]
            if (3*minStrength < DF) {
              DF <- 3*minStrength
              warning("DF was reduced according to the 3:1 rule (before adding hardness and home advantage).")
            }
            if (3*minStrength < MF) {
              MF <- 3*minStrength
              warning("MF was reduced according to the 3:1 rule (before adding hardness and home advantage).")
            }
            if (3*minStrength < ST) {
              ST <- 3*minStrength
              warning("ST was reduced according to the 3:1 rule (before adding hardness and home advantage).")
            }
            
            #3:1 rule after adding hardness and home advantage
            GKbefore <- GK
            SWbefore <- SW
            GK <- GK + hard[1]/2 + homeAdv[1]/2
            SW <- SW + hard[2]/2 + homeAdv[2]/2
            DF <- DF + hard[3] + homeAdv[3]
            MF <- MF + hard[4] + homeAdv[4]
            ST <- ST + hard[5] + homeAdv[5]
            
            if (GK  > 10 && GK > GKbefore) {
              GK <- max(GKbefore, 10)
              warning(paste("GK was reduced to the strength ", GK ," (after adding hardness and home advantage).",
                      sep = ""))
            }
            
            if (SW > 10 && SW > SWbefore) {
              SW <- max(SWbefore, 10)
              warning("SW was reduced to the strength ", SW," (after adding hardness and home advantage).",
                      sep = "")
            }
            
            minStrength <- min(DF, MF, ST)[1]
            if (3*minStrength < DF) {
              DF <- 3*minStrength
              warning("DF was reduced according to the 3:1 rule (after adding hardness and home advantage).")
            }
            if (3*minStrength < MF) {
              MF <- 3*minStrength
              warning("MF was reduced according to the 3:1 rule (after adding hardness and home advantage).")
            }
            if (3*minStrength < ST) {
              ST <- 3*minStrength
              warning("ST was reduced according to the 3:1 rule (after adding hardness and home advantage).")
            }
            
            return(c(GK, SW, DF, MF, ST))
          }
)

#' @rdname simRedCard
setMethod("simRedCard", 
          signature(obj = "formation", lineup = "numeric"),
          function(obj, lineup) {
            Hard <- matrix(c(90,10,0,0,0,0,0,0,70,30,0,0,0,0,0,0,50,40,10,
                             0,0,0,0,0,30,50,20,0,0,0,0,0,20,40,30,10,0,0,
                             0,0, 10,30,40,20,0,0,0,0,0,20,40,30,10,0,0,0,0,
                             10,30,40,20,0,0,0,0,0,20,40,30,10,0,0,0,0,10,20,
                             40,20,10,0,0,0,0,10,40,20,20,10), nrow = 8)
            sumHard <- sum(obj@hardness)
            greaterZero <- which(Hard[ ,sumHard+1] > 0)
            numberCards <- sample(greaterZero, 1, prob = Hard[ ,sumHard+1][greaterZero]) - 1
            if (numberCards <= 1) {
              return(lineup)
            }
            hard <- obj@hardness
            # allocation of the cards on the players
            # doubled probability for rows with hardness
            # gk has probability of zero, if there was no hardness
            if (hard[1] == 0) {
              propPos <- c(0,1,1,1,1)
            } else {
              propPos <- c(1,1,1,1,1)
            }  
              
            propPos[which(hard > 0)] <- propPos[which(hard > 0)] + 1
            # now transfer this to the singular players
            numberGK <- 1
            numberSW <- sum(!is.na(obj@SW))
            numberDF <- length(obj@DF)
            numberMF <- length(obj@MF)
            numberST <- length(obj@ST)
            strengthPlayers <- c(obj@GK, obj@SW, obj@DF, obj@MF, obj@ST)
            
            probPlayers <- rep(propPos, c(numberGK, numberSW, numberDF, 
                                          numberMF, numberST))
            players <- c("GK 1", rep("SW 1", numberSW), paste("DF", 1:numberDF), 
                         paste("MF", 1:numberMF), paste("ST", 1:numberST))
            givenCards <- sample(players, numberCards, replace = TRUE, 
                                 prob = probPlayers)
            if (!any(duplicated(givenCards))) {
            # scenario of no red card
                return(lineup)
            } else {
              # scenario of at least one red card  
              punishedPlayers <- unique(givenCards[duplicated(givenCards)])
              # punished players and their affect on the lineup
              for (i in punishedPlayers) {
                posNumb <- unlist(strsplit(i, split = " "))
                pos <- posNumb[1]
                numb <- as.numeric(posNumb[2])
                if (pos == "ST") {
                  strengthPunished <- obj@ST[numb]
                  strengthNew <- round(strengthPunished * sample(90, 1)/90 )
                  strengthLoss <- strengthPunished - strengthNew
                  lineup[5] <- max(0, lineup[5] - strengthLoss)
                }
                if (pos == "MF") {
                  strengthPunished <- obj@MF[numb]
                  strengthNew <- round(strengthPunished * sample(90, 1)/90 )
                  strengthLoss <- strengthPunished - strengthNew
                  lineup[4] <- max(0, lineup[4] - strengthLoss)
                }
                if (pos == "DF") {
                  strengthPunished <- obj@DF[numb]
                  strengthNew <- round(strengthPunished * sample(90, 1)/90 )
                  strengthLoss <- strengthPunished - strengthNew
                  lineup[3] <- max(0, lineup[3] - strengthLoss)
                }
                if (pos == "SW") {
                  strengthPunished <- obj@SW[numb]
                  strengthNew <- round(strengthPunished * sample(90, 1)/90 )
                  strengthLoss <- strengthPunished - strengthNew
                  lineup[2] <- max(0, lineup[2] - strengthLoss)
                }
                if (pos == "GK") {
                  strengthPunished <- obj@GK[numb]
                  strengthNew <- round(strengthPunished * sample(90, 1)/90 )
                  strengthLoss <- strengthPunished - strengthNew
                  lineup[1] <- max(0, lineup[1] - strengthLoss)
                }
                  
              }
              return(lineup)
            }
          }
)  
          
# --------------------------------------------
# Show function for formation
# --------------------------------------------

setMethod("show", "formation", function(object) {
    # headline
    cat("\nObject of class \"", class(object)[1],"\"\n\n", sep = "")
    cat("Your selected lineup is:\n")
    lineup <- toString(getLineup(object))
    lineup <- gsub(", ", "-", lineup)
    cat("\t", lineup)
    # iterate through all slots of the object
    #names <- slotNames(object)
    #for (name in names) {
    #  cat(name, "=", slot(object, name), "\n")
    #}
    cat("\n") 
  }  
)




