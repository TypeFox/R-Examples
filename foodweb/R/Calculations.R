calculate.metrics <- function(foodweb, maxlevels, name) {

S <<- ncol(foodweb)

#Set level=0 for all resources
  for (sp in 1:S) {
      if (sum(foodweb[1:S, sp])==0) {
        foodweb[S+1, sp] <- 0
      }
  }


#Assign integer trophic levels      
    for (level in 1:maxlevels) { #for each trophic level above basal resource, stopping at maxlevels (default = 8)
              for (sp in 1:S) { #for each species
                    #check if it eats any of the species one level below it  
                    for (resource in which(foodweb[S+1,]==level-1)) {
                          if (foodweb[resource, sp]==1 & foodweb[sp,resource]==0) { # if it does,
                              foodweb[S+1,sp] <- as.numeric(level) #assign the current level
                              break
                          }
                    }
              }
    }

#Omnivory
## Create a matrix where the trophic _level_ of the prey is represented
  ## i.e. multiply each column by a vector with the trophic _level_ of each species
  ## note: in order to differentiate between "feeds on level zero" and "doesn't feed on this species"
  ## I multiply by (trophic level + 1), turn remaining zeros into NA, then subtract 1
      by.level <- sweep(foodweb[1:S,], MARGIN=1, as.vector((foodweb[S+1,]+1), mode="numeric"), FUN="*")
      is.na(by.level[by.level==0]) <- TRUE #sets all values that are zero to NA
      by.level <- sweep(by.level[1:S,], MARGIN=1, 1, FUN="-") #substracts one in order to get back the original level            

#Record trophic position in row S+2
      for (sp in 1:S) {
        if (sum(foodweb[1:S,sp])==0) {
          foodweb[S+2,sp] <- 0
        } else {
          foodweb[S+2,sp] <- (sum(by.level[,sp], na.rm=TRUE)/sum(foodweb[1:S,sp]))+1
          }
      }


## Create a table specifying the species at each trophic position 
  levels.l <- data.frame(TrophicLevel = character(0), S=integer(0), sp= character(0))         
    for (level in unique(t((foodweb[S+2,])))) {
      levels.l <- rbind(levels.l, cbind(level, length(which(foodweb[S+2,]==level)), paste(which(foodweb[S+2,]==level), collapse=", ")))    
    }
  names(levels.l) <- c("Trophic position", "S", "Species")          

  levels.l <<- levels.l

## Create a table specifying which trophic position each species feeds on
  #Create a matrix where the trophic _positon_ of the prey is represented
  # i.e. multiply each column by a vector with the trophic position of each species
      by.position <- sweep(foodweb[1:S,], MARGIN=1, as.vector((foodweb[S+2,]+1), mode="numeric"), FUN="*")
      is.na(by.position[by.position==0]) <- TRUE #sets all values that are zero to NA
      by.position <- sweep(by.position[1:S,], MARGIN=1, 1, FUN="-") #substracts one in order to get back the original level            

 omn.l <- data.frame(Predator=integer(0), PredTrophPos=integer(0), Prey=integer(0))        
    for (sp in which(foodweb[S+1,]!=foodweb[S+2,])) {
      omn.l <- rbind(omn.l, cbind(sp, formatC(foodweb[S+2,sp],digits=2, format="f", drop0trailing=TRUE) , paste(formatC(as.numeric(na.omit(unique(by.position[,sp]))), digits=2, format="f", drop0trailing=TRUE), collapse=",  ")))
    }
 names(omn.l) <- c("Species", "Trophic position", "Prey(s)' position")
 omn.l <<- omn.l
 omn <- nrow(omn.l)

## Intraguild predation
# Tag cannibalistic species by placing a 1 in row S+3 


  for (sp in 1:S) {
      if (sum(foodweb[which(foodweb[S+1,-(sp)]==foodweb[S+1,sp]),sp])>0) {
        foodweb[S+3,sp] <- 1
      } else {foodweb[S+3,sp] <- 0}
  }
 
foodweb[S+3,sp]

cann <- sum(foodweb[S+3,])

#Assign the trophic level of the cannibal species to row S+3
temp <- sweep(foodweb[S+3,], MARGIN=2, as.vector((foodweb[S+1,]), mode="numeric"), FUN="*")
 is.na(temp[temp==0]) <- TRUE #sets all values that are zero to NA
foodweb[S+3,] <- temp


#Create a table specifying the trophic level at which cannibalism occurs

cann.l <- data.frame(TrophicLevel = character(0), S.cann=integer(0), sp= character(0))         
    for (level in na.omit(unique(t((foodweb[S+3,]))))) {
      if (level != 0) {
      cann.l <- rbind(cann.l, cbind(level, length(which(foodweb[S+3,]==level)), paste(which(foodweb[S+3,]==level), collapse=", ")))    
      }
    }  


names(cann.l) <- c("Trophic level", "# cannibalistic sp.", "Cannibal species")

cann.l <<- cann.l

#Fractions of species

  # Number of basal species
  Basal = sum(colSums(foodweb[1:S,])==0)

  # Number of top species
  Top = length(which(foodweb[S+1,]==max(foodweb[S+1,])))
  Top.l = which(foodweb[S+1,]==max(foodweb[S+1,]))

  # Number of intermediate species
  Int = S - Basal - Top

  # Number of herbivores
  Herb = length(which(foodweb[S+1,]==1))

#Other metrics
  # Total number of links
  L <- sum(foodweb[1:S,])

  # Connectance
  C <- L/(S^2)

  # Average number of links ("linkage density")
  AL = L/S

# Put indices in a vector

indices <<- vector(length=12)
indices <<- cbind(name,S, L, C, AL, omn/S, cann/S, length(unique(as.numeric(foodweb[S+1,]))),
                  Basal/S, Int/S, Top/S, Herb/S,(Basal+Int)/(Top+Int))
colnames(indices) <<- c("Web ID","Species richness", "Total # Links", "Connectance", "Link density", "Frac omniv", "Frac cannib", 
                       "Total # trophic positions", "Frac basal", "Frac intermediate", "Frac top", "Frac herbiv", "Prey:Predator")

foodweb <<- foodweb
}