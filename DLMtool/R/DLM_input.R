# DLM_input MPs

matlenlim<-function(x,DLM_data, ...){ # Length at maturity is knife-edge vulnerability
  dependencies="DLM_data@LFC, DLM_data@LFS"
  Allocate<-1
  Effort<-1
  Spatial<-c(1,1)
  
  newLFC <- DLM_data@L50[x] * 0.95
  newLFS <- DLM_data@L50[x] 
  Vuln <-c(newLFC, newLFS)
  c(Allocate, Effort, Spatial, Vuln)
}
class(matlenlim)<-"DLM_input"

matlenlim2 <-function(x,DLM_data, ...){ # Knife-edge vulnerability slightly higher than length at maturity 
  dependencies="DLM_data@LFC, DLM_data@LFS"
  Allocate<-1
  Effort<-1
  Spatial<-c(1,1)
  newLFS <- 14 + 1.1 * DLM_data@L50[x] 
  newLFC <- 0.95 * newLFS
  Vuln <-c(newLFC, newLFS)
  c(Allocate, Effort, Spatial, Vuln)
}
class(matlenlim2)<-"DLM_input"


# matagelim<-function(x,DLM_data, ...){ # Age at maturity is knife-edge vulnerability
  # dependencies="DLM_data@AM, DLM_data@MaxAge"
  # Allocate<-1
  # Effort<-1
  # Spatial<-c(1,1)
  # Vuln<-1/(1+exp((DLM_data@AM[x]-(1:DLM_data@MaxAge))/(DLM_data@AM[x]*DLM_data@AM[x]*0.05)))
  # c(Allocate, Effort, Spatial, Vuln)
# }
# class(matagelim)<-"DLM_input"

MRreal<-function(x,DLM_data, ...){ # A Marine reserve in area 1 with spatial reallocation of effort
  dependencies="DLM_data@MaxAge"
  Allocate<-1
  Effort<-1
  Spatial<-c(0,1)
  # Vuln<-rep(NA,DLM_data@MaxAge)
  Vuln<-rep(NA,2)
  c(Allocate, Effort, Spatial, Vuln)
}
class(MRreal)<-"DLM_input"

MRnoreal<-function(x,DLM_data, ...){ # A Marine reserve in area 1 with no spatial reallocation of effort
  dependencies="DLM_data@MaxAge"
  Allocate<-0
  Effort<-1
  Spatial<-c(0,1)
  # Vuln<-rep(NA,DLM_data@MaxAge)
  Vuln<-rep(NA,2)
  c(Allocate, Effort, Spatial, Vuln)
}
class(MRnoreal)<-"DLM_input"

curE<-function(x,DLM_data, ...){ # current effort
  dependencies="DLM_data@MaxAge"
  Allocate<-1
  Effort<-1
  Spatial<-c(1,1)
  # Vuln<-rep(NA,DLM_data@MaxAge)
  Vuln<-rep(NA,2)  
  c(Allocate, Effort, Spatial, Vuln)
}
class(curE)<-"DLM_input"

curE75<-function(x,DLM_data, ...){ #75% current effort
  dependencies="DLM_data@MaxAge"
  Allocate<-1
  Effort<-0.75
  Spatial<-c(1,1)
  # Vuln<-rep(NA,DLM_data@MaxAge)
  Vuln<-rep(NA,2)  
  c(Allocate, Effort, Spatial, Vuln)
}
class(curE75)<-"DLM_input"

LBSPR_ItEff <- function(x, DLM_data, yrsmth=1, reps=reps) {
 dependencies="DLM_data@CAL, DLM_data@CAL_bins, DLM_data@vbLinf, 
	DLM_data@vbK, DLM_data@Mort, LM_data@vbK, DLM_data@L50, DLM_data@L95, 
	DLM_data@wlb"
  MiscList <- LBSPR(x, DLM_data, yrsmth=yrsmth, reps=reps)
  
  XX <- 1:4 
  YY <- MiscList[[2]][(length(MiscList[[2]]) - (max(XX)-1)):length(MiscList[[2]])]
  
  EstSPR <- YY[length(YY)]
  
  TgSPR <- 0.4
  h <- DLM_data@steep[x]
  SPRLim <- -(2*(h-1))/(3*h+1) # SPR that results in 0.5 R0
  
  phi1 <- 6
  phi2 <- 1
  
  MaxDw <- -0.3
  MaxUp <- 0.3
  
  minSlope <- 0.01
  
  Slope <- coef(lm(YY~XX))[2]  
  # if (abs(Slope) < minSlope) Slope <- 0 
  Dist <- EstSPR - TgSPR 
  
  # Control Rule #
  Mod <- 0 
  Buff <- 0.1
  Buffer <- c(TgSPR - Buff,  TgSPR + Buff)
  inBuff <- FALSE
  belowTG <- FALSE 
  aboveTG <- FALSE
  slopeUp <- FALSE
  slopeDw <- FALSE 
  belowLim <- FALSE
  if (Dist < 0) belowTG <- TRUE 
  if (Dist > 0) aboveTG <- TRUE 
  if (EstSPR > min(Buffer) & EstSPR < max(Buffer)) inBuff <- TRUE
  if (Slope <= 0) slopeDw <- TRUE
  if (Slope > 0) slopeUp <- TRUE
  if (EstSPR < SPRLim) belowLim <- TRUE
   
  # If within buffer zone - only slope
  if (inBuff) Mod <- phi1 * Slope
  if (slopeUp & aboveTG) Mod <- phi1 * Slope +  phi2 * Dist
  if (slopeUp & belowTG) Mod <- phi1 * Slope 
  
  if (slopeDw & aboveTG) Mod <- phi1 * Slope 
  if (slopeDw & belowTG) Mod <- phi1 * Slope +  phi2 * Dist
  
  if (belowLim) Mod <- MaxDw
  
  Mod[Mod > MaxUp] <- MaxUp
  Mod[Mod < MaxDw] <- MaxDw
  Mod <- Mod + 1 
  
  Allocate <- 1
  Effort <- DLM_data@MPrec[x] * Mod
  MiscList[[5]] <- append(MiscList[[5]], Effort)
  Spatial <- c(1,1)
  Vuln <- rep(NA,2)
  out <- c(Allocate, Effort, Spatial, Vuln)
   
  Out <- list()
  Out[[1]] <- out 
  Out[[2]] <- MiscList
 
  return(Out) 
}
class(LBSPR_ItEff)<-"DLM_input"

LBSPR_ItSel <- function(x, DLM_data, yrsmth=1, reps=reps) {
 dependencies="DLM_data@CAL, DLM_data@CAL_bins, DLM_data@vbLinf, 
	DLM_data@vbK, DLM_data@Mort, LM_data@vbK, DLM_data@L50, DLM_data@L95, 
	DLM_data@wlb"
  MiscList <- LBSPR(x, DLM_data, yrsmth=yrsmth,reps=reps)
  
  XX <- 1:4 
  YY <- MiscList[[2]][(length(MiscList[[2]]) - (max(XX)-1)):length(MiscList[[2]])]
  
  EstSPR <- YY[length(YY)]
  
  TgSPR <- 0.4
  h <- DLM_data@steep[x]
  SPRLim <- -(2*(h-1))/(3*h+1) # SPR that results in 0.5 R0
 
  Allocate <- 1
  Effort <- 1
  Spatial <- c(1,1)

  if (EstSPR < TgSPR) {
    newLFC <- DLM_data@L50[x] * 1.05
    newLFS <- DLM_data@L50[x] * 1.1
    Vuln <-c(newLFC, newLFS)
  }
  if (EstSPR < SPRLim) {
    newLFC <- DLM_data@L50[x] * 1.2
    newLFS <- DLM_data@L50[x] * 1.25
    Vuln <-c(newLFC, newLFS)  
  }
  if (EstSPR >= TgSPR) {
    newLFC <- DLM_data@L50[x] * 0.85
    newLFS <- DLM_data@L50[x] * 0.9
    Vuln <-c(newLFC, newLFS)  
  }
   
 
  out <- c(Allocate, Effort, Spatial, Vuln)
   
  Out <- list()
  Out[[1]] <- out 
  Out[[2]] <- MiscList
 
  return(Out) 
}
class(LBSPR_ItSel)<-"DLM_input"

