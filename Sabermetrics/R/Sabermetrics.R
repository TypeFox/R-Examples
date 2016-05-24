## Title: sabermetrics
## Authors: Peter Xenopoulos, Fernando Crema
## Email: peter.xenopoulos@pomona.edu, fernando.crema@ciens.ucv.ve
## Website: www.peterxeno.com
## Description: Sabermetrics contains a series of useful sabermetrics functions



# The XML package is necessary for loading of the "weights" data frame, which contains pertinent data for seasonal constants useful to a vareity of stats

linearWeights <- XML::readHTMLTable("http://www.fangraphs.com/guts.aspx?type=cn",
                         which = 15,
                         stringsAsFactors = FALSE)

names(linearWeights) <- c("Season","wOBA","wOBAScale","wBB","wHBP","w1B","w2B","w3B","wHR","runSB","runCS","RPerPA","RPerW","cFIP")
linearWeights <- linearWeights[-1,]
linearWeights = as.data.frame(apply(linearWeights, 2, as.numeric))


# This loads the "parkfactors" data frame, which contains park constants
# Returns dataframe
parkfactors <- function(year,team) {
  defaulturl <- "http://www.fangraphs.com/guts.aspx?type=pf&teamid=0&season=2014"
  if(missing(year) && missing(team)) {
    parkfactors <- XML::readHTMLTable(defaulturl,which=18)
    parkfactors = as.data.frame(apply(parkfactors, 2, as.numeric))
    return(parkfactors)
  }
  if(missing(team)) {
    url.1 <- "http://www.fangraphs.com/guts.aspx?type=pf&teamid=0&season="
    url.2 <- as.character(year)
    url <- paste(url.1,url.2,sep="")
    parkfactors <- XML::readHTMLTable(defaulturl,which=18)
    return(parkfactors)
  }
  if(missing(year)) {
    parkfactors <- XML::readHTMLTable(defaulturl,which=18)
    return(parkfactors)
  }
  else {
    url.base <- "http://www.fangraphs.com/guts.aspx?type=pf&teamid="
    url.team <- as.character(team)
    url.year <- as.character(year)
    url.incomplete <- paste(url.base,url.team,sep="")
    url.incomplete2 <- paste(url.incomplete,"&season=",sep="")
    url.complete <- paste(url.incomplete2,url.year,sep="")
    parkfactors <- XML::readHTMLTable(url.complete,which=18)
    return(parkfactors)
  }
}

# This loads the "handedparkfactors" data frame, which contains handedness park constants
# Returns dataframe
handedparkfactors <- function(year,team) {
  defaulturl <- "http://www.fangraphs.com/guts.aspx?type=pfh&teamid=0&season=2014"
  if(missing(year) && missing(team)) {
    parkfactors <- XML::readHTMLTable(defaulturl,which=18)
    return(parkfactors)
  }
  if(missing(team)) {
    url.1 <- "http://www.fangraphs.com/guts.aspx?type=pfh&teamid=0&season="
    url.2 <- as.character(year)
    url <- paste(url.1,url.2,sep="")
    parkfactors <- XML::readHTMLTable(defaulturl,which=18)
    return(parkfactors)
  }
  if(missing(year)) {
    parkfactors <- XML::readHTMLTable(defaulturl,which=18)
    return(parkfactors)
  }
  else {
    url.base <- "http://www.fangraphs.com/guts.aspx?type=pfh&teamid="
    url.team <- as.character(team)
    url.year <- as.character(year)
    url.incomplete <- paste(url.base,url.team,sep="")
    url.incomplete2 <- paste(url.incomplete,"&season=",sep="")
    url.complete <- paste(url.incomplete2,url.year,sep="")
    parkfactors <- XML::readHTMLTable(url.complete,which=18)
    return(parkfactors)
  }
}

# Slugging Percentage
slg <- function(TB, AB) {
  slugging <- (TB/AB)
  return(slugging)
}

# On-base Percentage
obp <- function(H, BB, HBP, AB, SF) {
  onbase <- ((H+BB+HBP)/(AB+BB+SF+HBP))
  return(onbase)
}

# On-base plus Slugging
ops <- function(slg, obp) {
  ops <- slg+obp
  return(ops)
}

# OPS+
opsplus <- function(obp,slg,lgOBP,lgSLG) {
  opsplus <- ((obp/lgOBP)+(slg/lgSLG)-1)*100
  return(opsplus)
}

# Isolated Power
iso <- function(slg, avg) {
  iso <- slg-avg
  return(iso)
}

# Pythagorean Expectation
pyth <- function(rs, ra, alpha = 2) {
  pyth <- (rs^alpha)/(rs^alpha + ra^alpha)
  return(pyth)
}

# Batting Average on Balls in Play
babip <- function(h, hr, ab, k, sf){
  babip <- (h-hr)/(ab-hr-k-sf)
  return(babip)
}

# Walks plus Hits per Inning Pitched
whip <- function(h, bb, ip){
  whip <- (h+bb)/ip
  return(whip)
}

# Earned Run Average
era <- function(er, ip){
  era <- (9*er)/ip
  return(era)
}

# Fielding Percentage
fp <- function(po, a, e){
  fp <-(po+a)/(po+a+e)
  return(fp)
}

# Secondary Average
secavg <- function(ab, bb, tb, h, sb, cs){
  SecA <- (bb+tb-h+sb-cs)/ab
  return(SecA)
}

# Defense-Independent Componente ERA
dice <- function(hr, bb, hbp, k, ip){
  dice = (13*hr+3*bb+3*hbp-2*k)/ip + 3.0
  return(dice)
}

# Log 5
log5 <- function(pA, pB, order=0){
  if(order){
    aux = pB
    pB = pA
    pA = aux
  }
  log5 <- (pA-pA*pB)/(pA+pB-2*pA*pB)
  return(log5)
}

# Raw Equivalent Average
raweqa <- function(AB,H,TB,BB,IBB,HBP,CS,SB,SH,SF) {
  raweqa <- (H+TB+1.5*(BB+HBP+SB)+SH+SF-IBB/2)/(AB+BB+HBP+SH+SF+CS+SB)
  return(raweqa)
}

# FIP
fip <- function(HR,BB,HBP,K,IP,year) {
  constant <- linearWeights$cFIP[which(linearWeights$Season == year)]
  fip <- ((13*HR)+(3*(BB+HBP))-(2*K))/IP + constant
  return(fip)
}

# xFIP
xfip <- function(flyballs,lgHRFB,BB,HBP,K,ip,year) {
  constant <- linearWeights$cFIP[which(linearWeights$Season == year)]
  xfip <- (13*(flyballs*lgHRFB)+(3*(BB+HBP))-(2*K))/(ip)+constant
  return(xfip)
}

#wOBA
woba <- function(year,AB,BB,IBB,HBP,single,double,triple,HR,SF) {
  wBB <- linearWeights$wBB[which(linearWeights$Season == year)]
  wHBP <- linearWeights$wHBP[which(linearWeights$Season == year)]
  w1B <- linearWeights$w1B[which(linearWeights$Season == year)]
  w2B <- linearWeights$w2B[which(linearWeights$Season == year)]
  w3B <- linearWeights$w3B[which(linearWeights$Season == year)]
  wHR <- linearWeights$wHR[which(linearWeights$Season == year)]
  woba <- ((wBB*BB)+(wHBP*HBP)+(w1B*single)+(w2B*double)+(w3B*triple)+(wHR*HR))/(AB+BB-IBB+SF+HBP)
  return(woba)
}

#wRAA
wraa <- function(woba,year,PA) {
  wraa <- ((woba-linearWeights$wOBA[which(linearWeights$Season == year)])/linearWeights$wOBAScale[which(linearWeights$Season == year)])*PA
  return(wraa)
}

#ERA-
eraminus <- function(ERA,ParkFactor,LeagueERA) {
  # Park Factors for ERA must be divided by 100
  eraminus <- (ERA+(ERA-(ERA*(ParkFactor/100))))/(LeagueERA)*100
  return(eraminus)
}

#FIP-
fipminus <- function(FIP,ParkFactor,LeagueFIP) {
  fipminus <- (FIP+(FIP-(FIP*(ParkFactor/100))))/(LeagueFIP)*100
  return(fipminus)
}

#xFIP-
xfipminus <- function(xFIP,ParkFactor,LeaguexFIP) {
  xfipminus <- (xFIP+(xFIP-(xFIP*(ParkFactor/100))))/(LeaguexFIP)*100
  return(xfipminus)
}

#wRC
wrc <- function(wOBA,PA,year) {
  wrc <- (((wOBA-linearWeights$wOBA[which(linearWeights$Season == year)])/linearWeights$wOBAScale[which(linearWeights$Season == year)])+linearWeights$RPerPA[which(linearWeights$Season == year)])*PA
  return(wrc)
}

#wRC+
wrcplus <- function(wRAA,PA,year,parkfactor,leaguewRC) {
  leaguerpa <- linearWeights$RPerPA[which(linearWeights$Season == year)]
  wrcplus <- (((wRAA/PA + leaguerpa)+(leaguerpa-parkfactor*leaguerpa))/leaguewRC)*100
}

#wSB
wsb <- function(year,SB,CS,lgwSB,single,bb,hbp,ibb) {
  runSB <- linearWeights$runSB[which(linearWeights$Season == year)]
  runCS <- linearWeights$runCS[which(linearWeights$Season == year)]
  wsb <- (SB * runSB) + (CS * runCS) - (lgwSB * (single + bb + hbp - ibb))
}

#lgwSB
lgwSB <- function(year,SB,CS,single,bb,hbp,ibb) {
  runSB <- linearWeights$runSB[which(linearWeights$Season == year)]
  runCS <- linearWeights$runCS[which(linearWeights$Season == year)]
  lgwSB <- (SB * runSB + CS * runCS) / (single + bb + hbp - ibb)
}