## text book example from "Configuraional comparative Methods"
library(QCA3)
## csQCA
conditions <- c("GNPCAP", "URBANIZA", "LITERACY", "INDLAB", "GOVSTAB")
ans1a <- reduce(Lipset_cs,"SURVIVAL",conditions,explain="positive",remainder="exclude",case="CASEID")
QCA3:::prettyPI(ans1a)
## Formula 1 in Rihoux and De Meur(2009:57)
ans0a <- reduce(Lipset_cs,"SURVIVAL",conditions,explain="negative",remainder="exclude",case="CASEID")
QCA3:::prettyPI(ans0a)
## Formula 3 in Rihoux and De Meur(2009:59)
ans1 <- reduce(Lipset_cs,"SURVIVAL",conditions,explain="positive",remainder="include",case="CASEID")
QCA3:::prettyPI(ans1)
## Formula 4 in Rihoux and De Meur(2009:60)
QCA3:::prettyPI(SA(ans1)) ## 5 simplifying assumptions in p61
ans0 <- reduce(Lipset_cs,"SURVIVAL",conditions,explain="negative",remainder="include",case="CASEID")
QCA3:::prettyPI(ans0)
## Formula 5 in Rihoux and De Meur(2009:61)
QCA3:::prettyPI(SA(ans0)) ## 18 simplifying assumptions


## mvQCA
conditions <- c("GNPCAP", "URBANIZA", "LITERACY", "INDLAB")
if (packageDescription('QCA3')$Version <= "0.0-2") {
  mvTT <- cs_truthTable(Lipset_mv,"SURVIVAL",conditions,case="CASEID",nlevels=c(3,2,2,2))
} else mvTT <- mv_truthTable(Lipset_mv,"SURVIVAL",conditions,case="CASEID")

ans1a <- reduce(mvTT,explain="positive",remainder="exclude",case="CASEID")
QCA3:::prettyPI(ans1a)
## formula 1 Cronqvist and Berg-Schlosser(2009:80)
ans1 <- reduce(mvTT,explain="positive",remainder="include",case="CASEID")
QCA3:::prettyPI(ans1)
## formula 2 in Cronqvist and Berg-Schlosser(2009:81)
QCA3:::prettyPI(SA(ans1))
## 9 SAs (see end note 7)
ans0a <- reduce(mvTT,explain="negative",remainder="exclude",case="CASEID")
QCA3:::prettyPI(ans0a)
## formula 3 in Cronqvist and Berg-Schlosser(2009:81)
ans0 <- reduce(mvTT,explain="negative",remainder="include",contrad="positive",case="CASEID")
QCA3:::prettyPI(ans0)
## formula 4 in Cronqvist and Berg-Schlosser(2009:81)
QCA3:::prettyPI(SA(ans0))
## 7 SAs (see end note 9)


## fsQCA
conditions <- c("Developed.FZ","Urban.FZ","Literate.FZ","Industrial.FZ", "Stable.FZ")
ans1a <- reduce(Lipset_fs,"Survived.FZ",conditions,explain="positive",remaind="exclude",prepro="fs",consistency=0.7)
QCA3:::prettyPI(ans1a)
## Formula 1 in Ragin (2009:112)
ans0a <- reduce(Lipset_fs,"Survived.FZ",conditions,explain="positive",remaind="include",prepro="fs",consistency=0.7)
QCA3:::prettyPI(ans0a)
## Formula 2 in Ragin (2009:114)
ans0 <- reduce(Lipset_fs,"Survived.FZ",conditions,explain="negative",remaind="exclude",prepro="fs",consistency=0.7)
QCA3:::prettyPI(ans0)
## Formula 5 in Ragin (2009:115)
ans1 <- reduce(Lipset_fs,"Survived.FZ",conditions,explain="negative",remaind="include",prepro="fs",consistency=0.7)
QCA3:::prettyPI(ans1)
## Formula 6 in Ragin (2009:117)
