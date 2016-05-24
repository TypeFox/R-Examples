#
# 13_12_31 14_01_03 14_09_01 14_09_04
#
library(rbsa);
source8file(Sys.glob("*.code.r"));
#
ofif <- Sys.glob("*.code.r");
nfif <- paste0("../resu/",ofif);
fifi <- cbind(ofif,nfif);
#
anob <- c("analyse8description","code7objects4text6tags","parse8code","make8rd","write8lyse","component9","extract8object","object","text");
noob <- c(                "a8d",               "c7o4t6t",       "p8c",    "m8r",       "w8l",        "c9",           "e8o",   "NON","FAUX");
nana <- cbind(anob,noob);
#
change8names(fifi,nana);
pause("change8names");
#
documair7dir <- ".";
destination7dir <- "../resu";
napkg <- "documair";
resu <- build8pkg(napkg,documair7dir,destination7dir,what="lz");
if (length(resu)>0) {
  print(resu);
  stop("an error occurred!");
}
pause("build8pkg");


