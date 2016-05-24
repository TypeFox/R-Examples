#
# 13_10_30 13_11_04 13_12_27 13_12_28 13_12_30
# 13_12_31 14_01_03
#

library(rbsa);
#
fifi <- file2char("basic.code.r",
                  clean=FALSE,comme=character(0));
anat <- code7objects4text6tags(fifi,addbra=TRUE);
anan <- code7objects4text6tags(fifi);
print(anat);
pause("code7objects4text6tags");

code <- anan[[2]];
lyse <- parse8code(code);
print(lyse);
#
#pause("parse8code");

rdres <- make8rd(lyse);
print(rdres);

fifi <- file2char("objects.code.r",
                  clean=FALSE,comme=character(0));
anan <- code7objects4text6tags(fifi);
for (ii in bf(anan)) {
  cat("ii:",ii,"  ");
  lyse <- parse8code(anan[[ii]]);
  cat(lyse$name,"\n");
  rdf <- make8rd(lyse);
}
pause("documented objects");
