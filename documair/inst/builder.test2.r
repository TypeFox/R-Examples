#
library(rbsa);
# parameters
ouca <- "./perso/";
oume <- "../";
# sourcing everything
fico <- list.files(path = ouca, pattern = "*.code.r$");
for (fifi in paste0(ouca,fico)) {
  print(fifi);
  source(fifi);
}
#
# source doubtful /rbsa/ function
source("perso/rbsa.code.r");
# test of alias of parse8code
#
uu <- file2char("perso/builder.test2.txt",comme="#####");
vv <- parse8code(uu,
                 otags=documair0$tags$v,
                 tags=names(documair0$tags$v));
print(vv$ali);
print(vv$det);
cat("Test is finished\n");
