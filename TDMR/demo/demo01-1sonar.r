#*# --------- demo/demo00sonar.r ---------
#*# This demo shows a simple data mining process (level 1 of TDMR) for the classification task
#*# SONAR (from UCI repository, 
#*# http://archive.ics.uci.edu/ml/datasets/Connectionist+Bench+%28Sonar,+Mines+vs.+Rocks%29).
#*# The data mining process is in main_sonar.r, which calls tdmClassifyLoop and tdmClassify.
#*# with Random Forest as the prediction model. 
#*# main_sonar.r is sourced and executed from start_sonar.r.

path <- paste(find.package("TDMR"), "demo02sonar",sep="/");
#path <- paste("../inst", "demo02sonar",sep="/");

source(paste(path,"sonar_00.apd",sep="/"),local=TRUE);   # set opts 
source(paste(path,"start_sonar.r",sep="/"),chdir=TRUE,print.eval=TRUE);  # contains: result=main_sonar(opts);

