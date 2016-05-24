#*# This demo shows a level-2 example (SPOT tuning on task SONAR)

## load package and set working directory (dir with .apd, .conf and main_*.r file):
path <- paste(find.package("TDMR"), "demo02sonar",sep="/");
#path <- paste("../inst", "demo02sonar",sep="/");

tdm=list(mainFile="main_sonar.r"
        ,runList="sonar_01.conf"
        );
spotStep = "auto";    
source(paste(path,tdm$mainFile,sep="/"));    
source(paste(path,"start_bigLoop.r",sep="/"),chdir=TRUE);    # change dir to 'path' while sourcing
