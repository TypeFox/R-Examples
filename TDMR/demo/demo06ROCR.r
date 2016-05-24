#*# --------- demo/demo06ROCR.r ---------
#*# Run task SONAR with "area under ROC curve" as performance measure (rgain.type="arROC").
#*# Other settings are similar to demo01-1sonar.r (level 1 of TDMR).
#*# Finally, plot ROC curve for validataion data set and 
#+#          plot lift chart for training data set
#*#
    path <- paste(find.package("TDMR"), "demo02sonar",sep="/");
    #path <- paste("../inst", "demo02sonar",sep="/");

    ## set the elements of list opts. See ?tdmOptsDefaultsSet
    ## for a complete list of all default settings and explanatory comments
    opts = tdmOptsDefaultsSet();    
    opts$filename = "sonar.txt"
    opts$READ.CMD = "readCmdSonar(filename,opts)"    # def'd in main_sonar.r
    opts$data.title <- "Sonar Data";
    opts$rgain.type <- "arROC";

    source(paste(path,"start_sonar.r",sep="/"),chdir=TRUE);  # contains: result=main_sonar(opts);
    
    tdmGraphicNewWin(opts);
    cat("Area under ROC-curve for validation data set: ",
        tdmROCRbase(result),"\n");                                # side effect: plot ROC-curve
    tdmGraphicNewWin(opts);
    cat("Area under lift curve for  training data set: ",
        tdmROCRbase(result,dataset="training",typ="lift"),"\n");  # side effect: plot lift chart
    #setwd(oldwd);
