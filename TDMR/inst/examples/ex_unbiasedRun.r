    ## Load the best results obtained in a prior tuning for the configuration "sonar_04.conf" with tuning method "spot".
    ## The result envT from a prior run of tdmBigLoop with this .conf is read from demo02sonar/demoSonar.RData.
    ## Run task main_sonar again with these best parameters, using the default settings from tdmDefaultsFill
    ## umode="RSUB", tdm$nrun=5  and tdm$TST.testFrac=0.2.
    oldwd <- getwd();
    ## The best results are read from demo02sonar/demoSonar.RData relative to the TDMR package directory.
    setwd(paste(find.package("TDMR"), "demo02sonar",sep="/"));
    load("demoSonar.RData");
    tdm=envT$tdm;
    source("main_sonar.r");
    finals <- unbiasedRun("sonar_04.conf",envT,tdm=tdm);
    print(finals);
    setwd(oldwd);
