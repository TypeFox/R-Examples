    oldwd <- getwd(); setwd(paste(find.package("TDMR"), "demo02sonar",sep="/"));
    envT = tdmEnvTLoad("demoSonar.RData");    # loads envT
    source("main_sonar.r");
    envT$tdm$nrun=0;       # =0: no unbiasedRun, >0: perform unbiasedRun with opts$NRUN=envT$tdm$nrun
    finals = tdmEnvTSensi(envT,1);
    if (!is.null(finals)) print(finals);
    setwd(oldwd);
