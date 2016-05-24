
envT = tdmEnvTLoad("demoSonar.RData");    # load envT
source("main_sonar.r");
envT$tdm$nrun=2;       # =0: no, =2: two unbiased runs
finals = tdmEnvTSensi(envT,1);
if (!is.null(finals)) print(finals); 
