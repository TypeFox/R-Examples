# construct an initial environment envT from the TDMR settings given in tdm:
envT <- tdmEnvTMakeNew(tdm);
# start the big tuning loop:
envT <- tdmBigLoop(envT,spotStep);

