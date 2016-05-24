driveP <- function(){

if(BBB2$acc.ratio >  BBB2$prior.odds){
BBB2$tempting <- BBB2$tempting + BBB2$h.step.width
}

if(BBB2$acc.ratio < BBB2$prior.odds){
BBB2$tempting <- BBB2$tempting - BBB2$h.step.width 
}

BBB2$acc.ratio <- 0

}


