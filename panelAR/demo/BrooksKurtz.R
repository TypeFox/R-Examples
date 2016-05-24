BrooksKurtz <- function()
{
    message("Replication of Table 1 of Brooks and Kurtz (2012)'...")
    
    cat("\n")
    cat("Loading data...\n")
    data(BrooksKurtz)
    cat("> data(BrooksKurtz)\n")	
	
	# Table 1, Model 1: 
	message("Table 1, Model 1")
    cat("> out1.1 <- panelAR(kaopen ~ ldiffpeer + ldiffisi + ldiffgrowth + ldiffinflation  + ldiffneg + ldiffembi + limf + isi_objective + partisan + checks +  lusffr + linflation + lbankra + lcab + lgrowth +  ltradebalance + lngdpcap + lngdp + brk + timetrend + y1995, data=BrooksKurtz, panelVar='country', timeVar='year', autoCorr='psar1', panelCorrMethod='parks',rho.na.rm=TRUE, panel.weight='t', seq.times=TRUE)")
    out1.1 <- panelAR(kaopen ~ ldiffpeer + ldiffisi + ldiffgrowth + ldiffinflation  + ldiffneg + ldiffembi + limf + isi_objective + partisan + checks +  lusffr + linflation + lbankra + lcab + lgrowth +  ltradebalance + lngdpcap + lngdp + brk + timetrend + y1995, data=BrooksKurtz, panelVar='country', timeVar='year', autoCorr='psar1', panelCorrMethod='pwls',rho.na.rm=TRUE, panel.weight='t', seq.times=TRUE)
    cat("> summary(out1.1)")
    print(summary(out1.1))

	user.prompt <- function () {
 		ANSWER <- readline("\nType 'y' to continue to Model 2 or 'n' to quit: ")
      
	if (substr(ANSWER, 1, 1) == "n")
         {stop("Function terminated by user.")}
     
	}

	user.prompt()
	
	# Table 1, Model 2: 
	message("Table 1, Model 2")
    cat("> out1.2 <- panelAR(kaopen ~ ldiffisi + ldiffgrowth + ldiffinflation  + ldiffneg + ldiffembi + limf + isi_objective + partisan + checks +  lusffr + linflation + lbankra + lcab + lgrowth +  ltradebalance + lngdpcap + lngdp + brk + timetrend + y1995, data=BrooksKurtz, panelVar='country', timeVar='year', autoCorr='psar1', panelCorrMethod='parks',rho.na.rm=TRUE, panel.weight='t', seq.times=TRUE)")
    out1.2 <- panelAR(kaopen ~ ldiffisi + ldiffgrowth + ldiffinflation  + ldiffneg + ldiffembi + limf + isi_objective + partisan + checks +  lusffr + linflation + lbankra + lcab + lgrowth +  ltradebalance + lngdpcap + lngdp + brk + timetrend + y1995, data=BrooksKurtz, panelVar='country', timeVar='year', autoCorr='psar1', panelCorrMethod='pwls',rho.na.rm=TRUE, panel.weight='t', seq.times=TRUE)
    cat("> summary(out1.2)")
    print(summary(out1.2))
}

BrooksKurtz()