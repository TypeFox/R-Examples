## This demo creates figure 5
data(deff.rho)
data(dpp.rho)
par(mfrow=c(2,2))
plotresp.mix(deff.rho,TITLE="(a)")
plotresp.equiv(deff.rho,TITLE="(b)")
plotresp.mix(dpp.rho,XYLIM=c(0,100),RLAB="% Protected by",TITLE="(c)")
plotresp.equiv(dpp.rho,XLIM=c(0,100),RLAB="% Protected by",bounds=c(.01,99.9),TITLE="(d)")
par(mfrow=c(1,1))
## Each of the following lines are ways to output the graph in different formats
## change the filename as is appropriate, and remove the # to run
#dev.print(postscript,"H:\\main\\malaria\\saul\\combination\\figures\\figure5.ps")
#dev.print(pdf,"H:\\main\\malaria\\saul\\combination\\figures\\figure5.pdf")
#dev.print(win.metafile,"H:\\main\\malaria\\saul\\combination\\figures\\figure5.wmf")
 

