## This demo creates figure 3
data(deff.sigma)
data(dpp.sigma)
par(mfrow=c(2,2))
plotresp.mix(deff.sigma,TITLE="(a)")
plotresp.equiv(deff.sigma,TITLE="(b)")
plotresp.mix(dpp.sigma,XYLIM=c(0,100),RLAB="% Protected by",TITLE="(c)")
plotresp.equiv(dpp.sigma,XLIM=c(0,100),RLAB="% Protected by",bounds=c(.01,99.9),TITLE="(d)")
par(mfrow=c(1,1))
## Each of the following lines are ways to output the graph in different formats
## change the filename as is appropriate, and remove the # to run
#dev.print(postscript,"H:\\main\\malaria\\saul\\combination\\figures\\figure3.ps")
#dev.print(pdf,"H:\\main\\malaria\\saul\\combination\\figures\\figure3.pdf")
#dev.print(win.metafile,"H:\\main\\malaria\\saul\\combination\\figures\\figure3.wmf")
 

