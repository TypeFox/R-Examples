gplot.SSrat <- function(calcedgroup, maintitle=c('Sympathy','Antipathy','Preference','Impact'), 
                        SS.alpha="SS.05", thresh=0){
  maintitle=match.arg(maintitle)
  switch(maintitle,
         Sympathy = {plotdata=calcedgroup$S},
         Antipathy = {plotdata=calcedgroup$A},
         Preference = {plotdata=calcedgroup$P},
         Impact = {plotdata=calcedgroup$I}
  )
  if (thresh > 0) {maintitle=paste(maintitle,", threshold = ", thresh,sep='')}
  nrratees = nrow(calcedgroup$dataframe)
  status= as.factor(calcedgroup$dataframe[, SS.alpha])
  x=gplot(plotdata, usecurve = F, main = maintitle, label=calcedgroup$dataframe$respid, 
          vertex.col=as.numeric(status), vertex.sides=4, thresh=thresh)
  opar <- par(fig=c(0, 1, 0, 1), oma=c(0, 0, 0, 0), 
              mar=c(0, 0, 0, 0), new=TRUE)
  on.exit(par(opar))
  plot(0, 0, type='n', bty='n', xaxt='n', yaxt='n')
  legend("topright", legend=paste(calcedgroup$dataframe$respid, calcedgroup$dataframe$resplabel))
  legend("topleft", legend=levels(status), fill=1:length(levels(status)))
  invisible(x)
}


