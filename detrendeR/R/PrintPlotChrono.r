PrintPlotChrono = 
function(crn, rwl=NULL, file.name = "crono", crono.type="RAW", saveCronoJpg=TRUE ){

FILE_NAME = paste(file.name, crono.type, sep="")
CRONO_NAME =paste(file.name, crono.type, sep="-")
FILE_NAME_CSV = paste(file.name, crono.type, ".CSV", sep="")
CURRENT_DIRECTORY = getwd();

        if (!is.null(rwl)) rwl = rwl[crn[,ncol(crn)] > 0,]
        crn = crn[crn[,ncol(crn)] > 0,]

CRN = crn[,c(1,ncol(crn)-1)]

round(crn[,1:(ncol(crn)-2)],3) ->crn[,1:(ncol(crn)-2)]
  WriteMatrix(crn, File=FILE_NAME_CSV, sep = ";", row.names =TRUE,
                                row.name ="Year", ID=FALSE, col.width=4)        # Save CSV

  write.crn(CRN, paste(FILE_NAME, ".CRN", sep=""))                              # Save CRN

  jpeg(paste(CRONO_NAME, ".jpg", sep=""), width = 600, height = 600, quality = 100)                                       # Save Chono JPG
  PlotChron(crn)
  dev.off()

 if (!is.null(rwl)) {
  jpeg(paste(CRONO_NAME, "-series.jpg", sep=""), width = 600, height = 600, quality = 100)                                # Save Chono JPG
  yr.vec = as.numeric(rownames(rwl))
  matplot(yr.vec,rwl,ylim=c(0,max(rwl, na.rm=TRUE)), type="l", xlab = "Years", ylab = "", main = file.name, las = 1, col="grey20", lty =1, lwd=1)
  lines(yr.vec, crn[,1], col="red", lwd=2)
  dev.off()
  }

  setwd(dirname(getwd()))
  dir.create(paste("Cronos", crono.type, sep="-") , showWarnings = FALSE);
  setwd(paste("Cronos", crono.type, sep="-"));
  write.crn(CRN, paste(file.name, crono.type, ".CRN", sep=""));
  setwd(CURRENT_DIRECTORY);

}
