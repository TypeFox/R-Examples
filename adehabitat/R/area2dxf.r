"area2dxf" <- function(x, file, lay=1:nlevels(factor(x[,1])))
{
    ## Verifications of file format
    if (!inherits(x, "area"))
        stop("x should be of class area")
    if (substr(file, nchar(file)-3, nchar(file))!=".dxf")
        file<-paste(file, ".dxf", sep="")

    ## Verifications that the polygons are closed (identical first
    ## and last point for all polygons)
    lipol<-split(x, x[,1])
    for (i in 1:length(lipol)) {
        j<-lipol[[i]]
        if (!all(j[1,]==j[nrow(j),]))
            lipol[[i]]<-rbind.data.frame(lipol[[i]], lipol[[i]][1,])
    }
    x<-do.call("rbind.data.frame",lipol)

    ## header of file
    text<-"  0\nSECTION\n  2\nHEADER\n  9\n$EXTMIN\n 10\n"
    text<-paste(text, min(x[,2]),"\n", sep="")
    text<-paste(text, " 20\n", sep="")
    text<-paste(text, min(x[,3]),"\n", sep="")
    text<-paste(text, "  9\n$EXTMAX\n 10\n", sep="")
    text<-paste(text, max(x[,2]),"\n", sep="")
    text<-paste(text, " 20\n", sep="")
    text<-paste(text, max(x[,3]),"\n", sep="")
    text<-paste(text, "  0\nENDSEC\n  0\nSECTION\n", sep="")
    text<-paste(text, "2\nTABLES\n  0\nENDSEC\n  0\n", sep="")
    text<-paste(text, "SECTION\n  2\nBLOCKS\n  0\n", sep="")
    text<-paste(text, "ENDSEC\n  0\nSECTION\n  2\nENTITIES\n", sep="")

    ## The main part of the file
    lp<-split(x[,2:3], x[,1])
    for (i in 1:length(lp)) {
      text<-paste(text, "  0\nPOLYLINE\n  8\n", sep="")
      text<-paste(text, "  ",lay[i],"\n 66\n      1\n", sep="")
      for (j in 1:nrow(lp[[i]]))
        text<-paste(text, "  0\nVERTEX\n  8\n",
                    lay[[i]], "\n 10\n", lp[[i]][j,1],
                    "\n 20\n", lp[[i]][j,2], "\n", sep="")
      text<-paste(text, "  0\nSEQEND\n")
    }
    text<-paste(text, "      0\nENDSEC\n  0\nEOF\n")

    ## write the file
    cat(text, file=file)
  }

