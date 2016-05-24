visAFCM <- function (data, scannf=FALSE, nf=2, xax = 1, yax = 2, clab.row = FALSE, clab.col = 1,permute = FALSE
                     , posieig = "top", sub = NULL, graphstyle = "unique", graphrow = 1, graphcol = 3, cpoint=1, clabel=2, csub=2)
{
  Data<-as.data.frame.matrix(data)
  for (i in 1:ncol(Data)) {
      Data[[i]]<-as.factor(Data[[i]])
  }
    
  # Graphs
  if (scannf == "TRUE") {dev.new()} else {}
  acmData <- dudi.acm(Data, scannf=scannf, nf=nf)
    
  if (graphstyle == "unique") {
    acmDatatmp <- acmData
    class(acmDatatmp) <- "dudi"
    scatter(acmDatatmp, xax = xax, yax = yax, clab.row = clab.row, clab.col = clab.col, permute = permute, posieig = posieig, sub = sub)
  } else {
    }

  if (graphstyle == "multiple a") {
      par(mfrow=c(graphrow,graphcol))
      for (v in 1:ncol(Data)) {
        s.class (acmData$li[,c(xax,yax)], fac=Data [,v], col=1:nlevels(Data [,v]),label=levels(Data [,v]),sub=colnames(Data)[v], cpoint=cpoint, clabel=clabel, csub=csub)
      }
      par(mfrow=c(1,1))
  } else {
    }

  if (graphstyle == "multiple b") {
      par(mfrow=c(graphrow,graphcol))
      for (v in 1:ncol(Data)) {
        s.chull (acmData$li[,c(xax,yax)], fac=Data [,v], col=1:nlevels(Data [,v]),label=levels(Data [,v]),sub=colnames(Data)[v], optchull=1, cpoint=cpoint, clabel=clabel, csub=csub)
      }
  par(mfrow=c(1,1))
  } else {
    }
}

