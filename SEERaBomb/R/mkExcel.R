mkExcel=function(seerSet,tsdn,outDir="~/Results",outName=NULL,flip=FALSE) {
  if (length(tsdn)>1) {
    cat("collapsing brks vector to a tsdn string\n") 
    tsdn=paste0("b",paste(tsdn,collapse="_"))
  }
  if (is.null(seerSet$L)) stop("seerSet L field is empty. Please run tsd on your seerSet object!") else 
    L=seerSet$L[[tsdn]]
  if (!dir.exists(outDir)) dir.create(outDir,recursive=T)
  unlink(f<-paste0(outDir,"/",ifelse(is.null(outName),paste0(seerSet$bfn,tsdn),outName),ifelse(flip,"Flipped",""),".xlsx"))
  wb <- loadWorkbook(f,create=T) 
  OL=NULL
  intvs=names(L[[L$trtS[1]]][["Obs"]]) 
  #   picks=rownames(L[["noRad"]][["Obs"]][[intvs[1]]])
  sheetS=L$firstS
  if (flip) sheetS=seerSet$secondS
  for (icanc in sheetS) {
    #    icanc="prostate"
    createSheet(wb, name = icanc)
    M=NULL
    for (R in L$trtS) {
      D=NULL
      for (intv in intvs) {
        if (flip) {
          O=L[[R]][["Obs"]][[intv]][,icanc,drop=FALSE]
          E=L[[R]][["Exp"]][[intv]][,icanc,drop=FALSE]
        } else {
          O=L[[R]][["Obs"]][[intv]][icanc,,drop=FALSE]
          E=L[[R]][["Exp"]][[intv]][icanc,,drop=FALSE]
        }
        #         print(O)
        #         print(E)
        RR=O/E
        #         print(RR)
        LL=qchisq(.025,2*O) /(2*E)
        UL=qchisq(.975,2*O+2)/(2*E)
        if (flip) col=as.data.frame(round(cbind(RR,LL,UL,O=O,E=E),2)) else
          col=as.data.frame(round(cbind(t(RR),t(LL),t(UL),O=t(O),E=t(E)),2))
        names(col)=c("RR","LL","UL","O","E")
        D=cbind(D,paste0(col$RR,"(",col$LL,",",col$UL,") O=",O))
        #     print(R)
      }
      #   D
      colnames(D)=paste(intvs,"after",R)
      rownames(D)=rownames(col)
      M=cbind(M,D)
    } #rad
    #     writeWorksheet(wb, data.frame("second cancer"=picks,M), sheet = icanc,rownames=1)
    if (flip) writeWorksheet(wb, cbind("1st cancer"=L$firstS,M), sheet = icanc,rownames=1) else
      writeWorksheet(wb, cbind("2nd cancer"=seerSet$secondS,M), sheet = icanc,rownames=1)
    OL[[icanc]]=M
    setColumnWidth(wb,sheet = icanc, column = 1, width = 2500)
    for (j in 2:(dim(M)[2]+1)) setColumnWidth(wb,sheet = icanc, column = j, width = 4700)
    #     for (j in 1:(dim(M)[2]+1)) setColumnWidth(wb,sheet = icanc, column = j, width = 5600)
    createFreezePane(wb,sheet = icanc,2,2) 
  } #icanc
  saveWorkbook(wb)
  cat("Workbook was written to",f,"\n")
  invisible(OL)
}
