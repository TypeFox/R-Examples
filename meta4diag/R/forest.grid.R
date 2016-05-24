forest.grid <- function(x, ...) UseMethod("forest.grid")

forest.grid.meta4diag = function(x, accuracy.type="sens", est.type="mean", nameShow=T, dataShow=F, ciShow=T, graphwidth=1, main, xlab="",...){
  
  if(length(accuracy.type)!=1){stop("Argument \"accuracy.type\" could only be one character string.")}
  if(!is.character(accuracy.type)){stop("Argument \"accuracy.type\" could only be one character string.")}
  accuracy.type = tolower(accuracy.type)
  suitable.set = c("sens", "TPR", "spec", "TNR", "FPR", "FNR", "LRpos", "LRneg", "DOR")
  if(!(accuracy.type %in% tolower(suitable.set))){
    stop(paste("Please give the correct accuracy.type type, which could be ",paste(suitable.set, collapse=", "),".",sep=""))
  }
  if(!x$misc$sample.flag){
    if(accuracy.type %in% tolower(c("LRpos", "LRneg", "DOR"))){
      stop("The statistics is not the default return. Please let \"nsample=TRUE\" in the \"meta4diag()\" function.")
    }
  }
  if(!(est.type %in% c("mean","median"))){
    stop("Argument \"est.type\" could only be either \"mean\" or \"median\".")
  }
  ####### check nameShow
  if(!is.logical(nameShow)){
    if(!is.character(nameShow)){
      nameFlag = FALSE
      stop("Argument \"nameShow\" could only be FALSE, TRUE, \"left\", \"right\" or \"center\".")
    }else{
      if(tolower(nameShow) %in% c("left","right","center")){
        nameShow = tolower(nameShow)
        nameFlag = TRUE
      }else{
        nameFlag = FALSE
        stop("Argument \"nameShow\" could only be FALSE, TRUE, \"left\", \"right\" or \"center\".")
      }
    }
  }else{
    if(nameShow){
      nameShow = "right"
      nameFlag = TRUE
    }else{nameFlag = FALSE}
  }
  
  ####### check dataShow
  if(!is.logical(dataShow)){
    if(!is.character(dataShow)){
      dataFlag = FALSE
      stop("Argument \"dataShow\" could only be FALSE, TRUE, \"left\", \"right\" or \"center\".")
    }else{
      if(tolower(dataShow) %in% c("left","right","center")){
        dataShow = tolower(dataShow)
        dataFlag = TRUE
      }else{
        dataFlag = FALSE
        stop("Argument \"dataShow\" could only be FALSE, TRUE, \"left\", \"right\" or \"center\".")
      }
    }
  }else{
    if(dataShow){
      dataShow = "right"
      dataFlag = TRUE
    }else{dataFlag = FALSE}
  }
  
  ####### check ciShow
  if(!is.logical(ciShow)){
    if(!is.character(ciShow)){
      ciFlag = FALSE
      stop("Argument \"ciShow\" could only be FALSE, TRUE, \"left\", \"right\" or \"center\".")
    }else{
      if(tolower(ciShow) %in% c("left","right","center")){
        ciShow = tolower(ciShow)
        ciFlag = TRUE
      }else{
        ciFlag = FALSE
        stop("Argument \"ciShow\" could only be FALSE, TRUE, \"left\", \"right\" or \"center\".")
      }
    }
  }else{
    if(ciShow){
      ciShow = "left"
      ciFlag = TRUE
    }else{ciFlag = FALSE}
  }
  
  ######################### main estimates
  if(accuracy.type=="sens" || accuracy.type=="tpr"){
    fitname = "True positive rate (Sensitivity)"
    fit = x[["summary.fitted.(Se)"]]
    sfit = x[["summary.summarized.fitted"]]["mean(Se)",]
    xmin = min(fit[,"0.025quant"])
    xmax = max(fit[,"0.975quant"])
  }
  if(accuracy.type=="spec" || accuracy.type=="tnr"){
    fitname = "True negative rate (Specificity)"
    fit = x[["summary.fitted.(Sp)"]]
    sfit = x[["summary.summarized.fitted"]]["mean(Sp)",]
    xmin = min(fit[,"0.025quant"])
    xmax = max(fit[,"0.975quant"])
  }
  if(accuracy.type=="fpr"){
    fitname = "False positive rate (1-Specificity)"
    fit = x[["summary.fitted.(1-Sp)"]]
    sfit = x[["summary.summarized.fitted"]]["mean(1-Sp)",]
    xmin = min(fit[,"0.025quant"])
    xmax = max(fit[,"0.975quant"])
    
  }
  if(accuracy.type=="fnr"){
    fitname = "False negative rate (1-Sensitivity)"
    fit = x[["summary.fitted.(1-Se)"]]
    sfit = x[["summary.summarized.fitted"]]["mean(1-Se)",]
    xmin = min(fit[,"0.025quant"])
    xmax = max(fit[,"0.975quant"])
    
  }
  if(accuracy.type=="lrpos"){
    fitname = "Positive likelihood ratio (LR+)"
    fit = x[["summary.fitted.LRpos"]]
    sfit = x[["summary.summarized.statistics"]]["mean(LRpos)",]
    scale = 1.5*(sfit["0.975quant"]-sfit["0.025quant"])
    xmin = max(0, sfit["0.025quant"]-scale) 
    xmax = sfit["0.975quant"]+scale
    
  }
  if(accuracy.type=="lrneg"){
    fitname = "Negative likelihood ratio (LR-)"
    fit = x[["summary.fitted.LRneg"]]
    sfit = x[["summary.summarized.statistics"]]["mean(LRneg)",]
    scale = 1.5*(sfit["0.975quant"]-sfit["0.025quant"])
    xmin = max(0, sfit["0.025quant"]-scale) 
    xmax = sfit["0.975quant"]+scale
    
  }
  if(accuracy.type=="dor"){
    fitname = "Diagnostic odds ratio (DOR)"
    fit = x[["summary.fitted.DOR"]]
    sfit = x[["summary.summarized.statistics"]]["mean(DOR)",]
    scale = 1.5*(sfit["0.975quant"]-sfit["0.025quant"])
    xmin = max(0, sfit["0.025quant"]-scale) 
    xmax = sfit["0.975quant"]+scale
  }
  
  
  nr = dim(fit)[1]
  
  if(missing(main)){
    main = paste("Forest plot for ",fitname,sep="")
  }
  
  
  ###### grid use
  if(est.type=="mean"){
    est = round(fit[,1],2)
    sest = round(sfit[1],2)
    ml.est = max(nchar(as.character(est)))
    ci.est = format(est, nsmall = 2L, width=ml.est, justify="right")
    ci.sest = format(sest, nsmall = 2L)
  }
  if(est.type=="median"){
    est = round(fit[,"0.5quant"],2)
    sest = round(sfit["0.5quant"],2)
    ml.est = max(nchar(as.character(est)))
    ci.est = format(est, nsmall = 2L, width=ml.est, justify="right")
    ci.sest = format(sest, nsmall = 2L)
  }
  
  lb = round(fit[,"0.025quant"],2)
  slb = round(sfit["0.025quant"],2)
  ml.lb = max(nchar(as.character(lb)))
  ci.lb = format(lb, nsmall = 2L, width=ml.lb, justify="right")
  
  ub = round(fit[,"0.975quant"],2)
  sub = round(sfit["0.975quant"],2)
  ml.ub = max(nchar(as.character(ub)))
  ci.ub = format(ub, nsmall = 2L, width=ml.ub, justify="right")
  
  ci = paste(ci.est," [ ",ci.lb,", ",ci.ub," ]",sep="")
  sci = paste(ci.sest, " [ ", format(slb,nsmall=2L), ", ", format(sub,nsmall=2L), " ]", sep="")
  
  exlb = which(lb < xmin)
  exub = which(ub > xmax)
  
  xrange<-c(max(min(lb),xmin), min(max(ub),xmax))
  
  cwidth<-(ub-lb)
  
  info<-1/cwidth
  info<-info/max(info)
  ###### count the width first
  ##### remember to check the is.na(studynames)
  if(is.null(graphwidth)){graphwidth=2}
  graphwidths=unit(graphwidth,"inches")
  
  flags = c(nameFlag, dataFlag, ciFlag)*1
  ncFlag = sum(flags)
  
  nc = ncFlag + 1 
  
  if(nameFlag){
    namex = switch(nameShow,left=0,right=1,center=0.5)
    namejust = switch(nameShow,left="left",right="right",center="center")
    enamelabel = lapply(rownames(fit), function(fit_i) textGrob(fit_i, x=namex,just=namejust, gp=gpar(fontface="plain", col="black")))
    snamelabel = list(textGrob("Summary", x=namex,just=namejust, gp=gpar(fontface="bold", col="black")))
    tnamelabel = list(textGrob("Study", x=namex,just=namejust, gp=gpar(fontface="bold", col="black")))
    namelabel = c(tnamelabel, enamelabel, snamelabel)
    namewidths = unit.c(max(unit(rep(1,nr+2),"grobwidth",namelabel)),unit(3,"mm"))
  }
  
  if(dataFlag){
    ######### dataShow labels
    PP = x$data$TP + x$data$FN
    TPlab = format(x$data$TP, width=max(nchar(as.character(x$data$TP))))
    PPlab = format(PP, width=max(nchar(as.character(PP))))
    datalab1 = paste(TPlab,"/",PPlab,sep="")
    NN = x$data$TN + x$data$FP
    TNlab = format(x$data$TN, width=max(nchar(as.character(x$data$TN))))
    NNlab = format(NN, width=max(nchar(as.character(NN))))
    datalab2 = paste(TNlab,"/",NNlab,sep="")
    
    datax = switch(dataShow, left=0, right=1, center=0.5)
    datajust = switch(dataShow, left="left", right="right", center="center")
    edatalabel1 = lapply(datalab1, function(data_i) textGrob(data_i, x=datax,just=datajust, gp=gpar(fontface="plain", col="black")))
    tdatalabel1 = list(textGrob("TP/(TP+FN)", x=datax,just=datajust, gp=gpar(fontface="bold", col="black")))
    sdatalabel1 = list(textGrob(" ", x=datax,just=datajust, gp=gpar(fontface="bold", col="black")))
    datalabel1 = c(tdatalabel1, edatalabel1, sdatalabel1)
    datawidths1 = unit.c(max(unit(rep(1,nr+2),"grobwidth",datalabel1)),unit(6,"mm"))
    edatalabel2 = lapply(datalab2, function(data_i) textGrob(data_i, x=datax,just=datajust, gp=gpar(fontface="plain", col="black")))
    tdatalabel2 = list(textGrob("TN/(TN+FP)", x=datax,just=datajust, gp=gpar(fontface="bold", col="black")))
    sdatalabel2 = list(textGrob(" ", x=datax,just=datajust, gp=gpar(fontface="bold", col="black")))
    datalabel2 = c(tdatalabel2, edatalabel2, sdatalabel2)
    datawidths2 = unit.c(max(unit(rep(1,nr+2),"grobwidth",datalabel2)),unit(5,"mm"))
  }
  
  if(ciFlag){
    cix = switch(ciShow, left=0, right=1, center=0.5)
    cijust = switch(ciShow, left="left", right="right", center="center")
    ecilabel = lapply(ci, function(ci_i) textGrob(ci_i, x=cix,just=cijust, gp=gpar(fontface="plain", col="black")))
    scilabel = list(textGrob(sci, x=cix,just=cijust, gp=gpar(fontface="bold", col="black")))
    tcilabel = list(textGrob("Estimate", x=cix,just=cijust, gp=gpar(fontface="bold", col="black")))
    cilabel = c(tcilabel, ecilabel, scilabel)
    ciwidths = unit.c(unit(3,"mm"),max(unit(rep(1,nr+2),"grobwidth",cilabel)))
  }
  
  
  ############# grid plot
  plot.new()
  if(ncFlag==3){
    columnwidths = unit.c(namewidths, datawidths1, datawidths2, graphwidths, ciwidths)
    total_number = 9
    name_number =  1
    data_number = c(3,5)
    graph_number = 7
    ci_number = 9
  }
  if(ncFlag==2){
    if(!nameFlag){
      columnwidths = unit.c(datawidths1, datawidths2, graphwidths, ciwidths)
      total_number = 7
      data_number = c(1,3)
      graph_number = 5
      ci_number = 7
    }
    if(!dataFlag){
      columnwidths = unit.c(namewidths, graphwidths, ciwidths)
      total_number = 5
      name_number = 1
      graph_number = 3
      ci_number = 5
    }
    if(!ciFlag){
      columnwidths = unit.c(namewidths,datawidths1, datawidths2, graphwidths)
      total_number = 7
      name_number = 1
      data_number = c(3,5)
      graph_number = 7
    }
  }
  if(ncFlag==1){
    if(nameFlag){
      columnwidths = unit.c(namewidths, graphwidths)
      total_number = 3
      name_number =  1
      graph_number = 3
    }
    if(dataFlag){
      columnwidths = unit.c(datawidths1, datawidths2, graphwidths)
      total_number = 5
      data_number = c(1,3)
      graph_number = 5
    }
    if(ciFlag){
      columnwidths = unit.c(graphwidths, ciwidths)
      total_number = 3
      graph_number = 1
      ci_number = 3
    }
  }
  if(ncFlag==0){
    columnwidths = graphwidths
    total_number = 1
    graph_number = 1
  }
  pushViewport(viewport(layout=grid.layout(nr+2, total_number, widths=columnwidths, heights=unit(c(rep(1, nr+2),0.5), "lines"))))
  if(nameFlag){
    for(i in 1:(nr+2)){
      pushViewport(viewport(layout.pos.row=i,layout.pos.col=name_number))
      grid.draw(namelabel[[i]])
      popViewport()
    }
  }
  if(dataFlag){
    for(i in 1:(nr+2)){
      pushViewport(viewport(layout.pos.row=i,layout.pos.col=data_number[1]))
      grid.draw(datalabel1[[i]])
      popViewport()
    }
    for(i in 1:(nr+2)){
      pushViewport(viewport(layout.pos.row=i,layout.pos.col=data_number[2]))
      grid.draw(datalabel2[[i]])
      popViewport()
    }
  }
  pushViewport(viewport(layout.pos.col=graph_number, xscale=xrange))
  ticks<-pretty(xrange)
  xax<-xaxisGrob(gp=gpar(cex=0.6,col="black"),at=ticks,name="xax")
  xax1<-editGrob(xax, gPath("labels"), label=format(ticks,digits=2))
  grid.draw(xax1)
  grid.xaxis(gp=gpar(cex=0.6,col="black"))
  grid.text(xlab, y=unit(-3, "lines"),gp=gpar(col="black"))
  popViewport()
  
  for(i in 1:nr){
    pushViewport(viewport(layout.pos.row=i+1, layout.pos.col=graph_number, xscale=xrange))
    if(lb[i]<xmin && ub[i]<xmax){
      grid.lines(x=unit(c(xmin, est[i]), c("native","native")), y=0.5,
                 arrow=arrow(ends="first", angle=30, length=unit(0.05, "inches")),
                 gp=gpar(col="black"))
      grid.lines(x=unit(c(est[i], ub[i]), c("native","native")), y=0.5,
                 arrow=arrow(ends="last", angle=90, length=unit(0.05, "inches")),
                 gp=gpar(col="black"))
    }
    if(lb[i]>xmin && ub[i]>xmax){
      grid.lines(x=unit(c(lb[i], est[i]), c("native","native")), y=0.5,
                 arrow=arrow(ends="first", angle=90, length=unit(0.05, "inches")),
                 gp=gpar(col="black"))
      grid.lines(x=unit(c(est[i], xmax), c("native","native")), y=0.5,
                 arrow=arrow(ends="last", angle=30, length=unit(0.05, "inches")),
                 gp=gpar(col="black"))
    }
    if(lb[i]<xmin && ub[i]>xmax){
      grid.lines(x=unit(c(xmin, xmax), c("native","native")), y=0.5,
                 arrow=arrow(ends="both", angle=30, length=unit(0.05, "inches")),
                 gp=gpar(col="black"))
    }
    if(lb[i]>xmin && ub[i]<xmax){
      grid.lines(x=unit(c(lb[i], ub[i]), c("native","native")), y=0.5,
                 arrow=arrow(ends="both", angle=90, length=unit(0.05, "inches")),
                 gp=gpar(col="black"))
    }
    
    grid.rect(x=unit(est[i], "native"),
              width=unit(0.75*info[i], "snpc"), height=unit(0.75*info[i], "snpc"),
              gp=gpar(fill="black",col="black"))
    popViewport()
  }
  pushViewport(viewport(layout.pos.row=nr+2, layout.pos.col=graph_number, xscale=xrange))
  grid.lines(x=unit(c(slb, sub), c("native","native")), y=0.5,
             arrow=arrow(ends="both", angle=90, length=unit(0.05, "inches")),
             gp=gpar(col="black"))
  grid.polygon(x=unit(c(0.5*(sest+slb), sest, 0.5*(sub+sest), sest), "native"),
               y=unit(0.5 + c(0, 0.2, 0, -0.2), "npc"),gp=gpar(fill="black",col="black"))
  popViewport()
  
  pushViewport(viewport(layout.pos.col=graph_number, xscale=xrange))
  grid.lines(x=unit(slb, "native"), y=0:1, gp=gpar(col="darkgray", lty=2))
  grid.lines(x=unit(sub, "native"), y=0:1, gp=gpar(col="darkgray", lty=2))
  grid.polygon(x=unit(c(slb,sub,sub,slb), "native"),
               y=unit(c(0, 0, 1, 1), "npc"),gp=gpar(fill="gray",col="gray",alpha=0.3))
  popViewport()
  if(ciFlag){
    for(i in 1:(nr+2)){
      pushViewport(viewport(layout.pos.row=i,layout.pos.col=ci_number))
      grid.draw(cilabel[[i]])
      popViewport()
    }
  }
}