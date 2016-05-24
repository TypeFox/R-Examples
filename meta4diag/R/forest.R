forest <- function(x, ...) UseMethod("forest")

forest.meta4diag = function(x, accuracy.type="sens", est.type="mean", p.cex="scaled", p.pch=15, p.col="black",
                            nameShow="right", dataShow="center", ciShow="left", text.cex=1,
                            shade.col="gray", arrow.col="black", arrow.lty=1, arrow.lwd=1,
                            cut=c(0,1), intervals=c(0.025,0.975),
                            main="Forest plot", main.cex=1.5, axis.cex=1,...){
  
  op <- par(no.readonly = TRUE)
  if(length(accuracy.type)!=1){stop("Argument \"accuracy.type\" could only be one character string.")}
  if(!is.character(accuracy.type)){stop("Argument \"accuracy.type\" could only be one character string.")}
  accuracy.type = tolower(accuracy.type)
  suitable.set = c("sens", "TPR", "spec", "TNR", "FPR", "FNR", "LRpos", "LRneg", "RD", "LLRpos", "LLRneg", "LDOR", "DOR")
  if(!(accuracy.type %in% tolower(suitable.set))){
    stop(paste("Please give the correct accuracy.type type, which could be ",paste(suitable.set, collapse=", "),".",sep=""))
  }
  if(!x$misc$sample.flag){
    if(accuracy.type %in% tolower(c("LRpos", "LRneg", "RD", "LLRpos", "LLRneg", "LDOR", "DOR"))){
      stop("The statistics is not the default return. Please let \"nsample=TRUE\" in the \"meta4diag()\" function.")
    }
  }
  est.type = tolower(est.type)
  if(!(est.type %in% c("mean","median"))){
    stop("Argument \"est.type\" could only be either \"mean\" or \"median\".")
  }
  
  if(!is.numeric(intervals)){
    stop("Argument \"intervals\" has to be a numerical vector with length 2.")
  }
  if(length(intervals)!=2){
    stop("Argument \"intervals\" has to be a numerical vector with length 2.")
  }
  if(est.type=="median"){est.type="0.5quant"}
  
  ###### check intervals
  if(!all(intervals %in% x$misc$quantiles)){
    stop(paste("Argument \"intervals\" has to the values of quantiles. The options are ",paste(x$misc$quantiles,collapse=", "),sep=""))
  }
  if(intervals[1]>=0.5){
    stop("The first element of argument \"intervals\" has to be smaller than 0.5.")
  }
  if(intervals[2]<=0.5){
    stop("The first element of argument \"intervals\" has to be larger than 0.5.")
  }
  
  intervals = paste(intervals,"quant",sep="")
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
      ciShow = "right"
      ciFlag = TRUE
    }else{ciFlag = FALSE}
  }
  if(tolower(p.cex)!="scaled"){
    if(!is.numeric(p.cex)){
      stop("Argument \"p.cex\" could only be scaled or fixed to a positive numerical value.")
    }
  }
  
  if(!is.logical(cut)){
    if(!is.numeric(cut)){
      stop("cut has to be TRUE, FALSE or a numerical vector with length 2!")
    }else{
      if(length(cut)!=2){
        stop("cut has to be TRUE, FALSE or a numerical vector with length 2!")
      }
    }
  }

  ###### data
  datalab1 = format(x$data$tp, width=max(nchar(as.character(x$data$tp),type="width")))
  datalab2 = format(x$data$fp, width=max(nchar(as.character(x$data$fp),,type="width")))
  datalab3 = format(x$data$tn, width=max(nchar(as.character(x$data$tn),type="width")))
  datalab4 = format(x$data$fn, width=max(nchar(as.character(x$data$fn),type="width")))

  
  ######################### main estimates
  if(accuracy.type=="sens" || accuracy.type=="tpr"){
    fitname = "True positive rate (Sensitivity)"
    fit = x[["summary.fitted.(Se)"]]
  }
  if(accuracy.type=="spec" || accuracy.type=="tnr"){
    fitname = "True negative rate (Specificity)"
    fit = x[["summary.fitted.(Sp)"]]
  }
  if(accuracy.type=="fpr"){
    fitname = "False positive rate (1-Specificity)"
    fit = x[["summary.fitted.(1-Sp)"]]
  }
  if(accuracy.type=="fnr"){
    fitname = "False negative rate (1-Sensitivity)"
    fit = x[["summary.fitted.(1-Se)"]]
  }
  if(accuracy.type=="lrpos"){
    fitname = "Positive likelihood ratio (LR+)"
    fit = x[["summary.fitted.LRpos"]]
  }
  if(accuracy.type=="lrneg"){
    fitname = "Negative likelihood ratio (LR-)"
    fit = x[["summary.fitted.LRneg"]]
  }
  if(accuracy.type=="dor"){
    fitname = "Diagnostic odds ratio (DOR)"
    fit = x[["summary.fitted.DOR"]]
  }
  if(accuracy.type=="ldor"){
    fitname = "Log Diagnostic odds ratio (LDOR)"
    fit = x[["summary.fitted.LDOR"]]
  }
  if(accuracy.type=="rd"){
    fitname = "Risk difference (RD)"
    fit = x[["summary.fitted.RD"]]
  }
  if(accuracy.type=="llrpos"){
    fitname = "Log Positive likelihood ratio (LLR+)"
    fit = x[["summary.fitted.LLRpos"]]
  }
  if(accuracy.type=="llrneg"){
    fitname = "Log Negative likelihood ratio (LLR-)"
    fit = x[["summary.fitted.LLRneg"]]
  }
  est.fit = fit[,est.type]  
  lb.fit = fit[,intervals[1]]
  ub.fit = fit[,intervals[2]]
  
  ml.est.fit = max(nchar(as.character(round(est.fit,2))))
  ci.est.fit = format(round(est.fit,2), width=ml.est.fit, nsmall = 2L)
  
  ml.lb.fit = max(nchar(as.character(round(lb.fit,2)),type="width"))
  ci.lb.fit = format(round(lb.fit,2), nsmall = 2L, width=ml.lb.fit, justify="right")
  
  ml.ub.fit = max(nchar(as.character(round(ub.fit,2)),type="width"))
  ci.ub.fit = format(round(ub.fit,2), nsmall = 2L, width=ml.ub.fit, justify="right")
  
  ci.fit = paste(ci.est.fit," [ ",ci.lb.fit,", ",ci.ub.fit," ]",sep="")
  
  nr = dim(fit)[1]
  
  ###### different situations
  if(x$misc$covariates.flag){
    if(is.logical(cut)){
      if(cut){
        if(accuracy.type %in% c("lrpos","lrneg","dor")){
          scale = 0.3*(max(ub.fit)-min(lb.fit))
          xmax = min(ub.fit+scale)
          xmin = max(0, lb.fit-scale)
        }else{
          xmin = min(fit[,intervals[1]])
          xmax = max(fit[,intervals[2]])
        }
      }else{
        xmin = min(fit[,intervals[1]])
        xmax = max(fit[,intervals[2]])
      }
    }else{
      xmin = cut[1]
      xmax = cut[2]
    }
    exlb = lb.fit < xmin
    exub = ub.fit > xmax
  }else{ # no covariates
    if(x$misc$modality.flag){
      mod.level = x$misc$modality.level
      mod.name = x$misc$modality.name
      level.name = as.character(unique(x$data[,mod.name]))
      ind = lapply(1:mod.level, function(i) which(x$data[,mod.name]==level.name[i]))
      level.length = unlist(lapply(1:mod.level, function(i) length(ind[[i]])))
      if(accuracy.type=="sens" || accuracy.type=="tpr"){sfit = do.call(rbind, lapply(1:mod.level, function(i) x[["summary.expected.accuracy"]][paste("mean(Se.",level.name[i],")",sep=""),]))}
      if(accuracy.type=="spec" || accuracy.type=="tnr"){sfit = do.call(rbind, lapply(1:mod.level, function(i) x[["summary.expected.accuracy"]][paste("mean(Sp.",level.name[i],")",sep=""),]))}
      if(accuracy.type=="fpr"){sfit = do.call(rbind, lapply(1:mod.level, function(i) x[["summary.expected.accuracy"]][paste("mean(1-Sp.",level.name[i],")",sep=""),]))}
      if(accuracy.type=="fnr"){sfit = do.call(rbind, lapply(1:mod.level, function(i) x[["summary.expected.accuracy"]][paste("mean(1-Se.",level.name[i],")",sep=""),]))}
      if(accuracy.type=="lrpos"){sfit = do.call(rbind, lapply(1:mod.level, function(i) x[["summary.summarized.statistics"]][[i]]["mean(LRpos)",]))}
      if(accuracy.type=="lrneg"){sfit = do.call(rbind, lapply(1:mod.level, function(i) x[["summary.summarized.statistics"]][[i]]["mean(LRneg)",]))}
      if(accuracy.type=="dor"){sfit = do.call(rbind, lapply(1:mod.level, function(i) x[["summary.summarized.statistics"]][[i]]["mean(DOR)",]))}
      if(accuracy.type=="ldor"){sfit = do.call(rbind, lapply(1:mod.level, function(i) x[["summary.summarized.statistics"]][[i]]["mean(LDOR)",]))}
      if(accuracy.type=="rd"){sfit = do.call(rbind, lapply(1:mod.level, function(i) x[["summary.summarized.statistics"]][[i]]["mean(RD)",]))}
      if(accuracy.type=="llrpos"){sfit = do.call(rbind, lapply(1:mod.level, function(i) x[["summary.summarized.statistics"]][[i]]["mean(LLRpos)",]))}
      if(accuracy.type=="llrneg"){sfit = do.call(rbind, lapply(1:mod.level, function(i) x[["summary.summarized.statistics"]][[i]]["mean(LLRneg)",]))}
      
      est.sfit = sfit[,est.type]
      lb.sfit = sfit[,intervals[1]]
      ub.sfit = sfit[,intervals[2]]

    }else{ # no modality
      if(accuracy.type=="sens" || accuracy.type=="tpr"){sfit = x[["summary.expected.accuracy"]]["mean(Se)",]}
      if(accuracy.type=="spec" || accuracy.type=="tnr"){sfit = x[["summary.expected.accuracy"]]["mean(Sp)",]}
      if(accuracy.type=="fpr"){sfit = x[["summary.expected.accuracy"]]["mean(1-Sp)",]}
      if(accuracy.type=="fnr"){sfit = x[["summary.expected.accuracy"]]["mean(1-Se)",]}
      if(accuracy.type=="lrpos"){sfit = x[["summary.summarized.statistics"]]["mean(LRpos)",]}
      if(accuracy.type=="lrneg"){sfit = x[["summary.summarized.statistics"]]["mean(LRneg)",]}
      if(accuracy.type=="dor"){sfit = x[["summary.summarized.statistics"]]["mean(DOR)",]}
      if(accuracy.type=="ldor"){sfit = x[["summary.summarized.statistics"]]["mean(LDOR)",]}
      if(accuracy.type=="rd"){sfit = x[["summary.summarized.statistics"]]["mean(RD)",]}
      if(accuracy.type=="llrpos"){sfit = x[["summary.summarized.statistics"]]["mean(LLRpos)",]}
      if(accuracy.type=="llrneg"){sfit = x[["summary.summarized.statistics"]]["mean(LLRneg)",]}
      
      est.sfit = sfit[est.type]
      lb.sfit = sfit[intervals[1]]
      ub.sfit = sfit[intervals[2]]  
    }
    
    ml.est.sfit = max(nchar(as.character(round(est.sfit,2))))
    ci.est.sfit = format(round(est.sfit,2), width=ml.est.sfit, nsmall = 2L)
    
    ml.lb.sfit = max(nchar(as.character(round(lb.sfit,2)),type="width"))
    ci.lb.sfit = format(round(lb.sfit,2), nsmall = 2L, width=ml.lb.sfit, justify="right")
    
    ml.ub.sfit = max(nchar(as.character(round(ub.sfit,2)),type="width"))
    ci.ub.sfit = format(round(ub.sfit,2), nsmall = 2L, width=ml.ub.sfit, justify="right")
    
    ci.sfit = paste(ci.est.sfit," [ ",ci.lb.sfit,", ",ci.ub.sfit," ]",sep="")
    
    if(is.logical(cut)){
      if(cut){
        if(accuracy.type %in% c("lrpos","lrneg","dor")){
          scale = 1.5*(max(ub.sfit)-min(lb.sfit))
          xmin = max(0, lb.sfit-scale) 
          xmax = min(ub.sfit+scale)
        }else{
          xmin = min(fit[,intervals[1]])
          xmax = max(fit[,intervals[2]])
        }
      }else{
        xmin = min(fit[,intervals[1]])
        xmax = max(fit[,intervals[2]])
      }
    }else{
      xmin = cut[1]
      xmax = cut[2]
    }
    
    exlb = lb.fit < xmin
    exub = ub.fit > xmax
    
  } # end of covariates
  
  studynames = rownames(fit)
  ###### caculate width
  strwidth_names_main = strwidth(studynames,units="in",cex=text.cex,font=1,family="sans")
  strwidth_names_title = strwidth("Study",units="in",cex=text.cex,font=2,family="sans")
  strwidth_names_sum = strwidth("Summary",units="in",cex=text.cex,font=2,family="sans")
  
  name_width = max(strwidth_names_main,strwidth_names_title,strwidth_names_sum)
  
  strwidth_data_main1 = strwidth(datalab1,units="in",cex=text.cex,font=1,family="sans")
  strwidth_data_main2 = strwidth(datalab2,units="in",cex=text.cex,font=1,family="sans")
  strwidth_data_main3 = strwidth(datalab3,units="in",cex=text.cex,font=1,family="sans")
  strwidth_data_main4 = strwidth(datalab4,units="in",cex=text.cex,font=1,family="sans")
  
  strwidth_data_title1 = strwidth("TP",units="in",cex=text.cex,font=2,family="sans")
  strwidth_data_title2 = strwidth("FP",units="in",cex=text.cex,font=2,family="sans")
  strwidth_data_title3 = strwidth("TN",units="in",cex=text.cex,font=2,family="sans")
  strwidth_data_title4 = strwidth("FN",units="in",cex=text.cex,font=2,family="sans")
  
  data1_width = max(strwidth_data_main1,strwidth_data_title1)
  data2_width = max(strwidth_data_main2,strwidth_data_title2)
  data3_width = max(strwidth_data_main3,strwidth_data_title3)
  data4_width = max(strwidth_data_main4,strwidth_data_title4)
  
  if(accuracy.type %in% c("dor","lrneg","lrpos")){
    family = "mono"
  }else{
    family = "sans"
  }
  
  if(accuracy.type %in% c("dor","lrneg","lrpos")){
    strwidth_ci_main = strwidth(ci.fit,units="in",cex=text.cex,font=1,family=family)
    if(!x$misc$covariates.flag){
      strwidth_ci_sum = strwidth(ci.sfit,units="in",cex=text.cex,font=2,family=family)
    } 
  }else{
    strwidth_ci_main = strwidth(ci.fit,units="in",cex=text.cex,font=1,family=family)
    if(!x$misc$covariates.flag){
      strwidth_ci_sum = strwidth(ci.sfit,units="in",cex=text.cex,font=2,family=family)
    }
  }
  
  strwidth_ci_title = strwidth("Estimate",units="in",cex=text.cex,font=2,family="sans")
  
  if(!x$misc$covariates.flag){
    ci_width = max(strwidth_ci_main,strwidth_ci_title,strwidth_ci_sum)
  }else{
    ci_width = max(strwidth_ci_main,strwidth_ci_title)
  }
  
  figure_width = max(name_width,data1_width,data2_width,ci_width)*2.3
  
  if(missing(main)){
    main = paste("Forest plot for ",tolower(fitname),sep="")
  }
  
  flags = c(nameFlag, dataFlag, ciFlag)*1
  ncFlag = sum(flags)
  
  nc = ncFlag + 1
  
  if(tolower(p.cex)=="scaled"){
    xrange = ub.fit-lb.fit
    info = 1/xrange
    info = info/max(info)
    info = info[1:nr]*text.cex
  }else{
    info = rep(p.cex,nr)
  }
  
  if(ncFlag==3){
    layout(matrix(c(1,2,3,4,5,6,7),1,7),c(name_width,data1_width,data2_width,data3_width,data4_width,figure_width,ci_width),1) # name TP FP TN FN figure ci
  }
  if(ncFlag==2){
    if(!nameFlag){
      #figure_size = data1_width + data2_width + graph_width + ci_width
      layout(matrix(c(1,2,3,4,5,6),1,6),c(data1_width,data2_width,data3_width,data4_width,figure_width,ci_width),1) # TP FP TN FN figure ci
    }
    if(!dataFlag){
      #figure_size = name_width + graph_width + ci_width
      layout(matrix(c(1,2,3),1,3),c(name_width,figure_width,ci_width),1) # name figure ci
    }
    if(!ciFlag){
      #figure_size = name_width + data1_width + data2_width + graph_width 
      layout(matrix(c(1,2,3,4,5,6),1,6),c(name_width,data1_width,data2_width,data3_width,data4_width,figure_width),1) # name TP FP TN FN figure
    }
  }
  if(ncFlag==1){
    if(nameFlag){
      #figure_size = name_width + graph_width 
      layout(matrix(c(1,2),1,2),c(name_width,figure_width),1) # name figure
    }
    if(dataFlag){
      #figure_size = data1_width + data2_width + graph_width 
      layout(matrix(c(1,2,3,4,5),1,5),c(data1_width,data2_width,data3_width,data4_width,figure_width),1) # data figure
    }
    if(ciFlag){
      #figure_size = graph_width + ci_width
      layout(matrix(c(1,2),1,2),c(figure_width,ci_width),1) # figure ci
    }
  }
  if(ncFlag==0){
    layout(matrix(c(1),1,1),1,1) # figure
  }
  
  
  par(mar=c(4.5, 0, 3.5, 0),xaxs="i",family="sans",xaxt="n",yaxt="n",bty="n")
  
  
  if(nameFlag){
    plot.new()
    strwidth_names_main = strwidth(studynames,units="user",cex=text.cex,font=1,family="sans")
    strwidth_names_title = strwidth("Study",units="user",cex=text.cex,font=2,family="sans")
    strwidth_names_sum = strwidth("Summary",units="user",cex=text.cex,font=2,family="sans")
    
    name_width = max(strwidth_names_main,strwidth_names_title,strwidth_names_sum)
    
    name_adj = switch(nameShow,left=0,right=1,center=0.5)
    xlim = switch(nameShow,left=c(0,name_width),right=c(-name_width,0),center=c(-0.5,0.5)*name_width)
    if(x$misc$covariates.flag){
      plot.window(xlim=xlim,ylim=c(0,nr+1),mar=c(4.5, 1.5, 3.5, 1.5),xlab="",ylab="",xaxs="i",family="sans",xaxt="n",yaxt="n",bty="n")
      text(rep(0,nr),c(nr:1), studynames, adj=c(name_adj,0.5),family="sans",font=1,cex=text.cex)
      text(0,nr+1, "Study", adj=c(name_adj,0.5),family="sans",font=2,cex=text.cex)
    }else{
      if(x$misc$modality.flag){
        nlines = nr + 2*(mod.level-1) + 1
        cuml = 0
        plot.window(xlim=xlim,ylim=c(0,nlines),mar=c(4.5, 1.5, 3.5, 1.5),xlab="",ylab="",xaxs="i",family="sans",xaxt="n",yaxt="n",bty="n")
        for(i in 1:mod.level){
          cuml = cuml + level.length[i]
          text(rep(0,level.length[i]),(nlines-cuml+level.length[i]-1):(nlines-cuml), studynames[ind[[i]]], adj=c(name_adj,0.5),family="sans",font=1,cex=text.cex)
          cuml = cuml + 1
          text(0,nlines-cuml, "Summary", adj=c(name_adj,0.5),family="sans",font=2,cex=text.cex)
          cuml = cuml + 1
          if(i<mod.level){
            text(0,nlines-cuml, paste("Study ",level.name[i+1],sep=""), adj=c(name_adj,0.5),family="sans",font=2,cex=text.cex)
          } 
        }
        text(0,nlines, paste("Study ",level.name[1],sep=""), adj=c(name_adj,0.5),family="sans",font=2,cex=text.cex)
      }else{
        plot.window(xlim=xlim,ylim=c(0,nr+1),mar=c(4.5, 1.5, 3.5, 1.5),xlab="",ylab="",xaxs="i",family="sans",xaxt="n",yaxt="n",bty="n")
        text(rep(0,nr),c(nr:1), studynames, adj=c(name_adj,0.5),family="sans",font=1,cex=text.cex)
        text(0,0,"Summary", adj=c(name_adj,0.5),family="sans",font=2,cex=text.cex)
        text(0,nr+1, "Study", adj=c(name_adj,0.5),family="sans",font=2,cex=text.cex)
      }
    }
  }
  if(dataFlag){
    plot.new()
    
    strwidth_data_main1 = strwidth(datalab1,units="user",cex=text.cex,font=1,family="sans")
    strwidth_data_main2 = strwidth(datalab2,units="user",cex=text.cex,font=1,family="sans")
    strwidth_data_main3 = strwidth(datalab3,units="user",cex=text.cex,font=1,family="sans")
    strwidth_data_main4 = strwidth(datalab4,units="user",cex=text.cex,font=1,family="sans")
    strwidth_data_title1 = strwidth(" TP ",units="user",cex=text.cex,font=2,family="sans")
    strwidth_data_title2 = strwidth(" FP ",units="user",cex=text.cex,font=2,family="sans")
    strwidth_data_title3 = strwidth(" TN ",units="user",cex=text.cex,font=2,family="sans")
    strwidth_data_title4 = strwidth(" FN ",units="user",cex=text.cex,font=2,family="sans")
    
    data1_width = max(strwidth_data_main1,strwidth_data_title1)
    data2_width = max(strwidth_data_main2,strwidth_data_title2)
    data3_width = max(strwidth_data_main3,strwidth_data_title3)
    data4_width = max(strwidth_data_main4,strwidth_data_title4)
    
    data_adj = switch(dataShow,left=0,right=1,center=0.5)
    xlim1 = switch(dataShow,left=c(0,data1_width),right=c(-data1_width,0),center=data1_width*c(-0.5,0.5))
    xlim2 = switch(dataShow,left=c(0,data2_width),right=c(-data2_width,0),center=data2_width*c(-0.5,0.5))
    xlim3 = switch(dataShow,left=c(0,data3_width),right=c(-data3_width,0),center=data3_width*c(-0.5,0.5))
    xlim4 = switch(dataShow,left=c(0,data4_width),right=c(-data4_width,0),center=data4_width*c(-0.5,0.5))
    
    if(x$misc$covariates.flag){
      plot.window(xlim=xlim1,ylim=c(0,nr+1),oma=c(0,0.01,0,0.01),xlab="",ylab="",xaxs="i",family="sans",xaxt="n",yaxt="n",bty="n")
      text(rep(0,nr),c(nr:1), datalab1, adj=c(data_adj,0.5),family="sans",font=1,cex=text.cex)
      text(0,nr+1,"TP", adj=c(data_adj,0.5),family="sans",font=2, cex=text.cex)
      plot(-10,-10,xlim=xlim2,ylim=c(0,nr+1),xlab="",ylab="",xaxs="i",family="sans",xaxt="n",yaxt="n",bty="n")
      text(rep(0,nr),c(nr:1), datalab2, adj=c(data_adj,0.5),family="sans",cex=text.cex)
      text(0,nr+1,"FP", adj=c(data_adj,0.5),family="sans",font=2, cex=text.cex)
      plot(-10,-10,xlim=xlim3,ylim=c(0,nr+1),xlab="",ylab="",xaxs="i",family="sans",xaxt="n",yaxt="n",bty="n")
      text(rep(0,nr),c(nr:1), datalab3, adj=c(data_adj,0.5),family="sans",cex=text.cex)
      text(0,nr+1,"TN", adj=c(data_adj,0.5),family="sans",font=2, cex=text.cex)
      plot(-10,-10,xlim=xlim4,ylim=c(0,nr+1),xlab="",ylab="",xaxs="i",family="sans",xaxt="n",yaxt="n",bty="n")
      text(rep(0,nr),c(nr:1), datalab4, adj=c(data_adj,0.5),family="sans",cex=text.cex)
      text(0,nr+1,"FN", adj=c(data_adj,0.5),family="sans",font=2, cex=text.cex)
    }else{
      if(x$misc$modality.flag){
        nlines = nr + 2*(mod.level-1) + 1
        cuml = 0
        plot.window(xlim=xlim1,ylim=c(0,nlines),xlab="",ylab="",xaxs="i",family="sans",xaxt="n",yaxt="n",bty="n")
        for(i in 1:mod.level){
          cuml = cuml + level.length[i]
          text(rep(0,level.length[i]),(nlines-cuml+level.length[i]-1):(nlines-cuml), datalab1[ind[[i]]], adj=c(data_adj,0.5),family="sans",font=1,cex=text.cex)
          cuml = cuml + 1
          text(0,nlines-cuml, "", adj=c(data_adj,0.5),family="sans",font=2,cex=text.cex)
          cuml = cuml + 1
          text(0,nlines-cuml, "", adj=c(data_adj,0.5),family="sans",font=1,cex=text.cex)
        }
        text(0,nlines,"TP", adj=c(data_adj,0.5),family="sans",font=2, cex=text.cex)
        
        plot.new()
        plot.window(xlim=xlim2,ylim=c(0,nlines),xlab="",ylab="",xaxs="i",family="sans",xaxt="n",yaxt="n",bty="n")
        cuml = 0
        for(i in 1:mod.level){
          cuml = cuml + level.length[i]
          text(rep(0,level.length[i]),(nlines-cuml+level.length[i]-1):(nlines-cuml), datalab2[ind[[i]]], adj=c(data_adj,0.5),family="sans",font=1,cex=text.cex)
          cuml = cuml + 1
          text(0,nlines-cuml, "", adj=c(data_adj,0.5),family="sans",font=2,cex=text.cex)
          cuml = cuml + 1
          text(0,nlines-cuml, "", adj=c(data_adj,0.5),family="sans",font=1,cex=text.cex)
        }
        text(0,nlines,"FP", adj=c(data_adj,0.5),family="sans",font=2, cex=text.cex)
        
        plot.new()
        plot.window(xlim=xlim3,ylim=c(0,nlines),xlab="",ylab="",xaxs="i",family="sans",xaxt="n",yaxt="n",bty="n")
        cuml = 0
        for(i in 1:mod.level){
          cuml = cuml + level.length[i]
          text(rep(0,level.length[i]),(nlines-cuml+level.length[i]-1):(nlines-cuml), datalab3[ind[[i]]], adj=c(data_adj,0.5),family="sans",font=1,cex=text.cex)
          cuml = cuml + 1
          text(0,nlines-cuml, "", adj=c(data_adj,0.5),family="sans",font=2,cex=text.cex)
          cuml = cuml + 1
          text(0,nlines-cuml, "", adj=c(data_adj,0.5),family="sans",font=1,cex=text.cex)
        }
        text(0,nlines,"TN", adj=c(data_adj,0.5),family="sans",font=2, cex=text.cex)
        
        plot.new()
        plot.window(xlim=xlim4,ylim=c(0,nlines),xlab="",ylab="",xaxs="i",family="sans",xaxt="n",yaxt="n",bty="n")
        cuml = 0
        for(i in 1:mod.level){
          cuml = cuml + level.length[i]
          text(rep(0,level.length[i]),(nlines-cuml+level.length[i]-1):(nlines-cuml), datalab4[ind[[i]]], adj=c(data_adj,0.5),family="sans",font=1,cex=text.cex)
          cuml = cuml + 1
          text(0,nlines-cuml, "", adj=c(data_adj,0.5),family="sans",font=2,cex=text.cex)
          cuml = cuml + 1
          text(0,nlines-cuml, "", adj=c(data_adj,0.5),family="sans",font=1,cex=text.cex)
        }
        text(0,nlines,"FN", adj=c(data_adj,0.5),family="sans",font=2, cex=text.cex)
      }else{
        plot.window(xlim=xlim1,ylim=c(0,nr+1),oma=c(0,0.01,0,0.01),xlab="",ylab="",xaxs="i",family="sans",xaxt="n",yaxt="n",bty="n")
        text(rep(0,nr),c(nr:1), datalab1, adj=c(data_adj,0.5),family="sans",font=1,cex=text.cex)
        text(0,nr+1,"TP", adj=c(data_adj,0.5),family="sans",font=2, cex=text.cex)
        plot(-10,-10,xlim=xlim2,ylim=c(0,nr+1),xlab="",ylab="",xaxs="i",family="sans",xaxt="n",yaxt="n",bty="n")
        text(rep(0,nr),c(nr:1), datalab2, adj=c(data_adj,0.5),family="sans",cex=text.cex)
        text(0,nr+1,"FP", adj=c(data_adj,0.5),family="sans",font=2, cex=text.cex)
        plot(-10,-10,xlim=xlim3,ylim=c(0,nr+1),xlab="",ylab="",xaxs="i",family="sans",xaxt="n",yaxt="n",bty="n")
        text(rep(0,nr),c(nr:1), datalab3, adj=c(data_adj,0.5),family="sans",cex=text.cex)
        text(0,nr+1,"TN", adj=c(data_adj,0.5),family="sans",font=2, cex=text.cex)
        plot(-10,-10,xlim=xlim4,ylim=c(0,nr+1),xlab="",ylab="",xaxs="i",family="sans",xaxt="n",yaxt="n",bty="n")
        text(rep(0,nr),c(nr:1), datalab4, adj=c(data_adj,0.5),family="sans",cex=text.cex)
        text(0,nr+1,"FN", adj=c(data_adj,0.5),family="sans",font=2, cex=text.cex)
      }
    }
  }
  
  if(accuracy.type %in% c("lrpos", "lrneg", "rd", "llrpos", "llrneg", "ldor", "dor")){
    xlim = c(xmin,xmax)
  }else{
    xlim = c(max(xmin-0.1,0),min(xmax+0.1,1))
  }
  if(x$misc$covariates.flag){
    plot(-10,-10,xlim=xlim,ylim=c(0,nr+1),xlab="",ylab="",xaxs="i",family="sans",xaxt="s",yaxt="n",bty="n")
    arrows(lb.fit,c(nr:1),ub.fit,c(nr:1),col=arrow.col,angle = 90,length=0.03, code=3, lwd=arrow.lwd, lty=arrow.lty)
    points(est.fit,c(nr:1),pch=p.pch,cex=info,col=p.col)
    if(any(exlb)){
      ind.exlb = which(exlb)
      arrows(rep(xmin,length(ind.exlb)),nr-ind.exlb+1,rep(xmin+0.01*xmin,length(ind.exlb)),nr-ind.exlb+1,angle=12, length=0.1,code=1,col=arrow.col, lwd=arrow.lwd, lty=arrow.lty)
    }
    if(any(exub)){
      ind.exub = which(exub)
      arrows(rep(xmax-0.01*xmax,length(ind.exub)),nr-ind.exub+1,rep(xmax,length(ind.exub)),nr-ind.exub+1,angle=12, length=0.1,code=2,col=arrow.col, lwd=arrow.lwd, lty=arrow.lty)
    }
    axis(1,at=round(c(xmin,xmax),2),labels=round(c(xmin,xmax),2),cex.axis=axis.cex,...)
    
  }else{
    if(x$misc$modality.flag){
      nlines = nr + 2*(mod.level-1) + 1
      cuml = 0
      plot(-10,-10,xlim=xlim,ylim=c(0,nlines),xlab="",ylab="",xaxs="i",family="sans",xaxt="s",yaxt="n",bty="n")
      for(i in 1:mod.level){
        cuml = cuml + level.length[i]
        polygon(sfit[i,c(intervals[1],intervals[2],intervals[2],intervals[1])],c(nlines-cuml-1,nlines-cuml-1,nlines-cuml+level.length[i]-1,nlines-cuml+level.length[i]-1),col=shade.col,angle=45,density=10,border = NA)
        lines(rep(est.sfit[i],2),c(nlines-cuml-1,nlines-cuml+level.length[i]-1),col="darkgray",lwd=2,lty=1)
        lines(rep(lb.sfit[i],2),c(nlines-cuml-1,nlines-cuml+level.length[i]-1),col=shade.col,lwd=2,lty=2)
        lines(rep(ub.sfit[i],2),c(nlines-cuml-1,nlines-cuml+level.length[i]-1),col=shade.col,lwd=2,lty=2)
        
        arrows(lb.fit[ind[[i]]],c((nlines-cuml+level.length[i]-1):(nlines-cuml)),ub.fit[ind[[i]]],c((nlines-cuml+level.length[i]-1):(nlines-cuml)),col=arrow.col,angle = 90,length=0.03, code=3, lwd=arrow.lwd, lty=arrow.lty)
        points(est.fit[ind[[i]]],c((nlines-cuml+level.length[i]-1):(nlines-cuml)),pch=p.pch,cex=info[ind[[i]]],col=p.col)
        
        if(any(exlb[ind[[i]]])){
          ind.exlb = which(exlb[ind[[i]]])
          arrows(rep(xmin,length(ind.exlb)),nlines-cuml+level.length[i]-ind.exlb,rep(xmin+0.01*xmin,length(ind.exlb)),nlines-cuml+level.length[i]-ind.exlb,angle=12, length=0.1,code=1,col=arrow.col, lwd=arrow.lwd, lty=arrow.lty)
        }
        if(any(exub[ind[[i]]])){
          ind.exub = which(exub[ind[[i]]])
          arrows(rep(xmax-0.01*xmax,length(ind.exub)),nlines-cuml+level.length[i]-ind.exub,rep(xmax,length(ind.exub)),nlines-cuml+level.length[i]-ind.exub,angle=12, length=0.1,code=2,col=arrow.col, lwd=arrow.lwd, lty=arrow.lty)
        }
        cuml = cuml + 1
        arrows(lb.sfit[i],nlines-cuml,ub.sfit[i],nlines-cuml,col=arrow.col,angle = 90,length=0.03, code=3, lwd=arrow.lwd, lty=arrow.lty)
        polygon(c(0.5*(est.sfit[i]+lb.sfit[i]),est.sfit[i],0.5*(est.sfit[i]+ub.sfit[i]),est.sfit[i]),c(nlines-cuml,nlines-cuml-0.3,nlines-cuml,nlines-cuml+0.3),col=p.col,border = NA)
        cuml = cuml + 1
        text(0,nlines-cuml, "", family="sans",font=1,cex=text.cex)
      }
      axis(1,at=round(c(xmin,xmax),2),labels=round(c(xmin,xmax),2),cex.axis=axis.cex,...)
    }else{
      plot(-10,-10,xlim=xlim,ylim=c(0,nr+1),xlab="",ylab="",xaxs="i",family="sans",xaxt="s",yaxt="n",bty="n")
      polygon(sfit[c(intervals[1],intervals[2],intervals[2],intervals[1])],c(-1,-1,nr+1,nr+1),col=shade.col,angle=45,density=10,border = NA)
      arrows(lb.fit,c(nr:1),ub.fit,c(nr:1),col=arrow.col,angle = 90,length=0.03, code=3, lwd=arrow.lwd, lty=arrow.lty)
      arrows(lb.sfit,0,ub.sfit,0,col=arrow.col,angle = 90,length=0.03, code=3, lwd=arrow.lwd, lty=arrow.lty)
      points(est.fit,c(nr:1),pch=p.pch,cex=info,col=p.col)
      polygon(c(0.5*(est.sfit+lb.sfit),est.sfit,0.5*(est.sfit+ub.sfit),est.sfit),c(0,-0.3,0,0.3),col=p.col,border = NA)
      if(any(exlb)){
        ind.exlb = which(exlb)
        arrows(rep(xmin,length(ind.exlb)),nr-ind.exlb+1,rep(xmin+0.01*xmin,length(ind.exlb)),nr-ind.exlb+1,angle=12, length=0.1,code=1,col=arrow.col, lwd=arrow.lwd, lty=arrow.lty)
      }
      if(any(exub)){
        ind.exub = which(exub)
        arrows(rep(xmax-0.01*xmax,length(ind.exub)),nr-ind.exub+1,rep(xmax,length(ind.exub)),nr-ind.exub+1,angle=12, length=0.1,code=2,col=arrow.col, lwd=arrow.lwd, lty=arrow.lty)
      }
      abline(v=est.sfit,col="darkgray")
      abline(v=sfit[c(intervals[1],intervals[2])],col=shade.col,lty=2)
      axis(1,at=round(c(xmin,est.sfit,xmax),2),labels=round(c(xmin,est.sfit,xmax),2),cex.axis=axis.cex,...)
    }
  }
  
  if(ciFlag){
    plot.new()
    strwidth_ci_main = strwidth(ci.fit,units="user",cex=text.cex,font=1,family=family)
    if(!x$misc$covariates.flag){
      strwidth_ci_sum = strwidth(ci.sfit,units="user",cex=text.cex,font=2,family=family)
    }
    strwidth_ci_title = strwidth("Estimate",units="user",cex=text.cex,font=2,family="sans")
    if(!x$misc$covariates.flag){
      ci_width = max(strwidth_ci_main,strwidth_ci_title,strwidth_ci_sum)
    }else{
      ci_width = max(strwidth_ci_main,strwidth_ci_title)
    }
    
    
    ci_adj = switch(ciShow,left=0,right=1,center=0.5)
    xlim = switch(ciShow,left=c(0,ci_width),right=c(-ci_width,0),center=ci_width*c(-0.5,0.5))
    if(x$misc$covariates.flag){
      plot.window(xlim=xlim,ylim=c(0,nr+1),xlab="",ylab="",xaxs="i",family="sans",xaxt="n",xaxs="i",yaxt="n",bty="n")
      text(rep(0,nr),c(nr:1),ci.fit,adj=c(ci_adj,0.5),family=family,font=1,cex=text.cex)
      text(0,nr+1,"Estimates",adj=c(ci_adj,0.5),family="sans",font=2,cex=text.cex)
    }else{
      if(x$misc$modality.flag){
        nlines = nr + 2*(mod.level-1) + 1
        cuml = 0
        plot.window(xlim=xlim,ylim=c(0,nlines),xlab="",ylab="",xaxs="i",family="sans",xaxt="n",xaxs="i",yaxt="n",bty="n")
        for(i in 1:mod.level){
          cuml = cuml + level.length[i]
          text(rep(0,level.length[i]),(nlines-cuml+level.length[i]-1):(nlines-cuml), ci.fit[ind[[i]]], adj=c(ci_adj,0.5),family=family,font=1,cex=text.cex)
          cuml = cuml + 1
          text(0,nlines-cuml, ci.sfit[i], adj=c(ci_adj,0.5),family=family,font=2,cex=text.cex)
          cuml = cuml + 1
          text(0,nlines-cuml, "", adj=c(name_adj,0.5),family="sans",font=1,cex=text.cex)
        }
        text(0,nlines,"Estimates",adj=c(ci_adj,0.5),family="sans",font=2,cex=text.cex)
      }else{
        plot.window(xlim=xlim,ylim=c(0,nr+1),xlab="",ylab="",xaxs="i",family="sans",xaxt="n",xaxs="i",yaxt="n",bty="n")
        text(rep(0,nr),c(nr:1),ci.fit,adj=c(ci_adj,0.5),family=family,font=1,cex=text.cex)
        text(0,0,ci.sfit,adj=c(ci_adj,0.5),family=family,font=2,cex=text.cex)
        text(0,nr+1,"Estimates",adj=c(ci_adj,0.5),family="sans",font=2,cex=text.cex)
      }
    }
  }
  mtext(main, adj=0.5, line=-2, cex = main.cex, outer=TRUE) 
  par(op) 
  
  return(invisible())
}



