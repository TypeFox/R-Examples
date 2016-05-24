# =====================================================
# Finding the best performing workflows
# =====================================================
# Luis Torgo, Nov 2013
#
topPerformers <- function(compRes,maxs=rep(FALSE,dim(compRes[[1]][[1]]@iterationsScores)[2]),digs=3) {
  if (!inherits(compRes,'ComparisonResults')) stop(compRes,' needs to be of class "ComparisonResults".\n')
  if (length(maxs) == 1) maxs <- rep(maxs,dim(compRes[[1]][[1]]@iterationsScores)[2])
  else if (length(maxs) != dim(compRes[[1]][[1]]@iterationsScores)[2]) stop('"maxs" needs to have the same size as the number of evaluation statistics.\n')

  avgScs <- .statScores(compRes,'avg')
  bs <- lapply(avgScs,function(t) {
      r <- data.frame(Workflow=rep('',NROW(t)),Estimate=rep(0,NROW(t)),stringsAsFactors=F,row.names=rownames(t))
      for(i in 1:NROW(t))
          r[i,] <- if (maxs[i]) c(colnames(t)[which.max(t[i,])],round(max(t[i,],na.rm=TRUE),digs)) else c(colnames(t)[which.min(t[i,])],round(min(t[i,],na.rm=TRUE),digs))
      r
  })
  bs
}

# =====================================================
# Finding the best performing workflow for a given task and metric
# =====================================================
# Luis Torgo, Nov 2013
#
topPerformer <- function(compRes,metric,task,max=FALSE) {
  if (!inherits(compRes,'ComparisonResults'))
      stop(compRes,' needs to be of class "ComparisonResults".\n')
  if (!(metric %in% metricNames(compRes)))
      stop(metric,' not estimated in this "ComparisonResults" object.\n')
  if (!(task %in% taskNames(compRes)))
      stop(task,' not used in this "ComparisonResults" object.\n')

  avgScs <- .statScores(compRes,'avg')
  scs <- avgScs[[task]]
  w <- if (max) colnames(scs)[which.max(scs[metric,])] else colnames(scs)[which.min(scs[metric,])]
  getWorkflow(w,compRes)
}

# =====================================================
# Obtaining a ranking of the workflows
# =====================================================
# Luis Torgo, Nov 2013
#
rankWorkflows <- function(compRes,top=min(5,length(workflowNames(compRes))),maxs=rep(FALSE,dim(compRes[[1]][[1]]@iterationsScores)[2])) {
  if (!inherits(compRes,'ComparisonResults')) stop(compRes,' needs to be of class "ComparisonResults".\n')
  if (length(maxs) == 1) maxs <- rep(maxs,dim(compRes[[1]][[1]]@iterationsScores)[2])
  else if (length(maxs) != dim(compRes[[1]][[1]]@iterationsScores)[2]) stop('"maxs" needs to have the same size as the number of evaluation statistics.\n')

  avgScs <- .statScores(compRes,'avg')
  lapply(avgScs,function(t) {
      l <- vector("list",nrow(t))
      for(i in 1:nrow(t))
          l[[i]] <- data.frame(Workflow=colnames(t)[order(t[i,],decreasing=maxs[i])],Estimate=sort(t[i,],decreasing=maxs[i]),stringsAsFactors=F,row.names=1:ncol(t))[1:top,]
      names(l) <- rownames(t)
      l
  })
}
      
# =====================================================
# Apply a summary function over the individual iteration scores
# =====================================================
# Luis Torgo, Nov 2013
# =====================================================
# 
metricsSummary <- function(compRes,summary='mean',...) {
  if (!inherits(compRes,'ComparisonResults'))
    stop(compRes,' needs to be of class "ComparisonResults".\n')
  lapply(compRes,function(t)
         sapply(t,function(w)
                apply(w@iterationsScores,2,function(x) do.call(summary,list(x,...)))
                ))
}



# ======================================================================
# Construction of comparative analysis tables based on the results of 
# comparative experiment obtained with experimentalComparison() (i.e. based
# on a compExp object).
# The first argument is the compExp object that resulted from the 
# experiments. Then we have the system against which all remaning are
# compared to (defaults to first in the structure). Finally we can
# provide a vector of the names of the statistics we which to get a
# table (defaults to all).
# =====================================================
# Luis Torgo, Jan-Aug 2009, 2014
# =====================================================
pairedComparisons <-  function(obj,baseline,
                               maxs=rep(FALSE,length(metricNames(obj))),
                               p.value=0.05) {
    
    if (!inherits(obj,'ComparisonResults')) stop(obj,' is not of class "ComparisonResults".\n')

    ## basic information on the expriment
    ts <- taskNames(obj);     nts <- length(ts)
    ws <- workflowNames(obj); nws <- length(ws)
    ms <- metricNames(obj);   nms <- length(ms)

    if (nws < 2) stop("Paired comparisons only make sense with more than one workflow!")

    if (!missing(baseline)) {
        other <- setdiff(ws,baseline)
        pb <- which(ws==baseline)
    }

    ## this list will hold the results of the comparisons, one for each metric
    compResults <- vector("list",nms)
    names(compResults) <- ms

    ## Constructing the tests information for each metric
    for(p in 1:nms) {
        compResults[[p]] <- list()
        compResults[[p]]$setup <- list(nTasks=nts,nWorkflows=nws)
        compResults[[p]]$avgScores <- t(sapply(obj,function(m) sapply(m,function(i) mean(i@iterationsScores[,p],na.rm=TRUE))))
        compResults[[p]]$medScores <- t(sapply(obj,function(m) sapply(m,function(i) median(i@iterationsScores[,p],na.rm=TRUE))))
        compResults[[p]]$rks <- t(apply(if (maxs[p]) -compResults[[p]]$avgScores else compResults[[p]]$avgScores,1,rank))
        compResults[[p]]$avgRksWFs <- apply(compResults[[p]]$rks,2,mean)
        

        if (missing(baseline)) { # using the first workflow as baseline if none indicated
            base <- compResults[[p]]$baseline <- names(which.min(compResults[[p]]$avgRksWFs))
            other <- setdiff(ws,base)
            pb <- which(ws==base)
        } else base <- compResults[[p]]$baseline <- baseline
        
        ## Wilcoxon Signed Rank and t-Student tests
        compResults[[p]]$t.test <- array(NA,dim=c(nws,3,nts),
                 dimnames=list(c(base,other),c("AvgScore","DiffAvgScores","p.value"),
                     ts)
                         )
        compResults[[p]]$WilcoxonSignedRank.test <- array(NA,dim=c(nws,3,nts),
                 dimnames=list(c(base,other),c("MedScore","DiffMedScores","p.value"),
                     ts)
                         )
        for(t in ts) {
            compResults[[p]]$WilcoxonSignedRank.test[base,,t] <- NA
            compResults[[p]]$WilcoxonSignedRank.test[base,"MedScore",t] <- compResults[[p]]$medScores[t,base]
            compResults[[p]]$t.test[base,,t] <- NA
            compResults[[p]]$t.test[base,"AvgScore",t] <- compResults[[p]]$avgScores[t,base]
            for(o in other) {
                tst <- try(wilcox.test(obj[[t]][[o]]@iterationsScores[,p],
                                       obj[[t]][[base]]@iterationsScores[,p],
                                       paired=T))
                compResults[[p]]$WilcoxonSignedRank.test[o,"DiffMedScores",t] <- compResults[[p]]$medScores[t,base] - compResults[[p]]$medScores[t,o]
                compResults[[p]]$WilcoxonSignedRank.test[o,"MedScore",t] <- compResults[[p]]$medScores[t,o]
                compResults[[p]]$WilcoxonSignedRank.test[o,"p.value",t] <-
                    if (inherits(tst,"try-error"))  NA else tst$p.value
                tst <- try(t.test(obj[[t]][[o]]@iterationsScores[,p],
                                  obj[[t]][[base]]@iterationsScores[,p],
                                  paired=T))
                compResults[[p]]$t.test[o,"DiffAvgScores",t] <- compResults[[p]]$avgScores[t,base] - compResults[[p]]$avgScores[t,o]
                compResults[[p]]$t.test[o,"AvgScore",t] <- compResults[[p]]$avgScores[t,o]
                compResults[[p]]$t.test[o,"p.value",t] <-
                    if (inherits(tst,"try-error"))  NA else tst$p.value
            }
        }

        if (nts > 1) {
            ## Testing the null hypothesis that all WFs are equivalent
            chi <- 12*nts/(nws*(nws+1)) * (sum(compResults[[p]]$avgRksWFs^2) - (nws*(nws+1)^2)/4)
            FF <- (nts-1)*chi / (nts*(nws-1) - chi)
            critVal <- df(1-p.value,nws-1,(nws-1)*(nts-1))
            rejNull <- FF > critVal
            compResults[[p]]$F.test <- list(chi=chi,FF=FF,critVal=critVal,rejNull=rejNull)
            
            compResults[[p]]$Nemenyi.test <- NULL
            compResults[[p]]$BonferroniDunn.test <- NULL
            if (rejNull) {
                ## Nemenyi critical difference
                CD.n <- qtukey(1-p.value,nws,1e06)/sqrt(2)*sqrt(nws*(nws+1)/(6*nts))
                allRkDifs <- outer(compResults[[p]]$avgRksWFs,compResults[[p]]$avgRksWFs,
                                   function(x,y) abs(x-y))
                signifDifs <- allRkDifs >= CD.n
                compResults[[p]]$Nemenyi.test <- list(critDif=CD.n,
                                                      rkDifs=allRkDifs,
                                                      signifDifs=signifDifs)
                
                ## Bonferroni-Dunn test against the base
                
                ## Bonferroni-Dunn critical difference
                CD.bd <- qtukey(1-(p.value/(nws-1)),2,1e06)/sqrt(2)*sqrt(nws*(nws+1)/(6*nts))
                diffs2base <- abs(compResults[[p]]$avgRksWFs[-pb]-compResults[[p]]$avgRksWFs[pb])
                signifDifs <- diffs2base >= CD.bd
                compResults[[p]]$BonferroniDunn.test <- list(critDif=CD.bd,
                                                             baseline=base,
                                                             rkDifs=diffs2base,
                                                             signifDifs=signifDifs)
                
            }
        } else {
            compResults[[p]]$F.test <- compResults[[p]]$Nemenyi.test <- compResults[[p]]$BonferroniDunn.test <- NULL
            warning("With less 2 tasks the Friedman, Nemenyi and Bonferroni-Dunn tests are not calculated.")
        }
    }
    compResults
}



# ======================================================================
# Filtering the results of a call to pairedComparisons by a minimum p.value
# =====================================================
# Luis Torgo, Jan-Aug 2009, 2014
# =====================================================
signifDiffs <- function(ps,p.limit=0.05,metrics=names(ps),tasks=rownames(ps[[1]]$avgScores)) {
    res <- vector("list",length(metrics))
    names(res) <- metrics
    for(p in metrics) {
        res[[p]] <- list()
        res[[p]]$WilcoxonSignedRank.test <- vector("list",length(tasks))
        res[[p]]$t.test <- vector("list",length(tasks))
        names(res[[p]]$WilcoxonSignedRank.test) <- names(res[[p]]$t.test) <- tasks
        for(t in tasks) {
            res[[p]]$WilcoxonSignedRank.test[[t]] <- ps[[p]]$WilcoxonSignedRank.test[which(ps[[p]]$WilcoxonSignedRank.test[,"p.value",t] < p.limit | is.na(ps[[p]]$WilcoxonSignedRank.test[,"p.value",t]) & !is.nan(ps[[p]]$WilcoxonSignedRank.test[,"p.value",t])),,t]
            res[[p]]$t.test[[t]] <- ps[[p]]$t.test[which(ps[[p]]$t.test[,"p.value",t] < p.limit | is.na(ps[[p]]$t.test[,"p.value",t]) & !is.nan(ps[[p]]$t.test[,"p.value",t])),,t]
        }
    }
    res
}


# ======================================================================
# A CD diagram for the Nemenyi test
# =====================================================
# Luis Torgo, Nov, 2014
# =====================================================
CDdiagram.Nemenyi <- function(r,metric=names(r)[1]) {
    if (is.null(r[[metric]]$F.test) | is.null(r[[metric]]$Nemenyi.test)) stop("Results of both the F and Nemenyi tests are required for these diagrams.")
    o <- rank(r[[metric]]$avgRksWFs,ties.method="first")
    mxl <- ceiling(length(r[[metric]]$avgRksWFs)/2)
    data <- data.frame(avgRk=r[[metric]]$avgRksWFs,
                       invRk=length(r[[metric]]$avgRksWFs)+1-r[[metric]]$avgRksWFs,
                       sys=names(r[[metric]]$avgRksWFs),
                       line=o%%mxl + ifelse(o%%mxl==0,mxl,0) ,
#                   line=o%%mxl + 1,
                       side=ifelse(o <= mxl,-1,1)
                       )
    data$line <- ifelse(data$side==1,mxl+1-data$line,data$line)
    len <- length(r[[metric]]$avgRksWFs)
    cd <- r[[metric]]$Nemenyi.test$critDif
    g <- ggplot2::ggplot(data,ggplot2::aes_string(x="invRk",y="line")) + #geom_point() +
            ggplot2::geom_segment(ggplot2::aes_string(x="invRk",y=0,xend="invRk",yend="line",col="sys")) +
            ggplot2::geom_text(data=data[data$side==-1,],
                      ggplot2::aes_string(label = "sys", x = +Inf, y = "line",col="sys"),
                               hjust = 0,size=4) +
            ggplot2::geom_text(data=data[data$side==1,],
                      ggplot2::aes_string(label = "sys", x = -Inf, y = "line",col="sys"),
                      hjust = 1,size=4) +
            ggplot2::geom_segment(data=data[data$side==-1,],ggplot2::aes_string(x=max(o),y="line",xend="invRk",yend="line",col="sys")) +
            ggplot2::geom_segment(data=data[data$side==1,],ggplot2::aes_string(x="invRk",y="line",xend=0,yend="line",col="sys")) +
            ggplot2::scale_x_continuous(limits=c(0,len),
                               breaks=0:(len+1),
                               labels=paste((len+1):0)) +
            ggplot2::scale_y_continuous(limits=c(0,mxl+1)) +
            ggplot2::xlab("Average Rank") + ggplot2::ylab("") +
            ggplot2::theme(axis.text.y=ggplot2::element_blank(),
                  axis.ticks.y=ggplot2::element_blank(),
                  axis.title.y=ggplot2::element_blank(),
                  legend.position="none",
                  panel.background=ggplot2::element_blank(),
                  panel.border=ggplot2::element_blank(),
                  axis.line=ggplot2::element_line(size=1),
                  axis.line.y=ggplot2::element_blank(),
                  panel.grid.major=ggplot2::element_blank(),
                  plot.background=ggplot2::element_blank(),
                  plot.margin = grid::unit(c(3,6,1,5), "lines")
                  ) +
            ggplot2::coord_fixed(ratio=0.5) + 
            ggplot2::annotate("segment",x=0,xend=cd,
                     y=mxl+1,yend=mxl+1,size=1.5) +
            ggplot2::annotate("text",x=0,y=mxl+1,label=paste("Critical Difference",round(cd,1),sep=" = "),vjust=-0.5,hjust=0,size=3)

    wfsOrd <- names(r[[metric]]$avgRksWFs[order(r[[metric]]$avgRksWFs)])
    ss <- r[[metric]]$Nemenyi.test$signifDifs[wfsOrd,wfsOrd]
    mx <- ncol(ss)
    pos <- rep(1:mxl,each=mxl) - seq(0,1,by=1/(mxl+1))[-c(1,mxl+2)]
    frees <- rep(1,mxl)
    from <- 1
    currTill <- 0
    for(i in 1:nrow(ss)) {
        till <- which(ss[i,i:mx])
        till <- if (length(till)) till-1 else mx
        if (till > currTill) {
            theLine <- min(data[wfsOrd[i],"line"],data[wfsOrd[till],"line"])
            ypos <- pos[(theLine-1)*mxl+frees[theLine]]
            frees[theLine] <- frees[theLine]+1
            g <- g + ggplot2::annotate("segment",
                              x=data[wfsOrd[till],"invRk"]-.1,
                              xend=data[wfsOrd[i],"invRk"]+.1,
                              y=ypos, yend=ypos, size=1.2)
            currTill <- till
        }
        if (till == mx) break
    }
    grid::grid.newpage()
    gt <- ggplot2::ggplot_gtable(ggplot2::ggplot_build(g))
    gt$layout$clip[gt$layout$name == "panel"] <- "off"
    grid::grid.draw(gt)
}



# ======================================================================
# A CD diagram for the Boferroni-Dunn test
# =====================================================
# Luis Torgo, Nov, 2014
# =====================================================
CDdiagram.BD <- function(r,metric=names(r)[1]) {
    if (is.null(r[[metric]]$F.test) | is.null(r[[metric]]$BonferroniDunn.test)) stop("Results of both the F and Bonferroni-Dunn tests are required for these diagrams.")
    o <- rank(r[[metric]]$avgRksWFs,ties.method="first")
    mxl <- ceiling(length(r[[metric]]$avgRksWFs)/2)
    data <- data.frame(avgRk=r[[metric]]$avgRksWFs,
                       invRk=length(r[[metric]]$avgRksWFs)+1-r[[metric]]$avgRksWFs,
                       sys=names(r[[metric]]$avgRksWFs),
                       line=o%%mxl + ifelse(o%%mxl==0,mxl,0) ,
                       side=ifelse(o <= mxl,-1,1)
                       )
    data$line <- ifelse(data$side==1,mxl+1-data$line,data$line)
#    data$color <- rep("black",nrow(data))
    data$face <- rep(1,nrow(data))
    data[r[[metric]]$BonferroniDunn.test$baseline,"face"] <- 2
    data[names(which(r[[metric]]$BonferroniDunn.test$signifDifs)),"face"] <- 3
#    data[r[[metric]]$BonferroniDunn.test$baseline,"color"] <- "green"
#    data[names(which(r[[metric]]$BonferroniDunn.test$signifDifs)),"color"] <- "red"
    len <- length(r[[metric]]$avgRksWFs)
    cd <- r[[metric]]$BonferroniDunn.test$critDif
    baseScore <- data[r[[metric]]$BonferroniDunn.test$baseline,"invRk"]
    g <- ggplot2::ggplot(data,ggplot2::aes_string(x="invRk",y="line")) + #geom_point() +
            ggplot2::geom_segment(ggplot2::aes_string(x="invRk",y=0,xend="invRk",yend="line",col="sys")) +
            ggplot2::geom_text(data=data[data$side==-1,],
                      ggplot2::aes_string(label = "sys", x = +Inf, y = "line",col="sys",
 #                                fontface="face",colour="color"),
                                 fontface="face"),
                      hjust = 0,size=4) +
            ggplot2::geom_text(data=data[data$side==1,],
                      ggplot2::aes_string(label = "sys", x = -Inf, y = "line",col="sys",
#                          fontface=face,colour=color),
                          fontface="face"),
                      hjust = 1,size=4) +
            ggplot2::geom_segment(data=data[data$side==-1,],ggplot2::aes_string(x=max(o),y="line",xend="invRk",yend="line",col="sys")) +
            ggplot2::geom_segment(data=data[data$side==1,],ggplot2::aes_string(x="invRk",y="line",xend=0,yend="line",col="sys")) +
            ggplot2::scale_x_continuous(limits=c(0,len),
                               breaks=0:(len+1),
                               labels=paste((len+1):0)) +
            ggplot2::scale_y_continuous(limits=c(0,mxl+1)) +
            ggplot2::xlab("Average Rank") + ggplot2::ylab("") +
            ggplot2::theme(axis.text.y=ggplot2::element_blank(),
                  axis.ticks.y=ggplot2::element_blank(),
                  axis.title.y=ggplot2::element_blank(),
                  legend.position="none",
                  panel.background=ggplot2::element_blank(),
                  panel.border=ggplot2::element_blank(),
                  axis.line=ggplot2::element_line(size=1),
                  axis.line.y=ggplot2::element_blank(),
                  panel.grid.major=ggplot2::element_blank(),
                  plot.background=ggplot2::element_blank(),
                  plot.margin = grid::unit(c(1,6,1,5), "lines")
                  ) +
            ggplot2::coord_fixed(ratio=0.5) + 
#            ggplot2::annotate("segment",x=max(data$invRk),xend=max(data$invRk)-cd,
#                     y=0,yend=0,size=2) 
            ggplot2::annotate("segment",x=max(0,baseScore-cd),xend=min(baseScore+cd),
                     y=0,yend=0,size=2)  +
            ggplot2::annotate("text",x=-Inf,y=+Inf,label=paste0("Critical Difference = ",round(cd,1),"; Baseline = ",r[[metric]]$BonferroniDunn.test$baseline),vjust=0,hjust=0,size=3)


    grid::grid.newpage()
    gt <- ggplot2::ggplot_gtable(ggplot2::ggplot_build(g))
    gt$layout$clip[gt$layout$name == "panel"] <- "off"
    grid::grid.draw(gt)
}
