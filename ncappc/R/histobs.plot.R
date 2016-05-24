# Plot histogram of four NCA metrics

# roxygen comments
#' Plots histogram of selected set of NCA metrics.
#'
#' \pkg{histobs.plot} plots histogram of selected set of NCA metrics (e.g.
#' AUClast, AUCINF_obs, Cmax and Tmax).
#'
#' \pkg{histobs.plot} plots histogram of selected set of NCA metrics. The 
#' allowed NCA metrics for this histograms are "AUClast", "AUClower_upper", 
#' "AUCINF_obs", "AUCINF_pred", "AUMClast", "Cmax", "Tmax" and "HL_Lambda_z". By
#' default, this function produces histogram of AUClast, AUCINF_obs, Cmax and 
#' Tmax.
#' 
#' @param plotData A data frame with the estimated NCA metrics
#' @param figlbl Figure label based on dose identifier and/or population
#'   stratifier (\strong{NULL})
#' @param param A character array of the NCA metrics. The allowed NCA metrics
#'   for this histograms are "AUClast", "AUClower_upper", "AUCINF_obs",
#'   "AUCINF_pred", "AUMClast", "Cmax", "Tmax" and "HL_Lambda_z".
#'   (\strong{c("AUClast", "AUCINF_obs", "Cmax", "Tmax")})
#' @param cunit Unit for concentration (\strong{"[M].[L]^-3"})
#' @param tunit Unit for time (\strong{"[T]"})
#' @param spread Measure of the spread of simulated data (ppi (95\% parametric 
#'   prediction interval) or npi (95\% nonparametric prediction interval)) 
#'   (\strong{"npi"})
#'
#' @return returns a graphical object created by arrangeGrob function
#' @export
#'

histobs.plot <- function(plotData,figlbl=NULL,param=c("AUClast","AUCINF_obs","Cmax","Tmax"),cunit="[M].[L]^-3",tunit="[T]",spread="npi"){
  
  "..density.." <- "TYPE" <- "Obs" <- "arrangeGrob" <- "scale_linetype_manual" <- "scale_color_manual" <- "xlab" <- "ylab" <- "guides" <- "guide_legend" <- "theme" <- "element_text" <- "unit" <- "element_rect" <- "geom_histogram" <- "aes" <- "geom_vline" <- "melt" <- "ggplot" <- "coord_cartesian" <- "facet_grid" <- "labs" <- "gtable_filter" <- "ggplot_gtable" <- "ggplot_build" <- "textGrob" <- "gpar" <- "..count.." <- "..PANEL.." <- "scale_y_continuous" <- "percent" <- "sd" <- "quantile" <- "na.omit" <- "packageVersion" <- NULL
  rm(list=c("..density..","TYPE","Obs","arrangeGrob","scale_linetype_manual","scale_color_manual","xlab","ylab","guides","guide_legend","theme","element_text","unit","element_rect","geom_histogram","aes","geom_vline","melt","ggplot","coord_cartesian","facet_grid","labs","gtable_filter","ggplot_gtable","ggplot_build","textGrob","gpar","..count..","..PANEL..","scale_y_continuous","percent","sd","quantile","na.omit","packageVersion"))
  
  ggOpt_obs <- list(scale_linetype_manual(name="",values=c("mean(obs)"="solid","+/-spread"="dashed")),
                    scale_color_manual(name = "", values=c("mean(obs)"="blue","+/-spread"="blue")),
                    xlab(""), ylab(""),
                    scale_y_continuous(labels = percent),
                    guides(fill = guide_legend(override.aes = list(linetype = 0 )), shape = guide_legend(override.aes = list(linetype = 0))),
                    theme(plot.title = element_text(size=9, face="bold"),
                          plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm"),
                          axis.title.x = element_text(size=9,face="bold"),
                          axis.title.y = element_text(size=9,face="bold"),
                          axis.text.x  = element_text(size=9,face="bold",color="black",angle=45,vjust=1,hjust=1),
                          axis.text.y  = element_text(size=9,face="bold",color="black",hjust=0),
                          legend.position = "bottom", legend.direction = "horizontal",
                          legend.background = element_rect(),
                          legend.key.size = unit(0.8, "cm"),
                          legend.text  = element_text(size=8,face="bold"),
                          strip.text.x = element_text(size=8, face="bold")),
                    geom_histogram(aes(y=(..count..)/tapply(..count..,..PANEL..,sum)[..PANEL..]), size=0.6, color="black", fill="white"),
                    geom_vline(aes(xintercept=as.numeric(meanObs), color="mean(obs)", linetype="mean(obs)"), size=1, show_guide=T),
                    geom_vline(aes(xintercept=as.numeric(sprlow), color="+/-spread", linetype="+/-spread"), size=1),
                    geom_vline(aes(xintercept=as.numeric(sprhgh), color="+/-spread", linetype="+/-spread"), size=1))
  
  alwprm <- c("AUClast","AUClower_upper","AUCINF_obs","AUCINF_pred","AUMClast","Cmax","Tmax","HL_Lambda_z")
  npr    <- length(param)
  fctNm  <- data.frame()
  nc <- ifelse(npr<2, 1, ifelse(npr>=2 & npr<=6, 2, 3))
  if (!all(param%in%alwprm)){stop("Incorrect NCA metrics. Please select NCA metrics from \"AUClast\", \"AUClower_upper\", \"AUCINF_obs\", \"AUCINF_pred\", \"AUMClast\", \"Cmax\", \"Tmax\", \"HL_Lambda_z\".")}
  for (p in 1:npr){
    if (param[p] == "AUClast" | param[p] == "AUClower_upper" | param[p] == "AUCINF_obs" | param[p] == "AUCINF_pred"){
      fctNm <- rbind(fctNm, data.frame(prmNm=param[p],prmUnit=paste(param[p]," (",cunit,"*",tunit,")",sep="")))
    }else if (param[p] == "AUMClast"){
      fctNm <- rbind(fctNm, data.frame(prmNm=param[p],prmUnit=paste(param[p]," (",cunit,"*",tunit,"^2)",sep="")))
    }else if (param[p] == "Cmax"){
      fctNm <- rbind(fctNm, data.frame(prmNm=param[p],prmUnit=paste(param[p]," (",cunit,")",sep="")))
    }else if (param[p] == "Tmax"){
      fctNm <- rbind(fctNm, data.frame(prmNm=param[p],prmUnit=paste(param[p]," (",tunit,")",sep="")))
    }else if (param[p] == "HL_Lambda_z"){
      fctNm <- rbind(fctNm, data.frame(prmNm=param[p],prmUnit=paste(param[p]," (",tunit,")",sep="")))
    }
  }
  
  devtag <- "2.5th and 97.5th percentile boundaries"
  
  meanObs  <- sapply(plotData, FUN=function(x) mean(as.numeric(x), na.rm=T))
  sdObs    <- sapply(plotData, FUN=function(x) sd(as.numeric(x), na.rm=T))
  xlow     <- sapply(plotData, FUN=function(x) unname(quantile(as.numeric(x),0.01, na.rm=T)))
  xhgh     <- sapply(plotData, FUN=function(x) unname(quantile(as.numeric(x),0.99, na.rm=T)))
  if (spread=="ppi"){
    sprlow <- meanObs-1.96*sdObs
    sprhgh <- meanObs+1.96*sdObs
  }else if (spread=="npi"){
    sprlow <- sapply(plotData, FUN=function(x) unname(quantile(as.numeric(x),0.025, na.rm=T)))
    sprhgh <- sapply(plotData, FUN=function(x) unname(quantile(as.numeric(x),0.975, na.rm=T)))
  }
  longData <- melt(plotData, measure.vars = param)
  names(longData) <- c("TYPE","Obs")
  longData <- cbind(longData,meanObs=0,sdObs=0,sprlow=0,sprhgh=0,xlow=0,xhgh=0)
  for (p in 1:length(param)){
    longData[longData$TYPE==param[p],"meanObs"] <- meanObs[param[p]]
    longData[longData$TYPE==param[p],"sdObs"]   <- sdObs[param[p]]
    longData[longData$TYPE==param[p],"sprlow"]  <- sprlow[param[p]]
    longData[longData$TYPE==param[p],"sprhgh"]  <- sprhgh[param[p]]
    longData[longData$TYPE==param[p],"xlow"]    <- xlow[param[p]]
    longData[longData$TYPE==param[p],"xhgh"]    <- xhgh[param[p]]
  }
  longData <- na.omit(longData)
  longData$TYPE <- factor(longData$TYPE, levels=unique(longData$TYPE), labels=unique(longData$TYPE))
  param <- as.character(unique(longData$TYPE))
  fctNm <- fctNm[fctNm$prmNm%in%param,]
  
  gplt <- list()
  for (p in 1:length(param)){
    df <- subset(longData, TYPE==param[p])
    df$TYPE <- factor(df$TYPE, levels=param[p], labels=fctNm[fctNm$prmNm==param[p],"prmUnit"])
    xl <- df$xlow[1]; xu <- df$xhgh[1]
    gplt[[p]] <- ggplot(df,aes(x=as.numeric(Obs))) + ggOpt_obs + coord_cartesian(xlim=c(xl,xu)) + facet_grid(~TYPE, scales="free") +
      labs(title=paste("mean(obs)=",format(df$meanObs[1],digits=3),"\n+/-spread=(",format(df$sprlow[1],digits=3),",",format(df$sprhgh[1],digits=3),")\n",sep=""))
  }
  mylegend <- suppressMessages(suppressWarnings(gtable_filter(ggplot_gtable(ggplot_build(gplt[[1]])), "guide-box", trim=T)))
  lheight  <- sum(mylegend$heights)
  for (p in 1:length(param)){gplt[[p]] <- gplt[[p]] + theme(legend.position="none")}
  
  if(is.null(figlbl)){
    Label <- paste("Histogram of NCA metrics estimated from the observed data\n(spread = ",devtag,")\n\n",sep="")
  }else{
    Label <- paste("Histogram of NCA metrics estimated from the observed data (",figlbl,")\n(spread = ",devtag,")\n\n",sep="")
  }
  
  plot_args <- list(top = textGrob(Label,vjust=1,gp=gpar(fontface="bold",cex = 0.7)),
                    bottom = textGrob("Value\n\n",vjust=1,gp=gpar(fontface="bold",cex = 0.7)),
                    ncol=nc)
  if(packageVersion("gridExtra") < "0.9.2"){
    arg_names <- names(plot_args)
    arg_names <- sub("top","main",arg_names)
    arg_names <- sub("bottom","sub",arg_names)
    names(plot_args) <- arg_names
  }  
  gdr <- suppressMessages(suppressWarnings(do.call(arrangeGrob,c(gplt,plot_args))))
  #grid.arrange(gdr)
  histobsgrob <- list(gdr=gdr,legend=mylegend,lheight=lheight)
  return(histobsgrob)
}

