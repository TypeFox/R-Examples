# Population histogram

# roxygen comments
#' Plots population histogram of the NCA metrics selected for model diagnosis.
#'
#' \pkg{histpop.plot} plots population histogram of the NCA metrics selected 
#' for model diagnosis (e.g. AUClast, AUCINF_obs, Cmax and Tmax).
#'
#' \pkg{histpop.plot} plots histogram of the NCA metrics selected for the model 
#' diagnosis and compares with the corresponding metrics estimated from the 
#' observed data. The allowed NCA metrics for this histograms are "AUClast", 
#' "AUClower_upper", "AUCINF_obs", "AUCINF_pred", "AUMClast", "Cmax", "Tmax" and
#' "HL_Lambda_z". By default, this function produces histogram of AUClast and 
#' Cmax.
#' 
#' @param obsdata Data frame with the values of the NCA metrics estimated from
#'   the observed data
#' @param simdata Data frame with the values of the NCA metrics estimated from
#'   the simulated data
#' @param figlbl Figure label based on dose identifier and/or population
#'   stratifier (\strong{NULL})
#' @param param A character array of the NCA metrics. The allowed NCA metrics 
#'   for this histograms are "AUClast", "AUClower_upper", "AUCINF_obs", 
#'   "AUCINF_pred", "AUMClast", "Cmax", "Tmax" and "HL_Lambda_z". 
#'   (\strong{c("AUClast", "Cmax")})
#' @param cunit Unit for concentration (\strong{"[M].[L]^-3"})
#' @param tunit Unit for time (\strong{"[T]"})
#' @param spread Measure of the spread of simulated data (ppi (95\% parametric
#'   prediction interval) or npi (95\% nonparametric prediction interval))
#'   (\strong{"npi"})
#'
#' @return returns a graphical object created by arrangeGrob function
#' @export
#'

histpop.plot <- function(obsdata=outData,simdata=smeanData,figlbl=NULL,param=c("AUClast","Cmax"),cunit="[M].[L]^-3",tunit="[T]",spread="npi"){
  
  "..density.." <- "TYPE" <- "obs" <- "sim" <- "arrangeGrob" <- "scale_linetype_manual" <- "scale_color_manual" <- "xlab" <- "ylab" <- "guides" <- "guide_legend" <- "theme" <- "element_text" <- "unit" <- "element_rect" <- "geom_histogram" <- "aes" <- "geom_vline" <- "melt" <- "ggplot" <- "labs" <- "coord_cartesian" <- "facet_wrap" <- "gtable_filter" <- "ggplot_gtable" <- "ggplot_build" <- "textGrob" <- "gpar" <- "..count.." <- "..PANEL.." <- "scale_y_continuous" <- "percent" <- "sd" <- "quantile" <- "packageVersion" <- NULL
  rm(list=c("..density..","TYPE","obs","sim","arrangeGrob","scale_linetype_manual","scale_color_manual","xlab","ylab","guides","guide_legend","theme","element_text","unit","element_rect","geom_histogram","aes","geom_vline","melt","ggplot","labs","coord_cartesian","facet_wrap","gtable_filter","ggplot_gtable","ggplot_build","textGrob","gpar","..count..","..PANEL..","scale_y_continuous","percent","sd","quantile","packageVersion"))
  
  outData <- obsdata; smeanData <- simdata
  
  alwprm <- c("AUClast","AUClower_upper","AUCINF_obs","AUCINF_pred","AUMClast","Cmax","Tmax","HL_Lambda_z")
  npr    <- length(param)
  fctNm  <- data.frame()
  nc <- ifelse(npr<2, 1, ifelse(npr>=2 & npr<=6, 2, 3))
  if (!all(param%in%alwprm)){setwd("..");stop("Incorrect NCA metrics. Please select NCA metrics from \"AUClast\", \"AUClower_upper\", \"AUCINF_obs\", \"AUCINF_pred\", \"AUMClast\", \"Cmax\", \"Tmax\", \"HL_Lambda_z\".")}
  
  # ggplot variables
  ggOpt_pop <- list(scale_linetype_manual(name="",values=c("mean(obs)"="solid","mean(meanSim)"="solid","+/-spread"="dashed")),
                    scale_color_manual(name = "", values=c("mean(obs)"="red","mean(meanSim)"="blue","+/-spread"="blue")),
                    xlab(""), ylab(""),
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
                    geom_vline(aes(xintercept=as.numeric(obs), color="mean(obs)", linetype="mean(obs)"), size=1, show_guide=T),
                    geom_vline(aes(xintercept=as.numeric(mean), color="mean(meanSim)", linetype="mean(meanSim)"), size=1),
                    geom_vline(aes(xintercept=as.numeric(sprlow), color="+/-spread", linetype="+/-spread"), size=1),
                    geom_vline(aes(xintercept=as.numeric(sprhgh), color="+/-spread", linetype="+/-spread"), size=1),
                    scale_y_continuous(labels = percent))
  
  obsVal   <- sapply(obsdata, FUN=function(x) mean(as.numeric(x), na.rm=T))
  meanMean <- sapply(simdata, FUN=function(x) mean(as.numeric(x), na.rm=T))
  sdMean   <- sapply(simdata, FUN=function(x) sd(as.numeric(x), na.rm=T))
  xlow     <- sapply(simdata, FUN=function(x) unname(quantile(as.numeric(x),0.01,na.rm=T)))
  xhgh     <- sapply(simdata, FUN=function(x) unname(quantile(as.numeric(x),0.99,na.rm=T)))
  if (spread=="ppi"){
    sprlow <- meanMean-1.96*sdMean
    sprhgh <- meanMean+1.96*sdMean
  }else if (spread=="npi"){
    sprlow <- sapply(simdata, FUN=function(x) unname(quantile(as.numeric(x),0.025,na.rm=T)))
    sprhgh <- sapply(simdata, FUN=function(x) unname(quantile(as.numeric(x),0.975,na.rm=T)))
  }
  
  longData <- melt(simdata,measure=param)
  names(longData) <- c("TYPE","sim")
  longData <- cbind(longData,mean=0,sd=0,sprlow=0,sprhgh=0,obs=0,xlow=0,xhgh=0)
  
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
    longData[longData$TYPE==param[p],"mean"]   <- meanMean[param[p]]
    longData[longData$TYPE==param[p],"sd"]     <- sdMean[param[p]]
    longData[longData$TYPE==param[p],"sprlow"] <- sprlow[param[p]]
    longData[longData$TYPE==param[p],"sprhgh"] <- sprhgh[param[p]]
    longData[longData$TYPE==param[p],"obs"]    <- obsVal[param[p]]
    longData[longData$TYPE==param[p],"xlow"]   <- min(xlow[param[p]],sprlow[param[p]],obsVal[param[p]])
    longData[longData$TYPE==param[p],"xhgh"]   <- max(xhgh[param[p]],sprhgh[param[p]],obsVal[param[p]])
  }
  
  devtag <- ifelse (spread=="ppi","95% parametric prediction interval","95% nonparametric prediction interval")
  gplt <- list()
  for (p in 1:npr){
    df <- subset(longData, TYPE==param[p])
    df$TYPE <- factor(df$TYPE, levels=param[p], labels=fctNm[fctNm$prmNm==param[p],"prmUnit"])
    xl <- df$xlow[1]; xu <- df$xhgh[1]
    gplt[[p]] <- ggplot(df,aes(x=as.numeric(sim))) + ggOpt_pop +
      labs(title=paste("mean(obs)=",format(df$obs[1],digits=2),", mean(meanSim)=",format(df$mean[1],digits=3),"\n+/-spread=(",format(df$sprlow[1],digits=3),",",format(df$sprhgh[1],digits=3),")\n",sep="")) +
      coord_cartesian(xlim=c(xl,xu)) + facet_wrap(~TYPE, scales="free")
  }
  mylegend <- suppressMessages(suppressWarnings(gtable_filter(ggplot_gtable(ggplot_build(gplt[[1]])), "guide-box", trim=T)))
  lheight  <- sum(mylegend$heights)
  for (p in 1:npr){gplt[[p]] <- gplt[[p]] + theme(legend.position="none")}
  
  if(is.null(figlbl)){
    Label <- paste("Histogram of simulated population means\n(spread = ",devtag,")\n\n",sep="")
  }else{
    Label <- paste("Histogram of simulated population means (",figlbl,")\n(spread = ",devtag,")\n\n",sep="")
  }
  
  plot_args <- list(top = textGrob(Label,vjust=1,gp=gpar(fontface="bold",cex = 0.8)),
                    bottom = textGrob("Value\n\n",vjust=1,gp=gpar(fontface="bold",cex = 0.8)),
                    ncol=nc)
  if(packageVersion("gridExtra") < "0.9.2"){
    arg_names <- names(plot_args)
    arg_names <- sub("top","main",arg_names)
    arg_names <- sub("bottom","sub",arg_names)
    names(plot_args) <- arg_names
  }  
  gdr <- suppressMessages(suppressWarnings(do.call(arrangeGrob,c(gplt,plot_args))))
  #grid.arrange(gdr)
  
  histpopgrob <- list(gdr=gdr,legend=mylegend,lheight=lheight)
  return(histpopgrob)
}

