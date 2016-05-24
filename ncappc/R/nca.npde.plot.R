# NPDE plot
# roxygen comments
#' Plots population histogram of the NCA metrics selected for model diagnosis.
#'
#' \pkg{nca.npde.plot} plots individual NPDE values and histogram of the NPDE
#' values within a population group
#'
#' \pkg{nca.npde.plot} individual NPDE values and histogram of the NPDE
#' values of NCA metrics within a population group.
#' 
#' @param plotdata Data frame with the values of the NPDE values of each
#'   individual for the NCA metrics
#' @param xvar Name of the independent variable column against which NPDE values
#'   will be plotted
#' @param npdecol Column names or column numbers of containing the NPDE values
#' @param figlbl Figure label based on dose identifier and/or population
#'   stratifier (\strong{NULL})
#' @param cunit Unit for concentration (\strong{"M.L^-3"})
#' @param tunit Unit for time (\strong{"T"})
#'
#' @return returns a data frame with the mean and SD of population NPDE values
#'   of each NCA metric and two graphical objects created by arrangeGrob
#'   function for the individual and population histogram of the NPDE values
#' @export
#'

nca.npde.plot <- function(plotdata,xvar=NULL,npdecol=NULL,figlbl=NULL,cunit="[M].[L]^-3",tunit="[T]"){
  
  "npde" <- "type" <- "mcil" <- "mciu" <- "sdu" <- "sducil" <- "sduciu" <- "..density.." <- "sdl" <- "XVAR" <- "melt" <- "xlab" <- "ylab" <- "theme" <- "element_text" <- "unit" <- "geom_point" <- "facet_wrap" <- "scale_linetype_manual" <- "scale_color_manual" <- "guides" <- "guide_legend" <- "element_rect" <- "geom_histogram" <- "aes" <- "geom_vline" <- "ggplot" <- "labs" <- "..count.." <- "..PANEL.." <- "na.omit" <- "sd" <- "qt" <- "qchisq" <- "scale_y_continuous" <- "percent" <- "sdcil" <- "sdciu" <- NULL
  rm(list=c("npde","type","mcil","mciu","sdu","sducil","sduciu","..density..","sdl","XVAR","melt","xlab","ylab","theme","element_text","unit","geom_point","facet_wrap","scale_linetype_manual","scale_color_manual","guides","guide_legend","element_rect","geom_histogram","aes","geom_vline","ggplot","labs","..count..","..PANEL..","na.omit","sd","qt","qchisq","scale_y_continuous","percent","sdcil","sdciu"))
  
  if (!is.data.frame(plotdata)) stop("plotdata must be a data frame.")
  if (is.null(xvar)){
    stop("Missing x-variable against which the NPDE data will be plotted.")
  }else{
    names(plotdata)[which(names(plotdata)==xvar)] <- "XVAR"
  }
  if (is.null(npdecol)){
    stop("Missing column names/numbers of the data frame containing NPDE data.")
  }else if (!is.null(npdecol) && !is.numeric(npdecol)){
    if (!all(npdecol%in%names(plotdata))) stop("All column names given as NPDE data column must be present in plotdata")
  }else if (!is.null(npdecol) && is.numeric(npdecol)){
    if (any(npdecol <= 0) | (max(npdecol) > ncol(plotdata))) stop("Column number for NPDE data out of range of plotdata.")
    npdecol <- names(plotdata)[npdecol]
  }
  
  plotdata <- subset(plotdata, select=c("XVAR",npdecol))
  longdata <- melt(plotdata, measure.vars = npdecol)
  names(longdata)[c((ncol(longdata)-1):ncol(longdata))] <- c("type","npde")
  longdata <- na.omit(longdata)
  longdata <- subset(longdata, npde!="NaN")
  
  if (nrow(longdata)==0){
    return(list(forestdata=NULL))
  }else{
    npdecol         <- levels(factor(longdata$type))
    longdata$mean   <- numeric(nrow(longdata))
    longdata$sd     <- numeric(nrow(longdata))
    longdata$length <- numeric(nrow(longdata))
    longdata$ci     <- numeric(nrow(longdata))
    longdata$sd1    <- numeric(nrow(longdata))
    longdata$sd2    <- numeric(nrow(longdata))
    longdata$mcil   <- numeric(nrow(longdata))
    longdata$mciu   <- numeric(nrow(longdata))
    longdata$sdl    <- numeric(nrow(longdata))
    longdata$sdu    <- numeric(nrow(longdata))
    longdata$sdlcil <- numeric(nrow(longdata))
    longdata$sdlciu <- numeric(nrow(longdata))
    longdata$sducil <- numeric(nrow(longdata))
    longdata$sduciu <- numeric(nrow(longdata))
    
    for (i in 1:length(npdecol)){
      npdeval <- na.omit(as.numeric(longdata[longdata$type==npdecol[i] & longdata$npde!="NaN","npde"]))
      longdata[longdata$type==npdecol[i],"mean"]   <- mean(npdeval)
      longdata[longdata$type==npdecol[i],"sd"]     <- sd(npdeval)
      longdata[longdata$type==npdecol[i],"length"] <- length(npdeval)
      longdata[longdata$type==npdecol[i],"ci"]     <- abs(qt(0.025,length(npdeval)-1)*(sd(npdeval)/sqrt(length(npdeval))))
      longdata[longdata$type==npdecol[i],"sd1"]    <- abs(sd(npdeval)-(sd(npdeval)*sqrt((length(npdeval)-1)/qchisq(0.975,length(npdeval)-1))))
      longdata[longdata$type==npdecol[i],"sd2"]    <- abs(sd(npdeval)-(sd(npdeval)*sqrt((length(npdeval)-1)/qchisq(0.025,length(npdeval)-1))))
      longdata[longdata$type==npdecol[i],"sdcil"]  <- sd(npdeval)*sqrt((length(npdeval)-1)/qchisq(0.975,length(npdeval)-1))
      longdata[longdata$type==npdecol[i],"sdciu"]  <- sd(npdeval)*sqrt((length(npdeval)-1)/qchisq(0.025,length(npdeval)-1))
    }
    
    longdata$FCT    <- paste0(longdata$type,"\nmean=", signif(longdata$mean,3), "+/-", signif(longdata$ci,3), "\nSD=", signif(longdata$sd,3), "(-", signif(longdata$sd1,3), ",+", signif(longdata$sd2,3),")")
    
    longdata$mcil   <- longdata$mean - longdata$ci
    longdata$mciu   <- longdata$mean + longdata$ci
    longdata$sdl    <- longdata$mean - longdata$sd
    longdata$sdu    <- longdata$mean + longdata$sd
    longdata$sdlcil <- longdata$mean - longdata$sd - longdata$sd1
    longdata$sdlciu <- longdata$mean - longdata$sd + longdata$sd2
    longdata$sducil <- longdata$mean + longdata$sd - longdata$sd1
    longdata$sduciu <- longdata$mean + longdata$sd + longdata$sd2
    
    fctNm <- data.frame()
    npr   <- length(npdecol)
    nc    <- ifelse(npr<2, 1, ifelse(npr>=2 & npr<=6, 2, 3))
    for (p in 1:npr){
      if (npdecol[p] == "npdeAUClast" | npdecol[p] == "npdeAUClower_upper" | npdecol[p] == "npdeAUCINF_obs" | npdecol[p] == "npdeAUCINF_pred"){
        fctNm <- rbind(fctNm, data.frame(prmNm=npdecol[p],prmUnit=paste(gsub("^npde", "", npdecol)[p]," (",cunit,"*",tunit,")",sep="")))
      }else if (npdecol[p] == "npdeAUMClast"){
        fctNm <- rbind(fctNm, data.frame(prmNm=npdecol[p],prmUnit=paste(gsub("^npde", "", npdecol)[p]," (",cunit,"*",tunit,"^2)",sep="")))
      }else if (npdecol[p] == "npdeCmax"){
        fctNm <- rbind(fctNm, data.frame(prmNm=npdecol[p],prmUnit=paste(gsub("^npde", "", npdecol)[p]," (",cunit,")",sep="")))
      }else if (npdecol[p] == "npdeTmax"){
        fctNm <- rbind(fctNm, data.frame(prmNm=npdecol[p],prmUnit=paste(gsub("^npde", "", npdecol)[p]," (",tunit,")",sep="")))
      }else if (npdecol[p] == "npdeHL_Lambda_z"){
        fctNm <- rbind(fctNm, data.frame(prmNm=npdecol[p],prmUnit=paste(gsub("^npde", "", npdecol)[p]," (",tunit,")",sep="")))
      }else{
        fctNm <- rbind(fctNm, data.frame(prmNm=npdecol[p],prmUnit=npdecol[p]))
      }
    }
    longdata$type <- factor(longdata$type, levels=fctNm$prmNm, labels=fctNm$prmUnit)
    forestdata    <- subset(longdata[!duplicated(longdata$type),], select=c(type,mean,mcil,mciu,sd,sdcil,sdciu))
    
    # ggplot options for individual NPDE plots
    nc    <- ifelse(length(npdecol)<2, 1, ifelse(length(npdecol)>=2 & length(npdecol)<=6, 2, 3))
    ggOpt_npde <- list(xlab("\nID"), ylab("NPDE\n"),
                       theme(plot.title = element_text(size=10,face="bold"),
                             axis.title.x = element_text(size=10,face="bold"),
                             axis.title.y = element_text(size=10,face="bold"),
                             axis.text.x  = element_text(size=10,face="bold",color="black",angle=45,vjust=1,hjust=1),
                             axis.text.y  = element_text(size=10,face="bold",color="black",hjust=0),
                             panel.margin = unit(0.5, "cm"), plot.margin  = unit(c(0.5,0.5,0.5,0.5), "cm")),
                       geom_point(size=2), facet_wrap(~type, ncol=nc, scales="free"),
                       theme(strip.text.x = element_text(size=10, face="bold")))
    
    # ggplot options for histogram of NPDE
    ggOpt_histnpde <- list(scale_linetype_manual(name="",values=c("meanNormal"="solid","meanNPDE"="solid","SD"="dashed")),
                           scale_color_manual(name="",values=c("meanNormal"="red","meanNPDE"="blue","SD"="blue")),
                           xlab("\nNPDE"), ylab(""),
                           guides(fill = guide_legend(override.aes = list(linetype = 0 )), shape = guide_legend(override.aes = list(linetype = 0 ))),
                           theme(plot.title = element_text(size=10,face="bold"),
                                 axis.title.x = element_text(size=10,face="bold"),
                                 axis.title.y = element_text(size=10,face="bold"),
                                 axis.text.x  = element_text(size=10,face="bold",color="black",angle=45,vjust=1,hjust=1),
                                 axis.text.y  = element_text(size=10,face="bold",color="black",hjust=0),
                                 legend.text  = element_text(size=10,face="bold"),
                                 legend.background = element_rect(),
                                 legend.position = "bottom", legend.direction = "horizontal",
                                 legend.key.size = unit(0.8, "cm"),
                                 panel.margin = unit(0.5, "cm"),
                                 plot.margin  = unit(c(0.5,0.5,0.5,0.5), "cm")),
                           geom_histogram(aes(y=(..count..)/tapply(..count..,..PANEL..,sum)[..PANEL..]), size=0.6, color="black", fill="white"),
                           scale_y_continuous(labels = percent),
                           geom_vline(aes(xintercept=0, color="meanNormal", lty="meanNormal"), size=1, show_guide=T),
                           geom_vline(aes(xintercept=as.numeric(mean), color="meanNPDE", lty="meanNPDE"), size=1),
                           geom_vline(aes(xintercept=as.numeric(sdl), color="SD", lty="SD"), size=1),
                           geom_vline(aes(xintercept=as.numeric(sdu), color="SD", lty="SD"), size=1),
                           facet_wrap(~FCT, ncol=nc, scales="free"),
                           theme(strip.text.x = element_text(size=10, face="bold")))
    
    if (is.null(figlbl)){
      ttl   <- paste("NPDE vs. ",xvar,"\n\n",sep="")
    }else{
      ttl   <- paste("NPDE vs. ",xvar," (",figlbl,")\n\n",sep="")
    }
    ggnpde <- ggplot(longdata,aes(as.numeric(as.character(XVAR)),as.numeric(as.character(npde)))) + ggOpt_npde + labs(title = ttl)
    
    if (is.null(figlbl)){
      ttl   <- paste("Histogram of NPDE\n\n",sep="")
    }else{
      ttl   <- paste("Histogram of NPDE (",figlbl,")\n\n",sep="")
    }
    gghnpde <- ggplot(longdata,aes(as.numeric(as.character(npde)))) + ggOpt_histnpde + labs(title = ttl)
    return(list(forestdata=forestdata,ggnpde=ggnpde,gghnpde=gghnpde))
  }
}

