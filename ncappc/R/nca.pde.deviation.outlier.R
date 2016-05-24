# PDE and Deviation metrics
# roxygen comments
#' Calculates individual prediction distribution errors (PDE) and scaled
#' deviation of NCA metrics estimated from observed and simulated data.
#' Identifies outlier to population PK model.
#'
#' \pkg{nca.pde.deviation.outlier} calculates individual prediction distribution
#' errors (PDE) and scaled deviation of NCA metrics estimated from observed and
#' simulated data. Identifies outlier to population PK model.
#' 
#' \pkg{nca.pde.deviation.outlier} calculates individual prediction distribution
#' errors (PDE) and scaled deviation of NCA metrics estimated from observed and 
#' simulated data. The deviation of each estimated NCA metrics is scaled by the 
#' "spread" of the simulated values. The "spread" is measured either by the 95\%
#' parametric prediction interval or 95\% non-parametric prediction interval.
#' Any individual yielding an absolute value of the scaled deviation for any of
#' the selected NCA metrics greater than 1, is assigned as an outlier to the
#' corresponding population PK model. The allowed NCA metrics for this
#' diagnostic tests are "AUClast", "AUClower_upper", "AUCINF_obs",
#' "AUCINF_pred", "AUMClast", "Cmax", "Tmax" and "HL_Lambda_z". By default, this
#' function uses AUClast and Cmax metrics for the comparison.
#' 
#' @param obsdata A data frame containing the NCA metrics values estimated from
#'   the observed data
#' @param simdata A data frame containing the NCA metrics values estimated from
#'   the simulated data
#' @param idNm Column name for ID (\strong{"ID"})
#' @param id ID of the individual whose data is being evaluated
#' @param spread Measure of the spread of simulated data (ppi (95\% parametric
#'   prediction interval) or npi (95\% nonparametric prediction interval))
#'   (\strong{"npi"})
#' @param figlbl Figure label based on dose identifier and/or population
#'   stratifier, in addition to ID (\strong{NULL})
#' @param calcparam A character array of the NCA metrics used for calculations 
#'   of PDE and deviation. The allowed NCA metrics for this histograms are 
#'   "AUClast", "AUClower_upper", "AUCINF_obs", "AUCINF_pred", "AUMClast", 
#'   "Cmax", "Tmax" and "HL_Lambda_z". (\strong{c("AUClast", "Cmax")})
#' @param diagparam A character array of the NCA metrics used for diagnostic 
#'   test. The allowed NCA metrics for this histograms are "AUClast", 
#'   "AUClower_upper", "AUCINF_obs", "AUCINF_pred", "AUMClast", "Cmax", "Tmax" 
#'   and "HL_Lambda_z". (\strong{c("AUClast", "Cmax")})
#' @param cunit Unit for concentration (\strong{"[M].[L]^-3"})
#' @param tunit Unit for time (\strong{"[T]"})
#' @param noPlot Perform only NCA calculations without any plot generation
#'   (TRUE, FALSE) (\strong{FALSE})
#'
#' @return returns the observed data frame with added distance and simulation
#'   mean of the nCA metrics, and a data frame with the PDE values of the NCA
#'   metrics. If the individual is identified as an outlier for the PK model,
#'   histograms of the diagnostic NCA metrics are produced and a graphical
#'   object created by arrangeGrob function is returned.
#' @export
#'

nca.pde.deviation.outlier <- function(obsdata,simdata,idNm="ID",id=NULL,spread="npi",figlbl=NULL,calcparam=c("AUClast","Cmax"),diagparam=c("AUClast","Cmax"),cunit="[M].[L]^-3",tunit="[T]",noPlot=FALSE){
  
  "type" <- "..density.." <- "oval" <- "mval" <- "devl" <- "devu" <- "sval" <- "scale_color_manual" <- "scale_linetype_manual" <- "xlab" <- "ylab" <- "geom_histogram" <- "aes" <- "geom_vline" <- "facet_grid" <- "theme" <- "element_text" <- "unit" <- "element_rect" <- "ggplot" <- "labs" <- "coord_cartesian" <- "gtable_filter" <- "ggplot_gtable" <- "ggplot_build" <- "arrangeGrob" <- "textGrob" <- "gpar" <- "..count.." <- "..PANEL.." <- "sd" <- "quantile" <- "scale_y_continuous" <- "percent" <- "packageVersion" <- NULL
  rm(list=c("type","..density..","oval","mval","devl","devu","sval","scale_color_manual","scale_linetype_manual","xlab","ylab","geom_histogram","aes","geom_vline","facet_grid","theme","element_text","unit","element_rect","ggplot","labs","coord_cartesian","gtable_filter","ggplot_gtable","ggplot_build","arrangeGrob","textGrob","gpar","..count..","..PANEL..","sd","quantile","scale_y_continuous","percent","packageVersion"))
  
  # Check the mandatory arguments
  if(is.null(obsdata) | is.null(simdata)){stop("None of the obsdata and simdata arguments can be empty.")}
  alwprm <- c("AUClast","AUClower_upper","AUCINF_obs","AUCINF_pred","AUMClast","Cmax","Tmax","HL_Lambda_z")
  if (!all(calcparam%in%alwprm)){stop("Incorrect calcparam. Please select NCA metrics from \"AUClast\", \"AUClower_upper\", \"AUCINF_obs\", \"AUCINF_pred\", \"AUMClast\", \"Cmax\", \"Tmax\", \"HL_Lambda_z\".")}
  if (!all(diagparam%in%alwprm)){stop("Incorrect diagparam. Please select NCA metrics from \"AUClast\", \"AUClower_upper\", \"AUCINF_obs\", \"AUCINF_pred\", \"AUMClast\", \"Cmax\", \"Tmax\", \"HL_Lambda_z\".")}
  if (!all(calcparam%in%names(obsdata))){stop("obsdata data frame must contain columns with NCA mterics given in calcparam argument.")}
  if (!all(calcparam%in%names(simdata))){stop("simdata data frame must contain columns with NCA mterics given in calcparam argument.")}
  if (!all(diagparam%in%names(obsdata))){stop("obsdata data frame must contain columns with NCA mterics given in diagparam argument.")}
  if (!all(diagparam%in%names(simdata))){stop("simdata data frame must contain columns with NCA mterics given in diagparam argument.")}
  
  if (!(idNm%in%names(obsdata)) | !(idNm%in%names(simdata))){stop("Column name for ID is not present in observed and/or simulated data.")}
  if (is.null(id)){
    obsuid <- unique(obsdata[,idNm])
    simuid <- unique(simdata[,idNm])
    if (length(obsuid)>1 & length(simuid)>1){
      stop("Please provide the ID of the individual whose data is being evaluated.")
    }else{
      id <- ifelse(length(obsuid)==1, obsuid, simuid)
    }
  }else{
    if (length(id)!=1 | (!is.element(id, obsdata[,idNm])) | (!is.element(id, simdata[,idNm]))){stop("Either more than one values are provided to id argument, or provided id value is not present in observed and/or simulated data.")}
  }
  
  iobslst    <- as.list(apply(subset(obsdata, eval(parse(text=idNm))==id, select=calcparam), 2, FUN=function(x) as.numeric(as.character(x[x!="NaN" & !is.na(x) & x!=Inf & x!=-Inf]))))
  isimlst    <- as.matrix(apply(subset(simdata, eval(parse(text=idNm))==id, select=calcparam), 2, FUN=function(x) as.numeric(as.character(x[x!="NaN" & !is.na(x) & x!=Inf & x!=-Inf]))))
  if(is.null(colnames(isimlst))) isimlst <- t(isimlst)
  
  metric     <- ""    # NCA metric associated with the outlier
  diagNm     <- ""    # List of outliers for each NCA diagnostic metric
  pdata      <- data.frame(oval=numeric(0),sval=numeric(0),mval=numeric(0),devl=numeric(0),devu=numeric(0),xl=numeric(0),xu=numeric(0),type=character(0))
  pde        <- data.frame(matrix(ncol=length(calcparam),nrow=1))   # store PDE values
  names(pde) <- calcparam
  if (length(iobslst)==0 | length(isimlst)==0){
    pde[,calcparam]                         <- "NaN"
    obsdata[,paste("d",calcparam,sep="")]   <- "NaN"
    obsdata[,paste("sim",calcparam,sep="")] <- "NaN"
  }else{
    for (i in 1:length(iobslst)){
      pnm <- colnames(isimlst)[i]
      if (length(iobslst[[pnm]])==0 | length(unlist(isimlst[,pnm]))==0){
        pde[,pnm]                         <- "NaN"
        obsdata[,paste("d",pnm,sep="")]   <- "NaN"
        obsdata[,paste("sim",pnm,sep="")] <- "NaN"
      }else{
        obsval    <- iobslst[[pnm]]
        simval    <- unlist(isimlst[,pnm])
        msimval   <- mean(simval)
        sdsimval  <- sd(simval)
        sdsimmean <- sdsimval*(simval-msimval)
        sdobsmean <- sdsimval*(obsval-msimval)
        if (spread == "ppi"){
          distprm <- ifelse(sdsimval==0, (obsval-msimval), (obsval-msimval)/(2*sdsimval))
          lldist  <- msimval-1.96*sdsimval
          uldist  <- msimval+1.96*sdsimval
        }else if (spread == "npi"){
          distprm <- ifelse((obsval-msimval)>0, (obsval-msimval)/(unname(quantile(simval, 0.975))-msimval), (obsval-msimval)/(msimval-unname(quantile(simval, 0.025))))
          lldist  <- unname(quantile(simval, 0.025))
          uldist  <- unname(quantile(simval, 0.975))
        }
        pde[,pnm]                         <- ifelse (obsval<min(simval), 1/length(simval), ifelse (obsval>max(simval), 1-(1/length(simval)), sum(sdsimmean<sdobsmean)/length(simval)))
        obsdata[,paste("d",pnm,sep="")]   <- distprm
        obsdata[,paste("sim",pnm,sep="")] <- msimval
        pdata  <- rbind(pdata, data.frame(oval=obsval, sval=simval, mval=msimval, devl=lldist, devu=uldist, xl=min(lldist,obsval), xu=max(uldist,obsval), type=pnm, stringsAsFactors = F))
        if ((!is.na(distprm) & !is.nan(distprm)) && (pnm%in%diagparam) & (abs(distprm)>1)){
          metric <- paste(metric,paste("ID-",id,"_",pnm,sep=""),sep=", ")
        }
      }
    }
  }
  
  # If outlier exists plot the NCA metrics
  # Initiate grob
  gdr      <- NULL
  mylegend <- NULL
  lheight  <- NULL
  
  # ggplot options for outliers
  if (metric != ""){
    pdata     <- subset(pdata, type%in%diagparam)
    diagparam <- unique(pdata$type)
    npr       <- length(diagparam)
    metric    <- gsub("^, ", "", metric)
    fctNm     <- data.frame()
    nc        <- ifelse(npr<2, 1, ifelse(npr>=2 & npr<=6, 2, 3))
    for (p in 1:npr){
      if (diagparam[p] == "AUClast" | diagparam[p] == "AUClower_upper" | diagparam[p] == "AUCINF_obs" | diagparam[p] == "AUCINF_pred"){
        fctNm <- rbind(fctNm, data.frame(prmNm=diagparam[p],prmUnit=paste(diagparam[p]," (",cunit,"*",tunit,")",sep="")))
      }else if (diagparam[p] == "AUMClast"){
        fctNm <- rbind(fctNm, data.frame(prmNm=diagparam[p],prmUnit=paste(diagparam[p]," (",cunit,"*",tunit,"^2)",sep="")))
      }else if (diagparam[p] == "Cmax"){
        fctNm <- rbind(fctNm, data.frame(prmNm=diagparam[p],prmUnit=paste(diagparam[p]," (",cunit,")",sep="")))
      }else if (diagparam[p] == "Tmax"){
        fctNm <- rbind(fctNm, data.frame(prmNm=diagparam[p],prmUnit=paste(diagparam[p]," (",tunit,")",sep="")))
      }else if (diagparam[p] == "HL_Lambda_z"){
        fctNm <- rbind(fctNm, data.frame(prmNm=diagparam[p],prmUnit=paste(diagparam[p]," (",tunit,")",sep="")))
      }
    }
    
    if(noPlot==FALSE){
      ggOpt_otl <- list(scale_color_manual(name="",values=c("Obs"="red","meanSim"="blue","+/-spread"="blue")),
                        scale_linetype_manual(name="",values=c("Obs"="solid","meanSim"="solid","+/-spread"="dashed")),
                        xlab(""), ylab(""),
                        geom_histogram(aes(y=(..count..)/tapply(..count..,..PANEL..,sum)[..PANEL..]), size=0.6, color="black", fill="white"),
                        scale_y_continuous(labels = percent),
                        geom_vline(aes(xintercept=oval, color="Obs", linetype="Obs"), show_guide=T, size=1),
                        geom_vline(aes(xintercept=mval, color="meanSim", linetype="meanSim"), show_guide=T, size=1),
                        geom_vline(aes(xintercept=devl, color="+/-spread", linetype="+/-spread"), show_guide=T, size=1),
                        geom_vline(aes(xintercept=devu, color="+/-spread", linetype="+/-spread"), show_guide=T, size=1),
                        facet_grid(~type, scales="free"),
                        theme(plot.title = element_text(size=9, face="bold"),
                              axis.title.x = element_text(size=9,face="bold"),
                              axis.title.y = element_text(size=9,face="bold"),
                              axis.text.x  = element_text(size=9,face="bold",color="black",angle=45,vjust=1,hjust=1),
                              axis.text.y  = element_text(size=9,face="bold",color="black",hjust=0),
                              panel.margin = unit(0.5, "cm"), plot.margin  = unit(c(0.2,0.2,0.2,0.2), "cm"),
                              legend.position = "bottom", legend.direction = "horizontal",
                              legend.text  = element_text(size=8,face="bold"),
                              legend.background = element_rect(),
                              legend.key.size = unit(0.8, "cm"),
                              strip.text.x = element_text(size=8, face="bold")))
      devtag <- ifelse (spread=="ppi","95% parametric prediction interval","95% nonparametric prediction interval")
      gplt   <- list()
      
      if(is.null(figlbl)){
        figttl <- paste0("Outlier_ID-",id,"\n(spread = ",devtag,")\n\n")
      }else{
        figttl <- paste0("Outlier_ID-",id,"_",figlbl,"\n(spread = ",devtag,")\n\n")
      }
      
      for (p in 1:npr){
        df <- subset(pdata, type==diagparam[p])
        df$type <- factor(df$type, levels=diagparam[p], labels=fctNm[fctNm$prmNm==diagparam[p],"prmUnit"])
        xl <- df$xl[1]; xu <- df$xu[1]
        gplt[[p]] <- ggplot(df,aes(x=as.numeric(sval))) + ggOpt_otl +
          labs(title=paste("Obs=",format(df$oval[1],digits=2),", meanSim=",format(df$mval[1],digits=2),"\n+/-spread=(",format(df$devl[1],digits=2),",",format(df$devu[1],digits=2),")\n",sep="")) +
          coord_cartesian(xlim=c(xl,xu))
      }
      mylegend <- suppressMessages(suppressWarnings(gtable_filter(ggplot_gtable(ggplot_build(gplt[[1]])), "guide-box", trim=T)))
      lheight  <- sum(mylegend$heights)
      for (p in 1:npr){gplt[[p]] <- gplt[[p]] + theme(legend.position="none")}
      
      plot_args <- list(top = textGrob(figttl,vjust=1,gp=gpar(fontface="bold",cex = 0.8)),
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
    }
  }
  return(list(obsdata=obsdata,pde=pde,metric=metric,grob=gdr,legend=mylegend,lheight=lheight))
}
