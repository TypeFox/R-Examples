## Plot the standard curves of each analyte
## 
## Plot the standard curve for each analyte of the  
## \code{scluminex} object in a ggplot format so new parameters can be added 
## 
## @param x a \code{scluminex} object
## @param subset.plot list of analytes to be plotted. By default plot all 
## analytes.
## @param psize point size
## @param ncol number of columns to plot the analytes.
## @param nrow number of rows to plot the analytes.
## @param color.bkg character specifying the color of the background line. 
## NA for not showing background.
## @param size.legend size of the legend. NA for not showing.
## @param interval 'confidence' or 'prediction' character in order to plot 
## the fit and the corresponding bands. If NULL only points are plotted.
## @param level confidence for the curve interval. 
## @param ... other \code{ggplot} arguments
##
## @details Background information from  \code{scluminex} object is assumed 
## to be in natural scale. The plot transform it in log10 scale.
## 
## @return A \code{ggplot} object
## @import ggplot2
## @importFrom plyr ldply dlply rbind.fill . summarise
plotLum <- function(x, subset.plot = NULL , psize=1.8, ncol = NULL, 
                    nrow = NULL, color.bkg = "green",  
                    size.legend = 2.5,interval = "confidence", 
                    level = 0.95, ...){

    if(!inherits(x, "scluminex")) stop("'x' must be a scluminex class object")
    x.scluminex <- x
    if(!is.null(subset.plot)){
        if(any(subset.plot%nin%names(x))){
            stop("'subset.plot' does not match with 'scluminex' names object")
        } 
        x <- x[subset.plot]        
    }
    
    if(level>1 | level<0){
      stop("'level' must be a value between 0 and 1")  
    } 
    
    if(!is.null(interval)){
        tinterval <- charmatch(interval, c("confidence","prediction"))
        if(is.na(tinterval)){
          stop("'interval' argument must be 'confidence' or 'prediction'")  
        }
        interval <- switch(tinterval, "confidence", "prediction")
    }
    
    orix <- x
    fanalyte <- x[[1]]$fieldnames$fanalyte
    data_plot <-  ldply(lapply(x,function(y){y$data}),rbind)  
    data_plot$analyte <- data_plot[, fanalyte]
    data_curve <- data_plot
    data_plot$shape <- 1
    data_flag <- ldply(lapply(x,function(y){y$flag_data}),rbind)
    if (nrow(data_flag)> 0) {
        data_flag$shape <- 2
        data_plot <- rbind.fill(data_plot,data_flag)
    }

    data_plot$shape <- factor(data_plot$shape)  
    cvg <- ldply(lapply(x,function(y){y$convergence}),rbind)
    r <- ddply(data_plot, .(analyte), 
            summarise, x=min(log10_concentration), y=max(log10_mfi))
    names(cvg)[2] <- "cvgtf"
    cvg$convergence <- ifelse(cvg$cvgtf==1,"Converge","No converge")
    data_curve <- merge(data_curve,cvg)
    data_plot$analyte <- as.character(data_plot$analyte)
    ordanalyte <- order(as.character(data_plot$analyte),data_plot$shape )
    data_plot <- data_plot[ordanalyte,]
    nonconv <- cvg[cvg$cvgtf!=1,1]
    texts <- ldply(lapply(x,function(y){data.frame(
                rsquare=ifelse(is.null(y$rsquare),NA,y$rsquare),
                bkg_mean=log10(y$bkg_mean),
                bkg_method=y$bkg_method,fct=y$fct)}),rbind)
    texts$analyte <- texts$.id
    texts <- merge(texts,r)
    texts$rsq <- paste("FCT:", texts$fct, "\n",
                "BKG:", texts$bkg_method,"\n",
                "RS:" ,round(texts$rsquare, 3), sep="")
    texts$rsq <- sub("NA","", texts$rsq)
    texts$rsq <- with(texts, ifelse(analyte%in%nonconv, 
                sub("RS:","No-Converged", rsq), 
                as.character(rsq)))  

    if(!is.null(interval)){
        auxdata <- data_plot[,fanalyte]
        la <- unique(auxdata)
        xlimits <- lapply(la, function(y){  
        summary(data_plot[auxdata==y,"log10_concentration"])[c("Min.","Max.")]})
        suport <- lapply(xlimits, function(y) seq(y[1], y[2], length.out = 500))
        names(suport) <- la
        conf.info <- lapply(la, function(y) try(conf_bands(x.scluminex, y, 
                        xvalue = suport[[y]], 
                        level = level, 
                        interval = interval),
                        silent=TRUE))
        conf.info <- lapply(conf.info, function(x){ 
                        ans <- x
                        if(inherits(x,"try-error")){
                            ans <- data.frame(NA,NA,NA,NA,NA)
                        } 
                        ans})
        names(conf.info) <- la  
        conf.info <- ldply(conf.info)
        names(conf.info) <- c("analyte","est","estlo","estup","ase","xvalue")  
    }

    if(length(nonconv)>0)  conf.info <- subset(conf.info, analyte%nin%nonconv)
    analyte <- NULL
    log10_mfi <- NULL 
    log10_concentration <- NULL 
    shape <- NULL 
    x <- NULL 
    y <- NULL 
    rsq <- NULL 
    bkg_mean <- NULL
    ans <- ggplot(data_plot, aes(y=log10_mfi, x=log10_concentration) )
    ans <- ans + scale_colour_manual(values = c("blue","red","violet"))
    ans <- ans + geom_point(aes(shape=shape),size=psize)
    ans <- ans + scale_shape_manual(values=c(19,1))
    ans <- ans + geom_text(data=texts, 
            aes(x=x , y=y, label=rsq, hjust=0, vjust=1), 
            size=size.legend)
    ans <- ans + facet_wrap( ~ analyte, ncol=ncol, nrow=nrow, scales="free") 
    ans <- ans + theme(legend.position="none")
    if(orix[[1]]$alertbkg==0){
        ans <- ans + geom_hline(data=texts, aes(yintercept = bkg_mean), 
                color=color.bkg)
    }
    if(!is.null(interval)){
        xvalue <- est <- estlo <- estup <- NULL
        ans <- ans + geom_line(data=conf.info, aes(x=xvalue, y=est))
        ans <- ans + geom_line(data=conf.info, aes(x=xvalue, y=estlo),
                linetype="dotted")
        ans <- ans + geom_line(data=conf.info, aes(x=xvalue, y=estup),
                linetype="dotted")
    } 
    return(ans)
}


