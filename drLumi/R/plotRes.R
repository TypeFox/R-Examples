## Plot residuals of each analyte
## 
## Plot the residuals for each analyte of the  \code{scluminex} object in a 
## ggplot format so new parameters can be added. 
## 
## @param x a \code{scluminex} object
## @param psize point size
## @param out.limit value that defines an outlier. Must be positive value.
## @param size.text value that defines the size of the well into the plot.
## @param subset.plot list of analytes to be plotted. By default 
## plot all analytes.
## @param ncol number of columns to plot the analytes.
## @param nrow number of rows to plot the analytes.
## @param ... other \code{ggplot} arguments
##
## @return A \code{ggplot} object
## 
## @import ggplot2
plotRes <- function(x, psize=1.8, out.limit=2.5, size.text=1.5, 
                    subset.plot=NULL, ncol=NULL, nrow=NULL, ...) {

    if(!inherits(x,"scluminex")) stop("'x' must be a scluminex object")
    if(!is.null(subset.plot)){
        if(any(subset.plot%nin%names(x))){
            stop("'subset.plot' vector not match with 'scluminex' object")
        } 
        x <- x[subset.plot]
    }

    if(out.limit<0) stop("'out.limit' must be a positive value")

    fwell <- x[[1]]$fieldnames$fwell
    fanalyte <- x[[1]]$fieldnames$fanalyte

    nm <- names(x)
    nas <- lapply(nm, function(y) compair(x[[y]]$data))
    names(nas) <- nm
    nas <- do.call(rbind,nas)
    ana_var <- rownames(nas)
    nas <- data.frame(nas)
    nas[,fanalyte] <- rownames(nas)

    data_plot <-  ldply(lapply(x,function(y){y$data}), rbind.fill)
    data_plot <- merge(data_plot, nas, by.x = fanalyte,all.x=TRUE, sort=FALSE)
    data_plot$predicted.log10_mfi <- with(data_plot, ifelse(nas==TRUE,0,
                                    predicted.log10_mfi))
    data_plot$analyte <- data_plot[,fanalyte]
    data_plot$residuals <- with(data_plot, ifelse(nas==TRUE,0,residuals))
    data_plot$new_analyte <- ifelse(data_plot$nas==TRUE, 
                            paste(data_plot$analyte,"- No Convergence",sep=""),
                            as.character(data_plot$analyte))
    data_plot$new_analyte <- as.factor(data_plot$new_analyte)
    data_plot$col <- with(data_plot, ifelse(abs(residuals)>out.limit,
                            "red","black"))
    data_plot$well <- data_plot[, fwell]
    data_outlier <- subset(data_plot,abs(residuals)>out.limit)

    ordanalyte <- order(as.character(data_plot$new_analyte) )
    data_plot <- data_plot[ordanalyte,]
    data_plot$outliers <- out.limit

    myColors <- c("black", "red")
    names(myColors) <- c("black","red")
    colScale <- scale_colour_manual(name = "col",values = myColors)

    outliers <- NULL
    well <- NULL
    predicted.log10_mfi <- NULL
    ans <- ggplot(data_plot, 
            aes(x=predicted.log10_mfi, y=residuals, colour = col), ...) 
    ans <- ans + geom_hline(aes(yintercept= -outliers , linetype="dotted"), 
            alpha=0.5, size=0.5) 
    ans <- ans + geom_hline(aes(yintercept= 0 , linetype="dash"), 
            alpha=0.5, size=0.5) 
    ans <- ans + geom_hline(aes(yintercept= outliers , linetype="dotted"), 
            alpha=0.5, size=0.5) 
    if (nrow(data_outlier)>0) {
            vjust <- NULL
            data_outlier$vjust <- with(data_outlier, ifelse(residuals>0,1.3,-0.3))
            ans <- ans + geom_text(data=data_outlier,
                    aes(x=predicted.log10_mfi,y=residuals, colour = col, 
                    label=well, hjust= 0.5,vjust=vjust), 
                    size=size.text)
    }
    ans <- ans + geom_point(size=psize) 
    ans <- ans + facet_wrap(~ new_analyte,ncol=ncol, nrow=nrow, scales="free_x")
    ans <- ans + colScale + guides(colour=FALSE)    
    return(ans)
}

