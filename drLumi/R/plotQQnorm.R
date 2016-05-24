## Plot qqnorm of the residuals for each analyte
## 
## Plot qqnorm of the residuals for each analyte the  \code{scluminex} 
## object in a ggplot format so new parameters can be added
##  
## @param x a \code{scluminex} object
## @param subset.plot list of analytes to be plotted. By default 
## plot all analytes.
## @param psize point size
## @param ncol number of columns to plot the analytes.
## @param nrow number of rows to plot the analytes.
## @param ... other \code{ggplot} arguments
##
## @return A \code{ggplot} object
## 
## @importFrom reshape melt melt.data.frame
## @import ggplot2
plotQQnorm <- function(x, subset.plot = NULL, psize=1.8, 
                ncol=NULL, nrow = NULL, ...) {
    
    if(!inherits(x, "scluminex")) stop("'x' must be a scluminex object")
    if(!is.null(subset.plot)){
        if(any(subset.plot%nin%names(x))){
            stop("'subset.plot' vector not match with 'scluminex' object")
        } 
        x <- x[subset.plot]
    }

    fanalyte <- x[[1]]$fieldnames$fanalyte

    nm <- names(x)
    nas <- lapply(nm, function(y) compair(x[[y]]$data))
    names(nas) <- nm
    nas <- do.call(rbind,nas)
    ana_var <- rownames(nas)
    nas <- data.frame(nas)
    nas[, fanalyte] <- rownames(nas)

    data_plot <-  ldply(lapply(x,function(y){y$data}),rbind.fill)
    data_plot <- merge(data_plot, nas, by.x = fanalyte,all.x=TRUE, sort=FALSE)
    data_plot$analyte <- data_plot[,fanalyte]
    data_plot$residuals <- with(data_plot, ifelse(nas==TRUE,0,residuals))
    data_plot$new_analyte <- ifelse(data_plot$nas==TRUE,
                        paste(data_plot$analyte,"- No Convergence",sep=""),
                        as.character(data_plot$analyte))
    data_plot$new_analyte <- as.factor(data_plot$new_analyte)
    ordanalyte <- order(as.character(data_plot$new_analyte) )
    data_plot <- data_plot[ordanalyte,]
    coord_norm <- lapply(sort(nm), function(x){
            ans <- try(qqnorm(data_plot[data_plot$new_analyte==x,"residuals"], 
                    plot=FALSE),silent=TRUE)
                    ans})
    names(coord_norm) <- sort(nm)

    coord_norm <- lapply(coord_norm, function(x){ 
                ans <- x
                if(inherits(x,"try-error")){
                    ans <- list(x=0,y=0) 
                } 
                ans
                })


    Theoretical <- melt(ldply(lapply(coord_norm, 
                    function(y) y$x), "rbind"), id = ".id")
    Observed <- melt(ldply(lapply(coord_norm, 
                    function(x) x$y), "rbind"), id = ".id")
    coord_norm_xy <- as.data.frame(cbind(Theoretical, Observed[,3]))
    names(coord_norm_xy) <- c("analyte","pos","Theoretical","Observed")

    new_analyte <- unique(data_plot[,c("analyte","new_analyte")])
    coord_norm_xy <- merge(coord_norm_xy, new_analyte, 
                    by = "analyte", all.x=TRUE, sort=FALSE)

    ans <- ggplot(data = coord_norm_xy, aes(x=Theoretical, y=Observed)) 
    ans <- ans + geom_point(size=psize) + geom_smooth(method="lm", se=FALSE)
    ans <- ans + facet_wrap( ~ new_analyte, ncol=ncol, nrow=nrow, 
                            scales="free_x")
    return(ans)
}
