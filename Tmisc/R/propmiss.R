#' Missing stats
#' 
#' Returns the number of missing values, total length, and proportion
#' missing values for each variable in a data.frame
#' 
#' @author Stephen Turner
#' @keywords NA
#' 
#' @param df A data.frame.
#' 
#' @return A data.frame with missingness stats.
#' 
#' @examples
#' propmiss(data.frame(a=1:5, b=c(6,NA,NA,9,10)))
#' 
#' @export
propmiss <- function(df) {
    m <- sapply(df, function(x) {
        data.frame(
            nmiss=sum(is.na(x)), 
            n=length(x), 
            propmiss=sum(is.na(x))/length(x)
        )
    })
    d <- data.frame(t(m))
    d <- sapply(d, unlist)
    d <- as.data.frame(d)
    d$var <- row.names(d)
    row.names(d) <- NULL
    d <- cbind(d[ncol(d)],d[-ncol(d)])
    return(d[order(d$propmiss), ])
}