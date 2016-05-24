rm.dupl <- function(d1,d2,...){
         dd1 <- do.call("paste", as.data.frame(d1))
         dd2 <- do.call("paste", as.data.frame(d2))
         d <- d2[! dd2 %in% dd1, ]
         return(d)
}