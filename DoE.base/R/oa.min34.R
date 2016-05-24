oa.min34 <- function (ID, nlevels, variants=NULL, min3 = NULL, all = FALSE, rela = FALSE) 
{
    ## retrieve child array or array identified by character string
          ## gsub for case where ID is character string
    IDname <- gsub("\"", "", deparse(substitute(ID)))
    if (all(IDname %in% oacat$name)) {
        if (!exists(IDname)) 
            ID <- eval(parse(text = paste("oa.design(", IDname, 
                ")")))
        else if (is.character(ID)) 
            ID <- eval(parse(text = paste("oa.design(", IDname, 
                ")")))
    }
    if (is.null(min3)) 
        min3 <- oa.min3(ID, nlevels, all = TRUE, rela = rela, variants=variants)
    if (!is.list(min3)) 
        stop("min3 must be a list")
    if (!all(c("GWP3", "column.variants", "complete") %in% names(min3))) {
        stop("min3 is not of the appropriate form")
    }
    if (!min3$complete) 
        warning("The min3 object should not just contain one of the best designs but the complete set")
    variants <- min3$column.variants
    GWP3 <- min3$GWP3
    if (nrow(variants) == 1) 
        return(list(GWP = c(GWP3, `4` = length4(ID[, variants])), 
            column.variants = variants))
    else {
     ## initialize curMin
        curMin <- Inf
        MinVariants <- numeric(0)
        for (i in 1:nrow(variants)) {
            spalten <- variants[i, ]
            if (GWP3 == 0) 
                cur4 <- round(length4(ID[, spalten], rela = rela), 
                  4)
            else cur4 <- round(length4(ID[, spalten]), 4)
            if (cur4 == 0 & !all) 
                return(list(GWP = c(GWP3, `4` = curMin), column.variants = matrix(spalten, 
                  nrow = 1), complete = FALSE))
            if (cur4 == curMin) 
                MinVariants <- rbind(MinVariants, spalten)
            else if (cur4 < curMin) {
                curMin <- cur4
                MinVariants <- matrix(spalten, nrow = 1)
            }
        }
        rownames(MinVariants) <- 1:nrow(MinVariants)
        curMin <- c(GWP3, `4` = curMin)
        if (GWP3 == 0) 
            names(curMin)[2] <- "4.relative"
        list(GWP = curMin, column.variants = MinVariants, complete = TRUE)
    }
}
