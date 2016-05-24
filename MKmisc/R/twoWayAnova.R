## Modification of function Anova in package genefilter
twoWayAnova <- function(cov1, cov2, interaction = TRUE, na.rm = TRUE){
    function(x) {
        if (na.rm) {
            drop <- is.na(x)
            x <- x[!drop]
            cov1 <- cov1[!drop]
            cov2 <- cov2[!drop]
        }
        if(interaction){
          av <- anova(lm(x ~ cov1*cov2))
          res <- av[["Pr(>F)"]][1:3]
          names(res) <- c("cov1", "cov2", "interaction")
        }else{
          av <- anova(lm(x ~ cov1 + cov2))
          res <- av[["Pr(>F)"]][1:2]
          names(res) <- c("cov1", "cov2")
        }

        return(res)
    }
}
