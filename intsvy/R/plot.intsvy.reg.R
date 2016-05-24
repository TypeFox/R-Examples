plot.intsvy.reg <- function(x, ..., vars, se=TRUE, sort=FALSE) {
  x <- summary(x)
  x <- x[!is.na(x)]
  
  if (missing(vars)) {
    vars = row.names(x[[1]])
  } else {
  # Retrieve estimates and standard errors
  vars <- grep(paste(vars, collapse="|"), rownames(x[[1]]), value=T)
  }
  
  coefE <- lapply(x, function(x) x[row.names(x) %in% vars, "Estimate"])
  coefS <- lapply(x, function(x) x[row.names(x) %in% vars, "Std. Error"])
  estimates <- do.call('rbind', coefE)
  colnames(estimates) <- row.names(x[[1]][row.names(x[[1]]) %in% vars, ])
  sErrors <- do.call('rbind', coefS)
  colnames(sErrors) <- row.names(x[[1]][row.names(x[[1]]) %in% vars, ])
  
  # Sort argument removed
  ndf <- droplevels(na.omit(data.frame(reshape::melt.array(estimates), reshape::melt.array(sErrors)[,3])))
  colnames(ndf) <- c("Group", "Coefficient", "Value", "se")
  ndf$valueL <- ndf$Value - 1.96*ndf$se
  ndf$valueH <- ndf$Value + 1.96*ndf$se
  # To change order of variable labels in facet_wrap (R-squared appears always in the end)
  ndf$Coefficient <- factor(ndf$Coefficient, levels = vars)
  
  # Display plot in alphabetical order
  if (isTRUE(sort)) {
  ndf$Group <- factor(ndf$Group, levels=sort(levels(ndf$Group), decreasing=T))
  }
  
  # Plot
  pl <- ggplot(data=ndf, aes_string(x = "Value", y="Group", shape= "Coefficient", color="Coefficient")) +
    geom_point(size=5) +
    theme_bw() +
    facet_wrap(~Coefficient, scales="free_x") +
    theme(legend.position="top")
  if (se) {
    pl <- pl + geom_errorbarh(aes_string(xmin="valueL", xmax="valueH"))
  }
  pl
}