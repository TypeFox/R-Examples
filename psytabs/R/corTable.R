corTable <- function(data, use = "pairwise", method = "pearson", round = 2, significance = "stars", sd=FALSE, mean.sd.cols = FALSE) {
  
  corr.test.out <- psych::corr.test(plyr::colwise(makeNumeric)(data), use = use, method = method) # calculate correlation
  corr.test.out.r <- (data.frame(round(corr.test.out$r, round))) # get correlation coefs
  corr.test.out.p <- (data.frame(round(corr.test.out$p, 3))) # get p-values
  
  pasteStars <- function(x, y) paste(x, y, sep="")
  pastePval <- function(x, y) paste(x, " (", y, ")", sep="")
  
  if("stars" %in% significance && "p-values" %in% significance) {
    star.table <- data.frame(ifelse(corr.test.out.p < .001, "***", ifelse(corr.test.out.p < .01, "** ", ifelse(corr.test.out.p < .05, "* ", "")))) # generate significance stars  
    significance.table <- data.frame(mapply(pastePval, star.table, corr.test.out.p)) 
    cor.table <- mapply(pasteStars, corr.test.out.r, significance.table)
    cor.table <- gsub("(0)", "(< 0.001)", cor.table, fixed=T)
  } else if("stars" %in% significance) {
    significance.table <- data.frame(ifelse(corr.test.out.p < .001, "***", ifelse(corr.test.out.p < .01, "** ", ifelse(corr.test.out.p < .05, "* ", "")))) # generate significance stars  
    cor.table <- mapply(pasteStars, corr.test.out.r, significance.table)
  } else if("p-values" %in% significance){
    cor.table <- mapply(pastePval, corr.test.out.r, corr.test.out.p)
    cor.table <- gsub("(0)", "(< 0.001)", cor.table, fixed=T)
  } else {
    cor.table <- corr.test.out.r
  }
  
  rownames(cor.table) <- colnames(cor.table)
  rownames(cor.table) <- paste(1:length(rownames(cor.table)), ". ", rownames(cor.table), sep="")
  colnames(cor.table) <- paste("(", 1:length(rownames(cor.table)), ")", sep="")
  
  cor.table[upper.tri(cor.table)] <- ""
  cor.table <- as.matrix(cor.table)
  diag(cor.table) <- ""
  
  if(sd) {
    sd2 <- function(x) sd(x, na.rm = TRUE)
    diag(cor.table) <- paste("(", round(sapply(data, sd2), 2), ")", sep="")
  } else {
    diag(cor.table) <- "-"
  }
    
  
  if(mean.sd.cols) {
    cor.table <- cbind(cor.table, M = round(colMeans(data, na.rm=TRUE), 2))
    cor.table <- cbind(cor.table, SD = round(sapply(data, sd, na.rm=TRUE), 2))
    order <- c("M", "SD", paste("(", 1:length(rownames(cor.table)), ")", sep=""))
    cor.table <- cor.table[, order]
  }
  cor.table
}