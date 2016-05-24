ordinalizeFunc <-
function(dat, col.time, col.status, col.id, col.ctg, min.levels=3, details=FALSE){
  vnames <- colnames(dat)
  time <- dat[, col.time]; status <- dat[,col.status];  id <- dat[, col.id]; 
  n <- NROW(dat)	
  p <- length(col.ctg)
  OUT <- as.list(1:p); names(OUT) <- vnames[col.ctg]
  for (j in 1:p){
    col <- col.ctg[j]
    vname <- vnames[col]
    x <- as.character(dat[, col])
    x.level <- sort(unique(x))
    
    # SHOULD WE MERGE SOME LEVELS FIRST?
    tab <- table(x, status)
    if (min(tab) <= min.levels & max(dim(tab))>2) print(paste("You might want to merge some levels of", vname, "first before ordinalization.")) 
    if (details) {print(vname); print(x.level); print(tab)}
    
    # FITTING MARGINAL MODEL
    options(warn=-1)
    fit <- coxph(Surv(time, status) ~ x + cluster(id), data=dat)
    options(warn=0)
    betas <- c(baseline=0, coef(fit))
    
    # THE INFORMATION REGARDING "ORDINALZIATION"
    OUT[[j]] <- cbind(x.level, betas, rank=rank(betas))
    if (details) print(OUT[[j]])
    
    # REPLACE THE COLUMN IN DATA
    x.level.ordered <- x.level[order(betas)]
    x1 <- ordered(x, levels = x.level.ordered)
    dat[, col] <- as.numeric(x1)
  }
  return(list(dat=dat, info=OUT))
}
