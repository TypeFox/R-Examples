flimboot <-
function(fo, R, counter=FALSE) {
  org.data <- fo$clean
  clusters <- unique(org.data[, 1])
  n.obs <- table(org.data[, 1])
  flim.boot <- list()
  input <- fo$info$input
  id <- names(org.data)[1]
  obstime <- names(org.data)[2]
  tval <- fo$t.values
  met <- fo$method
  lam <- fo$lambda
  for(s in 1:R) {
    if(counter) print(s)
    samp <- sort(sample(clusters, length(clusters), replace=TRUE))
    samp.df <- data.frame(id=samp)
    names(samp.df)[1] <- names(org.data)[1]
    boot.data = merge(samp.df, org.data, all.x=TRUE)
    unique.samp.id <- unique(samp)
    unique.samp.obs <- as.vector(n.obs[match(unique.samp.id, names(n.obs))])
    unique.samp.rep <- as.vector(table(samp))    
    id.rep.obs <- cbind(unique.samp.id, unique.samp.rep, unique.samp.obs)
    NameFun <- function(info) {
      names <- as.numeric(paste0(rep(info[1], info[2]), ".", seq(1, info[2])))
      names.long <- sort(rep(names, info[3]))
    }
    boot.data[, 1] <- do.call("c", apply(id.rep.obs, 1, NameFun))
    boot.data <- boot.data[order(boot.data[, 1], boot.data[, 2]), ]
    boot.fit <- flim(input, boot.data, id, obstime,
                     t.values=tval, method=met, lambda=lam)
    flim.boot[[s]] <- boot.fit
  }
  returnme <- list(org = fo,
                   samples = flim.boot,
                   call = match.call())
  
  class(returnme) <- "flimboot"
  return(returnme)
}
