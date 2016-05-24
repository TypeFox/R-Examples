stoolsmlx <- function(identical.time){
  s <- "
gg_color_hue <- function(n) {
 hues = seq(15, 375, length=n+1)
 hcl(h=hues, l=65, c=100)[1:n]}

info_res <- function(f){
  if (!is.null(names(f))) 
    f <- list(f)
  nf <- length(f)
  info <- list()
  for (k in (1:nf)){
    fk <- f[[k]]
    nk <- length(fk$name)
    lk <- as.character(1:nk)
    valk <- gg_color_hue(nk)
    names(valk) <- lk
    labk <- fk$name
    names(labk) <- lk
    info[[k]] <- list(values=valk, labels=labk, colour=lk)
  }
  return(info)
}

merge_res <- function(r1,f){
  r2 <- list()
  nf <- length(f)
  for (j in (1:nf)){
    fj <- f[[j]]
    nj <- length(fj$name)
    r <- r1[[fj$name[1]]]
    if (nj>1){
      for (k in (2:nj)){
        rk <- r1[[fj$name[k]]]
        r <- merge(r,rk)
      }
    }
    r2[[j]] <- r
  }
"
  
  if (identical.time==TRUE){
    s <- paste0(s,"
  r3 <- r2[[1]]
  if (nf>1){
    for (k in (2:nf))
      r3 <- merge(r3,r2[[k]])
  }
  return(r3)
}
")
  }else{
    s <- paste0(s,"
  return(r2)
}
")
  }
  return(s)
}
