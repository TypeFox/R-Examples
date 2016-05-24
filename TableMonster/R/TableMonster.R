"%,%" <- function(x,y)paste(x,y,sep="")

"print.TableMonster" <-
function(x, ...)
{
  m <- match.call()
  ddd <- list()
  nmsddd <- names(m)[-(1:2)]
  n.ddd <- length(nmsddd)
  if(n.ddd>0) for(k in 1:n.ddd) ddd[[nmsddd[k]]] <- m[[2+k]]
    
  x.df <- as.data.frame(x)
  nr <- nrow(x.df)
  nc <- ncol(x.df)
  headings <- attr(x,"headings")
  ctypes <- attr(x,"ctypes")
  digits <- attr(x,"digits")
  n.h <- length(headings)
  depth <- rep(1, n.h)
  lngths <- NULL  
  for(k in 1:n.h)
  {
    ptr1 <- ptr0 <- headings[[k]] 
    if(!is.null(names(ptr1)))
    {
      ptr0 <- ptr1
      depth[k] <- depth[k] + 1
      ptr1 <- ptr0[[1]]
    }
    lnptr <- length(ptr0)
    lngths <- c(lngths, lnptr)
  }
  mxdpth <- max(depth)
  atmxdpth <- which(depth==mxdpth)

  for(k in 1:n.h)
  {
    j <- mxdpth - depth[k]
    out <- headings[[k]]
    while(j>0)
    {
      out <- list(` `=out)
      names(out) <- names(headings)[k]
      j <- j-1
    }
    headings[[k]] <- out
  }

  # extract headings and subheadings -- this solution for this step
  # currently is locked on depth 2 tables -- good to go for sometime.
  hdr <- list()
  hdr[[1]] <- names(headings[atmxdpth])
  hdr[[mxdpth]] <- substring(names(unlist(headings)),
                             unlist(regexec(".", names(unlist(headings)),
                                            fixed=TRUE))+1, nchar(names(unlist(headings))))
  h1 <- h1a <- NULL
  dpth2 <- any(depth>1)
  if(dpth2)
  {
    h1 <- h1a <- NULL
    h1[atmxdpth] <- "\\multicolumn{" %,% lngths[atmxdpth[1]] %,% "}{c}{" %,% hdr[[1]] %,% "}"

    h1[setdiff(1:n.h, atmxdpth)] <- ""
    h1 <- paste(h1, collapse="&") %,% "\\\\\n"

    nc1 <- length(hdr[[1]])
    tt <- cumsum(lngths)
    i0 <- tt[atmxdpth-1]+1
    i1 <- tt[atmxdpth]
    ni <- length(i0)
    prfx <- "\\cmidrule(r){" %,% i0[1] %,% "-" %,% i1[1] %,% "}"
    k.k <-  apply(cbind(i0,i1)[2:(ni-1),,drop=FALSE], 1, FUN=function(x)x[1]%,%"-"%,%x[2])
    bdy <- paste("\\cmidrule(lr){" %,% k.k, collapse="}")
    sfx <- "}\\cmidrule(l){" %,% i0[ni] %,% "-" %,% i1[ni] %,% "}\n"
    h1a <- prfx %,% bdy %,% sfx
  }
  
  h2 <- paste(hdr[[mxdpth]], collapse="&") %,% "\\\\\n"

  nc2 <- length(hdr[[mxdpth]])
  prfx <- "\\cmidrule(r){" %,% 1 %,% "-" %,% 1 %,% "}"
  k.k <-  sapply(2:(nc2-1), FUN=function(x)x%,%"-"%,%x)
  bdy <- paste("\\cmidrule(lr){" %,% k.k, collapse="}")
  sfx <- "}\\cmidrule(l){" %,% nc2 %,% "-" %,% nc2 %,% "}\n"
  h2a <- prfx %,% bdy %,% sfx
  
  add.to.row <- list()
  add.to.row[["command"]] <-
      c("\\toprule\n",
        h1,
        h1a,
        h2,
        h2a,
        "\\bottomrule\n")
  add.to.row[["pos"]] <- list()
  add.to.row[["pos"]][1:2] <- -1
  add.to.row[["pos"]][3:(3+2*dpth2)] <- 0
  add.to.row[["pos"]][4+2*dpth2] <- nr
  
  caption <- attr(x, "caption")
  
  xtbl.call <- as.call(expression(xtable, as.data.frame(x), digits=c(0,digits), align="ll" %,% paste(rep("r",nc-1), collapse=""), caption=caption))
  pr.xtbl.call <- as.call(expression(print, xtbl, hline.after=NULL, add.to.row=add.to.row, include.rownames=FALSE, include.colnames=FALSE, 
                                     type="latex"))

  if(n.ddd > 0)
  {
    lbl.idx <- grep("label", nmsddd)
    is.lbl <- (length(lbl.idx) > 0)
    if(is.lbl)
    {
      lbl.val <- ddd[[lbl.idx]]
      ddd <- ddd[-lbl.idx]
      n.ddd <- n.ddd - 1
      nmsddd <- names(ddd)
      xtbl.call[["label"]] <- lbl.val
    }
    is.ddd <- (n.ddd>0)
    if(is.ddd) for(k in 1:n.ddd) pr.xtbl.call[[nmsddd[k]]] <- ddd[[nmsddd[k]]]
  }
  
  xtbl <- eval(xtbl.call)
  eval(pr.xtbl.call)
}

"as.data.frame.TableMonster" <-
function(x, row.names = NULL, optional = FALSE, ...)
{
    attr(x, "headings") <- NULL
    attr(x, "ctypes") <- NULL
    attr(x, "digits") <- NULL
    class(x) <- "data.frame"
    x
}
