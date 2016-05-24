#' Plot expected and observed table from SitemFit
#'
#' @param sout output from SitemFit
#' @param itemName name of item to plot
#' @param ...  Not used.  Forces remaining arguments to be specified by name.
#' @param showSampleSize whether to show the sample size at the top of the plot
#' @import ggplot2
#' @importFrom reshape2 melt
#' @export
SitemPlot <- function(sout, itemName, ..., showSampleSize=TRUE) {
    garbageArguments <- list(...)
    if (length(garbageArguments) > 0) {
        stop("Values for the '...' argument are invalid; use named arguments")
    }

    s1 <- sout[[itemName]]
    obs <- s1$orig.observed
    ex <- s1$orig.expected
    rowTotal <- apply(obs, 1, sum)
    mask <- rowTotal > 0
    obs <- (obs / rowTotal)[mask,]
    ex <- (ex / rowTotal)[mask,]
    ss <- data.frame(sscore=as.numeric(names(rowTotal)), n=rowTotal)
    both <- rbind(cbind(type="expected", melt(ex)),
                  cbind(type="observed", melt(obs)))
    both$outcome <- factor(both$outcome, colnames(obs))
    plot <- ggplot(both, aes_string(x="sumScore", y="value")) + facet_wrap(~type) + ylim(0,1) +
        labs(y="probability", title=itemName)
    guide.style <- guide_legend(keywidth=.1, keyheight=.5, direction = "horizontal", title.position = "top",
                                label.position="bottom", label.hjust = 0.5, label.vjust = .5,
                                label.theme = element_text(angle = 90, size=8))
    plot <- plot + geom_line(aes_string(color="outcome")) + guides(color = guide.style)
    if (showSampleSize) {
        plot <- plot + geom_text(data=ss, aes_string(label="n", x="sscore"),
                                 y = 1, size=2, angle=90)
    }
    plot
}

#' Plot expected and observed table from SitemFit
#'
#' WARNING: This function is under development. The API may change in a future release.
#'
#' @param grp an IFA group
#' @param itemName name of item to plot
#' @param ...  Not used.  Forces remaining arguments to be specified by name.
#' @param width sets the x axis to [-width,width]
#' @param dataBins number of partitions for the latent scores
#' @param basis the basis vector in the latent space
#' @param factor the score to use (TODO: should be a function of the basis vector?)
#' @export
iccPlot <- function(grp, itemName, ..., width=3, dataBins=11, basis=c(1), factor=1) {
    garbageArguments <- list(...)
    if (length(garbageArguments) > 0) {
        stop("Values for the '...' argument are invalid; use named arguments")
    }

    basis <- basis / sqrt(sum(basis^2))
    
  ix <- match(itemName, colnames(grp$param))
  if (length(ix) != 1 || is.na(ix)) stop(paste("Can't find", itemName))
  
  labels <- levels(grp$data[[itemName]])
  spec1 <- grp$spec[[ix]]
  pm <- t(rpf.prob(spec1, grp$param[1:rpf.numParam(spec1),ix], basis %*% t(seq(-width, width, .1))))
  icc <- as.data.frame(melt(pm, varnames=c("theta",'category')))
  icc$theta <- seq(-width, width, .1)
  icc$category <- ordered(icc$category, labels=labels)
  icc$type <- 'expected'
  
  score <- grp$score[,factor]
  breaks <- seq(min(score, na.rm=TRUE),
                max(score, na.rm=TRUE),
                length.out=dataBins+1)
  bin <- unclass(cut(score, breaks, include.lowest = TRUE))
  
  eout <- array(dim=c(dataBins, spec1@outcomes+1))
  est <- numeric(dataBins)
  
  for (px in 1:dataBins) {
    t <- table(grp$data[[itemName]][bin==px])
    est[px] <- sum(t)
    eout[px,2:(spec1@outcomes+1)] <- t / sum(t)
  }
  eout[,1] <- ((c(breaks,0) + c(0,breaks))/2)[2:(dataBins+1)]
  bin.n <- data.frame(n=est, theta=eout[,1])
  
  edf <- melt(as.data.frame(eout), id.vars=c('V1'),
              variable.name="category")
  edf$category <- ordered(unclass(edf$category), labels=labels)
  edf$theta <- edf$V1
  edf$V1 <- NULL
  edf$type <- 'observed'
  
  both <- rbind(edf, icc)
  both$type <- mxFactor(both$type, levels=c('expected', 'observed'))
  
  plot <- ggplot(both, aes_string("theta", "value")) +
     facet_wrap(~type) +
    ylim(0,1) + xlim(-width,width) + labs(y="probability") +
    geom_text(data=bin.n, aes_string(label="n", x="theta"), y = 1, size=1.5, angle=90)
  guide.style <- guide_legend(keywidth=.1, keyheight=.5, direction = "horizontal", title.position = "top",
                              label.position="bottom", label.hjust = 0.5, label.vjust = .5,
                              label.theme = element_text(angle = 90, size=8))
  if (length(labels) <= 12) {
    plot <- plot + geom_line(aes_string(color="category", linetype="category")) +
      guides(color = guide.style, linetype = guide.style)
  } else {
    plot <- plot + geom_line(aes_string(color="category")) + 
      guides(color = guide.style)
  }
  plot + labs(title = itemName)
}

#' Create item response map table
#'
#' Categories are placed at the mean score of the examinees who picked
#' that category.
#'
#' @param grp an IFA group
#' @param ...  Not used.  Forces remaining arguments to be specified by name.
#' @param factor which factor to plot (defaults to 1)
#' @return A data.frame of the raw data backing the plot. Item outcomes
#' without any observations are omitted.
#' @export
itemResponseMap <- function(grp, ..., factor=1) {
    garbageArguments <- list(...)
    if (length(garbageArguments) > 0) {
        stop("Values for the '...' argument are invalid; use named arguments")
    }

  item.mask <- grp$param[factor,] > 0
  result <- NULL
  for (ix in rev(colnames(grp$param)[item.mask])) {
    lev <- levels(grp$data[,ix])
    for (ox in 1:length(lev)) {
      mask <- grp$data[,ix]==lev[ox]
      mask <- !is.na(mask) & mask
      if (all(!mask)) next
      result <- rbind(result, data.frame(item=ix,
                                         outcome=ox, outcome.name=lev[ox],
                                         score=mean(grp$score[mask, factor], na.rm=TRUE)))
    }
  }
  result[is.finite(result$score),]
}

#' Plot item information in the latent distribution
#'
#' For multidimensional items, you will need to supply a basis
#' vector. This vector is normalized to unit length.
#' 
#' @param grp an IFA group
#' @param ...  Not used.  Forces remaining arguments to be specified by name.
#' @param width the plot will span from -width to width
#' @param showTotal whether to plot the total item information
#' @param basis the basis vector (for multidimensional items)
#' @export
plotInformation <- function(grp, ..., width=3, showTotal=FALSE, basis=c(1)) {
    garbageArguments <- list(...)
    if (length(garbageArguments) > 0) {
        stop("Values for the '...' argument are invalid; use named arguments")
    }

    basis <- basis / sqrt(sum(basis^2))
  spec <- grp$spec
  param <- grp$param
  i.name <- colnames(grp$param)
  loc <- seq(-width, width, .1)
  grid <- basis %*% t(loc)
  df <- list(score=loc)
  total <- numeric(length(loc))
  for (ix in 1:length(spec)) {
    id <- i.name[ix]
    s <- spec[[ix]]
    df[[id]] <- rpf.info(s, param[1:rpf.numParam(s),ix], grid, basis)
    total <- total + df[[id]]
  }
  if (showTotal) df$total <- total
  df <- as.data.frame(df)
  long<- melt(df, id.vars=c('score'), variable.name="item")
  long$item <- factor(long$item)
  ggplot(long, aes_string("score", "value", group="item")) +
    geom_line(size=1.1,aes_string(color="item")) + ylab("information")
}
