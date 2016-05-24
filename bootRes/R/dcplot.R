dcplot <-
function(x, ci = TRUE, sig = TRUE, labels = NULL, vertical = FALSE) {

  ## TODO Schöne Monatsnamen erzeugen und als Achsenbeschriftung
  ## darstellen
  ## TODO Schöne Unterscheidungen für Variablen (horizontale Darstellung)
  
  op <- par(no.readonly = TRUE) 

  if(any(is.na(x$sig))) {
    sig <- FALSE
    ci <- FALSE
  }
  
  n <- dim(x)[1]
  
  if (is.null(labels)) {
    labels <- rownames(x)
    ## make nice labels
    labels2 <- labels
    for (i in 1:n) {
      if (strsplit(rownames(x)[i], "\\.")[[1]][2] == "prev") {
        upcase <- FALSE
      } else {
        upcase <- TRUE
      }
      if (upcase) {
        labels2[i] <- toupper(substr(strsplit(rownames(x)[i],
                                              "\\.")[[1]][3], 1, 1))
      } else {
        labels2[i] <- tolower(substr(strsplit(rownames(x)[i],
                                              "\\.")[[1]][3], 1, 1))
      }
    }
  } else {
    labels2 <- labels
  }

  vnames <- character(n)
  for (i in 1:n) {
    vnames[i] <- strsplit(rownames(x)[i], "\\.")[[1]][1]
  }
  vnames.u <- unique(vnames)
  nice.vnames <- paste(toupper(substring(vnames.u, 1, 1)),
                       substring(vnames.u, 2), sep = "")
  no.vars <- length(unique(vnames))

  if (no.vars == 1) {

    if (ci) {
      plot.range <- range(c(x$ci.upper, x$ci.lower))
    } else {
      plot.range <- range(x$coef)
    }     
    
    if (sig) {
      pb <- barplot(x$coef, ylim = plot.range, border = "#FFFFFF00", yaxt
                    = "n", col = ifelse(x$significant == 1, "grey60",
                             "grey85"), ylab = "Response coefficients") # draw plot
      axis(side = 2, at = round(seq(signif(min(x$ci.lower), 1), signif(max(x$ci.upper),
                       1),
                       by = 0.1), 2), las = 2)
      
    } else {
      pb <- barplot(x$coef, ylim = plot.range, border = "#FFFFFF00", yaxt
                    = "n", ylab = "Response coefficients") # draw plot
      axis(side = 2, at = round(seq(signif(min(x$coef), 1), signif(max(x$coef),
                       1),
                       by = 0.1), 2), las = 2)
    }
    
    if (ci) {
      lle <- (pb[2] - pb[1])/4                # lineend-length calculated
                                        # from bar-distance
      for (i in 1:n) {
        lines(c(pb[i], pb[i]), c(x$ci.lower[i], x$ci.upper[i]))
        lines(c(pb[i] -lle, pb[i] +lle), c(x$ci.lower[i],
                                           x$ci.lower[i]))
        lines(c(pb[i] -lle, pb[i] +lle), c(x$ci.upper[i],
                                           x$ci.upper[i]))
      }
    }
    axis(side = 1, at = pb, labels = labels2, line
      = 1)
    text(pb[2], ifelse(ci, max(x$ci.upper), max(x$coef))*1.2,
         nice.vnames[i], xpd = NA, text = 4)
    return(NULL)
  }

  if (ci) {
    plot.range <- range(c(x$ci.upper, x$ci.lower))
  } else {
    plot.range <- range(x$coef)
  }     

  ## NA padding/splitting positions
  
  repl <- n/no.vars
  na.pos <- c(0, c(repl * c(1:no.vars))[-no.vars])
  n.nas <- length(na.pos)

  
  if (!vertical) {                     # horizontal layout

    x2 <- NULL
    axis.pos <- NULL
    for (i in 1:n.nas) {
      x2 <- rbind(x2, NA, x[((na.pos[i]+1):(na.pos[i]+repl)),])
      axis.pos <- c(axis.pos, na.pos[i] + (i-1))
    }

    x2 <- x2[-1,]
    no.axis.pos <- axis.pos[-1]
    
    if (sig) {
      pb <- barplot(x2$coef, ylim = plot.range, border = "#FFFFFF00",
                    yaxt = "n", col = ifelse(!is.na(x2$significant),
                                  ifelse(x2$significant == 1, "grey60", "grey85"),
                                  "#FFFFFF00"), ylab = "Response coefficients") # draw plot
      axis(side = 2, at = round(seq(signif(min(x$ci.lower), 1), signif(max(x$ci.upper),
                       1),,
                       by = 0.1), 2), las = 2)
    } else {
      pb <- barplot(x2$coef, ylim = plot.range, border = "#FFFFFF00", yaxt
                    = "n", ylab = "Response coefficients") # draw plot
      axis(side = 2, at = round(seq(signif(min(x$coef), 1), signif(max(x$coef),
                       1),
                       by = 0.1), 2), las = 2)
    }
    
    if (ci) {
      lle <- (pb[2] - pb[1])/4                # lineend-length calculated
                                        # from bar-distance
      for (i in 1:(n+1)) {
        lines(c(pb[i], pb[i]), c(x2$ci.lower[i], x2$ci.upper[i]))
        lines(c(pb[i] -lle, pb[i] +lle), c(x2$ci.lower[i],
                                           x2$ci.lower[i]))
        lines(c(pb[i] -lle, pb[i] +lle), c(x2$ci.upper[i],
                                           x2$ci.upper[i]))
      }
    }
    axis(side = 1, at = pb[-no.axis.pos], labels = labels2, line
      = 1)
    for (i in 1:no.vars) {
      text(pb[axis.pos[i]+1], ifelse(ci, max(x$ci.upper), max(x$coef))*1.2,
           nice.vnames[i], xpd = NA, pos = 4)
    }

  } else {                              # vertical plot layout

    par(mfrow = c(no.vars, 1), oma = c(0.5, 1, 0, 0), mai = c(0.6, 0.5, 0.3, 0.1))

    ## loop through variables and create one barplot for each

    for (i in 1:no.vars) {

      xs <- x[((na.pos[i]+1):(na.pos[i]+repl)),]
      
      if (sig) {
      
        pb <- barplot(xs$coef, ylim = plot.range, border = "#FFFFFF00",
                      yaxt = "n", col = ifelse(!is.na(xs$significant),
                                    ifelse(xs$significant == 1, "grey60", "grey85"),
                                    "#FFFFFF00"), ylab = "Response coefficients")
        
        axis(side = 2, at = round(seq(signif(min(x$ci.lower), 1),
                         signif(max(x$ci.upper), 1), by = 0.1), 2), las
             = 2)
      
      } else {
        
        pb <- barplot(xs$coef, ylim = plot.range, border = "#FFFFFF00",
                      yaxt = "n", ylab = "Response coefficients") # draw
        
        axis(side = 2, at = round(seq(signif(min(x$coef), 1),
                         signif(max(x$coef), 1), by = 0.1), 2), las = 2)
      }
      
      if (ci) {
        lle <- (pb[2] - pb[1])/4                # lineend-length calculated
                                        # from bar-distance
        for (j in 1:(n/no.vars)) {
          
          lines(c(pb[j], pb[j]), c(xs$ci.lower[j], xs$ci.upper[j]))
          
          lines(c(pb[j] -lle, pb[j] +lle), c(xs$ci.lower[j],
                                             xs$ci.lower[j]))
          
          lines(c(pb[j] -lle, pb[j] +lle), c(xs$ci.upper[j],
                                             xs$ci.upper[j]))
          
        }
      }
      axis(side = 1, at = pb[-((na.pos[-1])+1)], labels = labels2[1:(n/no.vars)], line
           = 1)
      text(pb[1], ifelse(ci, max(x$ci.upper), max(x$coef))*1.2,
           nice.vnames[i], xpd = NA, pos = 4)
    }
  }
par(op)
}
