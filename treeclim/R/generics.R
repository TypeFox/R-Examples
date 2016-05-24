##' @export print tc_dcc
print.tc_dcc <- function(x, ...) {
  ci <- x$call$ci
  if (is.null(ci))
    ci <- 0.05
  cat("Coefficients (significance flags correspond to p < ",
      ci, "):\n", sep = "")
  print(x$coef, ...)
}

##' @export print tc_seascorr
print.tc_seascorr <- function(x, ...) {
  n <- length(x$coef)
  ## pretty season names
  if (is.null(x$call$season_lengths)) {
    season_lengths <- c(1, 3, 6)
  } else {
    season_lengths <- eval(x$call$season_lengths)
  }

  sl_levels <- ifelse(season_lengths == 1,
                      paste(season_lengths, "month"),
                      paste(season_lengths, "months"))
  ## what is the end month?
  if (is.null(x$call$complete)) {
    end_month <- 9
  } else {
    end_month <- x$call$complete
  }
  ## translate this and the preceeding 14 months in literals
  lmonths <- c(-1:-12, 1:12)
  end_month_l <- which(lmonths == end_month)
  used_months <- lmonths[(end_month_l - 13):end_month_l]
  fm <- function(x, single = TRUE) {
    if (single) {
      format_month(x)$single
    } else {
      format_month(x)$names
    }
  }
  used_months_ch <- fm(used_months, ...)
  for (i in 1:n) {
    cat("Results for a season length of ", sl_levels[i],
        ":\n", sep = "")
    print(data.frame(
      month = used_months_ch,
      type = rep(c("primary", "secondary"), each = 14),
      coef = c(rev(x$coef[[i]]$primary$coef),
        rev(x$coef[[i]]$secondary$coef)),
      significant = c(rev(x$coef[[i]]$primary$significant),
        rev(x$coef[[i]]$secondary$significant))
      ))
    cat("\n")
  }
}

##' @export print tc_coef
print.tc_coef <- function(x, ...) {
  rownames(x) <- abbrev_name(rownames(x))
  x$coef <- round(x$coef, 3)
  x$ci_lower <- round(x$ci_lower, 3)
  x$ci_upper <- round(x$ci_upper, 3)
  print.data.frame(x, ...)
}

##' @importFrom abind abind
##' @export print tc_mcoef
print.tc_mcoef <- function(x, ...) {
  mm <- abind(x$coef, x$significant, along = 3)
  ms <- apply(mm, c(1, 2), function(x) {
    if (!is.na(x[2])) {
      if (x[2]) {
        xc <- paste(round(x[1], 3), "*", sep = "")
      } else {
        xc <- round(x[1], 3)
      }
    } else {
      xc <- round(x[1], 3)
    }
    xc
  })
  rownames(ms) <- abbrev_name(rownames(ms))
  print(ms, ...)
}

##' @export print tc_design
print.tc_design <- function(x, ...) {
  pr <- x$aggregate
  names(pr) <- abbrev_name(x$names)
  years <- as.numeric(rownames(pr))
  cat(length(years), "observations of", dim(pr)[2],
      "variables.\nObserved years:", min(years), "-", max(years),
      "\nUsed (potentially aggregated) variables:\n")
  cat(paste(names(pr), collapse = "\n"), ...)
  cat("\n")
}

##' @export print tc_skills
print.tc_skills <- function(x, ...) {
  cat("Call:\n", paste(deparse(x$call), sep = "\n",
                       collapse = "\n"), "\n\n", sep = "")
  cat("Using climate target", x$target, "and calibration model with", x$cal.str,
      "as calibration period:\n\n")
  print.default(format(c(r = x$r.cal, "p-value" = x$p.cal)),
                print.gap = 2L, quote = FALSE)
  cat("\nCoefficients:\n")
  print.default(format(x$coef.cal), print.gap = 2L,
                quote = FALSE, ...)
  cat("\nVerification statistics:\n")
  verstat <- c("Red. of Error (RE)" = round(x$RE, 3),
               "Coeff. of Efficiency (CE)" = round(x$CE, 3))
  print.default(format(verstat), print.gap = 2L,
                quote = FALSE, ...)
  cat("\nDurbin-Watson Test (DW):\t", round(x$DW$statistic, 3), " (p = ",
      x$DW$p.value, ")\n\n", sep = "")
  cat("Model for whole period:\n\n")
  print.default(format(c(r = x$r.full, "p-value" = x$p.full)),
                print.gap = 2L, quote = FALSE)
  cat("\nCoefficients:\n")
  print.default(format(x$coef.full), print.gap = 2L,
                quote = FALSE, ...)
}

##' @export plot tc_skills
plot.tc_skills <- function(x, ...) {
  orig <- x
  y <- prediction <- NULL               # to keep R CMD check happy
  d <- data.frame(x = orig$years,
                  y = orig$full$x)
  v <- data.frame(x = c(orig$cal.years, orig$ver.years),
                  y = c(orig$pred.cal, orig$pred.ver),
                  prediction = c(rep("calibration",
                    length(orig$pred.cal)),
                    rep("verification", length(orig$pred.ver))))
  gg <- ggplot(data = d, aes(x = x, y = y))
  gg + geom_line(data = d, aes(x = x, y = y)) +
    geom_line(data = v, aes(x = x, y = y, color = prediction)) +
      theme_minimal() +
        xlab("years") + ylab("target")
}

##' @export coef tc_dcc
coef.tc_dcc <- function(object, ...) {
  coef(object$coef, ...)
}

##' @export coef tc_coef
coef.tc_coef <- function(object, ...) {
  print(data.frame(object), ...)
}

##' @export coef tc_mcoef
coef.tc_mcoef <- function(object, ...) {
  print(data.frame(object$coef), ...)
}

##' @export summary tc_dcc
summary.tc_dcc <- function(object, ...) {
  cat("Call:\n")
  print(object$call)
  cat("\nDesign matrix:\n")
  print(object$design)
  cat("\nCoefficients:\n")
  print(object$coef, ...)
}

##' @export + tc_paramlist
"+.tc_paramlist" <- function(p1, p2) {
  ## check if p1 already is _nested_ list
  if (is.list(p1[[1]])) {
    ## just write to this list
    param_list <- p1
    i <- length(param_list) + 1
    param_list[[i]] <- p2
  } else {
    ## make a wrapper list
    param_list <- list()
    param_list[[1]] <- p1
    param_list[[2]] <- p2
  }
  param_list
}

##' @import ggplot2
##' @import plyr
##' @export plot tc_dcc
plot.tc_dcc <- function(x, ...) {
  data <- x$coef
  if (is.null(x$call$boot)) {
    boot <- "std"
  } else {
    boot <- x$call$boot
  }

  y <- varname <- ci_lower <- ci_upper <- significant <- pid <- wid <-
    vid <- NULL                         # to keep R CMD check happy

  if (any(class(data) == "tc_coef")) {

    line0 <- data.frame(
      x = c(0.5, data$id, max(data$id) + 0.5),
      y = rep(0, dim(data)[1] + 2)
      )

    if (boot == "exact") {
      
      gg <- ggplot(data, aes(x = id, y = coef), ...) +
        geom_line(data = line0, aes(x, y), color = "grey") +
        geom_point(aes(color = varname, pch = significant, size = significant)) +
        scale_size_manual(values = c(2, 4)) +
        scale_x_continuous(breaks = data$id, labels = data$month) +
        ylab("Coefficients") +
        xlab("Months") +  
        theme_minimal() +
        theme(axis.title.x = element_blank())    
      
    } else {
      
      gg <- ggplot(data, aes(x = id, y = coef), ...) +
        geom_line(data = line0, aes(x, y), color = "grey") +
        scale_linetype_manual(values = c("dotted", "solid")) +
        geom_point(aes(color = varname), size = 3) +
        scale_x_continuous(breaks = data$id, labels = data$month) +
        ylab("Coefficients") +
        xlab("Months") +  
        theme_minimal() +
        theme(axis.title.x = element_blank()) +
        geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper, color =
                          varname, lty = significant), size = 1)
    }

    gg
    
  } else {

    ## mdcc case

    coef <- data$coef
    n <- dim(coef)[2]
    m <- dim(coef)[1]

    ## reformat into ggplot compatible data.frame

    pdata <- data.frame(
      varname = abbrev_name(rep(rownames(coef), n)),
      window = rep(names(coef), each = m),
      coef = as.vector(as.matrix(coef)),
      significant = as.vector(as.matrix(data$significant))
      )

    pdata$wid <- rep(1:n, each = m)
    pdata$vid <- rep(1:m, n)

    pdata$pid <- factor(paste(pdata$wid, pdata$vid, sep = "."))

    create_grid <- function(x) {
      w <- x$wid
      v <- x$vid
      data.frame(
        pid = factor(rep(paste(w, v, sep = "."), 4)),
        x = c(w - 1, w, w, w - 1),
        y = c(v - 1, v - 1, v, v)
        )
    }

    idgrid <- ddply(pdata, .variables = c("wid", "vid"), create_grid)

    idpgrid <- merge(pdata, idgrid, by = c("pid"))

    gg <- ggplot(idpgrid, aes(x = window, y = varname), ...) +
      geom_polygon(aes(x, y, fill = coef, group = pid)) +
      scale_fill_gradient2() +
      theme_minimal() +
      scale_x_continuous(breaks = seq(0.5, by = 1, length.out = n),
                         labels = names(coef)) +
      scale_y_continuous(breaks = seq(0.5, by = 1, length.out = m),
                         labels = abbrev_name(rownames(coef))) +
      theme(axis.text.x = element_text(angle = 90, vjust = 0),
            axis.title.x = element_blank(),
            axis.title.y = element_blank()) + 
      geom_point(data = subset(pdata, significant),
                 aes(x = wid - 0.5, y = vid - 0.5), pch = 8,
                 color = "grey")
    
    gg 
  }
}


plot.tc_seascorr <- function(x, ...) {

  end_month_order <- value <- significant <- y <-
    prediction <- NULL                  # to keep R CMD check happy

  ## how many season lengths do we consider?
  n <- length(x$coef)
  ## this *2 corresponds to the number of facets we need

  ## how long are the seasons?
  if (is.null(x$call$season_lengths)) {
    season_lengths <- c(1, 3, 6)
  } else {
    season_lengths <- eval(x$call$season_lengths)
  }

  sl_levels <- ifelse(season_lengths == 1,
                      paste(season_lengths, "month"),
                      paste(season_lengths, "months"))

  ## what is the end month?
  if (is.null(x$call$complete)) {
    end_month <- 9
  } else {
    end_month <- x$call$complete
  }

  ## translate this and the preceeding 14 months in literals
  lmonths <- c(-1:-12, 1:12)
  end_month_l <- which(lmonths == end_month)
  used_months <- lmonths[(end_month_l - 13):end_month_l]
  used_months_ch <- format_month(used_months)$single

  ## reassemble into single data.frame for plotting
  gd <- data.frame(
    season_length = rep(rep(sl_levels, each = 14), 2),
    end_month = rep(rev(used_months_ch), n * 2),
    end_month_order = rep(14:1, n * 2),
    type = c(rep("primary", 14*n), rep("secondary", 14*n)),
    significant = c(unlist(sapply(sapply(x$coef, "[", 1), "[", 2)),
      unlist(sapply(sapply(x$coef, "[", 2), "[", 2))),
    value = c(unlist(sapply(sapply(x$coef, "[", 1), "[", 1)),
      unlist(sapply(sapply(x$coef, "[", 2), "[", 1)))
    )

  
  
  ## draw plot
  gg <- ggplot(gd,
               aes(x = end_month_order, y = value,
                   fill = significant), ...)

  gg + facet_grid(type ~ season_length) +
    geom_bar(stat = "identity", position = "dodge") +
      scale_fill_manual(values = c("grey75", "grey50")) +
        theme_minimal() +
          scale_x_continuous(breaks = 1:14,
                             labels = substr(used_months_ch, 1, 1)) +
                               xlab("Ending month") +
                                 ylab("Correlation coefficient")
  
}
