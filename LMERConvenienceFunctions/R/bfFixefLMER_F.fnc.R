bfFixefLMER_F.fnc<-function (model, item = FALSE, method = c("F", "llrt", "AIC", 
                                          "BIC", "relLik.AIC", "relLik.BIC"), threshold = NULL, alpha = NULL, 
          alphaitem = NULL, prune.ranefs = TRUE, p.value = "upper", 
          set.REML.FALSE = TRUE, keep.single.factors = FALSE, reset.REML.TRUE = TRUE, 
          log.file = NULL) 
{
  if (length(item) == 0) {
    stop("please supply a value to argument \"item\".\n")
  }
  if (length(set.REML.FALSE) == 0) {
    stop("please supply a value to argumnet \"set.REML.FALSE\".\n")
  }
  if (length(reset.REML.TRUE) == 0) {
    stop("please supply a value to argument \"reset.REML.TRUE\".\n")
  }
  if (!method[1] %in% c("F", "llrt", "AIC", "BIC", "relLik.AIC", 
                        "relLik.BIC")) {
    stop("please supply a proper method name (F llrt AIC BIC, relLik.AIC, or relLik.BIC).\n")
  }
  if (as.vector(model@call[1]) == "glmer()") {
    stop("sorry, glmer models not yet supported")
  }
  if (is.null(threshold)) {
    threshold <- c(0.05, 0.05, 5, 5, 4, 4)[match(method[1], c("F", "llrt", "AIC", "BIC", "relLik.AIC", "relLik.BIC"))]
  }
  if (is.null(alpha)) {
    if (method[1] == "F") {
      alpha <- threshold
    }
    else {
      alpha <- 0
    }
  }
  if (is.null(alphaitem)) {
    if (method[1] == "F" || method[1] == "llrt") {
      alphaitem <- threshold
    }
    else {
      alphaitem <- 0.05
    }
  }
  if (is.null(log.file)) {
    log.file <- file.path(tempdir(), paste("bfFixefLMER_F_log_", gsub(":", "-", gsub(" ", "_", date())), ".txt", sep = ""))
  }
  data <- model@frame
  if (method[1] != "F" & set.REML.FALSE) {
    cat("setting REML to FALSE\n")
    model <- update(model, . ~ ., REML = FALSE)
  }
  options(warn = 1)
  temp.dir <- tempdir()
  tempdir()
  unlink(file.path(temp.dir, "temp.txt"))
  sink(file = NULL, type = "message")
  if (log.file != FALSE) 
    sink(file = log.file, split = TRUE)
  if (item != FALSE) {
    cat("checking", paste("by-", item, sep = ""), "random intercepts\n")
    model.updated <- NULL
    eval(parse(text = paste("model.updated<-update(model,.~.+(1|", item, "))", sep = "")))
    if (as.vector(anova(model, model.updated)[2, "Pr(>Chisq)"]) <= alphaitem) {
      cat("  log-likelihood ratio test p-value =", as.vector(anova(model, model.updated)[2, "Pr(>Chisq)"]), "\n")
      cat("  adding", paste("by-", item, sep = ""), "random intercepts to model\n")
      model <- model.updated
    }
    else {
      cat("  log-likelihood ratio test p-value =", as.vector(anova(model, model.updated)[2, "Pr(>Chisq)"]), "\n")
      cat("  not adding", paste("by-", item, sep = ""), 
          "random intercepts to model\n")
    }
  }
  coefs <- row.names(anova(model))
  smry <- pamer.fnc(model)
  smry.temp <- pamer.fnc(model)
  temp <- strsplit(coefs, ":")
  names(temp) <- coefs
  intr.order <- list()
  for (i in coefs) {
    intr.order[[i]] <- length(temp[[i]])
  }
  intr.order <- as.data.frame(unlist(intr.order))
  colnames(intr.order) <- "Order"
  intr.order$Coef <- row.names(intr.order)
  row.names(intr.order) <- 1:nrow(intr.order)
  orders <- sort(as.numeric(unique(intr.order$Order)), decreasing = TRUE)
  if (keep.single.factors) {
    orders <- orders[orders != 1]
  }
  smry.temp$Order <- intr.order$Order
  count <- 1
  order <- orders[1]
  for (order in orders) {
    cat("processing model terms of interaction level", order, "\n")
    keepers <- as.character(row.names(smry.temp[smry.temp$Order == order, ]))
    smry.temp2 <- smry.temp[keepers, ]
    if (smry.temp2[smry.temp2[, paste(p.value, ".p.val", sep = "")] == max(smry.temp2[, paste(p.value, ".p.val", sep = "")]), paste(p.value, ".p.val", sep = "")][1] < alpha) {
      cat("  all terms of interaction level", order, "significant\n")
    }
    else {
      while (smry.temp2[smry.temp2[, paste(p.value, ".p.val", sep = "")] == max(smry.temp2[, paste(p.value, ".p.val", sep = "")]), paste(p.value, ".p.val", sep = "")][1] >= alpha) {
        cat("  iteration", count, "\n")
        if (length(smry.temp2[, paste(p.value, ".p.val", sep = "")] == max(smry.temp2[, paste(p.value, ".p.val", sep = "")])) > 1) {
          evaluated.term <- as.character(row.names(smry.temp2[smry.temp2[, 3] == min(smry.temp2[, 3]), ]))
        }
        else {
          evaluated.term <- row.names(smry.temp2[smry.temp2[, paste(p.value, ".p.val", sep = "")] == max(smry.temp2[, paste(p.value, ".p.val", sep = "")]), ])[1]
        }
        cat("    p-value for term", paste("\"", evaluated.term, "\"", sep = ""), "=", smry.temp2[evaluated.term, paste(p.value, ".p.val", sep = "")], ">=", alpha, "\n") 
		hoi <- gsub("\\)", "", gsub("\\(", "", intr.order[intr.order$Order > order, "Coef"]))
        tt <- data.frame()
        for (i in 1:length(hoi)) {
          tt[i, 1] <- length(grep("FALSE", unique(unlist(strsplit(evaluated.term, ":")) %in% unlist(strsplit(hoi[i], ":"))))) < 1
        }
        if (length(grep("TRUE", tt[, 1])) != 0) {
          cat("    part of higher-order interaction\n")
          cat("    skipping term\n")
          keepers <- row.names(smry.temp2)
          keepers <- keepers[-grep(as.character(row.names(smry.temp2[smry.temp2[, 3] == min(smry.temp2[, 3]), ])), keepers)]
          smry.temp2 <- smry.temp2[keepers, ]
          smry.temp2 <- na.omit(smry.temp2)
        }
        else {
          cat("    not part of higher-order interaction\n")
          m.temp <- NULL
          eval(parse(text = paste("m.temp<-update(model,.~.-", evaluated.term, ")", sep = "")))
          if (method[1] == "llrt") {
            if (as.vector(anova(model, m.temp)[2, "Pr(>Chisq)"]) <= threshold) {
              cat("    log-likelihood ratio test p-value =", as.vector(anova(model, m.temp)[2, "Pr(>Chisq)"]), "<=", threshold, "\n")
              cat("    skipping term\n")
              keepers <- row.names(smry.temp2)
              cat("length =", length(keepers), "\n")
              keepers <- keepers[-which(keepers == as.character(row.names(smry.temp2[smry.temp2[, 3] == min(smry.temp2[, 3]), ])))]
              smry.temp2 <- smry.temp2[keepers, ]
              reduction <- FALSE
            }
            else {
              reduction <- TRUE
              cat("    log-likelihood ratio test p.value =", as.vector(anova(model, m.temp)[2, "Pr(>Chisq)"]), ">", threshold, "\n")
            }
          }
          if (method[1] %in% c("AIC", "BIC")) {
            m.temp.ic <- summary(m.temp)$AICtab[method[1]]
            model.ic <- summary(model)$AICtab[method[1]]
            ic.diff <- m.temp.ic - model.ic
            if (ic.diff >= threshold) {
              reduction <- FALSE
            }
            else {
              reduction <- TRUE
            }
            if (!reduction) {
              cat(paste("    ", method[1], " simple = ", round(m.temp.ic), "; ", method[1], " complex = ", round(model.ic), "; decrease = ", round(ic.diff), " >= ", threshold, "\n", sep = ""))
              cat("    skipping term\n")
              keepers <- row.names(smry.temp2)
              cat("length =", length(keepers), "\n")
              keepers <- keepers[-which(keepers == as.character(row.names(smry.temp2[smry.temp2[, 3] == min(smry.temp2[, 3]), ]))[1])]
              smry.temp2 <- smry.temp2[keepers, ]
            }
            else {
              cat(paste("    ", method[1], " simple = ", round(m.temp.ic), "; ", method[1], " complex = ", round(model.ic), "; decrease = ", round(ic.diff), " < ", threshold, "\n", sep = ""))
            }
          }
          if (method[1] %in% c("relLik.AIC", "relLik.BIC")) {
            ic <- gsub("relLik.(.*)", "\\1", method[1])
            m.temp.ic <- eval(parse(text = paste(ic, "(m.temp)", sep = "")))
            model.ic <- eval(parse(text = paste(ic, "(model)", sep = "")))
            my.relLik <- relLik(m.temp, model, ic)["relLik"]
            reduction <- FALSE
            cat(paste("    ", ic, " simple = ", round(m.temp.ic), "; ", ic, " complex = ", round(model.ic), "; relLik = ", round(my.relLik), "\n", sep = ""))
            if (m.temp.ic <= model.ic) {
              reduction <- TRUE
            }
            if (!reduction) {
              if (my.relLik >= threshold) {
                cat(paste("    ", "relLik of ", round(my.relLik), " >= ", threshold, "\n", sep = ""))
                cat("    skipping term\n")
                keepers <- row.names(smry.temp2)
                cat("length =", length(keepers), "\n")
                keepers <- keepers[-which(keepers == as.character(row.names(smry.temp2[smry.temp2[, 3] == min(smry.temp2[, 3]), ]))[1])]
                smry.temp2 <- smry.temp2[keepers, ]
              }
              else {
                cat(paste("    ", "relLik of ", round(my.relLik), " < ", threshold, "\n", sep = ""))
                reduction <- TRUE
              }
            }
            else {
              cat(paste("    simple model more likely than complex model\n"))
            }
          }
          if (method[1] == "F") {
            reduction <- TRUE
          }
          if (reduction) {
            cat("    removing term\n")
            model <- m.temp
            coefs <- row.names(anova(model))
            if (length(coefs) == 0) {
              break
            }
            smry <- pamer.fnc(model)
            smry.temp <- pamer.fnc(model)
            temp <- strsplit(coefs, ":")
            names(temp) <- coefs
            intr.order <- list()
            for (i in coefs) {
              intr.order[[i]] <- length(temp[[i]])
            }
            intr.order <- as.data.frame(unlist(intr.order))
            colnames(intr.order) <- "Order"
            intr.order$Coef <- row.names(intr.order)
            row.names(intr.order) <- 1:nrow(intr.order)
            smry.temp$Order <- intr.order$Order
            keepers <- as.character(row.names(smry.temp[smry.temp$Order == order, ]))
            smry.temp2 <- smry.temp[keepers, ]
            smry.temp2 <- na.omit(smry.temp2)
          }
        }
        count <- count + 1
        if (nrow(smry.temp2) == 0) {
          break
        }
      }
    }
  }
  if (reset.REML.TRUE) {
    cat("resetting REML to TRUE\n")
    model <- update(model, . ~ ., REML = TRUE)
  }
  if (prune.ranefs) {
    cat("pruning random effects structure ...\n")
    split1 <- gsub(" ", "", model@call)[2]
    split2 <- unlist(strsplit(split1, "\\~"))[2]
    split2 <- gsub("\\)\\+\\(", "\\)_____\\(", split2)
    split2 <- gsub("\\(0\\+", "\\(0\\&\\&\\&", split2)
    split2 <- gsub("\\(1\\+", "\\(0\\&\\&\\&", split2)
    split3 <- unlist(strsplit(split2, "\\+"))
    split4 <- grep("\\|", split3, value = TRUE)
    split4 <- gsub("\\&\\&\\&", "\\+", split4)
    split5 <- unlist(strsplit(split4, "_____"))
    m.ranefs <- vector("character")
    for (vb in 1:length(split5)) {
      if (length(which(unlist(strsplit(split5[vb], "")) == "|")) > 0) {
        m.ranefs <- c(m.ranefs, split5[vb])
      }
    }
    coefs <- row.names(anova(model))
    m.ranef.variables <- gsub("\\((.*)\\|.*", "\\1", m.ranefs)
    m.ranef.variables <- gsub(".\\+(.*)", "\\1", m.ranef.variables)
    m.ranef.variables <- m.ranef.variables[m.ranef.variables != 
                                             "1"]
    m.ranef.variables <- m.ranef.variables[m.ranef.variables != 
                                             "0"]
    ranef.to.remove <- vector("character")
    for (m.ranef.variable in m.ranef.variables) {
      if (!m.ranef.variable %in% coefs) {
        ranef.to.remove <- c(ranef.to.remove, m.ranef.variable)
      }
    }
    if (length(ranef.to.remove) > 0) {
      rtr <- vector("character")
      rtr2 <- vector("character")
      for (iii in ranef.to.remove) {
        rtr <- c(rtr, grep(iii, m.ranefs, value = TRUE))
        cat(paste("  ", iii, " in random effects structure, but not in fixed effects structure\n", sep = ""))
        cat("    removing", iii, "from model ...\n")
        rtr2 <- c(rtr2, sub(iii, 1, grep(iii, m.ranefs, 
                                         value = TRUE)))
      }
      eval(parse(text = paste("model<-update(model,.~.-", paste(rtr, collapse = "-"), "+", paste(rtr2, collapse = "+"), ",data=data)", sep = "")))
    }
    else {
      cat("  nothing to prune\n")
    }
  }
  options(warn = 0)
  sink(file = NULL, type = "message")
  if (log.file != FALSE) {
    cat("log file is", log.file, "\n")
    sink(file = NULL)
  }
  return(model = model)
}
