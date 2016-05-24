bfFixefLMER_t.fnc<-function (model, item = FALSE, method = c("t", "z", "llrt", "AIC", "BIC", "relLik.AIC", "relLik.BIC"), threshold = NULL, t.threshold = NULL, alphaitem = NULL, prune.ranefs = TRUE, set.REML.FALSE = TRUE, keep.single.factors = FALSE, reset.REML.TRUE = TRUE, log.file = NULL) {
  if (length(item) == 0) {
    stop("please supply a value to argument \"item\".\n")
  }
  if (length(set.REML.FALSE) == 0) {
    stop("please supply a value to argumnet \"set.REML.FALSE\".\n")
  }
  if (length(reset.REML.TRUE) == 0) {
    stop("please supply a value to argument \"reset.REML.TRUE\".\n")
  }
  if (!method[1] %in% c("t", "z", "llrt", "AIC", "BIC", "relLik.AIC", "relLik.BIC")) {
    stop("please supply a proper method name (t, z, llrt, AIC, BIC, relLik.AIC, or relLik.BIC).\n")
  }
  if ("factor" %in% attributes(attributes(model@frame)$terms)$dataClasses & 
        max(c(0, apply(as.data.frame(model@frame[, names(attributes(attributes(model@frame)$terms)$dataClasses[attributes(attributes(model@frame)$terms)$dataClasses == "factor"])]), MARGIN = 2, FUN = function(x) length(levels(factor(x)))))) > 2) {
	if(!as.vector(model@call[1]) == "glmer()"){
    	warning("factor variable with more than two levels in model terms, backfitting on t-values is not appropriate, please use function \"bfFixefLMER_F.fnc\" instead.\n")
	}
  }
  if (is.null(threshold)) {
    threshold <- c(2, 2, 0.05, 5, 5, 4, 4)[match(method[1], c("t", "z", "llrt", "AIC", "BIC", "relLik.AIC", "relLik.BIC"))]
  }
  if (is.null(t.threshold)) {
    if ((method[1] == "t") || (method[1] == "z")) {
      t.threshold <- 2
    }
    else {
      t.threshold <- Inf
    }
  }
  if (is.null(alphaitem)) {
    if (method[1] == "llrt") {
      alphaitem <- threshold
    }
    else {
      alphaitem <- 0.05
    }
  }
  if (is.null(log.file)) {
    log.file <- file.path(tempdir(), paste("bfFixefLMER_t_log_", gsub(":", "-", gsub("", "_", date())), ".txt", sep = ""))
  }
  if ((method[1] != "t" & set.REML.FALSE) || (method[1] != "z" & set.REML.FALSE)) {
	if(!as.vector(model@call[1]) == "glmer()"){
    	cat("setting REML to FALSE\n")
    	model <- update(model, . ~ ., REML = FALSE)
	}
  }
  data <- model@frame
  statistic <- "t-value"
  temp.dir <- tempdir()
  tempdir()
  options(warn = 1)
  unlink(file.path(temp.dir, "temp.txt"))
  sink(file = NULL, type = "message")
  if (log.file != FALSE) 
    sink(file = log.file, split = TRUE)
  if (item != FALSE) {
    cat("checking", paste("by-", item, sep = ""), "random intercepts\n")
    model.updated <- NULL
    eval(parse(text = paste("model.updated=update(model,.~.+(1|", item, "))", sep = "")))
    if (as.vector(anova(model, model.updated)[2, "Pr(>Chisq)"]) <= alphaitem) {
      cat("  log-likelihood ratio test p-value =", as.vector(anova(model, model.updated)[2, "Pr(>Chisq)"]), "\n")
      cat("  adding", paste("by-", item, sep = ""), "random intercepts to model\n")
      model <- model.updated
    }
    else {
      cat("  log-likelihood ratio test p-value =", as.vector(anova(model, model.updated)[2, "Pr(>Chisq)"]), "\n")
      cat("  not adding", paste("by-", item, sep = ""), "random intercepts to model\n")
    }
  }
  if (as.vector(model@call[1]) == "glmer()") {
    statistic <- "z-value"
    odv <- data[, as.character(unlist(as.list(model@call))$formula[2])]
    data[, as.character(unlist(as.list(model@call))$formula[2])] = rnorm(nrow(data), 0, 1)
	ow<-options()$warn
  	options(warn = -1)
    temp.lmer <- update(model, . ~ ., family = "gaussian", data = data)
	options(warn=ow)
    coefs <- row.names(anova(temp.lmer))
    data[, as.character(unlist(as.list(model@call))$formula[2])] <- odv
  }
  else {
    coefs <- row.names(anova(model))
  }
  if (is(model, "merMod")) {
    model <- asS4(model)
    smry <- as.data.frame(summary(model)$coefficients)
  }
  else {
    stop("the input model is not a mer object\n")
  }
  smry.temp <- smry[-c(1:nrow(smry)), ]
  smry <- as.data.frame(summary(model)$coefficients)
  for (coef in coefs) {
    intr.order <- length(unlist(strsplit(coef, ":")))
    orig.coef <- coef
    coef <- gsub(":", "\\.\\*:", coef)
    coef <- gsub("(^.*$)", "\\1\\.\\*", coef)
    coef <- gsub("\\(", "\\\\(", coef)
    coef <- gsub("\\)", "\\\\)", coef)
    smry.temp2 <- smry[grep(paste("^", coef, sep = ""), 
                            row.names(smry)), ]
    smry.temp3 <- smry.temp2
    smry.temp3 <- smry.temp3[-c(1:nrow(smry.temp3)), ]
    for (j in 1:nrow(smry.temp2)) {
      if (length(unlist(strsplit(row.names(smry.temp2)[j], ":"))) == intr.order) {
        smry.temp3 <- rbind(smry.temp3, smry.temp2[j, ])
      }
    }
    smry.temp2 <- smry.temp3
    smry.temp2[, 3] = abs(smry.temp2[, 3])
    smry.temp2 <- smry.temp2[smry.temp2[, 3] == max(smry.temp2[, 3]), ]
    row.names(smry.temp2) <- orig.coef
    smry.temp <- rbind(smry.temp, smry.temp2)
  }
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
  for (order in orders) {
    cat("processing model terms of interaction level", order, "\n")
    keepers <- as.character(row.names(smry.temp[smry.temp$Order == order, ]))
    smry.temp2 <- smry.temp[keepers, ]
    if (smry.temp2[smry.temp2[, 3] == min(smry.temp2[, 3]), 3][1] >= t.threshold) {
      cat("  all terms of interaction level", order, "significant\n")
    }
    while (smry.temp2[smry.temp2[, 3] == min(smry.temp2[, 3]), 3][1] < t.threshold) {
      cat("  iteration", count, "\n")
      cat("    ", statistic, "for term", paste("\"", row.names(smry.temp2[smry.temp2[, 3] == min(smry.temp2[, 3]), ])[1], "\"", sep = ""), 
          "=", smry.temp2[smry.temp2[, 3] == min(smry.temp2[, 3]), 3][1], "<", t.threshold, "\n")
      fac <- unlist(strsplit(gsub("\\)", "", gsub("\\(", "", as.character(row.names(smry.temp2[smry.temp2[, 3] == min(smry.temp2[, 3]), ]))[1])), ":"))
      hoi <- gsub("\\)", "", gsub("\\(", "", intr.order[intr.order$Order > order, "Coef"]))
      tt <- c()
      for (i in 1:length(hoi)) {
        tt[i] <- length(grep("FALSE", unique(fac %in% unlist(strsplit(hoi[i], ":"))))) < 1
      }
      if (length(grep("TRUE", tt)) != 0) {
        cat("     part of higher-order interaction\n")
        cat("     skipping term\n")
        keepers <- row.names(smry.temp2)
        keepers <- keepers[-grep(as.character(row.names(smry.temp2[smry.temp2[, 3] == min(smry.temp2[, 3]), ]))[1], keepers)]
        smry.temp2 <- smry.temp2[keepers, ]
      }
      else {
        cat("     not part of higher-order interaction\n")
        m.temp <- NULL
        eval(parse(text = paste("m.temp=update(model,.~.-", row.names(smry.temp2[smry.temp2[, 3] == min(smry.temp2[, 3]), ])[1], ")", sep = "")))
        if (method[1] == "llrt") {
          if (as.vector(anova(model, m.temp)[2, "Pr(>Chisq)"]) <= threshold) {
            cat("    log-likelihood ratio test p-value =", as.vector(anova(model, m.temp)[2, "Pr(>Chisq)"]), "<=", threshold, "\n")
            cat("    skipping term\n")
            keepers <- row.names(smry.temp2)
            cat("length =", length(keepers), "\n")
            keepers <- keepers[-which(keepers == as.character(row.names(smry.temp2[smry.temp2[, 3] == min(smry.temp2[, 3]), ]))[1])]
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
            cat(paste("     ", method[1], " simple = ", round(m.temp.ic), "; ", method[1], " complex = ", round(model.ic), "; decrease = ", round(ic.diff), " >= ", threshold, "\n", sep = ""))
            cat("     skipping term\n")
            keepers <- row.names(smry.temp2)
            cat("length =", length(keepers), "\n")
            keepers <- keepers[-which(keepers == as.character(row.names(smry.temp2[smry.temp2[, 3] == min(smry.temp2[, 3]), ]))[1])]
            smry.temp2 <- smry.temp2[keepers, ]
          }
          else {
            cat(paste("     ", method[1], " simple = ", round(m.temp.ic), "; ", method[1], " complex = ", round(model.ic), "; decrease = ", round(ic.diff), " < ", threshold, "\n", sep = ""))
          }
        }
        if (method[1] %in% c("relLik.AIC", "relLik.BIC")) {
          ic <- gsub("relLik.(.*)", "\\1", method[1])
          m.temp.ic <- eval(parse(text = paste(ic, "(m.temp)", sep = "")))
          model.ic <- eval(parse(text = paste(ic, "(model)", sep = "")))
          my.relLik <- relLik(m.temp, model, ic)["relLik"]
          reduction <- FALSE
          cat(paste("     ", ic, " simple = ", round(m.temp.ic), "; ", ic, " complex = ", round(model.ic), "; relLik = ", round(my.relLik), "\n", sep = ""))
          if (m.temp.ic <= model.ic) {
            reduction <- TRUE
          }
          if (!reduction) {
            if (my.relLik >= threshold) {
              cat(paste("     ", "relLik of ", round(my.relLik), " >= ", threshold, "\n", sep = ""))
              cat("     skipping term\n")
              keepers <- row.names(smry.temp2)
              cat("length =", length(keepers), "\n")
              keepers <- keepers[-which(keepers == as.character(row.names(smry.temp2[smry.temp2[, 3] == min(smry.temp2[, 3]), ]))[1])]
              smry.temp2 <- smry.temp2[keepers, ]
            }
            else {
              cat(paste("     ", "relLik of ", round(my.relLik), " < ", threshold, "\n", sep = ""))
              reduction <- TRUE
            }
          }
          else {
            cat(paste("     simple model more likely than complex model\n"))
          }
        }
        if (method[1] == "t") {
          reduction <- TRUE
        }
        if (reduction) {
          cat("     removing term\n")
          model <- m.temp
          if (as.vector(model@call[1]) == "glmer()") {
            odv <- data[, as.character(unlist(as.list(model@call))$formula[2])]
            data[, as.character(unlist(as.list(model@call))$formula[2])] <- rnorm(nrow(data), 0, 1)
			ow<-options()$warn
  			options(warn = -1)
            temp.lmer <- update(model, . ~ ., family = "gaussian", data = data)
			options(warn=ow)
            coefs <- row.names(anova(temp.lmer))
            data[, as.character(unlist(as.list(model@call))$formula[2])] <- odv
          }
          else {
            coefs <- row.names(anova(model))
            if (length(coefs) == 0) {
              break
            }
          }
          smry2 <- as.data.frame(summary(model)$coefficients)
          smry.temp2 <- smry2[-c(1:nrow(smry2)), ]
          for (coef in coefs) {
            intr.order <- length(unlist(strsplit(coef, 
                                                 ":")))
            orig.coef <- coef
            coef <- gsub(":", "\\.\\*:", coef)
            coef <- gsub("(^.*$)", "\\1\\.\\*", coef)
            coef <- gsub("\\(", "\\\\(", coef)
            coef <- gsub("\\)", "\\\\)", coef)
            smry.temp3 <- smry2[grep(paste("^", coef, 
                                           sep = ""), row.names(smry2)), ]
            smry.temp4 <- smry.temp3
            smry.temp4 <- smry.temp4[-c(1:nrow(smry.temp4)), 
                                     ]
            for (j in 1:nrow(smry.temp3)) {
              if (length(unlist(strsplit(row.names(smry.temp3)[j], ":"))) == intr.order) {
                smry.temp4 <- rbind(smry.temp4, smry.temp3[j, ])
              }
            }
            smry.temp3 <- smry.temp4
            smry.temp3[, 3] <- abs(smry.temp3[, 3])
            smry.temp3 <- smry.temp3[smry.temp3[, 3] == max(smry.temp3[, 3]), ]
            row.names(smry.temp3) <- orig.coef
            smry.temp2 <- rbind(smry.temp2, smry.temp3)
          }
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
          smry.temp2$Order <- intr.order$Order
          keepers <- as.character(row.names(smry.temp2[smry.temp2$Order == order, ]))
          smry.temp = smry.temp2
          smry.temp2 <- smry.temp2[keepers, ]
        }
      }
      count <- count + 1
      if (nrow(smry.temp2) == 0) {
        break
      }
    }
  }
  if (reset.REML.TRUE) {
	if(!as.vector(model@call[1]) == "glmer()"){
    	cat("resetting REML to TRUE\n")
    	model <- update(model, . ~ ., REML = TRUE)
	}
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
    m.ranefs <- split5
    if (as.vector(model@call[1]) == "glmer()") {
      statistic <- "z-value"
      odv <- data[, as.character(unlist(as.list(model@call))$formula[2])]
      data[, as.character(unlist(as.list(model@call))$formula[2])] = rnorm(nrow(data), 0, 1)
	  ow<-options()$warn
  	  options(warn = -1)
      temp.lmer <- update(model, . ~ ., family = "gaussian", data = data)
	  options(warn=ow)
      coefs <- row.names(anova(temp.lmer))
      data[, as.character(unlist(as.list(model@call))$formula[2])] <- odv
    }
    else {
      coefs <- row.names(anova(model))
    }
    m.ranef.variables <- gsub("\\((.*)\\|.*", "\\1", m.ranefs)
    m.ranef.variables <- gsub(".\\+(.*)", "\\1", m.ranef.variables)
    m.ranef.variables <- m.ranef.variables[m.ranef.variables != "1"]
    m.ranef.variables <- m.ranef.variables[m.ranef.variables != "0"]
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
        cat(paste("  ", iii, " in random effects structure but not in fixed effects structure\n", sep = ""))
        cat("    removing", iii, "from model ...\n")
        rtr2 <- c(rtr2, sub(iii, 1, grep(iii, m.ranefs, value = TRUE)))
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
