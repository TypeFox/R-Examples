testList<-function (list, type = "start") 
{
    argNames = ls()
    if (!is.list(list)) 
        res <- list(test = FALSE, message = paste(deparse(substitute(list)), 
            " is not a list", sep = ""))
    else {
        if (is.null(names(list))) 
            res <- list(test = FALSE, message = paste("list '", 
                deparse(substitute(list)), "' has no argument names", 
                sep = ""))
        else {
            elements <- switch(type, start = c("fixItems", "seed", 
                "nrItems", "theta", "D", "randomesque", "startSelect", 
                "nAvailable"), test = c("method", "priorDist", 
                "priorPar", "range", "D", "parInt", "itemSelect", 
                "infoType", "randomesque", "AP", "constantPatt"), 
                stop = c("rule", "thr", "alpha"), final = c("method", 
                  "priorDist", "priorPar", "range", "D", "alpha", 
                  "parInt"))
            if (is.null(elements)) 
                res <- list(test = FALSE, message = paste("invalid 'type' argument ('", 
                  type, "' is not allowed)", sep = ""))
            else {
                if (length(list) > length(elements)) 
                  res <- list(test = FALSE, mesage = paste("too many elements in ", 
                    deparse(substitute(list)), " for type '", 
                    type, "'", sep = ""))
                else {
                  res <- list(test = TRUE, message = "ok")
                  i <- 0
                  repeat {
                    i <- i + 1
                    if (sum(names(list)[i] == elements) == 0) {
                      res$test <- FALSE
                      break
                    }
                    else {
                      if (i == length(list)) 
                        break
                    }
                  }
                  if (!res$test) {
                    texte <- switch(i, `1` = "st", `2` = "nd", 
                      `3` = "rd")
                    if (is.null(texte)) 
                      texte <- "th"
                    res$message <- paste("invalid name '", names(list)[i], 
                      "' for ", i, texte, " element of '", deparse(substitute(list)), 
                      "'", sep = "")
                  }
                  else {
                    intNames <- c("fixItems")
                    seedNames <- c("seed")
                    singleIntNames <- c("nrItems")
                    numNames <- c("thr", 
                      "alpha", "D", "SETH", "AP")
                    metNames <- c("method")
                    priorNames <- c("priorDist")
                    parNames <- c("priorPar", "range")
                    ruleNames <- c("rule")
                    eapNames <- c("parInt")
                    itemSelectNames <- c("itemSelect")
                    infoTypeNames <- c("infoType")
                    startNames <- c("startSelect")
                    intOnlyNames <- c("randomesque")
                    boolNames <- c("nAvailable")
                    constantNames <- c("constantPatt")
                    severalNumNames<-c("theta")
                    i <- 0
                    repeat {
                      i <- i + 1
                      vect <- c(sum(names(list)[i] == intNames), 
                        sum(names(list)[i] == seedNames), sum(names(list)[i] == 
                          singleIntNames), sum(names(list)[i] == 
                          numNames), sum(names(list)[i] == metNames), 
                        sum(names(list)[i] == priorNames), sum(names(list)[i] == 
                          parNames), sum(names(list)[i] == ruleNames), 
                        sum(names(list)[i] == eapNames), sum(names(list)[i] == 
                          itemSelectNames), sum(names(list)[i] == 
                          infoTypeNames), sum(names(list)[i] == 
                          startNames), sum(names(list)[i] == 
                          intOnlyNames), sum(names(list)[i] == 
                          boolNames), sum(names(list)[i] == constantNames), sum(names(list)[i] == severalNumNames))
                      ind <- (1:16)[vect == 1]
                      prov <- switch(ind, `1` = ifelse(is.null(list[[i]]), 
                        TRUE, ifelse(is.numeric(list[[i]]), ifelse(max(abs(list[[i]] - 
                          round(list[[i]]))) <= 1e-04, TRUE, 
                          FALSE), FALSE)), `2` = ifelse(is.null(list[[i]]), 
                        TRUE, ifelse(is.numeric(list[[i]]), ifelse(length(list[[i]]) == 
                          1, TRUE, FALSE), ifelse(is.na(list[[i]]), 
                          TRUE, FALSE))), `3` = ifelse(is.numeric(list[[i]]) & 
                        length(list[[i]]) == 1, ifelse(abs(list[[i]] - 
                        round(list[[i]])) <= 1e-04, TRUE, FALSE), 
                        FALSE), `4` = (is.numeric(list[[i]]) & 
                        length(list[[i]]) == 1), `5` = (is.list(list[[i]]) == 
                        FALSE & length(list[[i]]) == 1 & sum(list[[i]] == 
                        c("ML", "BM", "WL", "EAP")) == 1), `6` = (is.list(list[[i]]) == 
                        FALSE & length(list[[i]]) == 1 & sum(list[[i]] == 
                        c("norm", "unif", "Jeffreys")) == 1), 
                        `7` = (is.numeric(list[[i]]) & length(list[[i]]) == 
                          2), `8` = (is.list(list[[i]]) == FALSE & 
                          length(list[[i]]) == 1 & sum(list[[i]] == 
                          c("length", "precision", "classification")) == 
                          1), `9` = (is.list(list[[i]]) == FALSE & 
                          length(list[[i]]) == 3 & is.numeric(list[[i]]) == 
                          TRUE & abs(list[[i]][3] - round(list[[i]][3])) <= 
                          1e-04), `10` = (is.list(list[[i]]) == 
                          FALSE & length(list[[i]]) == 1 & sum(list[[i]] == 
                          c("MFI", "bOpt", "thOpt", "MLWI", "MPWI", 
                            "MEI", "MEPV", "KL", "KLP", "GDI", 
                            "GDIP", "progressive", "proportional", 
                            "random")) == 1), `11` = (is.list(list[[i]]) == 
                          FALSE & length(list[[i]]) == 1 & sum(list[[i]] == 
                          c("observed", "Fisher")) == 1), `12` = (is.list(list[[i]]) == 
                          FALSE & length(list[[i]]) == 1 & sum(list[[i]] == 
                          c("bOpt", "thOpt", "MFI", "progressive", 
                            "proportional")) == 1), `13` = ifelse(is.numeric(list[[i]]), 
                          ifelse(max(abs(list[[i]] - round(list[[i]]))) <= 
                            1e-04, ifelse(list[[i]] > 0, TRUE, 
                            FALSE), FALSE), FALSE), `14` = ifelse(is.numeric(list[[i]]), 
                          ifelse(max(abs(list[[i]] - round(list[[i]]))) <= 
                            1e-04, ifelse((min(list[[i]]) == 
                            0 & max(list[[i]]) == 1), TRUE, FALSE), 
                            FALSE), FALSE), `15` = ifelse(is.null(list[[i]]), 
                          TRUE, ifelse(sum(list[[i]] == c("fixed4", 
                            "fixed7", "var", "BM", "EAP", "WL")) == 
                            1, TRUE, FALSE)), `16` = ifelse(is.numeric(list[[i]]), TRUE, FALSE))
                      if (!prov) {
                        res$test <- FALSE
                        res$message <- switch(ind, `1` = paste("element '", 
                          names(list)[i], "' of '", deparse(substitute(list)), 
                          "' must be a vector of integer values or NULL", 
                          sep = ""), `2` = paste("element '", 
                          names(list)[i], "' of '", deparse(substitute(list)), 
                          "' must be a single numeric value or NULL", 
                          sep = ""), `3` = paste("element '", 
                          names(list)[i], "' of '", deparse(substitute(list)), 
                          "' must be a single integer value", 
                          sep = ""), `4` = paste("element '", 
                          names(list)[i], "' of '", deparse(substitute(list)), 
                          "' must be a single numeric value", 
                          sep = ""), `5` = paste("element '", 
                          names(list)[i], "' of '", deparse(substitute(list)), 
                          "' must be either 'ML', 'BM', 'EAP' or 'WL'", 
                          sep = ""), `6` = paste("element '", 
                          names(list)[i], "' of '", deparse(substitute(list)), 
                          "' must be either 'norm', 'unif'or 'Jeffreys'", 
                          sep = ""), `7` = paste("element '", 
                          names(list)[i], "' of '", deparse(substitute(list)), 
                          "' must be a vector of two numeric values", 
                          sep = ""), `8` = paste("element '", 
                          names(list)[i], "' of '", deparse(substitute(list)), 
                          "' must be either 'length', 'precision'", 
                          "\n", " or 'classification'", sep = ""), 
                          `9` = paste("element '", names(list)[i], 
                            "' of '", deparse(substitute(list)), 
                            "' must be a vector of two numeric and", 
                            "\n", " one integer components", 
                            sep = ""), `10` = paste("element '", 
                            names(list)[i], "' of '", deparse(substitute(list)), 
                            "' must be either 'MFI', 'bOpt',", 
                            "\n", " 'MLWI', 'MPWI', 'MEI', 'MEPV', 'KL', 'KLP', 'GDI', 'GDIP'", 
                            "\n", " 'progressive', 'proportional' or 'random'", 
                            sep = ""), `11` = paste("element '", 
                            names(list)[i], "' of '", deparse(substitute(list)), 
                            "' must be either 'observed' or 'Fisher'", 
                            sep = ""), `12` = paste("element '", 
                            names(list)[i], "' of '", deparse(substitute(list)), 
                            "' must be either 'bOpt', 'thOpt', 'progressive', 'proportional' or 'MFI'", 
                            sep = ""), `13` = paste("element '", 
                            names(list)[i], "' of '", deparse(substitute(list)), 
                            "' must be a positive integer value", 
                            sep = ""), `14` = paste("element '", 
                            names(list)[i], "' of '", deparse(substitute(list)), 
                            "' must be a vector of boolean values", 
                            sep = ""), `15` = paste("element '", 
                            names(list)[i], "' of '", deparse(substitute(list)), 
                            "' must be either 'fixed4', 'fixed7', 'var' or NULL", 
                            sep = ""), `16` = paste("element '", 
                            names(list)[i], "' of '", deparse(substitute(list)), 
                            "' must be a vector of numeric values", 
                            sep = ""))
                        break
                      }
                      else {
                        if (i == length(list)) 
                          break
                      }
                    }
                  }
                }
            }
        }
    }
    return(res)
}
