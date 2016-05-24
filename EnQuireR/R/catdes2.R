
"catdes2"=function (donnee, num.var, proba = 0.05) 
{
    lab.sauv <- lab <- colnames(donnee)
    quali = NULL
    for (i in 1:length(lab)) {
        lab[i] = gsub(" ", ".", lab[i])
        if (is.factor(donnee[, i])) {
            if (any(is.na(donnee[, i]))) {
                levels(donnee[, i]) <- c(levels(donnee[, i]), 
                  "NA")
                donnee[, i][is.na(donnee[, i])] <- "NA"
            }
            if (levels(donnee[, i])[1] == "") 
                levels(donnee[, i])[1] = "NA"
            if (i != num.var) 
                quali = c(quali, i)
        }
    }
    quanti = (1:ncol(donnee))[-c(quali, num.var)]
    if (length(quanti) == 0) 
        quanti = NULL
    colnames(donnee) = lab
    res = list()
    nb.modalite <- length(levels(donnee[, num.var]))
    nb.quali = length(quali)
    old.warn = options("warn")
    if (nb.quali > 0) {
        options(warn = -1)
        Test.chi = matrix(NA, nrow = nb.quali, ncol = 2)
        marge.li = xtabs(~donnee[, num.var])
        nom = tri = structure(vector(mode = "list", length = nb.modalite), 
            names = levels(donnee[, num.var]))
        for (i in 1:nb.quali) {
            Table <- xtabs(~donnee[, num.var] + donnee[, quali[i]])
            marge.col = xtabs(~donnee[, quali[i]])
            Test <- chisq.test(donnee[, num.var], donnee[, quali[i]], correct = FALSE)
            Test.chi[i, 1] <- Test$p.value
            Test.chi[i, 2] <- Test$parameter
            for (j in 1:nlevels(donnee[, num.var])) {
                for (k in 1:nlevels(donnee[, quali[i]])) {
                  aux2 = Table[j, k]/marge.li[j]
                  aux3 = marge.col[k]/sum(marge.col)
                  if (aux2 > aux3) 
                    aux4 = phyper(Table[j, k] - 1, marge.li[j], 
                      sum(marge.li) - marge.li[j], marge.col[k], 
                      lower.tail = FALSE)
                  else aux4 = phyper(Table[j, k], marge.li[j], 
                    sum(marge.li) - marge.li[j], marge.col[k])
                  if (aux4 < proba) {
                    aux5 = (1 - 2 * as.integer(aux2 > aux3)) * 
                      qnorm(aux4)
                    aux1 = Table[j, k]/marge.col[k]
                    tri[[j]] = rbind(tri[[j]], c(aux1, aux2, 
                      aux3, aux4, aux5))
                    nom[[j]] = rbind(nom[[j]], c(levels(donnee[, 
                      quali[i]])[k], colnames(donnee)[quali[i]]))
                  }
                }
            }
            rownames(Test.chi) = colnames(donnee)[quali]
        }
        if (nrow(matrix(Test.chi, ncol = 2)) > 1) {
            if (sum(Test.chi[, 1] < proba) == 1) {
                nomaux = rownames(Test.chi[order(Test.chi[, 1]), 
                  ])[1]
                Test.chi = matrix(Test.chi[Test.chi[, 1] < proba, 
                  ], ncol = 2)
                rownames(Test.chi) = nomaux
            }
            else Test.chi = Test.chi[Test.chi[, 1] < proba, ]
        }
        else if (Test.chi[, 1] > proba) 
            Test.chi = NULL
        if (!is.null(Test.chi)) {
            if (nrow(matrix(Test.chi, ncol = 2)) > 1) {
                oo = order(Test.chi[, 1])
                Test.chi = Test.chi[oo, ]
            }
            colnames(Test.chi) = c("P.value", "df")
        }
        for (j in 1:nb.modalite) {
            if (!is.null(tri[[j]])) {
                oo = rev(order(tri[[j]][, 5]))
                tri[[j]] = tri[[j]][oo, ]
                nom[[j]] = nom[[j]][oo, ]
                if (nrow(matrix(tri[[j]], ncol = 5)) > 1) 
                  rownames(tri[[j]]) = paste(nom[[j]][, 2], nom[[j]][, 
                    1], sep = "=")
                else {
                  tri[[j]] = matrix(tri[[j]], ncol = 5)
                  rownames(tri[[j]]) = paste(nom[[j]][2], nom[[j]][1], 
                    sep = "=")
                }
                colnames(tri[[j]]) = c("Cla/Mod", "Mod/Cla", 
                  "Global", "p.value", "V-test")
            }
        }
        res$test.chi = Test.chi
        res$category = tri
    }
    if (!is.null(quanti)) {
        nom = result = structure(vector(mode = "list", length = nb.modalite), 
            names = levels(donnee[, num.var]))
        for (i in 1:length(quanti)) {
            moy.mod = by(donnee[, quanti[i]], donnee[, num.var], 
                mean, na.rm = TRUE)
            n.mod = summary(donnee[, num.var])
            sd.mod = by(donnee[, quanti[i]], donnee[, num.var], 
                sd, na.rm = TRUE)
            sd.mod = sd.mod * sqrt((n.mod - rep(1, nb.modalite))/n.mod)
            moy = mean(donnee[, quanti[i]], na.rm = TRUE)
            et = sd(donnee[, quanti[i]], na.rm = TRUE) * sqrt(1 - 
                1/sum(n.mod))
            for (j in 1:nb.modalite) {
                v.test = (moy.mod[j] - moy)/et * sqrt(n.mod[j])/sqrt((sum(n.mod) - 
                  n.mod[j])/(sum(n.mod) - 1))
                if (!is.na(v.test)) {
                  if (abs(v.test) > -qnorm(proba/2)) {
                    result[[j]] = rbind(result[[j]], c(v.test, 
                      moy.mod[j], moy, sd.mod[j], et))
                    nom[[j]] = c(nom[[j]], colnames(donnee)[quanti[i]])
                  }
                }
            }
        }
        for (j in 1:nb.modalite) {
            if (!is.null(result[[j]])) {
                oo = rev(order(result[[j]][, 1]))
                result[[j]] = result[[j]][oo, ]
                nom[[j]] = nom[[j]][oo]
                if (nrow(matrix(result[[j]], ncol = 5)) > 1) {
                  rownames(result[[j]]) = nom[[j]]
                  colnames(result[[j]]) = c("v.test", "Mean in category", 
                    "Overall mean", "sd in category", "Overall sd")
                }
                else {
                  result[[j]] = matrix(result[[j]], ncol = 5)
                  rownames(result[[j]]) = nom[[j]]
                  colnames(result[[j]]) = c("v.test", "Mean in category", 
                    "Overall mean", "sd in category", "Overall sd")
                }
            }
        }
        res$quanti = result
    }
    options(old.warn)
    return(res)
}
