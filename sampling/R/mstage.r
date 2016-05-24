mstage<-function (data, stage = c("stratified", "cluster", ""), varnames, 
    size, method = c("srswor", "srswr", "poisson", "systematic"), 
    pik, description = FALSE) 
{
    if (missing(size)) 
        stop("the size argument is missing")
    if (!missing(stage) & missing(varnames)) 
        stop("indicate the stage argument")
    if (!missing(stage)) {
        number = length(stage)
        for (i in 1:length(stage)) if (!(stage[i] %in% c("stratified", 
            "cluster", ""))) 
            stop("the stage argument is wrong")
    }     else number = length(size)
    if (number > 1) {
        if (!missing(varnames)) {
            if (!is.list(size)) 
                stop("the size must be a list")
            size = as.list(size)
            varnames = as.list(varnames)
            size1 = size[[1]]
            varnames1 = varnames[[1]]
            if (method[[1]] %in% c("systematic", "poisson")) 
                pik1 = pik[[1]]
        }
        else {
            size1 = size[[1]]
            varnames1 = NULL
            if (method[[1]] %in% c("systematic", "poisson")) 
                pik1 = pik[[1]]
        }
    }
    else {
        size1 = size
        if(missing(method)) method="srswor"
        else
        if (method %in% c("systematic", "poisson")) 
            pik1 = pik
    }
    if (description) 
        cat("STAGE 1", "\n")
    if (missing(stage)) {
        if (missing(varnames)) 
            if (missing(method)) 
                s = strata(data, stratanames = NULL, size = size1, description)
             else if (method[[1]] %in% c("systematic", "poisson")) 
                s = strata(data, stratanames = NULL, size = size1, 
                  method[[1]], pik = pik1, description)
            else s = strata(data, stratanames = NULL, size = size1, 
                method[[1]], description)
        else s = strata(data, stratanames = NULL, size1, method[[1]], 
            pik = pik1, description)
    }
    else if (stage[1] == "stratified") {
        s = strata(data, varnames1, size1, method="srswor",description)
        dimension_st = table(s$Stratum)
        if(description) cat("Number of strata:",length(dimension_st),"\n")  
                                       }
    else {
        s = cluster(data, varnames1, size1, method[[1]], pik1, description)
        if (is.null(s)) 
            stop("0 selected units in the first stage")
        m = match(varnames1, names(s))
        nl = nlevels(as.factor(s[, m]))
        lev = levels(as.factor(s[, m]))
        if (nl >= 1) {
            dimension_cl = NULL
            for (i in 1:nl) 
                   if(nrow(subset(s,s[, m] == unique(s[,m])[i]))>0)
                                dimension_cl = c(dimension_cl,nrow(subset(s,s[, m] == unique(s[,m])[i])))
                     }
        dimension = dimension_cl
          }
    if (is.null(s)) 
        stop("0 selected units in the first stage")
    if (number > 1) 
        if ((is.element("cluster", stage) & is.element("stratified", 
            stage))) 
            result = getdata(data, s)
        else result = s
    res = list()
    res[[1]] = s
    if (number >= 2) 
        for (j in 2:number) {
            if (description) 
                cat("STAGE ", j, "\n")
            if (!missing(varnames)) {
                if (stage[[j]] == "cluster") {
                  if (stage[[j - 1]] == "stratified") {
                    k = length(dimension_st)
                    s1 = NULL
                    limit = 0
                    dimension = list()
                    if (k >= 1) 
                      for (ii in 1:k) {
                        r = res[[j - 1]][(limit + 1):(limit + 
                          dimension_st[ii]), ]
                        r = getdata(data, r)
                        m = match(varnames[[j]], names(r))
                        if (method[[j]] %in% c("systematic", "poisson")) {
                          index = res[[j - 1]][(limit + 1):(limit + 
                            dimension_st[ii]), ]$ID_unit
                          pikk = pik[[j]][index]
                          if (!is.null(r)) 
                            s3 = cluster(r, clustername = varnames[[j]], 
                              size = size[[j]][ii], method = method[[j]], 
                              pik = pikk, description)
                          else s3 = NULL
                        }
                        else {

                          s3 = cluster(r, clustername = varnames[[j]], 
                            size = size[[j]][ii], method = method[[j]], 
                            description = description)
                        }
                        limit = limit + dimension_st[ii]
                        if (method[[j]] == "srswr") {
                          s3 = cbind.data.frame(r[s3$ID_unit, 
                            m], r[s3$ID_unit, ]$ID_unit, s3$Replicates, 
                            s3$Prob, r[s3$ID_unit, ]$Prob * s3$Prob)
                          colnames(s3) = c(varnames[[j]], "ID_unit", 
                            "Replicates", paste("Prob_", j, "_stage"), 
                            "Prob")
                        }
                        else if(!is.null(s3)) {
                          s3 = cbind.data.frame(r[s3$ID_unit, 
                            m], r[s3$ID_unit, ]$ID_unit, s3$Prob, 
                            result[s3$ID_unit, ]$Prob * s3$Prob)
                          colnames(s3) = c(varnames[[j]], "ID_unit", 
                            paste("Prob_", j, "_stage"), "Prob")
                        }
                        if (!is.null(s3)) {
                          m = match(varnames[[j]], names(s3))
                          for (l in 1:nlevels(as.factor(s3[, 
                            m]))) dimension = c(dimension, table(s3[, 
                            m])[l])
                          s1 = rbind(s1, s3)
                        }
                      }
                  }
                  else if (stage[[j - 1]] == "cluster") {
                 m_cl = match(varnames, names(res[[j-1]]),0)
                 mat=res[[j-1]][, m_cl]
                 nl = nlevels(as.factor(mat))
                 if (nl >= 1)     {
            dimension_cl =NULL
            for (i in 1:nl) 
               if(length(subset(mat, mat==unique(mat)[i]))>0)
                  dimension_cl = c(dimension_cl,length(subset(mat, mat==unique(mat)[i])))
            k = length(dimension_cl)                   
                  }
            else stop("error in the previous stage")
                    s1 = NULL
                    limit = 0
                    dimension = list()
                    if (k > length(size[[j]])) {
                      warning("the number of selected clusters in the previous stage is larger than the size argument")
                      warning("the size 1 is added")
                      size1 = size[[j]]
                      for (i in 1:(k - length(size[[j]]))) size1 = c(size1, 1)
                    }
                    else size1 = size[[j]]
                    if (k >= 1) 
                      for (ii in 1:k) {
                        r = res[[j - 1]][(limit + 1):(limit + 
                          dimension_cl[ii]), ]
                        r = getdata(data, r)
                        m = match(varnames[[j]], names(r))
                        if (method[[j]] %in% c("systematic", "poisson")) {
                          m1 = match(varnames[[j - 1]], names(r))
                          m2 = match(varnames[[j - 1]], names(data))
                          mm = match(r[1, m1], levels(factor(data[, 
                            m2])))
                          pikk = as.numeric(pik[[j]][[mm]])
                          if (!is.null(r)) 
                            s3 = cluster(r, clustername = varnames[[j]], 
                              size = size1[[ii]], method = method[[j]], 
                              pik = pikk, description)
                          else s3 = NULL
                        }
                        else s3 = cluster(r, clustername = varnames[[j]], 
                          size = size1[ii], method = method[[j]], 
                          pik, description)
                        limit = limit + dimension_cl[ii]
                        if (method[[j]] == "srswr") {
                          s3 = cbind.data.frame(r[s3$ID_unit, 
                            m], r[s3$ID_unit, ]$ID_unit, s3$Replicates, 
                            s3$Prob, r[s3$ID_unit, ]$Prob * s3$Prob)
                          colnames(s3) = c(varnames[[j]], "ID_unit", 
                            "Replicates", paste("Prob_", j, "_stage"), 
                            "Prob")
                        }
                        else if (!is.null(s3)) {
                          s3 = cbind.data.frame(r[s3$ID_unit, 
                            m], r[s3$ID_unit, ]$ID_unit, s3$Prob, 
                            result[s3$ID_unit, ]$Prob * s3$Prob)
                          colnames(s3) = c(varnames[[j]], "ID_unit", 
                            paste("Prob_", j, "_stage"), "Prob")
                        }
                        if (!is.null(s3)) {
                          m = match(varnames[[j]], names(s3))
                          for (l in 1:nlevels(as.factor(s3[, 
                            m]))) dimension = c(dimension, table(s3[, 
                            m])[l])
                          s1 = rbind(s1, s3)
                        }
                      }
                  }
                }
                else if (j > 1) {
                  k = length(dimension)
                  s1 = NULL
                  limit = 0
                  count = 0
                  if (k > length(size[[j]])) {
                    warning("the number of selected clusters at the previous stage is larger than the size argument")
                    warning("the size 1 is added")
                    size1 = size[[j]]
                    for (i in 1:(k - length(size[[j]]))) size1 = c(size1,1)
                  }
                  else size1 = size[[j]]
                  if (k >= 1) 
                    for (i in 1:k) for (ii in 1:length(dimension[[i]])) {
                      r = res[[j - 1]][(limit + 1):(limit + dimension[[i]][ii]), 
                        ]
                      count = count + 1
                      if (method[[j]] %in% c("systematic", "poisson")) {
                        index = res[[j - 1]][(limit + 1):(limit + 
                          dimension[[i]][ii]), ]$ID_unit
                        pikk = pik[[j]][index]
                        if (!is.null(r)) 
                          s2 = strata(r, NULL, size = size1[count], 
                            method = method[[j]], pik = pikk, description)
                        else s2 = NULL
                      }
                      else s2 = strata(r, NULL, size = size1[count], 
                        method = method[[j]], pik, description)
                      limit = limit + dimension[[i]][ii]
                      if (method[[j]] == "srswr") {
                        s2 = cbind.data.frame(r[s2$ID_unit, ]$ID_unit, 
                          s2$Replicates, s2$Prob, r[s2$ID_unit, 
                            ]$Prob * s2$Prob)
                        colnames(s2) = c("ID_unit", "Replicates", 
                          paste("Prob_", j, "_stage"), "Prob")
                      }
                      else if (!is.null(s2)) {
                        s2 = cbind.data.frame(r[s2$ID_unit, ]$ID_unit, 
                          s2$Prob, r[s2$ID_unit, ]$Prob * s2$Prob)
                        colnames(s2) = c("ID_unit", paste("Prob_", 
                          j, "_stage"), "Prob")
                      }
                      if (!is.null(s2)) 
                        s1 = rbind(s1, s2)
                    }
                }
            }
            else {
                if (missing(stage)) {
                  if (missing(method)) 
                    s1 = strata(result, stratanames = NULL, size = size[[j]], 
                      description = description)
                  else if (method[[j]] == "poisson" | method[[j]] == "systematic") 
                    s1 = strata(result, stratanames = NULL, size = size[[j]], 
                      method = method[[j]], pik = pik[[j]], description = description)
                  else s1 = strata(result, stratanames = NULL, 
                    size = size[[j]], method = method[[j]], description = description)
                  if (method[[j]] == "srswr") {
                    s1 = cbind.data.frame(result[s1$ID_unit, 
                      ]$ID_unit, s1$Replicates, s1$Prob, result[s1$ID_unit, 
                      ]$Prob * s1$Prob)
                    colnames(s1) = c("ID_unit", "Replicates", 
                      paste("Prob_", j, "_stage"), "Prob")
                  }
                  else {
                    s1 = cbind.data.frame(result[s1$ID_unit, 
                      ]$ID_unit, s1$Prob, result[s1$ID_unit, 
                      ]$Prob * s1$Prob)
                    colnames(s1) = c("ID_unit", paste("Prob_", 
                      j, "_stage"), "Prob")
                  }
                }
            }
            if (!is.null(s1)) {
                result = s1
                res[[j]] = result
            }
            else number = number - 1
        }
    if (!is.null(names(res[[1]]))) {
        m = match("Prob", names(res[[1]]))
        names(res[[1]])[m] = "Prob_ 1 _stage"
    }
    names(res) = c(1:number)
    res
}
