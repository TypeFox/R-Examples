cluster<-function (data, clustername, size, method = c("srswor", "srswr", "poisson", "systematic"), pik, description = FALSE) 
{
    if (size == 0) 
        stop("the size is zero")
    if (missing(method)) {
        warning("the method is not specified; by default, the method is srswor")
        method = "srswor"
    }
    if (!(method %in% c("srswor", "srswr", "poisson", "systematic"))) 
        stop("the name of the method is wrong")
    if (method %in% c("poisson", "systematic") & missing(pik)) 
        stop("the vector of probabilities is missing")
    data = data.frame(data)
    index = 1:nrow(data)
    if (missing(clustername)) {
        if (method == "srswor") 
            result = data.frame(index[srswor(size, nrow(data)) == 
                1], rep(size/nrow(data), size))
        if (method == "srswr") {
            s = srswr(size, nrow(data))
            st = s[s != 0]
            l = length(st)
            result = data.frame(index[s != 0])
            result = cbind.data.frame(result, st, prob = rep(1-(1-1/nrow(data))^size,l))
            colnames(result) = c("ID_unit", "Replicates", "Prob")
                               }
        if (method == "poisson") {
            pikk = inclusionprobabilities(pik, size)
            s = (UPpoisson(pikk) == 1)
            if (length(s) > 0) 
                result = data.frame(index[s], pikk[s])
            if (description) 
                cat("\nPopulation total and number of selected units:", 
                  nrow(data), length(s), "\n")
        }
        if (method == "systematic") {
            pikk = inclusionprobabilities(pik, size)
            s = (UPsystematic(pikk) == 1)
            result = data.frame(index[s], pikk[s])
        }
        if (method != "srswr") 
            colnames(result) = c("ID_unit", "Prob")
        if (description) 
            cat("\nPopulation total and number of selected units:", 
                nrow(data), sum(size), "\n")
    }
    else {
        data = data.frame(data)
        m = match(clustername, colnames(data))
        if (length(m) > 1) 
            stop("there are too many specified variables as clusters")
        if (is.na(m)) 
            stop("the cluster name is wrong")
        x1 = factor(data[, m]) 
        result = NULL
        if (nlevels(x1) == 0) 
            stop("the cluster variable has 0 modalities")
        else {
            nr_cluster = nlevels(x1)
            if (method == "srswor") {
                s = as.data.frame(levels(x1)[srswor(size, nr_cluster) == 1])
                names(s) = c("cluster")
                r = cbind.data.frame(index, data[, m])
                names(r) = c("index", "cluster")
                r = merge(r, s, by.x = "cluster", by.y = "cluster", sort = TRUE)
                result = cbind.data.frame(r, rep(size/nr_cluster, nrow(r)))
            }
            if (method == "srswr") {
                s = srswr(size, nr_cluster)
                st = cbind.data.frame(levels(x1)[s != 0], s[s != 0])
                names(st) = c("cluster", "repl")
                r = cbind.data.frame(index, data[, m])
                names(r) = c("index", "cluster")
                r = merge(r, st, by.x = "cluster", by.y = "cluster")
               result = cbind.data.frame(r, rep(1-(1-1/nr_cluster)^size, nrow(r)))
                            }
            if (method == "systematic") {
                pikk = inclusionprobabilities(pik, size)
                s = (UPsystematic(pikk) == 1)
                st = cbind.data.frame(levels(x1)[s], pikk[s])
                names(st) = c("cluster", "prob")
                r = cbind.data.frame(index, data[, m])
                names(r) = c("index", "cluster")
                result = merge(r, st, by.x = "cluster", by.y = "cluster")
            }
            if (method == "poisson") {
                pikk = inclusionprobabilities(pik, size)
                s = (UPpoisson(pikk) == 1)
                if (any(s)) {
                  st = cbind.data.frame(levels(x1)[s], pikk[s])
                  names(st) = c("cluster", "prob")
                  r = cbind.data.frame(index, data[, m])
                  names(r) = c("index", "cluster")
                  result = merge(r, st, by.x = "cluster", by.y = "cluster")
                  if (description) {
                    cat("Number of selected clusters:", sum(s), 
                      "\n")
                    cat("\nPopulation total and number of selected units:", 
                      nrow(data), nrow(result), "\n")
                  }
                }
                else {
                  if (description) {
                    cat("Number of selected clusters: 0\n")
                    cat("Population total and number of selected units:", 
                      nrow(data), 0, "\n")
                  }
                  result = NULL
                }
            }
            if (method == "srswr") {
                colnames(result) = c(clustername, "ID_unit", 
                  "Replicates", "Prob")
                if (description) {
                  cat("Number of selected clusters:", length(s[s != 
                    0]), "\n")
                  cat("Population total and number of selected units", 
                    nrow(data), nrow(result), "\n")
                }
            }
            else if (!is.null(result)) 
                colnames(result) = c(clustername, "ID_unit", 
                  "Prob")
            if (description & !(method %in% c("poisson", "srswr"))) {
                cat("Number of selected clusters:", size, "\n")
                cat("Population total and number of selected units", 
                  nrow(data), nrow(result), "\n")
            }
        }
    }
    result
}




