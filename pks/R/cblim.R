## Delineate knowledge structure by skill function
delineate <- function(skillfun, itemID = 1) {

  item.names <- as.character(unique(skillfun[, itemID]))
# mu <- t(as.matrix(skillfun[, colnames(skillfun) != itemID]))
  mu <- t(as.matrix(skillfun[, -itemID]))  # numeric part
  nskills <- nrow(mu)
  T <- as.matrix(expand.grid(rep(list(0:1), nskills)))
  colnames(T) <- rownames(mu)
  delineated.states <- matrix(0, nrow(T), length(item.names),
                              dimnames = list(NULL, item.names))
  for (i in seq_len(nrow(T))) {
    idx <- skillfun[, itemID][apply(mu * T[i, ] == mu, 2, all)]
    delineated.states[i, as.character(idx)] <- 1
  }
  K <- as.binmat(unique(pat.del.states <- as.pattern(delineated.states)))
  K <- K[order(rowSums(K)), ]
  rownames(K) <- as.pattern(K)
  colnames(K) <- item.names
  classes <- lapply(seq_len(nrow(K)),
                 function(i) rbind(T[grep(rownames(K)[i], pat.del.states), ]))
  names(classes) <- rownames(K)
  list(K = K, classes = classes)
}

