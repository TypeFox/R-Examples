# RG TM more memeory efficient
`rg_tm` <-
function (ni = 100, np = 100, ties = 300, weights = 1, seed = NULL) {
  if (!is.null(seed)) 
    set.seed(as.integer(seed))
  if (ties > 1) {
    net <- data.frame(i = sample(1:ni, (ties * 1.5), replace = TRUE), p = sample(1:np, (ties * 1.5), replace = TRUE))
    net <- net[!duplicated(net[, 1:2]), ]
    net <- net[1:ties, ]
  } else {
    net <- stats::runif(ni * np)
    net <- net <= ties
    net <- matrix(data = net, nrow = ni, ncol = np)
    net <- which(net, arr.ind = TRUE)
  }
  if(length(weights)==1 & weights[1] == 1) {
    return(as.tnet(net[, 1:2], type = "binary two-mode tnet"))
  } else {
    net <- data.frame(net, w = sample(weights, nrow(net), replace = TRUE))
    net <- as.tnet(net, type = "weighted two-mode tnet")
    return(net)
  }
}
