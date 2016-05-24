depths <-
function(depths, horizon.names) {
  ##remove spaces in the data
  tmp1 <- gsub(" ", "", depths)
  tmp_depths <- data.frame(upper1=rep(NA, length(tmp1)), upper2=rep(NA, length(tmp1)), lower1=rep(NA, length(tmp1)),
                           lower2=rep(NA, length(tmp1)))
  ##
  for (a in 1:length(tmp1)) {
    sep <- unlist(strsplit(tmp1[a], '-'))
    upper <- as.numeric(unlist(strsplit(sep[1], '/')))
    lower <- as.numeric(unlist(strsplit(sep[2], '/')))
    if (length(upper)==1) {tmp_depths[a, 1:2] <- as.numeric(c(upper, upper))
                         } else {tmp_depths[a, 1:2] <-upper}
    if (length(lower)==1) {tmp_depths[a, 3:4] <- as.numeric(c(lower, lower))
                         } else {tmp_depths[a, 3:4] <-lower}
  }
    tmp_depths <- -tmp_depths
    tmp_dp <- cbind(apply(tmp_depths[, 1:2], 1, mean), apply(tmp_depths[, 3:4], 1, mean))
    final_depths <- apply(tmp_dp, 1, mean)
  final <- data.frame(depths=final_depths, names=horizon.names)
  return(final)
}
