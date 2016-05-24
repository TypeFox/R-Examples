hintSudoku <- function(z, i,j) {
  tmp <- setdiff(1:9, z[i, ])
  if (length(tmp)) {tmp2 <- 1:9; tmp2[-tmp] <- " "} else tmp2 <- " "
  line1 <- paste("row:", paste(tmp2,collapse=" "))
  tmp <- setdiff(1:9, z[ ,j])
  if (length(tmp)) {tmp2 <- 1:9; tmp2[-tmp] <- " "} else tmp2 <- " "
  line2 <- paste("col:", paste(tmp2,collapse=" "))
  if (j < 4) {
    if (i < 4)      tmp <- setdiff(1:9, z[ 1:3, 1:3 ])
    else if (i < 7) tmp <- setdiff(1:9, z[ 4:6, 1:3 ])
    else            tmp <- setdiff(1:9, z[ 7:9, 1:3 ])
  } else if (j < 7) {
    if(i < 4)       tmp <- setdiff(1:9, z[ 1:3, 4:6 ])
    else if (i < 7) tmp <- setdiff(1:9, z[ 4:6, 4:6 ])
    else            tmp <- setdiff(1:9, z[ 7:9, 4:6 ])
  } else {
    if (i < 4)      tmp <- setdiff(1:9, z[ 1:3, 7:9 ])
    else if (i < 7) tmp <- setdiff(1:9, z[ 4:6, 7:9 ])
    else            tmp <- setdiff(1:9, z[ 7:9, 7:9 ])
  }
  if (length(tmp)) {tmp2 <- 1:9; tmp2[-tmp] <- " "} else tmp2 <- " "
  line3 <- paste("grp:", paste(tmp2,collapse=" "))
  paste(line1, line2, line3, "\n", sep="\n")
}
