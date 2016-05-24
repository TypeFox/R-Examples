GMask <- function(ES, Core, Y){
  Y <- as.factor(Y)
  levels(Y) = c(1,2)
  ES[, Y == 1] = (ES[, Y == 1] >= Core[,1] & ES[, Y == 1] <= Core[,2]) & !(ES[, Y == 1] >= Core[,3] & ES[, Y == 1] <= Core[,4])
  ES[, Y == 2] = (ES[, Y == 2] >= Core[,3] & ES[, Y == 2] <= Core[,4]) & !(ES[, Y == 2] >= Core[,1] & ES[, Y == 2] <= Core[,2])
  mask = ES
  mode(mask) = "integer"
  return(mask)
}