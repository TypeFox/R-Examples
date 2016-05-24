#rankscore on a given user.

getrankscore <- function(data, rec_items, fold_items_x_user, goodRating, alpha){
  
  #extract index of the hits
  match_TS <- which(rec_items %in% fold_items_x_user)
  
  if(length(match_TS) == 0 ) return(0)
  
  rankscoreMAX <- getrankscoreMAX(length(match_TS), alpha)
  
  rankscore_user <- (match_TS - 1)
  rankscore_user <- -rankscore_user/alpha
  rankscore_user <- 2^rankscore_user
  rankscore_user <- sum(rankscore_user)
  
  rankscore_user/rankscoreMAX
}

getrankscoreMAX<- function(n,alpha){
  
  rankscoreMAX <- 0
  
  rankscoreMAX <- 1/2^((c(1:n) - 1)/alpha)
  
  sum(rankscoreMAX)
}
