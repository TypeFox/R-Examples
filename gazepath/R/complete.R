## This function returns the row(s) with complete fixations or saccades
complete <- function(sim, x = c('f', 's')){
  if(x == 'f') y <- 's'
  if(x == 's') y <- 'f'
  if(x != 'f' & x != 's') stop('Indicate if the function should return complete fixations (f) or saccades (s)')
  
  fix <- which(sim[,1] == x)
  if(fix[1] == 1) fix <- fix[-1]
  row <- as.numeric(rownames(sim))
  
  return(fix[which(sim[fix - 1, 1] == y & sim[fix + 1, 1] == y & (row[fix + 1] - row[fix - 1]) == 2)])
}