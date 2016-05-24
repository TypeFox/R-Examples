# Automatically generated from all.nw using noweb

ped2df <- function(ped) {
  df <- data.frame(id=ped$id, findex=ped$findex, mindex=ped$mindex, sex=ped$sex)
  if(!is.null(ped$affected))
    df$affected = ped$affected
  if(!is.null(ped$status))
    df$status = ped$status
  return(df)
}

