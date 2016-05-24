findMatrix <- function(mat,type,Mats,Struc,Out,g)
{
  StrucUL <- unlist(Struc)
  
  # Define LISREL name:
  if (mat=="LX") lisName <- "LAMBDA-X"
  if (mat=="PH") lisName <- "PHI"
  if (mat=="TD") lisName <- "THETA-DELTA"
  if (mat=="GA") lisName <- "GAMMA"
  if (mat=="LY") lisName <- "LAMBDA-Y"
  if (mat=="PS") lisName <- "PSI"
  if (mat=="TE") lisName <- "THETA-EPS"
  if (mat=="BE") lisName <- "BETA"
  if (mat=="TX") lisName <- "TAU-X"
  if (mat=="TY") lisName <- "TAU-Y"
  if (mat=="AL") lisName <- "ALPHA"
  if (mat=="KA") lisName <- "KAPPA"
  
  Res <- NULL
  
  if (length(Struc[[type]]) > 0 & length(Mats[[mat]]) > 0 & any(Mats[[mat]]>Struc[[type]][g]))
  {
    IndStart <- min(Mats[[mat]][Mats[[mat]]>Struc[[type]][g]])
    if (!any(StrucUL[StrucUL > Struc[[type]][g]] < IndStart))
    {
      # Checn equal next group"
      if (grepl(paste(lisName,"EQUALS",lisName,"IN THE FOLLOWING GROUP"),Out[IndStart]))
      {
        Res <- "NextGroup"
      } else {
        Inds <- matRange(IndStart,lisName,Out)
        Res <- getMatrix(Out[Inds[1]:Inds[2]],lisName,mat %in% c("TD","TE","PS","PH"),mat %in% c("TD","TE","PS","PH"))
      }
    }
  }
  return(Res)
}
