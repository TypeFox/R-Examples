WelchSatter <- function(ufinal, usamp, df, alpha = 0.05)
{
  if (length(ufinal) != 1) stop("'ufinal' (propagated uncertainty) should be a single value!")
  if (length(usamp) != length(df)) stop("different number of values in 'usamp' and 'df'!")
    
  ws.df <- ufinal^4/sum(usamp^4/df)  
  k <- qt(1 - alpha/2, ws.df)
  u.exp <- k * ufinal
  
  return(list(ws.df = ws.df, k = k, u.exp))   
}