patternMatch.strong <- function(X.discret,unique.pat)
{

# filter full zero lines
  X.discret <- X.discret[which(apply(X.discret,1,sum) > 0),]
  X.string <- patternToString(X.discret)

  n.random <- array(0,length(unique.pat))
  for ( i in 1:length(unique.pat))
  {
   n.random[i] <- length(which(X.string %in% unique.pat[i]))
 }

return(n.random)

}
