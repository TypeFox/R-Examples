# Internal Function
cueCoding = function(cues=c("hello", "world"), maxn=1, adjacent=FALSE)
{
  if(adjacent & maxn!=2)
    { cat("'maxn=2' since 'adjacent=TRUE'.\n")
      maxn=2
    }

  if(length(cues)<maxn) maxn=length(cues)

  ngram.fnc = function(v, n) {
    len = length(v)
    ng = NULL
    combinations <- combn(len,n)
    for (i in 1:ncol(combinations)) {
       ng = c(ng, paste(v[combinations[,i]], collapse=""))
    }
    return(paste(ng, collapse="_"))
  }
  
  grams <- paste(cues, collapse="_")

  if(maxn==1)
    return(grams)
  else
    if(!adjacent)
      for(i in 2:maxn)
        { vv <- ngram.fnc(cues,i)
          grams = paste(grams, vv, sep="_")
        }
    else
      if(maxn==2)
        grams <- paste(c(grams,
                       paste(apply(cbind(cues[1:length(cues)-1],
                                   cues[2:length(cues)]),
                                   1, function(x) paste(x,collapse="")),
                       collapse="_")),
                collapse="_")

  return(grams)
}


