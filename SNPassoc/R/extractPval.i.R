`extractPval.i` <-
function (i, x, pos, models) 
{
    tt<-x[[i]]
    if (is.null(nrow(tt)))
      control<-1000
    else
      control<-nrow(tt)     

    if(length(models)*2>control)
     {
       tt <- tt[, pos]
       ans <- tt[!is.na(tt)][1]
       ans <- c(NA,ans, rep(NA, length(models) - 1))
     }
    else
     {
      if (!is.null(dim(tt))) {
        if (length(models) == 1) {
            tt <- tt[, pos]
            ans <- c(NA, tt[!is.na(tt)][1])
        }
        else {
            tt <- tt[, pos]
            ans <- c(NA, tt[!is.na(tt)])
            if ((length(ans) - 1) < length(models)) 
#                ans <- c(ans, rep(NA, length(models) - 1))
                ans <- c(ans, rep(NA, length(models) - (length(ans) - 1)))
          }
       }
      else if (!is.na(charmatch("Geno", tt))) 
        ans <- c(tt[1], rep(NA, length(models)))
      else ans <- c("Monomorphic", rep(NA, length(models)))
    }
    ans
}

