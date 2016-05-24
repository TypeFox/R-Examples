`TOMOinfo` <-
function(MOD, PLOT=FALSE)
{
if(missing(PLOT)) { PLOT=FALSE }
 newpar = NULL
 
 N = length(MOD$MOD)

## for(i in 1:N ) { print(any(is.na(MOD$MOD[[i]]))) }

 SLAY = NULL

 
  for( i in 1:N )
    {
      LAY = as.vector(MOD$MOD[[i]])

      if(all(is.na(LAY)))
        {
          next
        }
      ##  r1 = range(LAY[LAY>0.0], na.rm=TRUE)

      s1 = boxplot.stats(LAY, do.conf = FALSE, do.out = FALSE)
      m1 = mean(LAY, na.rm = TRUE)
      m2 = sqrt(var(LAY, na.rm = TRUE))
      
      
      SLAY = rbind(SLAY, c(i, s1$stats , m1, m2))

     ## A = format( c(SLAY[[i]]$bstats, SLAY[[i]]$mstats), digits=5 )
      
    ##   print(A)
     ##  newpar = rbind(newpar, jstats(LAY))
     ##  print( paste( sep=" ", i, r1[1],r1[2], 1/r1[1],  1/r1[2] ) )
    
    }

 rowname = paste(sep="-", "LAY", SLAY[,1])
 colname = c("I", "q1","q2","med", "q4","q5",  "mean", "sd")
 
 dimnames(SLAY)=list(rowname ,colname  )

 OLAY = data.frame(SLAY)
 print(OLAY)

if(PLOT) { boxplot(MOD$MOD) }
##  if(PLOT) { tufteboxplot(MOD$MOD) }
  
  invisible(OLAY)
}

