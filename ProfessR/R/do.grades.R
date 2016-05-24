`do.grades` <-
function(ggrades, divs=NULL, cut=0, tit="Exam Grades", breaks=length(ggrades)/3, ...)
{
  if(missing(divs)) { divs=NULL }
  if(missing(cut)) { cut=NULL }
  if(missing(tit)) { tit = "Exam Grades" }
  if(missing(breaks)) { breaks=length(ggrades)/3 }

  if(!is.null(cut))
    {
      agrades = ggrades[ggrades>cut]
    }
  else
    {
      agrades = ggrades
    }
  
  BIGN = length(ggrades)
  
  HA = hist(agrades,  main=tit, xlab="Scores", breaks=breaks, ...)

  points(agrades, rep(0, length(agrades)), col='purple')

  
   axis(1, at=seq(from=10*round(min(agrades[agrades>0])/10), to=max(agrades), by=10), labels=TRUE)
   axis(1, at=seq(from=10*round(min(agrades[agrades>0])/10), to=max(agrades), by=2), labels=FALSE)

  m = mean(agrades)
  s = sqrt(var(agrades))

  B = boxplot(agrades,  plot = FALSE)
  Bdivs = c(min(agrades), B$stats[1:4] + diff(B$stats)/2, max(agrades) )

  mstats = c( m-2*s,m-s, m,  m+s, m+2*s)
  abline(v=mstats, lty=2, col=rgb(1, .7, .7) )
  mtext(text="mean", side=3, at=m, line=0)

  abline(v=Bdivs, lty=2, col=rgb(.1, .7, .7) )
  
  u = par("usr")
###	print(u)
  if(is.null(divs))
    {	
      text( c(20, (m-2*s+m-s)/2,(m-s+ m)/2, (m+m+s)/2, (m+s+m+2*s)/2 ), rep(u[4], 5),
           labels=c("E", "D", "C", "B", "A"), pos=1 , col=rgb(1, .7, .7)  )
    }
  box()
  if(is.null(divs))
    {
      mtext(text="Click 4 divisions from LOW to HIGH", side=3, at = u[1], line=2, adj=0)
      K = list(x=NULL, y=NULL)
      
      for(LK in 1:4)
        {
          K1 = locator(type='p', col='blue', n=1)
          abline(v=K1$x, col='blue', lwd=2)
          K = list(x=c(K$x, K1$x) , y=c(K$y, K1$y))
        }

       divs =  c(min(ggrades), K$x, max(ggrades))
    }

  divs = sort(divs)

  abline(v=divs[2:(length(divs)-1)], col=rgb(0,0,1) )	

  ddivs = diff(divs)
  xgrad = divs[1:(length(divs)-1)] + ddivs/2

  xmin = HA$breaks[1]+(divs[2]-HA$breaks[1])/2
  xgrad[1]= xmin
  text(xgrad , rep(u[4], 5), labels=c("E", "D", "C", "B", "A"), pos=1   )
  
##### print(divs) 
  divs = sort(divs)

  cat("Grade divisions:", sep="\n")
  cat(divs, sep="\n")

#####  divisions are determined, now allocate grades.
###  this can be run independent of the divs determination

  KAPPA = getlet(ggrades, divs)	
###ggrades, lett=letts, scor=scores, divs=divs, LETS=LETS, SCRS=SCRS
  LETS = KAPPA$LETS
  letts = KAPPA$lett
  scores = KAPPA$scor	
  SCRS = KAPPA$SCRS
  cat("Letter Grade Distribution:", sep="\n")
  for(i in 1:length(LETS))
    {
      cat(paste(sep=' ', i, LETS[i], length(letts[letts==LETS[i]])), sep="\n")
    }

#####  cat("Numeric Grade Distribution:", sep="\n")
#####  for(i in 1:length(SCRS))
#####    {
 #####     cat(paste(sep=' ', i, SCRS[i], length(scores[scores==SCRS[i]])), sep="\n")
#####    }


  
  print(paste(sep=' ', "Mean Score=",mean(scores)))


  return(list(grades=ggrades, lett=letts, scor=scores, divs=divs, LETS=LETS, SCRS=SCRS, hist=HA, mstats=mstats, bstats=B$stats,  Bdivs=Bdivs))
}

