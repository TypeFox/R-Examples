summary.cpg <-
function(object,...) {
  print.cpg(object)
  original<-par()$mfrow
  cpg.window<-c(1,1)
  if (!is.factor(object$indep)) {
    cpg.window<-c(2,2)
      }
  if(nrow(object$results)>1) {
    cpg.window<-c(2,2)
      }
  else {
     cpg.window<-c(2,1)
      }
  par(mfrow=cpg.window)
  plot(object)
  if(sum(cpg.window)==3) scatterplot(object,tplot=T)
  if(sum(cpg.window)==4) {
     if(!is.factor(object$indep)) {
        plot(object,tplot=T)
        plot(object,classic=T)
        plot(object,tplot=T,classic=T)
        }
     else {
        plot(object,classic=T)
        scatterplot(object,c(1:2))
          }

        }
  par(mfrow=original)
  }
