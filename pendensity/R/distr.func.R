distr.func <- function(yi=NULL,obj,help.env=distr.func.help(obj)) {
  val <-0
  q <- obj$splines$q
  K <- obj$splines$K
  base.den <- obj$splines$base.den
  len.b <- length(base.den[,1])
  knots.val <- obj$splines$knots.val
  weight <- obj$results$ck
  x <- obj$values$x

  x.factor <- obj$values$covariate$x.factor
  func.env <- new.env()
  if(is.null(obj$values$x))  len.x.fac <- 1  else len.x.fac <- length(x.factor[,1])

  assign("allx",len.x.fac,envir=func.env)

  if(is.null(x)) assign("weight1",weight,envir=func.env)
  else {  
    for(i in 1:len.x.fac) {
      name <- paste("weight",i,sep="")
      assign(name,weight[i,],envir=func.env)
    }
  }
  assign("combi",obj$combi,envir=func.env)
  

  sum <- c(0,0)

  for(k in 1:len.x.fac) {
    weight <-  get(paste("weight",k,sep=""),envir=func.env)
     if(!is.null(yi)) {
      for(j in 1:(length(knots.val$help)-1)) {
        if(knots.val$help[j] <= yi & yi < knots.val$help[j+1]) val <- j
        if(val==length(base.den[,1])+1) val <- val-1
      }
      
      for(i in q:val) {
        track1 <- seq(1,q)+(i-q)
        track2 <- seq(q,1)
        if(i < val) {
          for(j in 1:q) {
            sum[k] <- sum[k] + poly.part(track1[j],track2[j],knots.val,help.env,q,yi=NULL)*weight[track1[j]]
          }
        }
        if(i==val) {
          for(l in 1:q) {
            if(!is.null(yi)) sum[k] <- sum[k] + poly.part(track1[l],track2[l],knots.val,help.env,q,yi,poly=FALSE)*weight[track1[l]]
          }
        }
      }
    }   
  
    else {

      val3 <- seq(q,length(knots.val$help)-q)
      for(i in 1:len.x.fac) {
        weight <- get(paste("weight",i,sep=""),envir=func.env)
        for(j in 1:length(val3)) {
          sum <- 0
          for(w in q:val3[j]) {
            track1 <- seq(1,q)+(w-q)
            track2 <- seq(q,1)
            if(w < val3[j]) {
              for(k in 1:q) {
                sum <- sum + poly.part(track1[k],track2[k],knots.val,help.env,q)*weight[track1[k]]
              }
            }
            if(w==(val3[j])) {
              term <- c()
              for(l in 1:q)  {
                term[l] <- paste("(",poly.part(track1[l],track2[l],knots.val,help.env,q,poly=TRUE),")",sep="")
                term[l] <- paste(term[l],weight[track1[l]],sep=" * ")
              }
              term[q+1] <- sum
              term <- paste(term,collapse="+")
              assign(paste("distr.func",i,".",j,sep=""),paste("obj <-function(x){",term,"}"),envir=func.env)
            }
          }
        }
      }
      return(func.env)
    }
  }
  if(!is.null(yi)) return(sum)
}
