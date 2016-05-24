#######################################################################
##
## Function:  anchors.order()
## Authors:   Jonathan Wand <wand(at)stanford.edu>
##            and Dan Hopkins 
## 
## DESCRIPTION: Calculate frequency of vignette ordering
##              Treatment of ties/intransitivies is a user option
##
## INPUT:
##   data  : object of class anchors.data
##   ties  : one of,
##             "set" 
##             "nominal"
##             "random"
##             "mset"
##
## MODIFIED: 
##   2006-10-02: JW
##   - WAS anchors.plot.R
##   - changed to formula entry, avoid attach()
##   - standardize plots, make optional
##
##   2008-04-28: JW
##   - WAS vignette.order()
##   - data must be a anchors.data object, created usually by anchors.data()
##
#######################################################################
anchors.order.calc <- function(data, ties = c("set","nominal","random","mset"))
{
  ties <- match.arg(ties)

  ## dimensions...
  N         <- NROW(data$z0)
  n.vign    <- NCOL(data$z0)
  
  ## compute relative ordering
  compare.mat <- matrix(0, n.vign,n.vign)
  rownames(compare.mat)  <- 1:n.vign
  namesvec <- c()
  for(j in 1:n.vign){
     namesvec <- c(namesvec, paste("<",j,sep="")) 
  }
  colnames(compare.mat)  <- namesvec
  for(i in 1:n.vign){
    for(j in 1:n.vign){ 
     compare.mat[i,j] <- sum(data$z0[,i] < data$z0[,j])/N
    }
  }
  compare.mat <- compare.mat
  diag(compare.mat) <- NA

  ## place holders
  ## only exist for ties=3
  ngroup <- nviol <- NULL
         
  if(ties=="nominal"){
    ## LEAVES TIED VIGNETTES IN ORDER OF ENTRY
    out <- NULL
    for (i in 1:N)
      out <- c(out, paste(order(data$z0[i,]),sep="",collapse="") )

    z   <- rev(sort(table(out)))
    main<- "Treatment of ties: using nomimnal order of vignettes\n"
  }
  
  if (ties=="random") {
    ##RANDOMLY BREAKS ALL TIES (Dan Hopkins + JW)
    out <- NULL
    for (i in 1:N){
      ## a random offest is added to all vignette responses
      ## JW: use runif NOT rnorm if we must do this
      data$z0[i,] <- data$z0[i,] + runif(n.vign,0,.001)	
      out <- c(out, paste(order(data$z0[i, ]), sep = "", 
        collapse = ""))
    }

    z    <- rev(sort(table(out)))
    main <- "Treatment of ties: randomly choose order"
  }
  
  if (ties=="mset") {
    ##LEAVES TIES IN PLACE (Dan Hopkins)
    out <- rep(NA,N)

    ##define a vector that identifies ties for consecutive pairings
    tied.neigh <- function(x){
      d <-length(x)
      tfvec <- rep(FALSE,(d-1))
      for(k in 1:(d-1)){
        if(x[k]==x[k+1]){
          tfvec[k]<- TRUE
        }
      }
      return(tfvec)
    }

    ##loop across all observations
    for (i in 1:N){
      tt1 <-sort(data$z0[i,])
      wh.ties <-tied.neigh(tt1)
      txtob <- ""
      opts <- sort(unique(data$z0[i,]))

      ##loop across each response category for observation
      for(l in sort(opts)){
        num <-which(data$z0[i,]==l)
        n2 <-length(num)
        txtob2 <- ""

        ##loop across vignette value
        for(m in 1:n2){
          if((! n2 == 1) && m==1)  txtob2<- paste(txtob2,"{",num[m],",",sep="")
          if((! n2 == 1) && m==n2) txtob2<- paste(txtob2,num[m],"}",sep="")
	  if(n2==1)                txtob2<- paste(txtob2,num[m],sep="")
          if((! n2==1) && (! m %in% c(1,n2)))
                                   txtob2<- paste(txtob2,num[m],",",sep="")
        }
        if(l==min(opts)){           txtob<- paste(txtob,txtob2,sep="")
        }else{                      txtob<- paste(txtob,txtob2,sep=",")  }
      }
      out[i] <- txtob
    }
    main <- "Treatment of ties: represent as sets"
  }

  if (ties=="set") {
    ##LEAVES TIES IN PLACE, plus additional calculations (JW)
    ngroup <-  out <- rep(NA,N)

    ss <- function(x) {
      if (length(x) < 2) return(x)
      paste("{", paste(x,collapse=","), "}",sep="")
    }

    for (i in 1:N){
      oo <- order(data$z0[i,])
      rr <- rank(data$z0[i,],ties.method="min")
      or <- rr[oo]

      a <- list()
      for (j in sort(unique(or)) ) {
        a[[ j ]] <- oo[ or == j]
      }

      b <- lapply(a,ss)
      out[i] <- paste(unlist(b),collapse=",")
      ngroup[i] <- length( unique(or))

      aa <- 0
      for (k in 1:(n.vign-1))
        aa <- aa + sum(rr[k] > rr[ -c(1:(k)) ])
      nviol[i] <- aa
      
    }
    main <- "Treatment of ties: represent as sets"
  }


  
  z <- rev(sort(table(out)))
#  z <- table(out)
  
  output <- list(labels=out,
                 ngroup=ngroup,
                 nviol=nviol,
                 vign = data$z0,
                 compare.matrix=compare.mat,
                 freq=z,
                 ties=ties,
                 main=main,
                 N=N)

  
  class(output) <- c("anchors.order")

  return(output)
}  


