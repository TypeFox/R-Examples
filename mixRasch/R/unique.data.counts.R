`unique.data.counts` <-
function(itemdata){

  n.i <- ncol(itemdata)
  
  #Makes sure categories start at zero
  makeZero <- min(itemdata,na.rm=TRUE)
  itemdata <- itemdata - makeZero

  temp <- itemdata
  temp <- cbind(temp,rowSums(itemdata))
    # used later to fill in missing table counts
    record.miss <- is.na(temp[,n.i + 1]) 
    #

  temp[is.na(temp)] <- "."
  temp   <- apply(temp,1, function(XXX) paste(XXX,sep="",collapse=""))

  people <- substr(temp,1,n.i) 
  n.x    <- table(temp)
  x.x    <- substr(names(n.x),1,n.i)
  r      <- suppressWarnings(as.numeric(substr(names(n.x),n.i+1,1000000000)))

  max.r <- sum( apply(itemdata,2,max,na.rm=T) )
  degen <- (r == max.r) | (r == 0)

  orderer <- order(r,decreasing = FALSE)
  n.x     <- as.vector(n.x[orderer])
  x.x     <- x.x[orderer]
  r       <- r[orderer] 
  degen   <- degen[orderer] 
 
  n.r <- sapply(r,function(XXX) sum(n.x[ ! is.na(r) & r==XXX ]) )

  #this section creates the table to be used with missing data
  if(sum(record.miss) > 0){
   #function for doing it
   missSTUFF <- function(XXX) {
                             use.items <- ! is.na(XXX)  
                             this.x    <- paste(ifelse(is.na(XXX),".",XXX),sep="",collapse="")
                             this.r    <- sum(XXX[use.items])
                             r.miss    <- rowSums( itemdata[,use.items] )
                             max.r    <- sum( apply(itemdata[,use.items],2,max,na.rm=T) )
                             degen    <- (this.r == max.r) | (this.r == 0)
                             n.r.miss  <- length(r.miss[r.miss==this.r & ! is.na(r.miss)])
                             c(this.x, this.r, n.r.miss, degen) }

   if(sum(record.miss)==1) fill.NA <- array(missSTUFF(itemdata[record.miss,]),dim=c(1,4))
   else                    fill.NA <- t(apply(itemdata[record.miss,], 1, missSTUFF))
 
   for(XXX in 1:nrow(fill.NA)){
    r[x.x==fill.NA[XXX,1]]   <- as.numeric(fill.NA[XXX,2])
    n.r[x.x==fill.NA[XXX,1]] <- as.numeric(fill.NA[XXX,3])
    degen[x.x==fill.NA[XXX,1]] <- fill.NA[XXX,4]
   }
  } # end conditional for any missing data

  # Makes a item vector from unique responses
  suppressWarnings(x.i <- matrix(as.numeric(unlist(strsplit(x.x,""))),ncol=n.i,byrow=T))

  colnames(x.i) <- colnames(itemdata)

  degen <- degen == "TRUE"

  degen.x.x <- x.x[degen]
    degen.r <- r[degen]
  degen.n.x <- n.x[degen]
  degen.n.r <- n.r[degen]
  degen.x.i <- x.i[degen,]

  x.x <- x.x[! degen]
    r <- r[! degen]
  n.x <- n.x[! degen]
  n.r <- n.r[! degen]
  x.i <- x.i[! degen,]
  resp <- (! is.na(x.i))

  theta.index <- apply(cbind(r,resp),1,paste,collapse="")
  theta.index <- index.array(theta.index)
  theta.use   <- theta.index[[2]]
  theta.index <- theta.index[[1]]

  rownames(x.i) <- r

  list(x.x = x.x, r = r, n.x = n.x, n.r = n.r, x.i = x.i, resp = resp, n.unique = length(x.x), 
       theta.index = theta.index, theta.use = theta.use,
       degen.x.x = degen.x.x, degen.r = degen.r, degen.n.x = degen.n.x, 
       degen.n.r = degen.n.r, degen.x.i = degen.x.i,
       people = people) 
} # end unique.dat.counts

