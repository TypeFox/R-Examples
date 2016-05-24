#
# reduce categories (ordinal to simple)
#
red.cat<-function(x,nrespcat){
  m<-median(range(x, na.rm=TRUE))
  if (nrespcat%%2 == 0) {     # 2 categories
     x <- ifelse(x < m,0,1)
  } else {                    # 3 categories
     x <- ifelse(x < m,0,ifelse(x > m,2,1))
  }
  x
}
