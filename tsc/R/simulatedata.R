simulatedata <-
function(x,y, num.mc){ 
  listvls <- c()
  for(i in 1:num.mc){
    x_null <- rnorm(length(x),0,1)
    y_null <- rnorm(length(y),0,1)
    test_stat<-.C("CWrapper1", n1=as.integer(length(x_null)),n2=as.integer(length(y_null)),y1=as.double(x_null),y2=as.double(y_null),test_stat=as.double(1))$test_stat
    listvls <- c(listvls,test_stat)

  }
 return(listvls)
}
