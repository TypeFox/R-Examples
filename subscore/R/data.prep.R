#' This function prepares the data that can be used by other functions 
#' @description This function generates a list of datasets based on the scored original dataset, 
#' which can be used as objects in subscore computing functions.
#' @param scored.data A dataset that contains scored responses only.
#' @param subtest.infor A numerical vector. The first number indicates the number of subtests,
#'                       followed by the number o items of each subtest.
#' @return A list that contains datasets of all subtests and the whole test. 
#' @examples 
#'         subtest.infor<-c(3,15,15,20) 
#'         # This test consists of 3 subtests, which have 15, 15 and 20 items respectively.
#'         test.datalist<-data.prep(scored.data,subtest.infor)
#' @export

data.prep<-function (scored.data, subtest.infor) {
  n.subtests<-subtest.infor[1]
  test.names<-rep(NA,n.subtests+1)
  test.names[n.subtests+1]<-paste('total.test')

    for (t in 1:n.subtests) { 
    test.names[t]<-paste ('subtest.',t,sep='')    
  } 
   
  test.data <- as.list(rep(NA, length(test.names)))
  names(test.data) <- test.names
  
  test.data[[1]]<-scored.data[,1:subtest.infor[2]]
  test.data[[n.subtests+1]]<-scored.data
  
  for (t in 2:n.subtests) {  
    test.data[[t]] <-scored.data[,(sum(subtest.infor[2:t],na.rm = TRUE)+1):sum(subtest.infor[2:(t+1)],na.rm = TRUE)]
    } 
  return(test.data)
} 