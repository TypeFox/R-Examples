#library(testthat)
library(RSentiment)

expect_list_equal<-function(t,s)
{
  if(length(setdiff(t,s))>0)
    return(1)
  else
    return(0)
}
#test_check("RSentiment")
test_that("calculate_score is score of sentences", {
  
  
  expect_equal(expect_list_equal(calculate_score(c("This is very good","This is bad")), c(2,-1)),0)
  
})
test_that("calculate_sentiment is sentiment of sentences", {
  
  df<-data.frame(text=c("This is very good","This is bad"),sentiment=c("Very Positive","Negative"))
  expect_equal(expect_list_equal(calculate_sentiment(c("This is very good","This is bad")), df),0)
  
})
test_that("calculate_sentiment is sentiment of sentences", {
  
  score_array<-array(0,dim=c(2,6))
 score_array[1,1]<-'Sarcasm'
   score_array[1,2]<-'Neutral'
  score_array[1,3]<-'Negative'
  score_array[1,4]<-'Positive'
  score_array[1,5]<-'Very Negative'
  score_array[1,6]<-'Very Positive'
  score_array[2,1]<-0
  score_array[2,2]<-0
  score_array[2,3]<-1
  score_array[2,4]<-0
  score_array[2,5]<-0
  score_array[2,6]<-1
  expect_equal(expect_list_equal(calculate_total_presence_sentiment(c("This is very good","This is bad")), score_array),0)
  
})
