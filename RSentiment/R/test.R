
#' Calculate the score of sentences
#'
#' This function loads text and calculates score of each sentence on basis 
#' of presence of words of positive and negative sentiment, presence of negation,
#' and checking for sarcasm. 0 indicates neutral sentiment. Positive value indicates
#' positive sentiment. Negative value indicates negative sentiment. 99 indicates
#' sarcasm.
#' @param text A vector of sentences or a sentence (English).
#' @return A vector containing polarity of each sentence.
#' @examples
#'calculate_score("This is good")
#'calculate_score(c("This is good","This is bad"))
 
#'@export 


calculate_score<-function(text){
  text<-as.character(text)
  text <- unlist(lapply(text, function(x) { stringr::str_split(x, "\n") })) 
  
  
#function to calculate number of words in each category within a sentence
getpolarity <- function(sentences, negative_words,positive_words){
  negation<-c("no","not")
  polaritys <- plyr::laply(sentences, function(sentence, negative_words,positive_words){
    
    
    if(is.na(sentence))
      return(-1)
    if(regexpr("[?]", sentence) > 0 )
      return(99)
    #remove unnecessary characters and split up by word
    sentence<- iconv(sentence,"WINDOWS-1252","UTF-8")
    sentence <- gsub('[[:punct:]]', '', sentence)
    sentence <- gsub('[[:cntrl:]]', '', sentence)
    sentence <- gsub('\\d+', '', sentence)
    sentence <- tolower(sentence)
    wordList <- stringr::str_split(sentence, '\\s+')
    words <- unlist(wordList)
    #build vector with matches between sentence and each category
    positive.matches <-match(words, positive_words)
    negative.matches <-match(words, negative_words)
    # get the position of the matched term or NA
    # we just want a TRUE/FALSE
    positive_matches <-!is.na(positive.matches)
    negative_matches <-!is.na(negative.matches)
    # final score
    score <-sum(positive_matches) - sum(negative_matches)
    very.matches<-match(words,c("very","most","more"))
    very_matches <- !is.na(very.matches)
    if(score>=0)
      score<-score+sum(very_matches)
    else
      score<-score-sum(very_matches)
    negation.matches<-match(words,negation)
    negation_matches <- !is.na(negation.matches)
    
    if(score>=0)
    {
      if(sum(negation_matches)>0)
        score<-(-score)
    }
    else
    {
      if(sum(negation_matches)<0)
        score<-(-score)
    }
    
    return(score)
    
    
  }, negative_words, positive_words)
  
  return(polaritys)
}    

#build tables of positive and negative sentences with polaritys
negative_words<- iconv(negative_words,"WINDOWS-1252","UTF-8")
positive_words<- iconv(positive_words,"WINDOWS-1252","UTF-8")
negative_words <- tolower(negative_words)
positive_words <- tolower(positive_words)

res<-getpolarity(text,negative_words,positive_words )
return (res)
}



#' Predicts the sentiment of sentences
#'
#' This function loads text and calculates sentiment of each sentence. It classifies
#' sentences into 6 categories: Positive, Negative, Very Positive, Very Negative 
#' Sarcasm and Neutral.
#'
#' @param text A vector of sentences or a sentence (English).
#' @return A vector containing sentiment of each sentence.

#' @examples
#'calculate_sentiment("This is good")
#'calculate_sentiment(c("This is good","This is bad"))
#'@export
calculate_sentiment<-function(text)
{
  
  res<-calculate_score(text)
  
sentiment<-c()


for(i in 1 : length(res))
{
  if(res[i]==99)
  {
    sentiment[i]<-'Sarcasm'
    
  }
  else
  {
    if(res[i]==0)
    {
      sentiment[i]<-'Neutral'
      
      
    }
    else if (res[i]==-1){
      sentiment[i]<-'Negative'
      
    }
    else if (res[i]==1){
      sentiment[i]<-'Positive'
      
    }
    else if (res[i]>1){
      sentiment[i]<-'Very Positive'
      
      
      
    }
    else{
      sentiment[i]<-'Very Negative'
            
    }
    
    
    
  }
}
results<-data.frame(text,sentiment)
return (results)
}
#' Calculate the number of sentences in each category of sentiment.
#'
#' This function loads text and calculates number of sentences which are positive,
#' negative, very positive, very negative, neutral and sarcasm.
#' 
#' @param text A vector of sentences or a sentence (English).
#' @return A 2-D matrix with two rows and 6 columns where first row contains the name of sentiment
#' category and the second row contains the number in each category in string format.
#' @examples
#'calculate_total_presence_sentiment(c("This is good","This is bad"))
#'@export
calculate_total_presence_sentiment<-function(text){
  
  
  res<-calculate_score(text)
  
  
  score_array<-array(0,dim=c(2,6))
  score_array[1,1]<-'Sarcasm'
  score_array[1,2]<-'Neutral'
  score_array[1,3]<-'Negative'
  score_array[1,4]<-'Positive'
  score_array[1,5]<-'Very Negative'
  score_array[1,6]<-'Very Positive'
  
  for(i in 1 : length(res))
  {
    if(res[i]==99)
    {
     
      score_array[2,1]<-as.numeric(score_array[2,1])+1
    }
    else
    {
      if(res[i]==0)
      {
        
        score_array[2,2]<-as.numeric(score_array[2,2])+1
        
      }
      else if (res[i]==-1){
        
        score_array[2,3]<-as.numeric(score_array[2,3])+1
      }
      else if (res[i]==1){
        
        score_array[2,4]<-as.numeric(score_array[2,4])+1
      }
      else if (res[i]>1){
        
        score_array[2,6]<-as.numeric(score_array[2,6])+1
        
        
      }
      else{
        
        score_array[2,5]<-as.numeric(score_array[2,5])+1
        
      }
      
      
      
    }
  }
  
  return (score_array)
}


