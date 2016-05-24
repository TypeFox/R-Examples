#' @import httr
#' @import RCurl
#' @import jsonlite

library(httr)
library(RCurl)
library(jsonlite)

#' Gets a list of translation directions supported by the service
#'
#' @param api_key yandex API key
#' @param lang If set, the response contains explanations of language codes. Language names are output in the language corresponding to the code in this parameter.
#' @return data frame giving supported translation direction
#' @export
#' @examples
#' data=get_translation_direction(api_key)

get_translation_direction=function(api_key,lang="")
{
  url="https://translate.yandex.net/api/v1.5/tr.json/getLangs?"
  
  url=paste(url,"key=",api_key,sep="")
  
  if(lang != "")
  {
    url=paste(url,"&ui=",lang,sep="")
  }
  
  d=getURL(url,ssl.verifyhost = 0L, ssl.verifypeer = 0L)
  d=fromJSON(d)
  d
}


#' Detects the language of the specified text.
#'
#' @param api_key yandex API key
#' @param text The text to detect the language for
#' @return data frame giving detected language
#' @export
#' @examples
#' data=detect_language(api_key,text="how are you?")


detect_language=function(api_key,text="")
{
  url="https://translate.yandex.net/api/v1.5/tr.json/detect?"
  
  url=paste(url,"key=",api_key,sep="")

  if(text != "")
  {
    url=paste(url,"&text=",text,sep="")
  }
  
  d=getURL(url,ssl.verifyhost = 0L, ssl.verifypeer = 0L)
  d=fromJSON(d)
  d$lang
}  

#' Translates text to the specified language
#'
#' @param api_key yandex API key
#' @param text The text to translate.The maximum size of the text being passed is 10000 characters.
#' @param lang The translation direction.You can use any of the following ways to set it:As a pair of language codes separated by a hyphen ("from"-"to"). For example, en-ru indicates translating from English to Russian.As the final language code (for example, ru). In this case, the service tries to detect the source language automatically.
#' @return data frame giving translated text
#' @examples
#' data=translate(api_key,text="how are you?",lang="hi")


translate=function(api_key,text="",lang="")
{
  url="https://translate.yandex.net/api/v1.5/tr.json/translate?"
  
  url=paste(url,"key=",api_key,sep="")
  
  if(text != "")
  {
    url=paste(url,"&text=",text,sep="")
  }
  
  if(lang != "")
  {
    url=paste(url,"&lang=",lang,sep="")
  }
  
  d=getURL(url,ssl.verifyhost = 0L, ssl.verifypeer = 0L)
  d=fromJSON(d)
  d$code=NULL
  d
}
