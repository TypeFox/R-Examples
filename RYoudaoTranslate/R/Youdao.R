youdaoUrl = function(word,api,keyfrom)
{
  if(any(is.null(c(word,api,keyfrom))))
    stop("word, api and keyfrom are all necessary.")
  word = gsub(pattern=" ",replacement="+",x=word)
  paste("http://fanyi.youdao.com/openapi.do?keyfrom=",keyfrom,"&key=",api,"&type=data&doctype=json&version=1.1&q=",
        word,sep="")
}

youdaoTranslate = function(word,api,keyfrom)
{
  url = getURL(youdaoUrl(word,api,keyfrom))
  obj = fromJSON(url)  
  return(obj$web)
}

youdaoDisplay = function(youdaoObj,word)
{
  data = NULL
  for(i in youdaoObj)
  {
    data = c(data,paste(i$key,paste(i$value,collapse=""),sep=","))
  }
  data = paste(data,collapse="; ")
  names(data) = word
  return(data)
}

youdaoLookUp = function(word,api,keyfrom)
{
  youdaoDisplay(youdaoObj=youdaoTranslate(word,
                                          api,
                                          keyfrom),word) 
}
