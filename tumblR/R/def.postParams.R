def.postParams <-
function(type="text",state="published",tags=NA, tweet=NA, 
                        date=date,format="html", slug=NA,
                        title_text=NA,body=NA,
                        caption_photo=NA,link=NA,source_photo=NA,data_photo=NA,
                        quote=NA,source_quote=NA,
                        url_link=NA,title_link=NA,description=NA,
                        title_chat=NA, conversation=NA,
                        external_url=NA,data_audio=NA,caption_audio=NA,
                        embed=NA,data_video=NA,caption_video=NA,
                        reblog=FALSE,editing=FALSE,
                        id=NA,reblog_key=NA,comment=NA){
  
  type_type<-c("text", "photo", "quote", "link", "chat", "audio", "video")
  
  if(!(type %in% type_type))
    stop("type must be a string type. Available values for type are: text, photo, quote, link, chat, audio, video") 
  
  state_type<-c("published", "draft", "queue", "private")
  
  if(!(state %in% state_type))
    stop("state must be a string type. Available values for type are: published, draft, queue, private")

  if(!is.character(tags))
    stop("tags must be a string type")
    
  
  if(!is.na(tweet)){
    if(!is.character(tweet))
      stop("tweet must be a string type")
  }
  
  format_type<-c("html", "markdown")
  
  if(!(format %in% format_type))
    stop("format must be a string type. Available values for type are: html, markdown")
  
  if(!is.na(slug)){
    if(!is.character(slug))
      stop("slug must be a string type")
  }
  
  if(editing==TRUE){
    if(is.na(id)){
      stop("For editing a post, id is a requested field")
      } else{
        if(!is.numeric(id))
          stop("id must be a numeric type")
      }
  }
  
  if(reblog==TRUE){
    if(is.na(id) || is.na(reblog_key)){
      stop("id and reblog_key are requested parameter for reblogging")
    } else{
      if(!is.numeric(id) && !is.numeric(reblog_key))
        stop("id and reblog_key must be numeric type")
    }
    
    if(!is.na(comment)){
      if(!is.character(comment))
        stop("comment must be a string type")
    }  
  }
  
  if(type=="text"){
    
    if(is.na(body)){
      stop("body is a requested parameter")
    } else{
      if(!is.character(body))
        stop("body must be a string type")
    }
    
    if(!is.na(title_text)){
      if(!is.character(title_text))
        stop("title_text must be a string type")
    }    
    
    Text_Params<-list(type=type,state=state,tags=tags, tweet=tweet, 
                            date=date,format=format, slug=slug,
                            title=title_text,body=body,
                            id=id,reblog_key=reblog_key,comment=comment)
    
    len<-length(Text_Params)
    s<-NULL
    for(i in 1:len){
      if(!is.na(Text_Params[[i]][1]))
        s<-c(s,i)
    }
    
    bodyParams<-Text_Params[s]
    
  }
  
  if(type=="photo"){
    
    if(is.na(source_photo) && is.na(data_photo)){
      stop("Either source_photo or data_photo is a requested parameter")
    } else{
      if(!is.character(source_photo) && !is.character(data_photo))
        stop("source_photo or data_photo must be a string type")
    }
    
    if(!is.na(caption_photo)){
      if(!is.character(caption_photo))
        stop("caption_photo must be a string type")
    }    
    
    if(!is.na(link)){
      if(!is.character(link))
        stop("link must be a string type")
    }    
    
    Photo_Params<-list(type=type,state=state,tags=tags, tweet=tweet, 
                             date=date,format=format, slug=slug,
                             caption=caption_photo,link=link,source=source_photo,data=data_photo,
                             id=id,reblog_key=reblog_key,comment=comment,
                             stringsAsFactors=FALSE)
    
    len<-length(Photo_Params)
    s<-NULL
    for(i in 1:len){
      if(!is.na(Photo_Params[[i]][1]))
        s<-c(s,i)
    }
    
    bodyParams<-Photo_Params[s]
    
  }
  
  if(type=="quote"){
    
    if(is.na(quote)){
      stop("quote is a requested parameter")
    } else{
      if(!is.character(quote))
        stop("quote must be a string type")
    }
    
    if(!is.na(source_quote)){
      if(!is.character(source_quote))
        stop("source_quote must be a string type")
    }    
    
    Quote_Params<-list(type=type,state=state,tags=tags, tweet=tweet, 
                             date=date,format=format, slug=slug,
                             quote=quote,source=source_quote,
                             id=id,reblog_key=reblog_key,comment=comment,
                             stringsAsFactors=FALSE)
    
    len<-length(Quote_Params)
    s<-NULL
    for(i in 1:len){
      if(!is.na(Quote_Params[[i]][1]))
        s<-c(s,i)
    }
    
    bodyParams<-Quote_Params[s]
    
  }
  
  if(type=="link"){
    
    if(is.na(url_link)){
      stop("url_link is a requested parameter")
    } else{
      if(!is.character(url_link))
        stop("url_link must be a string type")
    }
    
    if(!is.na(title_link)){
      if(!is.character(title_link))
        stop("title_link must be a string type")
    }  
    
    if(!is.na(description)){
      if(!is.character(description))
        stop("description must be a string type")
    }
    
    Link_Params<-list(type=type,state=state,tags=tags, tweet=tweet, 
                            date=date,format=format, slug=slug,
                            url=url_link,title=title_link,description=description,
                            id=id,reblog_key=reblog_key,comment=comment,
                            stringsAsFactors=FALSE)
    
    len<-length(Link_Params)
    s<-NULL
    for(i in 1:len){
      if(!is.na(Link_Params[[i]][1]))
        s<-c(s,i)
    }
    
    bodyParams<-Link_Params[s]
    
  }
  
  if(type=="chat"){
    
    if(is.na(conversation)){
      stop("conversation is a requested parameter")
    } else{
      if(!is.character(conversation))
        stop("conversation must be a string type")
    }
    
    if(!is.na(title_chat)){
      if(!is.character(title_chat))
        stop("title_chat must be a string type")
    }    
    
    Chat_Params<-list(type=type,state=state,tags=tags, tweet=tweet, 
                            date=date,format=format, slug=slug,
                            title_chat=title_chat, conversation=conversation,
                            id=id,reblog_key=reblog_key,comment=comment,
                            stringsAsFactors=FALSE)
    
    len<-length(Chat_Params)
    s<-NULL
    for(i in 1:len){
      if(!is.na(Chat_Params[[i]][1]))
        s<-c(s,i)
    }
    
    bodyParams<-Chat_Params[s]
    
  }
  
  if(type=="audio"){
    
    if(is.na(external_url) && is.na(data_audio)){
      stop("Either external_url or data_audio is a requested parameter")
    } else{
      if(!is.character(external_url) && !is.character(data_audio))
        stop("external_url or data_audio must be a string type")
    }
    
    if(!is.na(caption_audio)){
      if(!is.character(caption_audio))
        stop("caption_audio must be a string type")
    }    
    
    Audio_Params<-list(type=type,state=state,tags=tags, tweet=tweet, 
                             date=date,format=format, slug=slug,
                             external_url=external_url,data=data_audio,caption=caption_audio,
                             id=id,reblog_key=reblog_key,comment=comment,
                             stringsAsFactors=FALSE)
    
    len<-length(Audio_Params)
    s<-NULL
    for(i in 1:len){
      if(!is.na(Audio_Params[[i]][1]))
        s<-c(s,i)
    }
    
    bodyParams<-Audio_Params[s]
    
  }
  
  if(type=="video"){
    
    if(is.na(embed) && is.na(data_video)){
      stop("Either embed or data_video is a requested parameter")
    } else{
      if(!is.character(embed) && !is.character(data_video))
        stop("embed or data_video must be a string type")
    }
    
    if(!is.na(caption_video)){
      if(!is.character(caption_video))
        stop("caption_audio must be a string type")
    }  
    
    Video_Params<-list(type=type,state=state,tags=tags, tweet=tweet, 
                             date=date,format=format, slug=slug,
                             embed=embed,data=data_video,caption=caption_video,
                             id=id,reblog_key=reblog_key,comment=comment,
                             stringsAsFactors=FALSE)
    len<-length(Video_Params)
    s<-NULL
    for(i in 1:len){
      if(!is.na(Video_Params[[i]][1]))
        s<-c(s,i)
    }
    
    bodyParams<-Video_Params[s]
    
  }
  
  return(bodyParams)
  
}
