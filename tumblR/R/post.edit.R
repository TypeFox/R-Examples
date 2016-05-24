post.edit <-
function(base_hostname=NA, type="text",state="published",tags=NA, tweet=NA, 
                    date=as.character(Sys.time()),format="html", slug=NA,
                    title_text=NA,body=NA,
                    caption_photo=NA,link=NA,source_photo=NA,data_photo=NA,
                    quote=NA,source_quote=NA,
                    url_link=NA,title_link=NA,description=NA,
                    title_chat=NA, conversation=NA,
                    external_url=NA,data_audio=NA,caption_audio=NA,
                    embed=NA,data_video=NA,caption_video=NA,
                    id=NA,
                    token=NA,consumer_key=NA,consumer_secret=NA){
  
  if(!is.character(base_hostname))
    stop("base_hostname must be a string")
  
  if(class(token)[1]!="Token1.0")
    stop("token must be a Token1.0 type")
  
  if(!is.character(consumer_key))
    stop("consumer_key must be a string")
  
  if(!is.character(consumer_secret))
    stop("consumer_secret must be a string")
  
  bodyParams<-def.postParams(type=type,state=state,tags=tags, tweet=tweet, 
                             date=date,format=format, slug=slug,
                             title_text=title_text,body=body,
                            caption_photo=caption_photo,link=link,source_photo=source_photo,data_photo=data_photo,
                             quote=quote,source_quote=source_quote,
                             url_link=url_link,title_link=title_link,description=description,
                             title_chat=title_chat, conversation=conversation,
                             external_url=external_url,data_audio=data_audio,caption_audio=caption_audio,
                             embed=embed,data_video=data_video,caption_video=caption_video,
                             reblog=FALSE,editing=TRUE,
                             id=id,reblog_key=NA,comment=NA)
  
  url<-paste("http://api.tumblr.com/v2/blog/",base_hostname,"/post/edit",sep="")
  connection="POST"
  
  res<-fromJSON(http.connection(url,token,bodyParams,consumer_key,consumer_secret,connection))
  
  return(res)
}
