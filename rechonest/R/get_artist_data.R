
#' To get artist's data
#'
#' @param api_key Echo Nest API key
#' @param name artist's name
#' @param id artist's id
#' @param hotttnesss artist's hotttnesss
#' @param terms artist's terms
#' @param blogs blogs about artist
#' @param news news articles about artist
#' @param reviews reviews about artist
#' @param familiarity artist's familiarity
#' @param audio artist's audio details
#' @param images artist's images details
#' @param songs artist's songs details
#' @param discovery artist's discovery details
#' @param partner partner catalog
#' @param biographies artist's biographies
#' @param doc_counts artist's doc_counts
#' @param artist_location artist location
#' @param years_active years active
#' @param urls urls of artist websites
#' @return data frame giving artist's hotttnesss
#' @export
#' @examples
#' \dontrun{
#' data=get_artist_data(api_key,name="coldplay",terms=T,blogs=T)
#' }

get_artist_data=function(api_key,name=NA,id=NA,hotttnesss=T,
                         terms=F,blogs=F,news=F,familiarity=F,
                         audio=F,images=F,songs=F,reviews=F,
                         discovery=F,partner=NA,biographies=F,
                         doc_counts=F,artist_location=F,years_active=F,urls=F)
{
  url=paste("http://developer.echonest.com/api/v4/artist/profile?api_key=",api_key,"&format=json",sep="")
  final=""
  
  if(!is.na(name))
  {
    name=gsub(" ","+",name)
    url=paste(url,"&name=",name,sep="")
  }
  if(!is.na(id))
  {
    url=paste(url,"&id=",id,sep="")
  }
  
  if(hotttnesss)
  {
    url=paste(url,"&bucket=hotttnesss&bucket=hotttnesss_rank",sep="")
  }
  
  if(discovery)
  {
    url=paste(url,"&bucket=discovery&bucket=discovery_rank",sep="")
  }
  
  if(familiarity)
  {
    url=paste(url,"&bucket=familiarity&bucket=familiarity_rank",sep="")
  }
  
  if(blogs)
  {
    url=paste(url,"&bucket=blogs",sep="")
  }
  
  if(doc_counts)
  {
    url=paste(url,"&bucket=doc_counts",sep="")
  }
  
  if(years_active)
  {
    url=paste(url,"&bucket=years_active",sep="")
  }
  
  if(artist_location)
  {
    url=paste(url,"&bucket=artist_location",sep="")
  }
  
  if(urls)
  {
    url=paste(url,"&bucket=urls",sep="")
  }
  
  
  
  if(biographies)
  {
    url=paste(url,"&bucket=biographies",sep="")
  }
  
  if(news)
  {
    url=paste(url,"&bucket=news",sep="")
  }
  
  if(reviews)
  {
    url=paste(url,"&bucket=reviews",sep="")
  }
  
  if(songs)
  {
    url=paste(url,"&bucket=songs",sep="")
  }
  
  if(images)
  {
    url=paste(url,"&bucket=images",sep="")
  }
  
  if(audio)
  {
    url=paste(url,"&bucket=audio",sep="")
  }
  
  if(terms)
  {
    url=paste(url,"&bucket=terms",sep="")
  }
  
  if(!is.na(partner))
  {
    partner=gsub(" ","+",partner)
    url=paste(url,"&bucket=id:",partner,sep="")
  }
  
  rd=getURL(url)
  rd=fromJSON(rd)
  
  data=rd$response$artist
  
}
