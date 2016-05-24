get_site_ID <- function(lake_name){
  siteID <- get_metadata(lake_name, 'siteID')
  return(siteID)
}

temp_types = function(){
  return(c('Lake.Temp.Summer.Satellite','Lake.Temp.Summer.InSitu'))
}