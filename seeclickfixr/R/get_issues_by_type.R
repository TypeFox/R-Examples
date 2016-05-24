get_issues_by_type <- function(city, issue_type, status = "open,acknowledged,closed,archived", limit = 100) {
  total <- 0
  page <- 1
  url <- paste("https://seeclickfix.com/api/v2/issues?place_url=",city,"&search=",issue_type,"&status=",status, "&per_page=",limit,"&page=",page,sep = "")
  url <- gsub(" ","%20",x=url)
  rawdata <- RCurl::getURL(url)
  scf <- jsonlite::fromJSON(txt=rawdata,simplifyDataFrame = T,flatten=F)
  
  issue_id = scf$issues$id
  issue_status = scf$issues$status
  summary = scf$issues$summary
  description = scf$issues$description
  rating = scf$issues$rating
  lat = scf$issues$lat
  lng = scf$issues$lng
  issue_address = scf$issues$address
  created_at = scf$issues$created_at
  acknowledged_at = scf$issues$acknowledged_at
  closed_at = scf$issues$closed_at
  reopened_at = scf$issues$reopened_at
  updated_at = scf$issues$updated_at
  shortened_url = scf$issues$shortened_url
  video_url = scf$issues$media$video_url
  image_full = scf$issues$media$image_full
  image_square_100x100 = scf$issues$media$image_square_100x100
  representative_image_url = scf$issues$media$representative_image_url
  issue_types = scf$issues$point$type
  # scf$issues$point$coordinates # duplicate of lat/lng
  url = scf$issues$url
  html_url = scf$issues$html_url
  comment_url = scf$issues$comment_url
  flag_url = scf$issues$flag_url
  close_url = if(length(scf$issues$transitions$close_url)>0){scf$issues$transitions$close_url} else{NA}
  open_url = if(length(scf$issues$transitions$open_url)>0){scf$issues$transitions$open_url} else{NA}
  reporter_id = scf$issues$reporter$id
  reporter_name = scf$issues$reporter$name
  reporter_wittytitle = scf$issues$reporter$witty_title
  reporter_role = scf$issues$reporter$role
  reporter_civicpoints = scf$issues$reporter$civic_points
  reporter_avatar_full = scf$issues$reporter$avatar$full
  reporter_avatar_square = scf$issues$reporter$avatar$square_100x100
  
  allout <- data.frame(
    issue_id,
    issue_status,
    summary,
    description,
    rating,
    lat,
    lng,
    issue_address,
    created_at,
    acknowledged_at,
    closed_at,
    reopened_at,
    updated_at,
    shortened_url,
    video_url,
    image_full,
    image_square_100x100,
    representative_image_url,
    issue_types,
    url,
    html_url,
    comment_url,
    flag_url,
    close_url,
    open_url,
    reporter_id,
    reporter_name,
    reporter_wittytitle,
    reporter_role,
    reporter_civicpoints,
    reporter_avatar_full,
    reporter_avatar_square 
  )
  
  total <- nrow(allout)
  
  while(limit>total){
    page <- page+1
    if((limit-total)<100){limit <- (limit-total)}
    url <- paste("https://seeclickfix.com/api/v2/issues?place_url=",city,"&search=",issue_type,"&status=",status, "&per_page=",limit,"&page=",page,sep = "")
    url <- gsub(" ","%20",x=url)
    rawdata <- RCurl::getURL(url)
    scf <- jsonlite::fromJSON(txt=rawdata,simplifyDataFrame = T,flatten=F)
    
    issue_id = scf$issues$id
    issue_status = scf$issues$status
    summary = scf$issues$summary
    description = scf$issues$description
    rating = scf$issues$rating
    lat = scf$issues$lat
    lng = scf$issues$lng
    issue_address = scf$issues$address
    created_at = scf$issues$created_at
    acknowledged_at = scf$issues$acknowledged_at
    closed_at = scf$issues$closed_at
    reopened_at = scf$issues$reopened_at
    updated_at = scf$issues$updated_at
    shortened_url = scf$issues$shortened_url
    video_url = scf$issues$media$video_url
    image_full = scf$issues$media$image_full
    image_square_100x100 = scf$issues$media$image_square_100x100
    representative_image_url = scf$issues$media$representative_image_url
    issue_types = scf$issues$point$type
    # scf$issues$point$coordinates # duplicate of lat/lng
    url = scf$issues$url
    html_url = scf$issues$html_url
    comment_url = scf$issues$comment_url
    flag_url = scf$issues$flag_url
    close_url = if(length(scf$issues$transitions$close_url)>0){scf$issues$transitions$close_url} else{NA}
    open_url = if(length(scf$issues$transitions$open_url)>0){scf$issues$transitions$open_url} else{NA}
    reporter_id = scf$issues$reporter$id
    reporter_name = scf$issues$reporter$name
    reporter_wittytitle = scf$issues$reporter$witty_title
    reporter_role = scf$issues$reporter$role
    reporter_civicpoints = scf$issues$reporter$civic_points
    reporter_avatar_full = scf$issues$reporter$avatar$full
    reporter_avatar_square = scf$issues$reporter$avatar$square_100x100
    
    holder <- data.frame(
      issue_id,
      issue_status,
      summary,
      description,
      rating,
      lat,
      lng,
      issue_address,
      created_at,
      acknowledged_at,
      closed_at,
      reopened_at,
      updated_at,
      shortened_url,
      video_url,
      image_full,
      image_square_100x100,
      representative_image_url,
      issue_types,
      url,
      html_url,
      comment_url,
      flag_url,
      close_url,
      open_url,
      reporter_id,
      reporter_name,
      reporter_wittytitle,
      reporter_role,
      reporter_civicpoints,
      reporter_avatar_full,
      reporter_avatar_square 
    )
    allout <- rbind(allout,holder)
    total <- nrow(allout) 
  }
  
  return(allout)
}
