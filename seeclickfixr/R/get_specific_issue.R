get_specific_issue <- function(issue_id) {
  
  url <- paste("http://seeclickfix.com/api/v2/issues/", issue_id, sep = "")
  url <- gsub(" ","%20",x=url)
  rawdata <- RCurl::getURL(url)
  scf <- jsonlite::fromJSON(rawdata)
  
  allout <- data.frame(issue_id = scf$issues$id,
                       status = scf$issues$status,
                       summary = scf$issues$summary,
                       description = scf$issues$description,
                       rating = scf$issues$rating,
                       lat = scf$issues$lat,
                       lng = scf$issues$lng,
                       issue_address = scf$issues$address,
                       created_at = scf$issues$created_at,
                       acknowledged_at = scf$issues$acknowledged_at,
                       closed_at = scf$issues$closed_at,
                       reopened_at = scf$issues$reopened_at,
                       updated_at = scf$issues$updated_at,
                       shortened_url = scf$issues$shortened_url,
                       video_url = scf$issues$media$video_url,
                       image_full = scf$issues$media$image_full,
                       image_square_100x100 = scf$issues$media$image_square_100x100,
                       representative_image_url = scf$issues$media$representative_image_url,
                       issue_type = scf$issues$point$type,
                       # scf$issues$point$coordinates, # duplicate of lat/lng
                       url = scf$issues$url,
                       html_url = scf$issues$html_url,
                       comment_url = scf$issues$comment_url,
                       flag_url = scf$issues$flag_url,
                       close_url = scf$issues$transitions$close_url,
                       open_url = scf$issues$transitions$open_url,
                       reporter_id = scf$issues$reporter$id,
                       reporter_name = scf$issues$reporter$name,
                       reporter_wittytitle = scf$issues$reporter$witty_title,
                       reporter_role = scf$issues$reporter$role,
                       reporter_civicpoints = scf$issues$reporter$civic_points,
                       reporter_avatar_full = scf$issues$reporter$avatar$full,
                       reporter_avatar_square = scf$issues$reporter$avatar$square_100x100
                       )
  
  return(allout)
}
