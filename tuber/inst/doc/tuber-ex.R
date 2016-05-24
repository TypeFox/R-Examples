## ---- eval=FALSE, install------------------------------------------------
#  install.packages("tuber")

## ---- eval=FALSE, install_g----------------------------------------------
#  # install.packages('devtools')
#  devtools::install_github("soodoku/tuber", build_vignettes = TRUE)

## ---- eval=FALSE, loadlib------------------------------------------------
#  library(tuber)

## ---- eval=FALSE, auth---------------------------------------------------
#  yt_oauth("998136489867-5t3tq1g7hbovoj46dreqd6k5kd35ctjn.apps.googleusercontent.com", "MbOSt6cQhhFkwETXKur-L9rN")

## ---- eval=FALSE, prints-------------------------------------------------
#  ## Waiting for authentication in browser...
#  ## Press Esc/Ctrl + C to abort
#  ## Authentication complete.

## ---- eval=FALSE, getstats-----------------------------------------------
#  get_stats(video_id="N708P-A45D0")

## ---- eval=FALSE, getdetails---------------------------------------------
#  get_details(video_id="N708P-A45D0")

## ---- eval=FALSE, getcaptions--------------------------------------------
#  get_captions(video_id="yJXTXN4xrI8")

## ---- eval=FALSE, searchvids---------------------------------------------
#  yt_search("Barack Obama")

## ---- eval=FALSE, get_comments-------------------------------------------
#  res <- get_comments(video_id="N708P-A45D0")
#  # First comment
#  res$items[[1]]$snippet$topLevelComment$snippet$textDisplay

