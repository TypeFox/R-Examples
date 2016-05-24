get_reddit <-
function(links, sleep.time=0) {
    rdd.call <- paste0("http://buttons.reddit.com/button_info.json?url=",links)
    if(.Platform$OS.type == "windows") { if(!file.exists("cacert.perm")) download.file(url="https://curl.haxx.se/ca/cacert.pem", destfile="cacert.perm") }
    if(.Platform$OS.type == "windows") { api_scrapper <- function(x) try(RCurl::getURL(x, cainfo = "cacert.perm", timeout = 240, ssl.verifypeer = FALSE)) } else { 
        api_scrapper <- function(x) try(RCurl::getURL(x, timeout = 240, ssl.verifypeer = FALSE)) }
    Sys.sleep(sleep.time)
    temp <- try(jsonlite::fromJSON(api_scrapper(rdd.call)))
    rdd.response <- data.frame(
        domain = temp$data$children$data$domain,
        banned_by = temp$data$children$data$banned_by,
        media_embed = temp$data$children$data$media_embed,
        subreddit = temp$data$children$data$subreddit,
        selftext_html = temp$data$children$data$selftext_html,
        selftext = temp$data$children$data$selftext,
        likes = temp$data$children$data$likes,
        secure_media = temp$data$children$data$secure_media,
        link_flair_text = temp$data$children$data$link_flair_text,
        id = temp$data$children$data$id,
        gilded = temp$data$children$data$gilded,
        secure_media_embed = temp$data$children$data$secure_media_embed,
        clicked = temp$data$children$data$clicked,
        stickied = temp$data$children$data$stickied,
        author = temp$data$children$data$author,
        media = temp$data$children$data$media,
        score = temp$data$children$data$score,
        approved_by = temp$data$children$data$approved_by,
        over_18 = temp$data$children$data$over_18,
        hidden = temp$data$children$data$hidden,
        thumbnail = temp$data$children$data$thumbnail,
        subreddit_id = temp$data$children$data$subreddit_id,
        edited = temp$data$children$data$edited,
        link_flair_css_class = temp$data$children$data$link_flair_css_class,
        author_flair_css_class = temp$data$children$data$author_flair_css_class,
        downs = temp$data$children$data$downs,
        saved = temp$data$children$data$saved,
        is_self = temp$data$children$data$is_self,
        permalink = temp$data$children$data$permalink,
        name = temp$data$children$data$name,
        created = temp$data$children$data$created,
        url = temp$data$children$data$url,
        author_flair_text = temp$data$children$data$author_flair_text,
        title = temp$data$children$data$title,
        created_utc = temp$data$children$data$created_utc,
        distinguished = temp$data$children$data$distinguished,
        num_comments = temp$data$children$data$num_comments,
        visited = temp$data$children$data$visited,
        num_reports = temp$data$children$data$num_reports,
        ups = temp$data$children$data$ups)
    return(rdd.response)
}
