get_socialmedia <-
function(links, sleep.time=0) {
    for(i in 1:length(links)) {
        
        # extract first batch
        link.now <- links[i]
        
        # sleep time
        Sys.sleep(sleep.time)
        
        # percentage of work done
        if(length(links)>10) {
            is.even <- function(x) x %% 2 == 0
            if(is.even(i)){ print(paste(round(100 * i/length(links), digits=2), "%")) }
        }
        
        # download cert if not available
        if(.Platform$OS.type == "windows") { if(!file.exists("cacert.perm")) download.file(url="https://curl.haxx.se/ca/cacert.pem", destfile="cacert.perm") }
        
        # function for looping thru URLs
        if(.Platform$OS.type == "windows") { api_scrapper <- function(x) try(RCurl::getURL(x, cainfo = "cacert.perm", timeout = 240, ssl.verifypeer = FALSE)) } else { 
            api_scrapper <- function(x) try(RCurl::getURL(x, timeout = 240, ssl.verifypeer = FALSE)) }        
        
        # create call URLs
        fbk.call <- paste0("https://api.facebook.com/method/links.getStats?urls=",link.now,"&format=json")
        rdd.call <- paste0("http://buttons.reddit.com/button_info.json?url=",link.now)
        lkn.call <- paste0("https://www.linkedin.com/countserv/count/share?url=",link.now,"&format=json")
        stu.call <- paste0("http://www.stumbleupon.com/services/1.01/badge.getinfo?url=",link.now)
        pin.call <- paste0("http://api.pinterest.com/v1/urls/count.json?callback=receiveCount&url=",link.now,"&format=json")
        
        # prepare for response
        fbk.response <- data.frame()
        rdd.response <- data.frame()
        lkn.response <- data.frame()
        stu.response <- data.frame()
        pin.response <- data.frame()
        
        # collect responses
        fbk.response <- try(data.frame(jsonlite::fromJSON(api_scrapper(fbk.call))))
        rdd.response <- try(jsonlite::fromJSON(api_scrapper(rdd.call))$data$children$data)
        lkn.response <- try(data.frame(jsonlite::fromJSON(api_scrapper(lkn.call))))
        stu.response <- try(data.frame(jsonlite::fromJSON(api_scrapper(stu.call))))
        pin.response <- try(data.frame(jsonlite::fromJSON(gsub("receiveCount\\(|\\)", "", api_scrapper(pin.call)))))
        
        # prepare data frame for aggregation
        response <- data.frame(
            url=NA,
            normalized_url=NA,
            fbk_shares=NA,
            fbk_likes=NA,
            fbk_comments=NA,
            fbk_total=NA,
            fbk_clicks=NA,
            rdt_score=NA,
            rdt_downs=NA,
            rdt_ups=NA,
            rdt_comments=NA,
            lkn_shares=NA,
            stu_views=NA,
            pin_counts=NA
            )
        
        # prepare data frame for output
        if(i==1){
            response.temp <- data.frame()
            start.time <- Sys.time()
        }
                
        # aggregate results
        if(is.data.frame(fbk.response)) { if(length(fbk.response$url)>0) { try(response$url<-as.character(fbk.response$url)) } }
        if(is.data.frame(fbk.response)) { if(length(fbk.response$normalized_url)>0) { try(response$normalized_url<-as.character(fbk.response$normalized_url)) } }
        if(is.data.frame(fbk.response)) { if(length(fbk.response$share_count)>0) { try(response$fbk_shares<-as.numeric(as.character(fbk.response$share_count))) } }
        if(is.data.frame(fbk.response)) { if(length(fbk.response$like_count)>0) { try(response$fbk_likes<-as.numeric(as.character(fbk.response$like_count))) } }
        if(is.data.frame(fbk.response)) { if(length(fbk.response$comment_count)>0) { try(response$fbk_comments<-as.numeric(as.character(fbk.response$comment_count))) } }
        if(is.data.frame(fbk.response)) { if(length(fbk.response$total_count)>0) { try(response$fbk_total<-as.numeric(as.character(fbk.response$total_count))) } }
        if(is.data.frame(fbk.response)) { if(length(fbk.response$click_count)>0) { try(response$fbk_clicks<-as.numeric(as.character(fbk.response$click_count))) } }
        if(is.data.frame(rdd.response)) { if(length(rdd.response$score)>0) { try(response$rdt_score<-rdd.response$score) } }
        if(is.data.frame(rdd.response)) { if(length(rdd.response$downs)>0) { try(response$rdt_downs<-rdd.response$downs) } }
        if(is.data.frame(rdd.response)) { if(length(rdd.response$ups)>0) { try(response$rdt_ups<-rdd.response$ups) } }
        if(is.data.frame(rdd.response)) { if(length(rdd.response$num_comments)>0) { try(response$rdt_comments<-rdd.response$num_comments) } }
        if(is.data.frame(lkn.response)) { if(length(lkn.response$count)>0) { try(response$lkn_shares<-as.numeric(as.character(lkn.response$count))) } }
        if(is.data.frame(stu.response)) { if(length(stu.response$result.views)>0) { try(response$stu_views<-as.numeric(as.character(stu.response$result.views))) } }
        if(is.data.frame(pin.response)) { if(length(pin.response$count)>0) { try(response$pin_counts<-as.numeric(as.character(pin.response$count))) } }
            
        # collate results
        response.temp <- rbind(response.temp, response)
            }
    print(paste0("Query execution time: ", round(as.numeric(difftime(Sys.time(), start.time, units = "mins")), digits=2), " minutes"))
    return(response.temp)
}
