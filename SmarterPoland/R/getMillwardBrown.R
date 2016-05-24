get_millwardbrown <- function(url) {
  doc <- rjson::fromJSON(file=url, method = "C")
  return(doc)
}

get_millwardbrown_detail <- function(x) {
  detailed_url <- paste0('http://wybory.millwardbrown.com/partie-polityczne-parlament-krajowy/',x,'.json')
  doc <- rjson::fromJSON(file=detailed_url, method = "C")
#  download.file(detailed_url,destfile = basename(detailed_url))
#  doc <- jsonlite::fromJSON(basename(detailed_url)) #jsonlite::
  return(doc)
}

# getMillwardBrown <- function() {
#   main_url <- 'http://wybory.millwardbrown.com/partie-polityczne-parlament-krajowy.json'
#   
#   ## dane z wykresu
#   dane <- get_millwardbrown(main_url)
#   
#   ## dodatkowe informacje
#   ids <- sapply(dane$polls, function(x) x$id)
#   pojedyncze <- lapply(ids,get_millwardbrown_detail)
#   
#   # to the data frame
#   namess <- sapply(dane$data, function(x) x$name)
#   inds <- which(sapply(namess, length) > 0)
#   partie <- unlist(namess[inds])
#   lista <- lapply(inds, function(ind) {
# #    tmp <- dane$data[[id]]$data
#     tmp <- as.data.frame(t(as.data.frame(dane$data[[ind]]$data)))
#     rownames(tmp) <- NULL
#     colnames(tmp) <- c("data", "poparcie")
#     tmp$partia <- namess[ind]
#     tmp
#   })
#   do.call(rbind, lista)
# }

getMillwardBrown <- function() {
  main_url <- 'http://wybory.millwardbrown.com/partie-polityczne-parlament-krajowy.json'
  
  main_data <- get_millwardbrown(main_url)
  
  ids <- sapply(main_data$polls, function(x) x$id)
  poll_details <- lapply(ids, get_millwardbrown_detail)
  
  poll_desc <- do.call(rbind,
            lapply(poll_details,function(x) data.frame(ID = x$id,
                                    collected_on = as.character(x$collectedOn),
                                    sample_size = x$sampleSize,
                                    turnout = x$turnout)))
  
  poll_desc_2 <- do.call(rbind,
            lapply(poll_details,function(x) {
            data.frame(
              do.call(rbind,
                    lapply(x$results, function(z) {
                    data.frame(value = z$y, name = z$name)
                  })),
              ID = x$id) 
          }))
  
  poll_full_details <- merge(poll_desc, poll_desc_2, by="ID", all.y = TRUE)
  
  return(poll_full_details)
}

