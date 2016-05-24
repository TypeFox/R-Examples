index.CN <- function(exchange){


  # check the validity of the argument --------------------------------------
  if(exchange %in% c("sh", "sz") ==FALSE){
    stop("The exchange entered is invalid. The value of 'exchange' can only be 'sh', or 'sz'")
  }

  if(exchange == "sh"){

    # Shanghai Exchange Index -------------------------------------------------

    raw_content_sh <- scan("http://hq.sinajs.cn/list=s_sh000001",
                           what = "raw", encoding = "UTF-8", quiet = TRUE)

    content_sh <- strsplit(raw_content_sh[2], split = ",")[[1]]
    content_sh[1] <- "Shanghai Index"
    content_sh[5] <- as.numeric(content_sh[5])*100  # Jan 2016: to correct the unit. For "sz", no need to do this
    content_sh[6] <- strsplit(content_sh[6], split = "\"")[[1]][1]

    result_sh <-as.list(content_sh)

    names(result_sh) <- c("index", "index.value", "change", "change.percentage",
                          "volume.hand", "amount.10k")
    return(result_sh)
  }

    # Shenzhen Exchange Index -------------------------------------------------

  if(exchange == "sz"){
      raw_content_sz <- scan("http://hq.sinajs.cn/list=s_sz399001",
                             what = "raw", encoding = "UTF-8", quiet = TRUE)

      content_sz <- strsplit(raw_content_sz[2], split = ",")[[1]]
      content_sz[1] <- "Shenzhen Index"
      content_sz[6] <- strsplit(content_sz[6], split = "\"")[[1]][1]

      result_sz <-as.list(content_sz)

      names(result_sz) <- c("index", "index.value", "change", "change.percentage",
                            "volume.hand", "amount.10k")

      return(result_sz)
    }

}
