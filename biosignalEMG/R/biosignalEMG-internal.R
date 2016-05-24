.validchannel <- function(data, channel) {
    if (length(channel) == 1) {
        if (is.numeric(channel)) {
            if (channel < 1 | channel > dim(data$values)[2]) {
                valid <- (-1)
                ch <- numeric()
                m <- paste("'channel' must be between 0 and", dim(data$values)[2])
            } else {
                valid <- 1
                ch <- channel
                m <- ""
            }
        } else {
            if (is.character(channel)) {
                ch <- grep(channel, data$data.name)
                if (length(ch) > 0) {
                  if (length(ch) == 1) {
                    valid <- 1
                    m <- ""
                  } else {
                    valid <- (-5)
                    ch <- numeric()
                    m <- "Multiple possible channels"
                  }
                } else {
                  valid <- (-4)
                  ch <- numeric()
                  m <- paste(channel, "is not a valid parameter 'channel'")
                }
            } else {
                valid <- (-2)
                ch <- numeric()
                m <- paste(channel, "is not a valid parameter 'channel'")
            }
        }
    } else {
        valid <- (-3)
        ch <- numeric()
        m <- "Choose a single channel"
    }
    return(list(val = valid, message = m, channel = ch))
} 
