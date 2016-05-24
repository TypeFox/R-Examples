nickel.expand <- local({
    source("nickel.R")
    source("ewrates.R")
    trim <- function(start)
    {
        end   <- start + 5
        entry <- pmax(nickel$agein, start)
        exit  <- pmin(nickel$ageout, end)
        valid <- (entry < exit)
        cens  <- (nickel$ageout[valid] > end)
        result <- nickel[valid,]
        within(result, {
            icd[cens] <- 0
            agein <- entry[valid]
            ageout <- exit[valid]
            agr <- start
            ygr <- trunc((dob+agein-1)/5)*5+1
        })
    }
    nickel.expand <- do.call("rbind", lapply(seq(20,95,5), trim))
    merge(nickel.expand, ewrates,
          by.x=c("agr","ygr"), by.y=c("age","year"))
})
