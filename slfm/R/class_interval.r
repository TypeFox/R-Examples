class_interval_uni <- function(x, cv = 0.5) {
    if(x[2] < cv) {
        clas <- "N"
    } else if(x[1] > cv) {
        clas <- "S"
    } else {
        clas <- "I"
    }
    return(clas)
}

class_interval <- function(x) {
    apply(x, 1, class_interval_uni)
}

format_classification <- function(x) {
  x <- as.factor(x)
  lvls <- levels(x)
  lvls[lvls == "S"] <- "S"
  lvls[lvls == "I"] <- "I"
  lvls[lvls == "N"] <- "N"
  levels(x) <- lvls
  x
}

class_mean <- function(x, cv=0.5) {
  x[,1] > cv
}