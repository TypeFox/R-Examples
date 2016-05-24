print.ldknn.run <-
function(x,...){
    cat("Data:\n")
    print(summary(x$data))
    cat("\nReference Level: ")
    cat(x$reference.level)
    cat("\n\nOdds:\n")
    print(rbind(x$odds))
}

