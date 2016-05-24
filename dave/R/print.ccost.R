print.ccost <-
function(x,...) {
    o.ccost<- x
    cat("Confusion matrix:","\n")
    print(o.ccost$conf.matrix)
    cat("\n")
    cat("Weight matrix:","\n")
    print(o.ccost$weight.matrix)
     }
