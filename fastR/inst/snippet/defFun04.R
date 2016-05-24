favstats <- function(x) {
	result <- c(min(x),max(x),mean(x),median(x), sd(x))
    names(result) <- c("min","max","mean","median","sd")
    return(result)
}
favstats((1:20)^2)
summary(Sepal.Length~Species,data=iris,fun=favstats)
