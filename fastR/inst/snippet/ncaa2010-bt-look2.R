# nicer output this way
ratings <- 
    BTabilities(ncaa.model)
ratings[
    rev(order(ratings[,1]))[1:30],]
