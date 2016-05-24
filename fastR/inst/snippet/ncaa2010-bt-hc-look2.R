ratings <- 
    BTabilities(ncaa.model2)
ratings[
    rev(order(ratings[,1]))[1:30],]
