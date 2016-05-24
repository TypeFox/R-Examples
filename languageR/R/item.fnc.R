`item.fnc` <-
function(data) {
    itemMeans = tapply(data$RT, data$Item, mean) 
    itemMeans = itemMeans[sort(names(itemMeans))]
    item.dat = unique(data[ , c(2, 3, 4, 5)])
    item.dat = item.dat[order(item.dat$Item), ]
    item.dat$Means = itemMeans
    item.dat.lm = lm(Means ~ X + Y + Z, data = item.dat)
    return(item.dat.lm)
}

