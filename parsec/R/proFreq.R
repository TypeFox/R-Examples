# proFreq <- function(profiles, population, popwght=rep(1, nrow(population))) {
proFreq <- function(profiles, population) {
    popwght <- rep(1, nrow(population)) # parametro futuro per introdurre i pesi nella popolazione
    cnt <- popelem(profiles, population)
    nm <- names(profiles$freq)
    profiles$freq <- sapply(1:length(nm), function(i) sum((cnt == i)*popwght, na.rm = TRUE))
    names(profiles$freq) <- nm
    return(profiles)
}