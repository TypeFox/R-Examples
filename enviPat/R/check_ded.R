check_ded <-
function (formulas, deduct) 
{
    warn <- c()
    formel2 <- as.character(deduct)
    formel2 <- gsub("D", "[2]H", formel2)
    ende2 <- nchar(formel2)
    element2 <- c()
    number2 <- c()
    j <- c(1)
    while (j <= ende2) {
        if (substr(formel2, j, j) == c("[")) {
            b <- j
            while (any(substr(formel2, j, j) == c("]")) != TRUE) {
                j <- c(j + 1)
            }
            k <- j
            while (any(substr(formel2, j, j) == c("0", "1", "2", 
                "3", "4", "5", "6", "7", "8", "9")) != TRUE) {
                j <- c(j + 1)
            }
            m <- c(j - 1)
            element2 <- c(element2, substr(formel2, b, m))
        }
        if (any(substr(formel2, j, j) == c("0", "1", "2", "3", 
            "4", "5", "6", "7", "8", "9")) != TRUE) {
            k <- j
            while (any(substr(formel2, j, j) == c("0", "1", "2", 
                "3", "4", "5", "6", "7", "8", "9")) != TRUE) {
                j <- c(j + 1)
            }
            m <- c(j - 1)
            j <- c(j - 1)
            element2 <- c(element2, substr(formel2, k, m))
        }
        if (any(substr(formel2, j, j) == c("0", "1", "2", "3", 
            "4", "5", "6", "7", "8", "9")) == TRUE) {
            k <- j
            while (any(substr(formel2, j, j) == c("0", "1", "2", 
                "3", "4", "5", "6", "7", "8", "9")) == TRUE) {
                j <- c(j + 1)
            }
            m <- c(j - 1)
            j <- c(j - 1)
            number2 <- c(number2, as.numeric(substr(formel2, 
                k, m)))
        }
        j <- j + 1
    }
    element3 <- c()
    number3 <- c()
    getit <- as.character(levels(as.factor(element2)))
    for (j in 1:length(getit)) {
        element3 <- c(element3, getit[j])
        number3 <- c(number3, sum(number2[element2 == getit[j]]))
    }
    element2 <- element3
    number2 <- number3
    for (i in 1:length(formulas)) {
        warn <- c(warn, FALSE)
        element1 <- c()
        number1 <- c()
        formel1 <- as.character(formulas[i])
        formel1 <- gsub("D", "H", formel1)
        ende1 <- nchar(formel1)
        j <- c(1)
        while (j <= ende1) {
            if (substr(formel1, j, j) == c("[")) {
                b <- j
                while (any(substr(formel1, j, j) == c("]")) != 
                  TRUE) {
                  j <- c(j + 1)
                }
                k <- j
                while (any(substr(formel1, j, j) == c("0", "1", 
                  "2", "3", "4", "5", "6", "7", "8", "9")) != 
                  TRUE) {
                  j <- c(j + 1)
                }
                m <- c(j - 1)
                element1 <- c(element1, substr(formel1, b, m))
            }
            if (any(substr(formel1, j, j) == c("0", "1", "2", 
                "3", "4", "5", "6", "7", "8", "9")) != TRUE) {
                k <- j
                while (any(substr(formel1, j, j) == c("0", "1", 
                  "2", "3", "4", "5", "6", "7", "8", "9")) != 
                  TRUE) {
                  j <- c(j + 1)
                }
                m <- c(j - 1)
                j <- c(j - 1)
                element1 <- c(element1, substr(formel1, k, m))
            }
            if (any(substr(formel1, j, j) == c("0", "1", "2", 
                "3", "4", "5", "6", "7", "8", "9")) == TRUE) {
                k <- j
                while (any(substr(formel1, j, j) == c("0", "1", 
                  "2", "3", "4", "5", "6", "7", "8", "9")) == 
                  TRUE) {
                  j <- c(j + 1)
                }
                m <- c(j - 1)
                j <- c(j - 1)
                number1 <- c(number1, as.numeric(substr(formel1, 
                  k, m)))
            }
            j <- j + 1
        }
        element3 <- c()
        number3 <- c()
        getit <- as.character(levels(as.factor(element1)))
        for (j in 1:length(getit)) {
            element3 <- c(element3, getit[j])
            number3 <- c(number3, sum(number1[element1 == getit[j]]))
        }
        element1 <- element3
        number1 <- number3
        for (j in 1:length(element2)) {
            if (any(element2[j] == element1) == FALSE) {
                warn[i] <- c("TRUE")
            }
            else {
                if (as.numeric(number2[element2 == element2[j]]) > 
                  as.numeric(number1[element1 == element2[j]])) {
                  warn[i] <- c("TRUE")
                }
            }
        }
    }
    return(warn)
}


