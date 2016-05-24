# Imported from coxme package

formula1 <- function(x) {
    if (class(x)=='formula') {  #top level call
        n <- length(x)  # 2 if there is no left hand side, 3 otherwise
        temp <- formula1(x[[n]])
        if (is.null(temp$fixed)) x[[n]] <- 1  # only a random term!
        else x[[n]] <- temp$fixed
        return(list(fixed=x, random=temp$random))
        }
    
    if (class(x) == '(' ) {
        if (class(x[[2]])== 'call' && x[[2]][[1]] == as.name('|')) {
            return(list(random = list(x)))
            }
            
        temp <- formula1(x[[2]])  # look inside the parenthesised object
        if (is.null(temp$fixed)) return(temp) #doubly parenthesised random 
        else {
            # A random term was inside a set of parentheses, pluck it out
            #  An example would be (age + (1|group))
            if (length(temp$fixed) <= 2) x <- temp$fixed  #remove unneeded (
            else      x[[2]] <- temp$fixed
               return(list(fixed= x, random=temp$random))
            }
        }
    if (class(x) == 'call' && x[[1]] == as.name('+')) {
        temp1 <- formula1(x[[2]])
        temp2 <- formula1(x[[3]])

        if (is.null(temp1$fixed)) {
            # The left-hand side of the '+' had no fixed terms
            return(list(fixed=temp2$fixed, 
                        random=c(temp1$random, temp2$random)))
            }
        else if (is.null(temp2$fixed)) # right had no fixed terms
            return(list(fixed=temp1$fixed, 
                        random=c(temp1$random, temp2$random)))
        else {
            return(list(fixed= call('+', temp1$fixed, temp2$fixed),
                        random=c(temp1$random, temp2$random)))
            }
        }
    if (class(x)== 'call' && x[[1]] == as.name('-')) {
        temp1 <- formula1(x[[2]])
        temp2 <- formula1(x[[3]])
        if (!is.null(temp2$random))
            stop("You cannot have a random term after a - sign")

        if (is.null(temp1$fixed))  #no fixed terms to the left
            return(list(fixed=temp2$fixed, 
                        random= temp1$random))
        else {  #there must be fixed terms to the right
            return(list(fixed= call('-', temp1$fixed, temp2$fixed),
                        random= temp1$random))
            }
       }            
    if (class(x)== 'call' && (x[[1]] == '*' || x[[1]] == ':')) {
        temp1 <- formula1(x[[2]])
        temp2 <- formula1(x[[3]])

        if (is.null(temp1$random) && is.null(temp2$random))
            return(list(fixed=x))   # The simple case, no random terms

        if (!is.null(temp1$random) && !is.null(temp2$random))
                stop ("The interaction of two random terms is not defined")
        if (is.null(temp1$fixed) || is.null(temp2$fixed)) {
            if (x[[1]] == ':') fixed <- NULL
            else if (is.null(temp1$fixed)) fixed <- temp2$fixed
            else fixed <- temp1$fixed
            }
        else  fixed <- call(deparse(x[[1]]), temp1$fixed, temp2$fixed)
        if (is.null(temp2$random))  #left hand side was random
            random <- lapply(temp1$random, 
                             function(x,y) call(':', x, y), y=temp2$fixed)
        else  #right side was
            random = lapply(temp2$random,
                                 function(x,y) call(':', x, y), y=temp1$fixed)

        if (is.null(fixed)) return(list(random= random))
        else return(list(fixed=fixed, random=random))
        }
    return(list(fixed=x))
}
formula2 <- function(term) {
    if (is.call(term) && term[[1]] == as.name(':')) {
        interact <- term[[3]]
        term <- term[[2]]
        }
    else interact <- NULL
   
    if (class(term) != '(' || !is.call(term[[2]]) || 
                              term[[2]][[1]] != as.name('|')) 
        stop("Formula error: Expected a random term") 

    term <- term[[2]]  # move past the parenthesis
    out <- list(intercept=findIntercept(term[[2]]))
    out$group<- term[[3]]
    out$interaction <- interact
    out$fixed <- term[[2]]
    out
  }
findIntercept <- function(x) {
   if (is.call(x)) {
       if (x[[1]] == as.name('+')) findIntercept(x[[2]]) |findIntercept(x[[3]])
       else FALSE
       }
   else if (x==1) TRUE
        else FALSE
}
hasAbar <- function(x) {
  if (class(x)== 'call') {
        if (x[[1]]== as.name('|')) return(TRUE)
        else if (x[[1]]==as.name( '+') || x[[1]]== as.name('-') ||
                 x[[1]]==as.name( '*') || x[[1]]== as.name(':'))
            return(hasAbar(x[[2]]) || hasAbar(x[[3]]))
        else return(FALSE)
        }
    else if (class(x) == '(') return(hasAbar(x[[2]]))
    else return(FALSE)
    }

subbar <- function(x) {
    if (class(x)=='formula') x[[length(x)]] <- subbar(x[[length(x)]])
    if (class(x)== 'call') {
        if (x[[1]]==as.name( '+') || x[[1]]== as.name('-') ||
            x[[1]]==as.name( '*') || x[[1]]== as.name(':')) {
            x[[2]] <- subbar(x[[2]])
            x[[3]] <- subbar(x[[3]])
            }
        }
    else if (class(x)== '(') {
        if (class(x[[2]])== 'call' && x[[2]][[1]] == as.name('|')) 
            x[[2]][[1]] <- as.name('+')
        else x[[2]] <- subbar(x[[2]])
        }
    x
    }
   
