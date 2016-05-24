createStack <- function(...) {
    stack <- list(...)
    
    push <- function(x)
        stack <<- c(list(x), stack)
    
    pop <- function() {
        stopifnot(length(stack) > 0)
        
        tmp <- stack[[1]]
        stack <<- stack[- 1]
        tmp
    }
    
    view <- function() {
        cat("Stack of length", length(stack), "\n\n")
        print(stack)
    }
    
    convertToList <- function()
        stack
    
    getLength <- function()
        length(stack)
    
    structure(
        function(dispatch, ...)
            switch(dispatch,
                   push = push(...),
                   pop = pop(),
                   view = view(),
                   list = convertToList(),
                   length = getLength(),
                   stop("Unknown 'dispatch': ", dispatch)),
        class = "stack")
}
push <- function(stack, x)
    stack('push', x)
pop <- function(stack)
    stack('pop')
print.stack <- function(x, ...)
    x('view')
as.list.stack <- function(x, ...)
    x('list')
length.stack <- function(x)
    x('length')