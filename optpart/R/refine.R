refine <- function (x, clustering, ...) 
{
   UseMethod("refine")
}

refine.pco <- function (x, clustering, ax=1, ay=2, ...) 
{
    clustering <- as.integer(clustify(clustering))

    for (i in 1:max(clustering)) {
        plot(x, ax, ay)
        cat(paste("Refining cluster # ", i, "\n"))
        hilight(x, clustering, ax, ay)
        chullord(x, clustering == i, ax, ay, col = i + 1)
        new <- plotid(x)
        clustering[new] <- i
        points(x, clustering == i, ax, ay, col = i + 1)
    }
    plot(x, ax, ay)
    hilight(x, clustering, ax, ay)
    for (i in 1:max(clustering)) {
        chullord(x, clustering == i, ax, ay, col = i + 1)
    }
    out <- list(clustering=clustering)
    attr(out, "class") <- "clustering"
    return(out)
}

refine.nmds <- function (x, clustering, ax=1, ay=2, ...) 
{
    clustering <- as.integer(clustify(clustering))

    for (i in 1:max(clustering)) {
        plot(x, ax, ay)
        cat(paste("Refining cluster # ", i, "\n"))
        hilight(x, clustering, ax, ay)
        chullord(x, clustering == i, ax, ay, col = i + 1)
        new <- plotid(x)
        clustering[new] <- i
        points(x, clustering == i, ax, ay, col = i + 1)
    }
    plot(x, ax, ay)
    hilight(x, clustering, ax, ay)
    for (i in 1:max(clustering)) {
        chullord(x, clustering == i, ax, ay, col = i + 1)
    }
    out <- list(clustering=clustering)
    attr(out, "class") <- "clustering"
    return(out)
}

refine.default <- function (x,clustering, ...) 
{
    clustering <- as.integer(clustify(clustering))
 
    repeat {
        plots <- readline(' enter the plots    : ')
        if (plots == "") break
        new <- as.numeric(readline(' New cluster        : '))
        for (i in strsplit(plots,",")[[1]]){
            ord <- 1:nrow(x)
            y <- match(i,row.names(x))
            if (!is.na(y)) {
                clustering[y] <- new
            }
            else print('no such plot')
        }
    }
    out <- list(clustering=clustering)
    class(out) <- 'clustering'
    out 
}

