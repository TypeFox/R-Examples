readnet <- function(con=file("default.net")) {
    ## read from .net file and create a network object.
    ## note: not all info from the .net file is used, so information
    ## may be lost (!) if overwriting the .net file with savenet(nw)
    ## The function is not foolproof if the .net files do not have the
    ## same structure as the deal generated .net files or the hugin
    ## net files we have seen after manipulating a Deal net file.
#    fn <- filename
#    zz <- file(fn,"r")
    open(con,"r")
    l <- readLines(con)
    lno <- length(l)

    lcount <- 0
    nodes <- list()
    nodecount <- 1
    nnames <- c()
    ## look for line with 'node' in it
    while (lcount <= lno) {
        lcount <- lcount + 1
        nodeptr <- grep(" node ",l[lcount],value=TRUE)
        poteptr <- grep("potential ",l[lcount],value=TRUE)
        if (length(nodeptr)>0) {
            ## we have a node definition
            ss <- unlist(strsplit(l[lcount]," "))
            ss <- ss[ss!=""]
            nd         <- list()
            nd$idx     <- nodecount
            nd$type    <- ss[1]
            nd$name    <- ss[3]
            nnames <- c(nnames,ss[3])

            ## read position
            i <- 0
            slut <- FALSE
            while (!slut) {
                i <- i+1
                posstr <- grep("position",l[lcount+i],value=TRUE)
                if (length(posstr)>0) slut <- TRUE
            }

            ## extract coordinates from posstr
            c1 <- regexpr("[(]",posstr)
            x <- substr(posstr,c1+1,nchar(posstr)-2)
            y <- unlist(strsplit(x," "))
            y <- y[y!=""]
            
            nd$position<- as.numeric(y)

            ## read levels if discrete
            if (nd$type=="discrete") {

                i <- 0
                slut <- FALSE
                while (!slut) {
                    i <- i+1
                    statestr <- grep("states",l[lcount+i],value=TRUE)
                    if (length(statestr)>0) slut <- TRUE
            }
                
                ## extract states from statestr
                c1 <- regexpr("[(]",statestr)
                x <- substr(statestr,c1+1,nchar(statestr)-2)
                x <- gsub("\"","",x)
                y <- unlist(strsplit(x," "))
                y <- y[y!=""]
                nd$levelnames <- y
                nd$levels     <- length(nd$levelnames)
            }
            
            class(nd) <- "node"
            
            nodes[[nodecount]] <- nd
            nodecount <- nodecount + 1
        }
        if (length(poteptr)>0) {
            ## we have a potential definition
            str <- poteptr
            
            c1 <- regexpr("[(]",str)
            c2 <- regexpr("[)]",str)
            x <- substr(str,c1+1,c2-1)
            c3 <- regexpr("[|]",x)

            if (c3==-1) { ## no conditional
                x <- gsub(" ","",x)
                nodenumber <- match(x,nnames)
                nodes[[nodenumber]]$parents <- c()
            }
            else { ## potentials
                lhs <- gsub(" ","",substr(x,1,c3-1))
                nodenumber <- match(lhs,nnames)
                
                rhs <- substr(x,c3+1,nchar(x))
                rhsy <- unlist(strsplit(rhs," "))
                rhsy <- rhsy[rhsy!=""]
                parents <- match(rhsy,nnames)
                nodes[[nodenumber]]$parents <- parents
            }
        }
    }

    ## update network info
    nw <- list()
    names(nodes) <- nnames
    nw$nodes <- nodes
    nw$n <- length(nodes)
    ltype <- unlist(lapply(nw$nodes,function(x) x$type))
    nw$discrete <- (1:nw$n)[ltype=="discrete"]
    nw$continuous <- (1:nw$n)[ltype=="continuous"]
    nw$nd <- length(nw$discrete)
    nw$nc <- length(nw$continuous)
    class(nw) <- "network"

    close(con)
    nw
}
