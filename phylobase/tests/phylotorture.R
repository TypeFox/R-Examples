## torture-testing phylo4 objects.
library(phylobase)
library(ape)

set.seed(10101)
n <- 200
p1 <- vector("list", n)
## don't want to slow down R CMD check by doing this every time:
## n <- 10000
for (i in 1:n) {
    if (i <= n/2) {
        e <- matrix(sample(1:10, replace=TRUE, size=10), ncol=2)
    }
    else {
        e <- cbind(sample(rep(11:19, 2)), sample(1:19))
        e <- rbind(c(0, sample(11:19, 1)), e)
    }
    p1[[i]] <- try(phylo4(e), silent=TRUE)
}
OKvals <- sapply(p1, class) != "try-error"
## table(sapply(p1[!OKvals], as.character)) # I think this is causing issues with
##  R check because of different width of terminal/output, trying something simpler:
message(unique(sapply(p1[!OKvals], as.character)))
sort(unname(table(sapply(p1[!OKvals], as.character))))
if (sum(OKvals))     message("There are ", sum(OKvals), " valid trees...")

if (any(OKvals)) {
    p2 <- p1[OKvals]
    length(p2)
    has.poly <- sapply(p2, hasPoly)
    has.sing <- sapply(p2, hasSingle)
    has.retic <- sapply(p2, hasRetic)   
    message("number of trees with polytomies: ", sum(has.poly))
    message("number of trees with singletons: ", sum(has.sing))
    message("number of trees with reticulation: ", sum(has.retic))
    if (any(has.sing)) {
        p4 <- p2[has.sing]
        plot(p4[[1]])  ## gives descriptive error
        t2 <- try(plot(collapse.singles(as(p2[[1]],"phylo"))))
        ## "incorrect number of dimensions"
    }
    if (any(!has.sing)) {
        ## first tree without singles -- HANGS!
        ## don't try the plot in an R session you care about ...
        p3 <- p2[!has.sing]
        ## plot(p2[[13]])
    }
}

## elements 8 and 34 are 
## what SHOULD the rules for trees be?

## (a) reduce node numbers to 1 ... N ?
## (b) check: irreducible, non-cyclic, ... ?

## convert to matrix format for checking?

reduce_nodenums <- function(e) {
    matrix(as.numeric(factor(e)),ncol=2)
}

# make an illegal phylo4 object, does it pass checks?
# a disconnected node:

t1 <- read.tree (text="((a,b), (c,(d, e)));")
plot(t1)

broke1 <- t1
broke1$edge[broke1$edge[,2] ==9, 1] <- 9  # disconnect the node, two subtrees, ((a, b), c)  and (d,e)

try(as(broke1, "phylo4") -> tree, silent=TRUE)   # makes a phylo4  object with no warning
try(phylo4(broke1$edge), silent=TRUE)    # constructor makes a phylo4 object with no warning
## error message comes from ape, not phylo? -- AND
##   error is about singles, not disconnected nodes
## print(try(plot(tree), silent=TRUE ))  ## pdc couldn't get this to work, so temporarily commenting

# root node value != ntips + 1:

broke2 <- t1
broke2$edge[broke2$edge==6] <- 10

## warning, but no error
## plot(broke2)  ## seems to hang R CMD check??
## generates error, but it's about wrong number of tips, not wrong value at root.
message(try(as(broke2, "phylo4"), silent=TRUE))
## error regarding number of tip labels vs edges and nodes
message(try(phylo4(broke2$edge), silent=TRUE))

# switch root node value (6) with next internal node (7):

broke3 <- broke2
broke3$edge[broke3$edge==7] <- 6
broke3$edge[broke3$edge==10] <- 7

## both of the following now fail with
## "root node is not at position (nTips+1)
try(as(broke3,"phylo4") -> tree3)  # works with no error message
try(phylo4(broke3$edge))    # works with no error message
## plot(tree3)  # would work if we could create it?


# tips have larger numbers than root node:

broke4 <- t1
broke4$edge[broke4$edge==1] <- 11
broke4$edge[broke4$edge==2] <- 12
broke4$edge[broke4$edge==3] <- 13
broke4$edge[broke4$edge==4] <- 14
broke4$edge[broke4$edge==5] <- 15

message(try(as(broke4, "phylo4"), silent=TRUE))
message(try(phylo4(broke4$edge), silent=TRUE))
# print(try(plot(broke4), TRUE))   ## CAUSES R TO HANG!

###
foo <- new('phylo4')

foo@edge <- rcoal(10)$edge
message(try(plot(foo)))

foo@label <- c(rep('blah',10), rep("",9))

#####
## tree with only 2 tips: will fail under previous versions
## with "Error in if (which(nAncest == 0) != nTips + 1) { : 
##  argument is of length zero"

edge <- matrix(c(3, 1, 3, 2), byrow=TRUE, ncol=2)
try(p2 <- phylo4(edge), silent=TRUE)
