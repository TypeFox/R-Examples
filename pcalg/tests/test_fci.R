library(pcalg)

showProc.time <- local({
    pct <- proc.time()
    function() { ## CPU elapsed __since last called__
	ot <- pct ; pct <<- proc.time()
	cat('Time elapsed: ', (pct - ot)[1:3],'\n')
    }
})

########################################################
##
##       Example 1: Zhang (2008), Fig. 6, p.1882
##                  Paper with rules
##
########################################################
## Removed here: is already in ../man/fci.Rd
## -------                     -------------
showProc.time()

########################################################
##
##       Example 2: Zhang (2006), Fig. 5.2, p.198
##                  Dissertation
##
########################################################

## create the graph g
p <- 5; . <- 0
V2 <- LETTERS[1:6]
edL2 <- setNames(vector("list", length=length(V2)), V2)
edL2[[1]] <- list(edges=4,weights=1)
edL2[[2]] <- list(edges=c(4,5),weights=c(0.5,1))
edL2[[3]] <- list(edges=4,weights=1)
edL2[[4]] <- list(edges=5,weights=1)
edL2[[6]] <- list(edges=c(3,4),weights=c(1,1))
g2 <- new("graphNEL", nodes=V2, edgeL=edL2,edgemode="directed")
print.table(1*(as(g2, "matrix") != 0), zero.print=".")
##   A B C D E F
## A . . . 1 . .
## B . . . 1 1 .
## C . . . 1 . .
## D . . . . 1 .
## E . . . . . .
## F . . 1 1 . .

## hidden:
L2 <- 6

## compute the true covariance matrix of g
cov.mat2 <- trueCov(g2)
## delete rows and columns which belong to L
true.cov2 <- cov.mat2[-L2,-L2]
## transform it into a correlation matrix
true.corr2 <- cov2cor(true.cov2)

## PAG
suffStat2 <- list(C = true.corr2, n = 10^9)
true.pag2 <- fci(suffStat2, indepTest=gaussCItest, alpha = 0.99, p=p)

## define correct PAG
corr.pag2 <- rbind(c(.,.,.,2,.),
                   c(.,.,.,2,2),
                   c(.,.,.,2,.),
                   c(1,1,1,.,2),
                   c(.,3,.,3,.))

correctEst2 <- all(corr.pag2 == true.pag2@amat)
if (!correctEst2) stop("Test fci wrong: example 2!")
showProc.time()



########################################################
##
##             Example 3: random DAG
##
########################################################


set.seed(40)
##Random graph only R1-R10
g3 <- randomDAG(14,0.3)

## Define the latent variables
L3 <- c(8,10)

##pcAlgo.Perfect with true correlation matrix
##______________________________________________________
p <- 12
amat.g <- as(g3,"matrix")
colnames(amat.g) <- rownames(amat.g) <- graph::nodes(g3)
amat.g[amat.g!=0] <- 1
print.table(amat.g, zero.print=".")

##Compute the true covariance matrix of g
cov.mat3 <- trueCov(g3)

##Delete rows and columns which belong to L
true.cov3 <- cov.mat3[-L3,-L3]
##Transform it in a correlation matrix
true.corr3 <- cov2cor(true.cov3)

##PAG
suffStat3 <- list(C = true.corr3, n = 10^9)
true.pag3 <- fci(suffStat3, indepTest=gaussCItest, alpha = 0.99, p=p)

##define correct PAG
corr.pag3 <- rbind(c(.,.,2,.,.,2,.,.,.,2,2,2),
                   c(.,.,2,.,2,.,.,.,.,2,2,2),
                   c(1,1,.,.,2,.,2,2,.,.,.,.),
                   c(.,.,.,.,.,2,.,.,2,.,2,2),
                   c(.,3,3,.,.,2,2,.,2,2,2,.),
                   c(3,.,.,1,3,.,2,.,2,.,.,.),
                   c(.,.,3,.,3,3,.,.,.,.,.,.),
                   c(.,.,3,.,.,.,.,.,.,.,.,.),
                   c(.,.,.,3,3,3,.,.,.,.,.,.),
                   c(3,3,.,.,3,.,.,.,.,.,2,.),
                   c(3,3,.,1,3,.,.,.,.,1,.,2),
                   c(3,3,.,3,.,.,.,.,.,.,3,.))

correctEst3 <- all(corr.pag3 == true.pag3@amat)
if (!correctEst3) stop("Test fci wrong: example 3!")
showProc.time()


#########################################################################################
##
##      Example 4: Spirtes 1997 p.21 DAG with latent variables and p.24 PAG
##
#########################################################################################

p <- 5; . <- 0
amat4 <- rbind(c(.,.,.,.,1,1,.),
               c(.,.,.,1,.,.,1),
               c(.,.,.,1,.,1,.),
               c(.,.,.,.,1,.,.),
               c(.,.,.,.,.,.,.),
               c(.,.,.,.,.,.,1),
               c(.,.,.,.,.,.,.))
colnames(amat4) <- rownames(amat4) <- as.character(1:7)
L4 <- c(1,2)
V4 <- as.character(1:7)
edL4 <- vector("list",length=7)
names(edL4) <- V4
edL4[[1]] <- list(edges=c(5,6),weights=c(abs(rnorm(1)),abs(rnorm(1))))
edL4[[2]] <- list(edges=c(4,7),weights=c(abs(rnorm(1)),abs(rnorm(1))))
edL4[[3]] <- list(edges=c(4,6),weights=c(abs(rnorm(1)),abs(rnorm(1))))
edL4[[4]] <- list(edges=5,weights=c(abs(rnorm(1))))
edL4[[6]] <- list(edges=7,weights=c(abs(rnorm(1))))
g4 <- new("graphNEL", nodes=V4, edgeL=edL4,edgemode="directed")

## compute the true covariance matrix of g1
cov.mat4 <- trueCov(g4)

## delete rows and columns which belong to L1
true.cov4 <- cov.mat4[-L4,-L4]

## transform it into a correlation matrix
true.corr4 <- cov2cor(true.cov4)

##PAG
suffStat4 <- list(C = true.corr4, n = 10^9)
true.pag4 <- fci(suffStat4, indepTest=gaussCItest, alpha = 0.99, p=p)

##define correct PAG
corr.pag4 <- rbind(c(.,2,.,2,.),
                   c(1,.,2,.,2),
                   c(.,3,.,2,.),
                   c(1,.,2,.,2),
                   c(.,2,.,3,.))

correctEst4 <- all(corr.pag4 == true.pag4@amat)
if (!correctEst4) stop("Test fci wrong: example 4!")
showProc.time()

