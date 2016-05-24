

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# Requires: Initialized objects (as done by test-all.R header)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #


##                      1       2       3       4       5       6       7 ##
qry<-data.frame(id = 1:7, seqid = "1",
            lstart = c(10100L, 11800L, 12220L, 12220L, 12220L, 32000L, 40000L),
            lend =   c(10100L, 12000L, 12225L, 12227L, 12227L, 32100L, 40100L),
            rstart = c(10200L, 12200L, 12057L, 12613L, 12650L, 32200L, 40200L),
            rend =   c(10300L, 12250L, 12179L, 12620L, 12700L, 32300L, 40300L))
##                      1       2       3       4       5       6       7 ##


# C) 
res<-overlapJuncs(qry,junc)

if(! all(is.na(res$sod)==c(TRUE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE)) )
    stop("[test_overlap_juncs] Wrong res$sod NA's.")

if(sum(res$nref) != 27)
    stop("[test_overlap_juncs] Wrong sum of res$nref.")


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# Test overhang codes (No requirements)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
rid     <- 1L
rLstart <- 1000L
rLend   <- 2000L
rRstart <- 3000L
rRend   <- 4000L
rMaxEnd <- 4000L
#
qLstart <- as.integer(c(0700, 0800, 0900, 1000, 1000, 1000, 1000))
qLend   <- as.integer(c(0800, 0900, 2000, 2000, 2000, 2000, 5000))
qRstart <- as.integer(c(0900, 3000, 3000, 3000, 3000, 5000, 6000))
qRend   <- as.integer(c(4000, 4000, 4000, 4000, 5000, 7000, 7000))
qid     <- as.integer(1:length(qLstart))
dfr <- .Call("gap_overlap", qid, qLstart, qLend, qRstart, qRend,
             rid, rLstart, rLend, rRstart, rRend, rMaxEnd)

if(!all(dfr$ovhl==c("inp","int","ext",rep("no",4))))
    stop("[test_overlap_juncs] Wrong ovhl codes")
if(!all(dfr$ovhr==c(rep("no",4),"ext","int","inp")))
    stop("[test_overlap_juncs] Wrong ovhr codes")

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# END OF FILE
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
