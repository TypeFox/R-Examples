require(kinship2)

twindat <- c(1,3,4,2,
             2,0,0,1,
             3,8,7,1,
             4,6,5,2,
             5,0,0,2,
             6,0,0,1,
             7,0,0,2,
             8,0,0,1,
             100,3,4,1,
             101,3,4,2,
             102,3,4,2,
             103,3,4,2,
             104,3,4,2,
             105,3,4,2,
             106,3,4,2,
             107,0,0,1,
             108,0,0,1,
             201,2,1,1,
             202,2,1,1,
             203,2,1,1,
             204,2,1,1,
             205,107,102,1,
             206,108,103,2)
twindat <- matrix(twindat, ncol=4, byrow=T)
dimnames(twindat) <- list(NULL, c('id', 'dadid', 'momid', 'sex'))
twindat <- data.frame(twindat)
 
relate=data.frame(id1=c(101,102,101,104,203), id2=c(102,103,103,105,204), 
                  code=c(1,1,1,2,1))
#
# Renumber everyone as 1,2,....; makes the all.equal checks easier
indx <- sort(unique(unlist(twindat[,1:3])))
twindat$id <- match(twindat$id, indx) -1
twindat$dadid <- match(twindat$dadid, indx) -1
twindat$momid <- match(twindat$momid, indx) -1
relate$id1 <- match(relate$id1, indx) -1
relate$id2 <- match(relate$id2, indx) -1
 
# Build the pedigree and kinship
tped <- with(twindat, pedigree(id, dadid, momid, sex,
                               relation=relate))
kmat <- kinship(tped)

truth <- matrix(c(5,6, 0,
                  5,4, .25,    #parent child
                  10,11,.5,    # mz twins
                  22,12, .25,  # aunt, mz with mother
                  22, 13, .125, # aunt, dz
                  13, 14, .25,  # dz twins
                  20, 21, .5,   # mz twins
                  19, 16, 0 ,   # marry in uncle
                  19, 11, .125, # aunt who is a twin
                  19, 3,  .125), #grandmother
                byrow=T, ncol=3)

all.equal(kmat[truth[,1:2]], truth[,3])
