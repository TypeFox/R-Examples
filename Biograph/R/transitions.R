transitions <-
function (Bdata,newnamstates)
# Called in Parameters and Remove.intrastate
# namstates needed if the function is called from Parameters to 
# create the param attribute of the Biograph object (attr(Bdata,"param))

{
# --------- Determine OR and DE and generate flow table ------------------------------------
if (missing(Bdata)) stop("Function 'transitions': Data missing.")
if (missing(newnamstates)) namstates <- StateSpace (Bdata)$namstates else
       namstates <- newnamstates
# namstates <- newnamstates
#if (exists("namstates")==FALSE | is.null(namstates)) namstates <- newnamstates

#stop ("Function 'transitions': namstates missing. First run StateSpace function. ")
numstates <- length (namstates)
nsample <- nrow(Bdata)
ns <- nchar(Bdata$path)
maxns <- max(ns)
str_char <- array(" ",c(maxns))
nntrans <- array(0,c(numstates,numstates))    
dimnames (nntrans) <- list(Origin=namstates,Destination=namstates)
Bdata$path <- as.character(Bdata$path)
for (i in 1:nsample)
 { if (ns[i] > 1) 
    { str_char <- stringf(Bdata$path[i])
      for (k in 2:(ns[i]))
      { io <- grep(substr(Bdata$path[i],k-1,k-1),namstates)
        id <- grep(substr(Bdata$path[i],k,k),namstates)
        nntrans[io,id] <- nntrans[io,id] + 1  # nntrans = flow table
   }}}                                                  
  # Allocate to each possible transition a number (for survival package etc)
tmat <- matrix(NA,ncol=numstates,nrow=numstates)
dimnames(tmat) <- list(From=namstates,To=namstates)
or <- numeric(length=numstates*numstates)
des <- numeric(length=numstates*numstates) 
kk <- 0
for (i in 1:numstates) { for (j in 1:numstates) 
     {if (nntrans[i,j] > 0) {kk <- kk + 1
                            tmat[i,j] <- kk
                            or[kk] <- i
                            des[kk] <- j}}}
#                        or and des: see alsp Putter_GLHS.r
#for  (i in 1:numstates) msdata$from[msdata$OR==namstates[i]] <- i
#for  (i in 1:numstates) msdata$to[msdata$DES==namstates[i]] <- i
ntrans <-  length(nntrans[nntrans>0]) # number of transitions
    #ntrans <- sum(apply(!is.na(tmat),1,sum)) # number of transitions 
if (ntrans==0) stop("Function transitions stops because there are no transitions in data set.")
transitions <- data.frame(cbind(Trans=1:ntrans,OR=or[1:ntrans],DES=des[1:ntrans],ORN=namstates[or[1:ntrans]],DESN=namstates[des[1:ntrans]]))
transitions$ODN <- paste(transitions$ORN,transitions$DESN,sep="")

# Possible transitions
trans_possible <- array(TRUE,c(numstates,numstates))
dimnames(trans_possible) <- list(Origin=namstates,Destination=namstates)
for (i in 1:numstates)
     { for (j in 1:numstates)
       { if (nntrans[i,j] == 0) trans_possible[i,j] <-FALSE
     }}
return (list (nsample = nsample,
              namstates = namstates,
              ntrans = ntrans,
              trans_possible = trans_possible,
              tmat = tmat, 
              transitions = transitions,
              nntrans = nntrans))
}
