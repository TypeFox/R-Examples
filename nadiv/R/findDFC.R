findDFC <- function(pedigree, exact = FALSE, parallel = FALSE, ncores = getOption("mc.cores", 2L))
{
  numeric.pedigree <- numPed(pedigree)
  ped <- cbind(numeric.pedigree, genAssign(numeric.pedigree), rep(0, dim(numeric.pedigree)[1]))
  num.out <- ped[ped[,4] >= 2, ]
  ni <- dim(num.out)[1]
  maxid <- max(num.out[,1]) 

  i <- unlist(mapply(rep, num.out[-ni, 1], each = seq((ni-1), 1)))
  j <- unlist(lapply(seq(2,ni), FUN = function(x) num.out[x:ni, 1]))
  if(exact) exct <- 1 else exct <- 0

  if(parallel) {
     wrap_DFC <- function(x){
         i.tmp <- i[min(x):max(x)]
         j.tmp <- j[min(x):max(x)]
         Cout <- .C("dfc",
	    as.integer(numeric.pedigree[, 2] - 1),
            as.integer(numeric.pedigree[, 3] - 1),
	    as.integer(i.tmp - 1),
	    as.integer(j.tmp - 1),
	    as.integer(length(i.tmp)),
	    as.integer(exct))
        Cout[[3]]
     }
     dfcs.vec <- parallel::pvec(seq.int(length(i)), FUN = wrap_DFC, mc.set.seed = FALSE, mc.silent = TRUE, mc.cores = ncores, mc.cleanup = TRUE)
     } else{ 
          Cout <- .C("dfc",
	     as.integer(numeric.pedigree[, 2] - 1),
             as.integer(numeric.pedigree[, 3] - 1),
	     as.integer(i - 1),
	     as.integer(j - 1),
	     as.integer(length(i)),
	     as.integer(exct))
          dfcs.vec <- Cout[[3]]
       }


  yes.dfcs <- which(dfcs.vec == 1)

return(list(PedPositionList = data.frame(i = i[yes.dfcs], j = j[yes.dfcs]), DFC = data.frame(i = pedigree[i[yes.dfcs], 1], j = pedigree[j[yes.dfcs], 1]), FamilyCnt = dim(unique(cbind(pedigree[i[yes.dfcs], 2:3], pedigree[j[yes.dfcs], 2:3])))[1]))
}

