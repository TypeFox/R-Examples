blockpick.big <- function(k, gen, k.block, design=NULL, show=10, alias.block.2fis=FALSE, select.catlg=catlg){
  ## function that picks design right number of blocks (power of 2, 2^k.block) 
  ## it may be possible to find an appropriate design within 
  ##    the same base design with modified generator columns
  ##  This is currently not pursued (and may never be).

  ## k.block is the number of independent block factors needed for block construction
  ##   (log2 of number of blocks)
  ## usual assumption: their interactions with treatment factors can be neglected

  ## 2^k.block-1 first columns are used exclusively for blocks

  ## if requested (default), alias.block.2fis=FALSE stops confounding of treatment 2fis with blocks
  ## if this is not requested, no attempt is made to reduce the number of aliased 2fis, 
  ##       but they are reported

  if (!is.logical(alias.block.2fis)) stop("alias.block.2fis must be a logical.")
  if (!is.null(design)){ 
     if (!is.character(design)) stop("If specified, design must be a character string.")
     if (!design %in% names(select.catlg)) stop("If specified, design must be the name of a design in select.catlg.")
     if (!missing(gen)) warning("If design is specified, gen is ignored.")
     }
  if (is.null(design) & missing(gen)) stop("One of gen or design must be given.")
  if (!"catlg" %in% class(select.catlg)) stop("select.catlg must be a catalogue of designs of class catlg.")

  if (is.character(design)) {
       if (!design %in% names(select.catlg)) stop ("Design must be a valid design name.")
       if (select.catlg[[design]]$res == 3 & !alias.block.2fis)
           stop("Resolution III design. Package FrF2 does not allow to request all 2fis to be clear of blocks.")
       gen <- select.catlg[[design]]$gen
  }
  g <- length(gen)

  hilf <- c(k,gen,k.block,show)
  if (!is.numeric(hilf))
      stop ("k, gen, k.block and show must be numeric.")
  if (!all(hilf == floor(hilf) & hilf > 0))
      stop ("k, gen, k.block and show must contain positive integer numbers.")
  if (!k >= 3) stop ("blockpick requires k>=3.")
  if (!k.block < k) stop ("blockpick requires k.block < k.")
  if (any(gen %in% 2^(0:(k-1)))) 
        stop ("gen must not contain column numbers of base factors.")
  if (any(!gen %in% 3:(2^k-1))) 
        stop ("Column numbers in gen must be smaller than ", 2^k,".")
  perm <- permutations(k)

  hilf <- digitsBase(gen,ndigits=k) ## always k rows
         ## generating columns as binary numbers
  ergeb <- matrix(0,nrow=nrow(perm),ncol=g)
         ## will contain equivalent generating columns after permuting base factors
         ## block columns are the first k.block base factors
         ## are related to other factors in the same way they originally were
         ## in their position before permuting
  nfacs.block <- rep(NA, nrow=nrow(perm))
         ## will contain number of factors aliased with block generators
  n2fis.block <- rep(NA, nrow=nrow(perm))
         ## will contain check whether 2fis aliased with block generators
  if (g>1) hilfc <- combn(g,2)
         ## interactions among generated factors
  count <- 0 ## number of possibilities found
  
  ## loop through permutations of base factors
  for (i in 1:factorial(k)){
      ergeb[i,] <- as.integer(reord(hilf,perm[i,]))
      nfacs.block[i] <- sum(ergeb[i,] < 2^k.block) ## generated main effects confounded with blocks

      if (nfacs.block[i]==0 & !alias.block.2fis){
           ## if in principle eligible, check for aliasing of 2fis with blocks
           hilf1 <- digitsBase(perm[i,-(1:k.block)],ndigits=k)

           n2fis.block[i] <- sum(colSums(digitsBase(ergeb[i,])[-((k-k.block+1):k),,drop=FALSE])<=1)>0
                ## non-block rows of digitsBase(ergeb[i,]) 
                ## if these contain at most one "1", 
                ## interactions of blocks with 2fis between base factors and generated factors occur
                ## only for clean rows of ergeb, interactions AMONG generated factors are checked
       if (g > 1 & !n2fis.block[i]){
             ## check interactions among generated factors
             hilf2 <- digitsBase(ergeb[i,],ndigits=k)
             for (j in 1:choose(g,2)){
                    if (as.intBase(paste((hilf2[,hilfc[1,j]]+hilf2[,hilfc[2,j]])%%2,collapse="")) < 2^k.block) {
                           n2fis.block[i]<-TRUE
                           break
                    }
             }
             }
             if (!n2fis.block[i]) {
                 count <- count + 1
                 if (count>=show) break
             }
    }
  }  ## end of loop over permutations
  
  if (!alias.block.2fis) pick <- which(nfacs.block==0 & !n2fis.block) else pick <- which(nfacs.block==0)
      ### sometimes length of the two vectors in the first which seems to differ
      ### have not yet found the reason
  if (length(pick)<show) show <- length(pick) 
  if (show==0 & !alias.block.2fis) stop("no adequate block design found with 2fis unconfounded with blocks\n")
  if (show==0) stop("no adequate block design found\n")
  else {
    ntreat=k+length(gen)-k.block
    blocks=2^k.block
    blockcols <- 1:(blocks-1)
    if (alias.block.2fis){
        alias.2fis.block <- vector("list",show)
        nam2fis <- outer(Letters[1:ntreat],Letters[1:ntreat],function(X,Y) paste(X,Y,sep=":"))
        nam2fis <- nam2fis[upper.tri(nam2fis)]
        for (i in 1:show) {
                       treats <- c(2^(k.block:(k-1)),ergeb[pick[i],])
                       hilf <- outer(treats,treats,
                                      function(X,Y) as.integer((digitsBase(X,ndigits=k) + 
                                                                digitsBase(Y,ndigits=k)) %% 2)) 
                       hilf <- hilf[upper.tri(hilf)]
                       alias.2fis.block[[i]] <- nam2fis[hilf<blocks]
        }
    }
    else alias.2fis.block <- "none"
    gens <- ergeb[pick[1:show],,drop=FALSE]  
    list(orig=gen, basics = c(nruns=2^k,nblocks=blocks, ntreat=ntreat), 
           perms=perm[pick[1:show],,drop=FALSE],blockcols=blockcols,
           alias.2fis.block=alias.2fis.block, gen=gens)
  }
}

