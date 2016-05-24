blockpick <- function(k, gen, k.block, design=NULL, show=10, alias.block.2fis=FALSE, select.catlg=catlg){
  ## function that picks design with right number of blocks (power of 2, 2^k.block) 

  ## k.block is the number of independent block factors needed for block construction
  ##   (log2 of number of blocks)
  ## usual assumption: their interactions with treatment factors can be neglected

  ## 2^k.block-1 columns are used exclusively for blocks

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
       if (!select.catlg[[design]]$nruns == 2^k) stop("mismatch between k and design")
       gen <- select.catlg[[design]]$gen
       res <- select.catlg[[design]]$res
  }
  if (identical(gen,0) | length(gen)==0) {
      ## treat full factorials by catalogue (block.catlg)
      gen <- numeric(0)
      if (k <= max(block.catlg$k) & k.block <= ncol(block.catlg)-2){
        if (length(which(block.catlg$k==k & block.catlg$k.block==k.block))>0)
        bcols <- block.catlg[which(block.catlg$k==k & block.catlg$k.block==k.block)[1],
                         paste("b",1:k.block,sep="")]
        else stop("This full factorial is not in the catalogue.")
      }
      else bcols <- catlg[[paste(k+k.block,"-",k.block,".",1,sep="")]]$gen
      }
  g <- length(gen)
  minus <- 1
  if (g>0){
     minus <- which(gen<0)
     gen <- abs(gen)
  }

  hilf <- c(k,gen,k.block,show)
  if (!is.numeric(hilf))
      stop ("k, gen, k.block and show must be numeric.")
  if (!all(hilf == floor(hilf) & hilf > 0))
      stop ("k, gen, k.block and show must contain integer numbers.")
  if (!k >= 3) stop ("blockpick requires k>=3.")
  if (!k.block < k) stop ("blockpick requires k.block < k.")
  if (any(gen %in% 2^(0:(k-1)))) 
        stop ("gen must not contain column numbers of base factors.")
  if (any(!gen %in% 3:(2^k-1))) 
        stop ("Column numbers in gen must be smaller than ", 2^k,".")
        
  ## assignment of k.block factors to the remaining columns
  ## 
  hilf <- cols(k, gen)
  fi2s <- hilf$fi2s
        names(fi2s) <- apply(combn(k+g,2),2,function(obj) paste(Letters[obj],collapse=""))
        minus2fis <- which(apply(combn(k+g,2),2,function(obj) length(obj %in% minus)==1))
        names(fi2s)[minus2fis] <- paste("-",names(fi2s)[minus2fis],sep="")
  nfi2s <- table(fi2s)   
  fi2cols <- as.numeric(names(nfi2s))
  if (is.null(design)) {
      if (any(2^(0:(k-1)) %in% fi2cols)) res <- 3
      else if (any(nfi2s>1)) res <- 4
      else res="5+"
  }
  if (alias.block.2fis) {banned <- hilf$main; eligible <- c(fi2cols, hilf$freecols)}
  else {banned <- c(hilf$main,fi2cols); eligible <- hilf$freecols}
  if (length(eligible)<2^k.block-1 & !alias.block.2fis) 
          stop("no adequate block design found with 2fis unconfounded with blocks")
  if (length(eligible)<2^k.block-1) stop("no adequate block design found")
  if (g>0)
  perm <- t(combn(length(eligible),k.block))  ## rows are possible selections of k.block
  else {
      if (g==0 & !alias.block.2fis) if (!all(bcols %in% eligible))
             stop("no adequate block design found with 2fis unconfounded with blocks")
      perm <- matrix(which(eligible %in% bcols),nrow=1)
  }

  hilf <- digitsBase(eligible,ndigits=k) ## always k rows
         ## potential block columns as binary numbers
  ergeb <- matrix(0,nrow=nrow(perm),ncol=2^k.block-1)
         ## will contain generating and main effects columns for blocks
  banned.block <- 99
         ## will contain number of block main effects aliased with banned columns
  dependent.block <- TRUE
         ## will contain number of block main effects aliased with banned columns
  n2fis.block <- NA
         ## number of 2fis aliased with block main effects
  n2fis.clear <- NA
         ## number of 2fis aliased with block main effects
  hilfc <- Yates[1:(2^k.block-1)]
         ## combination patterns of block factors into all block main effects
  count <- 0 ## number of possibilities found
  i <- 0 ## current position in vectors
  stopp <- FALSE
  last <- numeric(0)
  
  ## loop through selections of block columns
  while (!stopp){
      perm <- perm[setdiff(1:nrow(perm),last),,drop=FALSE]
      i <- i+1
      bg <- eligible[perm[1,]]  ## block generators
     bcols <- sapply(lapply(hilfc, function(obj) mult.gen(Yates[bg[obj]])), 
                 function(obj) {h2 <- which(sapply(Yates[1:(2^k-1)],function(obj2){
                                     hh <- FALSE 
                                     if (length(obj2)==length(obj)) 
                                          if (all(obj2==obj)) hh <- TRUE
                                     hh}))
                              if (length(h2)==0) 0 else h2
                             }
               )

      last <- which(apply(perm,1,function(obj) all(eligible[obj] %in% bcols)))
           ### remove redundant rows from perm
      ergeb[i,] <- bcols  ## 2^k.block-1 potential block main effect columns
      
      if (any(bcols==0)) next
           ## dependent.block is TRUE per default and remains so for bcols with 0
      banned.block[i] <- sum(bcols %in% banned) ## generated block main effects confounded with banned columns
      dependent.block[i] <- length(table(bcols)) < length(bcols) ## k.block columns not independent
      n2fis.block[i] <- sum(nfi2s[which(as.numeric(names(nfi2s)) %in% bcols)])
      n2fis.clear[i] <- length(nfi2s[which(nfi2s==1 & !as.numeric(names(nfi2s)) %in% c(2^(0:(k-1)),gen,bcols))])
      if (banned.block[i] == 0 & !dependent.block[i]){
                 count <- count + 1
                 if (count>=show) break
             }
      if (nrow(perm) <= length(last)) stopp <- TRUE
  }  ## end of loop over permutations
  
  pick <- which(banned.block==0 & !dependent.block)

  if (length(pick)<show) show <- length(pick) 
  if (show==0 & !alias.block.2fis) stop("no adequate block design found with 2fis unconfounded with blocks")
  if (show==0) stop("no adequate block design found") else {
    ntreat <- k+length(gen)
    blocks <- 2^k.block
    blockcols <- ergeb[pick[1:show],,drop=FALSE]
    if (alias.block.2fis){
        alias.2fis.block <- vector("list",show)
        for (i in 1:show) 
                 alias.2fis.block[[i]] <- names(fi2s)[which(abs(fi2s) %in% blockcols[i,])]
    }
    else alias.2fis.block <- "none"
    nclear.2fis <- n2fis.clear[pick[1:show]]
    nblock.2fis <- n2fis.block[pick[1:show]]
    if (sum(nclear.2fis)>0) {
        clear.2fis <- vector("list",show)
        for (i in 1:show) 
        clear.2fis[[i]] <- names(fi2s)[which(fi2s %in% 
             setdiff(as.numeric(names(nfi2s))[which(nfi2s==1)], c(2^(0:(k-1)),gen,blockcols[i,])))]
        }
    else clear.2fis <- character(0)
    gen[minus] <- -gen[minus]
    list(gen=gen, 
           basics = c(nruns=2^k,nblocks=blocks, ntreat=ntreat, res.base=res), 
           blockcols=blockcols[,2^(0:(k.block-1))],
           alias.2fis.block=alias.2fis.block,
           nblock.2fis=nblock.2fis,
           nclear.2fis=nclear.2fis,
           clear.2fis=clear.2fis)
  }
}

