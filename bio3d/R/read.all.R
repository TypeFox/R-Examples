"read.all" <-
function(aln, prefix ="", pdbext="", sel=NULL, ...) {

  ## Usage:
  ## sel <- c("N", "CA", "C", "O", "CB", "*G", "*D",  "*E", "*Z")
  ## pdbs.all <- read.all(aln, sel=sel)
  
  files  <- paste(prefix, aln$id, pdbext,sep="")

  ##cat(files,sep="\n")
  toread <- file.exists(files)

  ## check for online files
  toread[ substr(files,1,4)=="http" ] <- TRUE


  if(all(!toread))
    stop("No corresponding PDB files found")


  coords <- NULL; res.nu <- NULL
  res.bf <- NULL; res.ch <- NULL
  res.id <- NULL
  blank <- rep(NA, ncol(aln$ali))
  ## all atom data
  coords.all <- NULL
  elety.all <- NULL; resid.all <- NULL; resno.all <- NULL
  
  for (i in 1:length(aln$id)) {

    cat(paste("pdb/seq:",i,"  name:", aln$id[i]),"\n")

    if(!toread[i]) {
      warning(paste("No PDB file found for seq", aln$id[i],
              ": (with filename) ",files[i]), call.=FALSE)
      coords <- rbind(coords, rep(blank,3))
      res.nu <- rbind(res.nu, blank)
      res.bf <- rbind(res.bf, blank)
      res.ch <- rbind(res.ch, blank)
      res.id <- rbind(res.id, blank)
      ##
      ##coords.all
      ##
    } else {
      pdb <- read.pdb( files[i], verbose=FALSE, ... )
      pdbseq  <- aa321(pdb$atom[pdb$calpha,"resid"])
      aliseq  <- toupper(aln$ali[i,])
      tomatch <- gsub("X","[A-Z]",aliseq[aliseq!="-"])
      
      ##-- Search for ali residues (1:15) in pdb
      start.num <- regexpr(pattern = paste(c(na.omit(tomatch[1:15])),collapse=""),
                           text = paste(pdbseq,collapse=""))[1]
      if (start.num == -1) {
        stop("Starting residues of sequence not found in PDB")
      }

      ##-- Numeric vec, 'nseq', for mapping aln to pdb
      nseq <- rep(NA,length(aliseq))
      ali.res.ind <- which(aliseq != "-")
      if( length(ali.res.ind) > length(pdbseq) ) {
        warning(paste(aln$id[i],
         ": sequence has more residues than PDB has Calpha's"))
        ali.res.ind <- ali.res.ind[1:length(pdbseq)] ## exclude extra
        tomatch <-  tomatch[1:length(pdbseq)]        ## terminal residues
      }
      nseq[ali.res.ind] = start.num:((start.num - 1) + length(tomatch))

      ##-- Check for miss-matchs
      match <- aliseq != pdbseq[nseq] 
      if ( sum(match, na.rm=TRUE) >= 1 ) {
        mismatch.ind <- which(match)
        mismatch <- cbind(aliseq, pdbseq[nseq])[mismatch.ind,]
        n.miss <- length(mismatch.ind)

        if(sum(mismatch=="X") != n.miss) { ## ignore masked X res        
          details <- rbind(aliseq, !match,
                           pdbseq[nseq],
                           pdb$atom[pdb$calpha,"resno"][nseq] )
                           #### calpha[,"resno"][nseq] ) ###- typo??
          rownames(details) = c("aliseq","match","pdbseq","pdbnum")
          msg <- paste("ERROR:", aln$id[i],
                       "alignment and pdb sequences do not match")
          cat(msg,"\n"); print(details); cat(msg,"\n")
          print( cbind(details[,mismatch.ind]) )
          stop(msg)
        }
      }
      
      ##-- Store nseq justified PDB data
      ca.ali <- pdb$atom[pdb$calpha,][nseq,]
      coords <- rbind(coords, as.numeric( t(ca.ali[,c("x","y","z")]) ))
      res.nu <- rbind(res.nu, ca.ali[, "resno"])
      res.bf <- rbind(res.bf, as.numeric( ca.ali[,"b"] ))
      res.ch <- rbind(res.ch, ca.ali[, "chain"])
      res.id <- rbind(res.id, ca.ali[, "resid"])

      raw <- store.atom(pdb)
      if(is.null(sel)) {
        coords.all <- rbind(coords.all, as.numeric( raw[c("x","y","z"),,nseq] ) )
        elety.all <- rbind(elety.all, c(raw[c("elety"),,nseq]) )
        resid.all <- rbind(resid.all, c(raw[c("resid"),,nseq]) )
        resno.all <- rbind(resno.all, c(raw[c("resno"),,nseq]) )


      } else {
        coords.all <- rbind(coords.all, as.numeric( raw[c("x","y","z"), sel, nseq] ) )
        elety.all <- rbind(elety.all, c(raw[c("elety"),sel,nseq]) )
        resid.all <- rbind(resid.all, c(raw[c("resid"),sel,nseq]) )
        resno.all <- rbind(resno.all, c(raw[c("resno"),sel,nseq]) )
      } 
##      raw <- store.main(pdb)
##      b <- cbind(b, raw[,,nseq])

    } # end for
  } # end else

  rownames(aln$ali) <- aln$id
##  out<-list(xyz=coords, resno=res.nu, b=res.bf,
##            chain = res.ch, id=aln$id, ali=aln$ali)
  out<-list(xyz=coords, all=coords.all, resno=res.nu, b=res.bf,
            chain = res.ch, id=aln$id, ali=aln$ali, resid=res.id,
            all.elety=elety.all, all.resid=resid.all, all.resno=resno.all)

  atm <- rep( rep(sel,each=3), ncol(aln$ali))
  colnames(out$all) = atm
  atm <- rep( sel, ncol(aln$ali))
  colnames(out$all.elety) = atm
  colnames(out$all.resid) = atm
  colnames(out$all.resno) = atm
  
  class(out)=c("pdbs", "fasta")
  return(out)
  
}

