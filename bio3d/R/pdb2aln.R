"pdb2aln" <-
function(aln, pdb, id="seq.pdb", aln.id=NULL, file="pdb2aln.fa", ...) {

   # check inputs
   if(!inherits(aln, "fasta") || !is.pdb(pdb))
      stop("Incorrect type of input object:
        Should be 'pdb2aln(aln, pdb, ...)'")

   cl <- match.call()

   # Mask the gaps in the first sequence to get the 
   # reference of original alignment positions 
   aln$ali[1, is.gap(aln$ali[1,])] <- "X"
   
   aa1 <- pdbseq(pdb)
   
   if(!is.null(aln.id)) findid <- grep(aln.id, aln$id)
   if(is.null(aln.id) || length(findid)==0) {
      # do sequence-profile alignment
      naln <- seq2aln(seq2add=aa1, aln=aln, id=id, file=tempfile(), ...)
      # check if the old alignment doesn't change
      if(!identical(aln$ali, naln$ali[1:(nrow(naln$ali)-1), !is.gap(naln$ali[1,]), drop = FALSE])) 
         warning("Alignment changed! Try aln.id with the closest sequence ID in the alignment")
   } else {
      # do pairwise sequence alignment     
      if(length(findid) > 1) {
         warning(paste("Multiple entities found in alignment for id=", 
            aln.id, ". Use the first one...", sep=""))
      }
      idhit <- findid[1]

      ##- Align seq to masked template from alignment
      tmp.msk <- aln$ali[idhit, ]
      tmp.msk[is.gap(tmp.msk)] <- "X"
    
      dots = list(...)
      dots$outfile = tempfile()
      args = c(list(aln = seqbind(tmp.msk, aa1)), dots)
      seq2tmp <- do.call(seqaln.pair, args)
      
      ##- check sequence identity
      ii <- seq2tmp$ali[1,]=='X'
      ide <- seqidentity(seq2tmp$ali[, !ii])[1,2]
        
      if(ide < 0.4) {
         warning(paste("Sequence identity is too low (<40%).",
             "You may want profile alignment (set aln.id=NULL)", sep=" "))
      }
 
      ##- Insert gaps to adjust alignment
      ins <- is.gap( seq2tmp$ali[1,] )
      ntmp <- matrix("-", nrow=nrow(aln$ali), ncol=(ncol(aln$ali)+sum(ins)))
      ntmp[,!ins] <- aln$ali
    
      ## Add seq to bottom of adjusted alignment
      naln <- seqbind(ntmp, seq2tmp$ali[2,])
      rownames(naln$ali) <- c(rownames(aln$ali), id)
      naln$id <- c(aln$id, id)
   }
   
   # original alignment positions (include gaps)
   # and CA indices of PDB
   ref <- matrix(NA, nrow=2, ncol=ncol(naln$ali))
   rownames(ref) <- c("ali.pos", "ca.inds")
   ref[1, !is.gap(naln$ali[1,])] <- 1:ncol(aln$ali)
   ref[2, !is.gap(naln$ali[id,])] <- atom.select(pdb, "calpha", verbose=FALSE)$atom
   
   # remove X
   naln$ali[1, naln$ali[1,]=="X"] <- "-"
 
   if(!is.null(file)) 
      write.fasta(naln, file=file)

   out <- list(id=naln$id, ali=naln$ali, ref=ref, call=cl)
   class(out) <- "fasta"
 
   return (out)
}
