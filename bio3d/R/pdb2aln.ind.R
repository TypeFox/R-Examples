"pdb2aln.ind" <-
function(aln, pdb, inds = NULL, ...) {

   # get the new alignment; also check arguments internally
   naln <- pdb2aln(aln=aln, pdb=pdb, ...)

   if(is.null(inds)) 
      inds <- gap.inspect(aln$ali)$f.inds

   ninds <- which(naln$ref["ali.pos",] %in% inds)
   ca.inds <- naln$ref["ca.inds", ninds]

   if(any(is.na(ca.inds))) {
      warning("Gaps are found in equivalent positions in PDB")
   }

   inds.a = inds[!is.na(ca.inds)]
   inds.b = ca.inds[!is.na(ca.inds)]

   a = list(atom=inds.a, xyz=atom2xyz(inds.a)) 
   class(a) = "select"  
   b = list(atom=inds.b, xyz=atom2xyz(inds.b)) 
   class(b) = "select"  

   out = list(a = a, b = b)

   return(out)
}
