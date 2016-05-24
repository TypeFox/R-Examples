#' Biological Units Construction
#'
#' Construct biological assemblies/units based on a 'pdb' object.
#'
#' @details
#' A valid structural/simulation study should be performed on the biological 
#' unit of a protein system. For example, the alpha2-beta2 tetramer form of
#' hemoglobin. However, canonical PDB files usually contain the asymmetric unit of 
#' the crystal cell, which can be:
#' \enumerate{
#'      \item One biological unit 
#'      \item A portion of a biological unit 
#'      \item Multiple biological units
#' }
#' The function performs symmetry operations to the coordinates based on the 
#' transformation matrices stored in a 'pdb' object returned by 
#' \code{\link{read.pdb}}, and returns biological units stored as a list of
#' \code{pdb} objects.
#'
#' @param pdb an object of class \code{pdb} as obtained from
#'    function \code{\link{read.pdb}}.
#' @param biomat a list object as returned by \code{read.pdb} 
#'    (pdb$remark$biomat), containing matrices for 
#'    symmetry operation on individual chains to build biological units.
#'    It will override the matrices stored in \code{pdb}.
#' @param multi logical, if TRUE the biological unit is returned as a 
#'    'multi-model' \code{pdb} object with each symmetric copy a distinct
#'    structural 'MODEL'. Otherwise, all copies are represented 
#'    as separated chains.
#' @param ncore number of CPU cores used to do the calculation. By default
#'          (\code{ncore=NULL}), use all available CPU cores.
#'
#' @return 
#'    a list of \code{pdb} objects with each representing an individual 
#'    biological unit.
#'
#' @seealso \code{\link{read.pdb}}
#'
#' @author Xin-Qiu Yao
#'
#' @examples
#' \donttest{
#'    pdb <- read.pdb("2dn1")
#'    biounit <- biounit(pdb)
#'    pdb
#'    biounit
#' }
#' \dontrun{
#'    biounit <- biounit(read.pdb("2bfu"), multi=TRUE)
#'    write.pdb(biounit[[1]], file="biounit.pdb")
#'    # open the pdb file in VMD to have a look on the biological unit
#' } 
biounit <- function(pdb, biomat = NULL, multi = FALSE, ncore = NULL) {

    if(!is.pdb(pdb)) 
       stop("Please provide a 'pdb' object as obtained from 'read.pdb()'")

    if(!is.null(biomat)) 
       remarks <- biomat
    else
       remarks <- pdb$remark$biomat
   
    if(is.null(remarks))
       stop("Can't find 'remark' records for building biological units")

    ncore = setup.ncore(ncore) 

    cl <- match.call()

    if(!is.null(remarks)) { 

       # check max number of copies
       ncopy.max <- max(sapply(remarks$mat, length))

       if(!multi && ncopy.max > 10)
          cat("It is slow to represent many symmetric copies as separated chains\n Try multi = TRUE\n")

       # Are chains treated differently?
       nn <- sapply(remarks$mat, function(x) length(unique(names(x))))

       if(any(nn > 1) && multi)
          stop("Can't store as multiple models as separated symmetry operations are performed on distinct chains within one biological unit")
 
       biounits <- lapply(1:remarks$num, function(i) {
          # the transformation matrices
          mats <- remarks$mat[[i]]

          # applied to the chains
          chain <- remarks$chain[[i]]
 
          # number of copies
          ncopy <- length(mats)
  
          if(!multi) {    
             ## save copies as individual chains 

             # The original copy stored as spearated chains
             biounit0 <- lapply(chain, function(x) trim.pdb(pdb, chain=x, verbose=FALSE) )
             # available chain ID repository
             chains0 <- setdiff(c(LETTERS, letters, 0:9), chain)
             
             jch <- 1
             used.chain <- NULL
             biounit <- NULL
             for(j in 1:ncopy) {
                mt <- mats[[j]]
                chs <- strsplit(names(mats)[j], split=" ")[[1]]
                for(k in chs) {
                   bio <- biounit0[[match(k, chain)]]
                   xyz <- rbind(matrix(bio$xyz, nrow=3), 1)
                   xyz <- matrix(mt %*% xyz, nrow = 1)
                   if(! k %in% used.chain) {
                      ch <- k
                      used.chain <- c(used.chain, k)
                   } else {
                      ch <- chains0[jch]
                      jch = jch + 1
                   }
                   bio$xyz <- xyz
                   bio$atom[, "chain"] <- ch
                   bio$atom[, c("x", "y", "z")] <- round(matrix(xyz, ncol=3, byrow=TRUE), digits=3)
                   
                   biounit <- c(biounit, list(bio))
                }
             } 
             biounit <- do.call(cat.pdb, c(biounit, list(rechain = FALSE)))

#             # temporarily write the pdb of biounit and re-read it
#             tmpf <- tempfile()
#             write.pdb(biounit, file=tmpf)
#             biounit = read.pdb(tmpf, verbose=FALSE)
          } 
          else {
             ## save copies as multi-models

             # The original copy
             biounit <- trim.pdb(pdb, chain=chain, verbose=FALSE)

             xyz = rbind(matrix(biounit$xyz, nrow=3), 1)
             ll <- mclapply(2:ncopy, function(j) {

                mt <- mats[[j]]
                xyz = matrix(mt %*% xyz, nrow=1)
                xyz
             }, mc.cores = ncore )
             biounit$xyz <- rbind(biounit$xyz, do.call(rbind, ll))
             class(biounit$xyz) <- "xyz"
          }
          
          biounit$call <- cl
          return(biounit)
       } ) # end of lapply(1:remarks$num)

       ## multimeric state
       nchs <- sapply(biounits, function(x) length(unique(x$atom[, "chain"])) * nrow(x$xyz))
       mer <- c("monomer", "dimer", "trimer", "tetramer", "multimer")
       names(biounits) <- paste(remarks$method, ".determined.", 
            mer[ifelse(nchs>5, 5, nchs)], " (",  nchs, " chains)", sep="") 
#       if(length(biounits) == 1) biounits = biounits[[1]]  
    }

    return(biounits)
}
