# rel2qtl: convert data from QTLRel format to R/qtl format
#
# Karl Broman
# first written: 24 Oct 2012
# last modified: 24 Oct 2012

rel2qtl <-
function(gdat, pdat, gmap)
{
  # simple checks of data
  id.gdat <- rownames(gdat)
  id.pdat <- rownames(pdat)
  if(any(!is.element(id.pdat, id.gdat))) {
    missingind <- id.pdat[!is.element(id.pdat, id.gdat)]
    stop("Individuals in pdat that are not in gdat: ",
               paste(missingind, collapse=":"))
  }
  if(any(!is.element(id.gdat, id.pdat))) {
    missingind <- id.gdat[!is.element(id.gdat, id.pdat)]
    stop("Individuals in gdat that are not in pdat: ",
               paste(missingind, collapse=":"))
  }

  # reorder pdat as in gdat
  if(!all(id.gdat == id.pdat))
    pdat <- pdat[id.gdat,,drop=FALSE]

  # marker names
  mn.gdat <- colnames(gdat)
  mn.gmap <- as.character(gmap[,1])

  if(length(mn.gdat) != length(mn.gmap)) # different numbers of markers
    stop("Different numbers of markers in gdat and gmap.")

  if(any(!is.element(mn.gdat, mn.gmap))) {
    missingmar <- mn.gdat[!is.element(mn.gdat, mn.gmap)]
    stop("Markers in gdat that are not in gmap: ",
         paste(missingmar, collapse=":"))
  }
  if(any(!is.element(mn.gmap, mn.gdat))) {
    missingmar <- mn.gmap[!is.element(mn.gmap, mn.gdat)]
    stop("Markers in gmap that are not in gdat: ",
         paste(missingmar, collapse=":"))
  }

  if(!all(mn.gmap == mn.gdat)) # reorder markers as in gmap
    gdat <- gdat[,mn.gmap,drop=FALSE]

  # determine chromosomes
  chr <- unique(as.character(gmap[,2]))
  if(any(chr=="X")) {
    warning("Cannot convert X chromosome data.")
    chr <- chr[chr != "X"]
  }

  # split up genotype data
  geno <- vector("list", length(chr))
  names(geno) <- chr
  for(i in chr) {
    geno[[i]] <- list(data = as.matrix(gdat[,gmap[,2]==i,drop=FALSE]),
                      map = gmap[gmap[,2]==i,3])
    names(geno[[i]]$map) <- as.character(gmap[gmap[,2]==i,1])

    # if chr=="X", then X chromosome, otherwise assume autosome
    class(geno[[i]]) <- "A"
  }

  # paste phenotypes
  cross <- list(geno=geno, pheno=pdat)

  # cross class
  class(cross) <- c("f2", "cross")

  # need to jitter map to move markers apart?
  if(any(sapply(cross$geno, function(a) any(diff(a$map)==0)))) {
    cross <- qtl::jittermap(cross)
  }

  # find sex column, convert to 2-level factor
  # (had a problem with the levels being "", "F", "M")
  sexcol <- grep("sex", names(cross$pheno), ignore.case=TRUE)
  if(length(sexcol) == 1) { # found sex column
    sex <- cross$pheno[,sexcol]
    if(is.factor(sex) && length(levels(sex)) > 2) {
      sex <- as.character(sex)
      cross$pheno[,sexcol] <- factor(sex, levels=sort(unique(sex)))
    }
  }

  cross
}

# end of rel2qtl.R
