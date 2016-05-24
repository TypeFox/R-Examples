# qtl2rel: convert data from R/qtl format to QTLRel format
#
# Karl Broman
# first written: 24 Oct 2012
# last modified: 24 Oct 2012

qtl2rel <-
function(cross)
{
  #require(qtl) # need the qtl package for this function

  if(class(cross)[1] != "f2" || class(cross)[2] != "cross")
    stop("Input must be an intercross.")

  n.ind <- qtl::nind(cross)

  # create IDs; I think they must be numeric
  id <- qtl::getid(cross)
  if(is.null(id))
    id <- 1:n.ind
  if(!is.numeric(id)) {
    warning("Creating arbitrary numeric IDs.")
    id <- 1:n.ind
  }

  # deal with sex info
  sexpgm <- qtl::getsex(cross)
  if(is.null(sexpgm$sex)) {
    warning("No sex information; treating all as females")
    sex <- factor(rep("F", n.ind), levels=c("F", "M"))
  }
  else { # convert to F/M
    sex <- factor(c("F","M")[sexpgm$sex+1], levels=c("F","M"))
  }

  # Convert X chromosome genotypes to simple 1/2/3 format
  chrtype <- sapply(cross$geno, class)
  if(any(chrtype =="X")) {
    # convert to simple 1/2/3 format
    for(i in which(chrtype=="X"))
      cross$geno[[i]]$data <- qtl::reviseXdata("f2", "simple", sexpgm, geno=cross$geno[[i]]$data)
  }

  # create pedigree info
  parents <- c(-1002, -1001, -102, -101)
  while(any(parents %in% id))
    parents <- parents - 1000

  ped <- data.frame(id=c(parents, id),
                    sex=factor(c("M","F","M","F",as.character(sex)), levels=c("F","M")),
                    generation=factor(c("F0","F0","F1","F1",rep("F2", n.ind)), levels=c("F0","F1","F2")),
                    sire=c(0, 0, parents[1], parents[1], rep(parents[3], n.ind)),
                    dam=c(0, 0, parents[2], parents[2], rep(parents[4], n.ind)),
                    family=factor(c("F0","F0","F1","F1",rep("F1-1", n.ind)), levels=c("F0","F1","F1-1")))

  # combine genotypes into gdat
  gdat <- qtl::pull.geno(cross)
  rownames(gdat) <- as.character(id)

  # pull out phenotypes
  pdat <- cross$pheno
  rownames(pdat) <- as.character(id)

  # construct genetic map info
  gmap <- data.frame(snp=factor(qtl::markernames(cross)),
                     chr=rep(names(cross$geno), qtl::nmar(cross)),
                     dist=unlist(lapply(cross$geno, function(a) as.numeric(a$map))))

  # return list
  list(ped=ped, gdat=gdat, pdat=pdat, gmap=gmap)
}

# end of qtl2rel.R
