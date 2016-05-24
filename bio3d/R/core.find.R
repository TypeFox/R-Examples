core.find <- function(...)
  UseMethod("core.find")

core.find.default <- function(xyz, ...)
  core.find.pdbs(xyz, ...)

core.find.pdb <- function(pdb, verbose=TRUE, ...) {
  if(nrow(pdb$xyz)<4)
    stop("provide a multi model PDB file with 4 or more frames")

  inds1 <- atom.select(pdb, "calpha", verbose=verbose)
  inds2 <- atom.select(pdb, "nucleic", elety="P", verbose=verbose)
  inds <- combine.select(inds1, inds2, operator="OR", verbose=verbose)

  tmp <- trim.pdb(pdb, inds)
  core <- core.find.pdbs(tmp$xyz, verbose=verbose, ...)

  ## map to pdb inds
  full.ids <- paste(pdb$atom$elety, pdb$atom$resno, pdb$atom$chain, sep="-")
  tmp.ids  <- paste(tmp$atom$elety, tmp$atom$resno, tmp$atom$chain, sep="-")

  core$step.inds  <- which(full.ids %in% tmp.ids[core$step.inds])[core$step.inds]
  core$atom       <- which(full.ids %in% tmp.ids[core$atom])
  core$c1A.atom   <- which(full.ids %in% tmp.ids[core$c1A.atom])
  core$c0.5A.atom <- which(full.ids %in% tmp.ids[core$c0.5A.atom])

  core$xyz        <- atom2xyz(core$atom)
  core$c1A.xyz    <- atom2xyz(core$c1A.atom)
  core$c0.5A.xyz  <- atom2xyz(core$c0.5A.atom)

  ##core$resno       <- which(full.ids %in% tmp.ids[core$resno])[core$resno]
  ##core$c1A.resno   <- which(full.ids %in% tmp.ids[core$c1A.resno])[core$c1A.resno]
  ##core$c0.5A.resno <- which(full.ids %in% tmp.ids[core$c0.5A.resno])[core$c0.5A.resno]

  core$resno       <- tmp$atom$resno[core$resno]
  core$c1A.resno   <- tmp$atom$resno[core$c1A.resno]
  core$c0.5A.resno <- tmp$atom$resno[core$c0.5A.resno]

  return(core)
}


"core.find.pdbs" <-
function(pdbs,
         shortcut  = FALSE,
         rm.island = FALSE,
         verbose   = TRUE,
         stop.at   = 15,
         stop.vol  = 0.5,
         write.pdbs = FALSE,
         outpath="core_pruned",
         ncore = 1,
         nseg.scale = 1, ...) {

  ##  Itterative core deffination for lsq fit optimisation
  ##  (core positions are those with low ellipsoid volume)

  # Parallelized by parallel package (Fri Apr 26 16:49:38 EDT 2013)
  ncore <- setup.ncore(ncore)
  if(ncore > 1) {
     # Issue of serialization problem
     # Maximal number of cells of a double-precision matrix
     # that each core can serialize: (2^31-1-61)/8
     R_NCELL_LIMIT_CORE = 2.68435448e8
     R_NCELL_LIMIT = ncore * R_NCELL_LIMIT_CORE

     if(nseg.scale < 1) {
        warning("nseg.scale should be 1 or a larger integer\n")
        nseg.scale=1
     }
  }

  error.ellipsoid<-function(pos.xyz) {
    S<-var(pos.xyz)
    prj  <- eigen(S, symmetric = TRUE)
    prj$values[prj$values < 0 & prj$values >= -1.0E-12]<-1.0E-12
    vol<-4/3*pi*prod( sqrt( prj$values ) )
    out<-list(vol=vol, U=prj$vectors, L=prj$values)
  }


  if(is.matrix(pdbs)) {
    xyz <- pdbs

    xyz.inds <- which(apply(is.na( xyz ), 2, sum)==0)
    res.inds<-xyz.inds[seq(3,length(xyz.inds),by=3)]/3

    pdbseq = rep("ALA",length(xyz.inds)/3)
    pdbnum = c(1:(length(xyz.inds)/3))

  } else {
    if( (is.list(pdbs)) && (class(pdbs)=="pdbs") ) {
      xyz=pdbs$xyz

      xyz.inds <- which(apply(is.na( xyz ), 2, sum)==0)
      res.inds <- which(apply(pdbs$ali=="-", 2, sum)==0)

      pdbseq = aa123(pdbs$ali[1,]); pdbnum = pdbs$resno[1,]
    } else {
      stop("input 'pdbs' should either be:
               a list object from 'read.fasta.pdb'
               or a numeric 'xyz' matrix of aligned coordinates")
    }
  }

  # First core = all non gap positions
  res.still.in <- res.inds # indices of core residues
  xyz.still.in <- xyz.inds # indices of core xyz's
  new.xyz.inds <- xyz.inds # indices of core xyz's
  xyz.moved    <- xyz      # core-fitted coords
  throwout.res <- NULL     # non-core res inds
  throwout.xyz <- NULL     # non-core xyz inds
  remain.vol   <- NULL
  core.length  <- NULL


  fit.to = rep(FALSE,ncol(xyz.moved))        # Preliminary fitting
  fit.to[ as.vector(xyz.still.in) ]<-TRUE    # on first structure
#  xyz.tmp <- t(apply(xyz.moved, 1,           # to find mean structure
#                       rot.lsq,              # for next fitting
#                       yy=xyz.moved[1,],
#                       xfit=fit.to))
  xyz.tmp <- fit.xyz(xyz.moved[1,], xyz.moved, which(fit.to), which(fit.to),
                      ncore=ncore, nseg.scale=nseg.scale)

  mean.xyz <- apply(xyz.tmp,2,mean)

  if(write.pdbs) { dir.create(outpath,FALSE)  }

  while(length(res.still.in) > stop.at) {

    # Core fitting, (core => pdbnum[ res.still.in ])
    fit.to = rep(FALSE,ncol(xyz.moved))
    fit.to[ as.vector(xyz.still.in) ]<-TRUE
#    xyz.moved <- t(apply(xyz.moved, 1,
#                         rot.lsq,
#                         #yy=xyz.moved[1,],
#                         yy=mean.xyz,
#                         xfit=fit.to))
    xyz.moved <- fit.xyz(mean.xyz, xyz.moved, which(fit.to), which(fit.to),
                    ncore=ncore, nseg.scale=nseg.scale)

    mean.xyz <- apply(xyz.moved,2,mean)

    i<-1; j<-3
    volume<-NULL # ellipsoid volume
    if(ncore > 1) {
       e <- mclapply(1:(length(new.xyz.inds)/3), function(j) {
          error.ellipsoid( xyz.moved[, new.xyz.inds[atom2xyz(j)]] )$vol
       })
       volume <- unlist(e)
    } else {
       while(j<=length( new.xyz.inds )) {
         e<-error.ellipsoid(xyz.moved[,new.xyz.inds[i:j]])
         volume<-c(volume,e$vol)
         i<-i+3;j<-j+3
       }
    }

    record <- cbind(res.still.in ,   # store indices and volumes
                    matrix(new.xyz.inds,ncol=3,byrow=3),
                    volume)

    # Find highest volume (most variable position)
    if (shortcut) {
      if (length(res.still.in) >= 35) {
        # remove four at a time
        highest.vol.ind <- rev(order(volume))[1:4]
      } else { highest.vol.ind <- which.max(volume) }
    } else {
      # no shortcut rm one at a time
      highest.vol.ind <- which.max(volume)
    }

    if (rm.island) {
      # Exclude length 4 residue islands
      check <- bounds( res.still.in )
      check.ind <- which(check[,"length"] < 4)
      if ( length(check.ind) > 0 ) {
        res.cut=NULL
        for (r in 1:length(check.ind)) {
          res.cut <- c(res.cut, check[check.ind[r],"start"]:
                       check[check.ind[r],"end"])
        }
        highest.vol.ind <- unique( c(highest.vol.ind,
                             which( is.element(res.still.in, res.cut)) ))
      }
    }

    # rm position from "new.xyz.inds"
    xyz.exclude <- record[highest.vol.ind,c(2:4)]
    inds.torm <- which(is.element( new.xyz.inds, as.vector(xyz.exclude) ))
    new.xyz.inds <- new.xyz.inds[ -inds.torm ]

    # Store details of the residue we excluded
    tmp.vol <- sum(record[-highest.vol.ind,5])
    throwout.res <- c( throwout.res, as.vector(record[highest.vol.ind,1]))
    throwout.xyz <- rbind( throwout.xyz, record[highest.vol.ind,2:4] )
    remain.vol   <- c(remain.vol, tmp.vol)
    res.still.in <- record[-highest.vol.ind,1]
    xyz.still.in <- record[-highest.vol.ind,2:4]
    core.length  <- c(core.length,length(res.still.in))

    if(verbose) {
      # Progress report
      cat( paste(" core size",length(res.still.in),"of",
                length(res.inds))," vol =",
          round(tmp.vol,3),"\n" )

      if(write.pdbs) {
        # Write current core structure
        write.pdb(file  = file.path(outpath, paste("core_",
                    sprintf("%04.0f", length(res.still.in)),".pdb",sep="")),
                  #xyz   = xyz[1, new.xyz.inds ],
                  xyz   = mean.xyz[ new.xyz.inds ],
                  resno = pdbnum[ res.still.in ],
                  resid = pdbseq[ res.still.in ],
                  b     = round((volume[-highest.vol.ind] /
                    max(volume[-highest.vol.ind]) * 1),2) )
      }
    }

    if(tmp.vol < stop.vol) {
      cat(paste(" FINISHED: Min vol (",stop.vol,") reached\n"))
      break
    }
  }

  # ordered thro-out lists
  ordered.res<-as.vector(c(throwout.res, res.still.in))
  ordered.xyz<-rbind(throwout.xyz, xyz.still.in)
  rownames(ordered.xyz)=NULL
  vol = c(remain.vol, rep(NA,stop.at))
  len = c(core.length,rep(NA,stop.at))
  blank<-rep(NA, len[1]); blank[na.omit(len)]=na.omit(vol)
  ordered.vol<-c(rev(blank),NA); blank[na.omit(len)]=na.omit(len)
  ordered.len<-c(rev(blank),NA)

  # sample cores (volume < 1 A^3, < 0.5 A^3, or the final core)
  if( min(ordered.vol,na.rm=TRUE) < 1) {
    a.atom <- sort(ordered.res[which(ordered.vol<1)[1]:length(ordered.vol)])
    a.xyz  <- sort(as.vector(ordered.xyz[which(ordered.vol<1)[1]:
                                         length(ordered.vol),]))
    a.resno <- as.numeric(pdbnum[a.atom])
  } else {
    a.atom  <- NULL
    a.xyz   <- NULL
    a.resno <- NULL
  }
  if( min(ordered.vol,na.rm=TRUE) < 0.5) {
    b.atom <- sort(ordered.res[which(ordered.vol<0.5)[1]:length(ordered.vol)])
    b.xyz  <- sort(as.vector(ordered.xyz[which(ordered.vol<0.5)[1]:
                                         length(ordered.vol),]))
    b.resno <- as.numeric(pdbnum[b.atom])
  } else {
    b.atom  <- NULL
    b.xyz   <- NULL
    b.resno <- NULL
  }
  tmp.inds <- which(!is.na(ordered.vol))
  c.atom <- sort(ordered.res[tmp.inds[length(tmp.inds)]:length(ordered.vol)])
  c.xyz  <- atom2xyz(c.atom)
  c.resno <- as.numeric(pdbnum[c.atom])

  output <- list(volume      = ordered.vol,
                 length      = ordered.len,
                 resno       = pdbnum[ ordered.res ],
                 step.inds   = ordered.res,
                 atom        = c.atom,
                 xyz         = c.xyz,
                 c1A.atom    = a.atom,
                 c1A.xyz     = a.xyz,
                 c1A.resno   = a.resno,
                 c0.5A.atom  = b.atom,
                 c0.5A.xyz   = b.xyz,
                 c0.5A.resno = b.resno )
  class(output)="core"; return(output)

}

