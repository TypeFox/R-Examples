"torsion.pdb" <-
function(pdb) {

  
  colpaste <- function(x,col.names=colnames(x)) { 
    apply(x, 1, function(row) paste(row[col.names], collapse="."))
  }
  getinds <- function(atoms,ref=atom.names) {
    sort(atom2xyz(charmatch(atoms, ref)))
  }
  repadd <- function(num, nrep=nres, toadd=nxyz) {
    c(num, rep(num, (nrep-1)) +
      rep(cumsum(rep(toadd, (nrep-1))), each=length(num)))
  }

  ##-- List atoms form each residue of each chain
  atom.data <- colpaste(pdb$atom,c("elety","resno","chain"))
  atom.list <- matrix(unlist(strsplit(atom.data,"\\.")), ncol=3, byrow=TRUE)
  res.data  <- colpaste(pdb$atom,c("resno","chain"))
  res.list  <- unique(res.data)

  atom.names <- c("N","CA","C","O","CB", "*G","*G1","*G2","*D","*D1",
                  "*D2","*E","*E1","*E2","*Z", "NH1","NH2")

  atom.greek <- c("N","CA","C","O","CB", "G","G1","G2","D","D1",
                   "D2", "E","E1","E2", "Z", "*","*")

  coords <- NULL; blank <- matrix(NA, nrow=3, ncol=length(atom.names))

  ##-- Store coords as a 3 x Natm x Nres = [xyz,atm,res] matrix
  for(i in 1:length(res.list)) {
    res.blank <- blank
    res.ind <- which(res.list[i]==res.data)

### --- Start Edit: Wed Dec  8 18:30:07 PST 2010
###    blank.ind <- charmatch(atom.list[res.ind,1], atom.names, nomatch=0) +
###      charmatch(substr(atom.list[res.ind,1],2,4), atom.greek, nomatch=0)
###
### Bug Fix for NMR structures: make sure we have no Hydrogens 
    atoms.noh <- atom.list[res.ind,1]
    atoms.noh[ grep("H",atoms.noh) ] = "H"
    
    blank.ind <- charmatch(atoms.noh, atom.names, nomatch=0) +
      charmatch(substr(atoms.noh,2,4), atom.greek, nomatch=0)
### --- End Edit
    
    res.blank[,blank.ind[blank.ind!=0]] <-
      matrix(pdb$xyz[atom2xyz(res.ind[blank.ind!=0])], nrow=3)
    coords <- cbind(coords,res.blank)
  }
  
  natm <- length(atom.names);
  nxyz <- 3*natm
  nres <- length(coords)/(nxyz)
  dim(coords) <- c(3, natm, nres)
  dimnames(coords)=list(xyz=c("x","y","z"), atm=atom.names, res=res.list)

  ##-- Torsions for selected atoms
  co <- c(coords)
  chi1  <- torsion.xyz( co[ repadd( getinds( c("N","CA","CB","*G") )) ] )
  chi11 <- torsion.xyz( co[ repadd( getinds( c("N","CA","CB","*G1") )) ])
  ###chi12 <- torsion.xyz( co[ repadd( getinds( c("N","CA","CB","*G2") )) ])

  chi2  <- torsion.xyz( co[ repadd( getinds( c("CA","CB","*G","*D") )) ])
  chi21 <- torsion.xyz( co[ repadd( getinds( c("CA","CB","*G","*D1") )) ])
  ###chi22 <- torsion.xyz( co[ repadd( getinds( c("CA","CB","*G","*D2") )) ])
  ## New catch for atom name CG1 of ILE residues 
  chi2.ILE <- torsion.xyz( co[ repadd( getinds( c("CA","CB","*G1","*D1") )) ])

  chi3  <- torsion.xyz( co[ repadd( getinds( c("CB","*G","*D","*E") )) ])
  chi31 <- torsion.xyz( co[ repadd( getinds( c("CB","*G","*D","*E1") )) ])
  ##chi32 <- torsion.xyz( co[ repadd( getinds( c("CB","*G","*D","*E2") )) ])

  chi4  <- torsion.xyz( co[ repadd( getinds( c("*G","*D","*E","*Z") )) ])

  chi51 <- torsion.xyz( co[ repadd( getinds( c("*D","*E","*Z", "NH1") )) ])
  ###chi52 <- torsion.xyz( co[ repadd( getinds( c("*D","*E","*Z", "NH2") )) ])
  
  omega <- torsion.xyz( co[ repadd( c(4:9, 52:57) ) ])
  alpha <- c(NA, torsion.xyz( co[ repadd( c(4:6,55:57,106:108,157:159) ) ]))
  phi   <- c(NA, torsion.xyz( co[ repadd( c(7:9,52:60) ) ]))
  psi   <- torsion.xyz( co[ repadd( c(1:9,52:54) ) ])
  ## alpha c("CA","CA","CA","CA")
  ## omega c("CA","C","N","CA")
  ## phi   c("C","N","CA","C")
  ## psi   c("N","CA","C","N")


##- Old Output with redundent angles (e.g. chi22 etc.).
###  out <- list(psi=psi, phi=phi[-(nres+1)], omega=omega,
###              chi1=chi1, chi11=chi11, chi12=chi12,
###              chi2=chi2, chi21=chi21, chi22=chi22,
###              chi3=chi3, chi31=chi31, chi32=chi32,
###              chi4=chi4,
###              chi51=chi51, chi52=chi52,
###              alpha=alpha[-(nres+1)], coords=coords)


  ##- New reduced output with only one chi per sidechain position  
  tor.collapse <- function(a1, a11) {
    a <- a1
    got.a11 <- !(is.na(a11))
    a[got.a11] <- a11[got.a11]
    return(a)
  }

  chi1.F <- tor.collapse(chi1, chi11)
  chi2.F <- tor.collapse(chi2, chi21)
  chi2.F <- tor.collapse(chi2.F, chi2.ILE)
  chi3.F <- tor.collapse(chi3, chi31)
  
  ## New table/matrix for output
  tbl=cbind(phi[-(nres+1)], psi, chi1.F, chi2.F, chi3.F, chi4, chi51)
  colnames(tbl) <- c("phi", "psi", "chi1", "chi2", "chi3", "chi4", "chi5")
  
  out <- list(psi=psi, phi=phi[-(nres+1)], omega=omega,
              chi1=chi1.F, ##chi11=chi11, chi12=chi12,
              chi2=chi2.F, ##chi21=chi21, chi22=chi22,
              chi3=chi3.F, ##chi31=chi31, chi32=chi32,
              chi4=chi4,
              chi5=chi51, ##chi51=chi51, chi52=chi52,
              alpha=alpha[-(nres+1)], coords=coords,
              tbl=tbl )
  
}

