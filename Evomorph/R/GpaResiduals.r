GpaResiduals <- function(lands,gpa_coords){
  if (missing(gpa_coords))
  {print("Gpa coordinates not provided. Using gpagen from geomorph...")
  if (missing(lands)){stop("Landmark data not provided")}
  if (!(is.array(lands)) | is.matrix(lands))
  {stop("Landmark data must be array type")}
  proc=gpagen(lands,PrinAxes = F)
  coords=proc$coords
  }
  else{if(!(is.array(gpa_coords)) | is.matrix(gpa_coords)){stop("Gpa coordinates must be array type")}
  coords=gpa_coords
  }
  coords2d=two.d.array(coords)
  consensus=apply(coords, c(1,2), mean);
  consensusvec=apply(coords2d,2,mean);
  resi=t(t(coords2d)-consensusvec);
  results=list(consens=consensus,cvectorized=consensusvec,resid=resi);
  return(results)
}

