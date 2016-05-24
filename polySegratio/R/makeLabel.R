`makeLabel` <-
function(x, columns=c(1,2),squash=TRUE,sep="") {

  ## Description: generates marker labels from two columns where the
  ## first column is a nucleotide sequence which is mainly blank in
  ## that it is the same as the previous one while the second column
  ## is increasing numbers (fragment size) for each nucleotide
  ## combination

  ## Arguments:
  ## x: data frame of markers including labels
  ## columns: the column numbers containing labels (default: c(1,2))
  ## squash: remove blanks (default:TRUE)
  ## sep: separator when combining two label columns (default: "")

  ## Value:
  ## returns vector of marker names
  
  return(paste(autoFill(x[,columns[1]],squash=squash),x[,columns[2]],sep=sep))
}

