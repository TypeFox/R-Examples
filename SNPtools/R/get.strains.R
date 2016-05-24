get.strains <-
function(file = "http://cgd.jax.org/tools/SNPtools/Build38/sanger.snps.NCBI38.txt.gz") {

  # Open the connection SNP file.
  con = TabixFile(file)
  open(con)

  # Get the column names from the last row of the header info.
  hdr = headerTabix(con)

  # Close the connection.
  close(con)

  strains = sub("^#", "", strsplit(hdr$header[[length(hdr$header)]],
                split = "\t")[[1]])
  strains = strains[-1:-5]
  qual.cols = grep("quality", strains)
  if(length(qual.cols) > 0) {
    strains = strains[-qual.cols]
  } # if(length(qual.cols) > 0)

  return(strains)

}
