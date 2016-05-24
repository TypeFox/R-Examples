strip.quality.columns = function(variants) {

  qual.cols = grep("qual", colnames(variants))

  if(length(qual.cols) > 0) {
    variants = variants[,-qual.cols]
  } # if(length(qual.cols) > 0)

  return(variants)

}
