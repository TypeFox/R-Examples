ConvDiffCoef = function(db, dt){
#  require(pracma)

  # Make db real if its imaginary part is small; otherwise error
  if(sum(abs(Im(db))) < sum(abs(Re(db))/1000)){
    db = Re(db)
  }else{
    stop('ConvDiffCoef: Im(db) is not negligible')
  }
  Pb = t(pascal(length(db), 1))
  b = as.vector(Pb %*% (db / dt^(1:length(db) - 1)))
  return(b)
}
