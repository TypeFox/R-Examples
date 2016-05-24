ridge <-
function(wavelet.spectrum, band = 5, scale.factor = 0.1){

  min.level = scale.factor * max(wavelet.spectrum)

  ridge.column = function(column.vec, band=band){
            
             nrows = length(column.vec)
             
             ind = seq(1,nrows)
             band.max.vec = column.vec
             
             for (i in (1:band)) {
                  
                  lower.ind  = ind - i
                  lower.ind[lower.ind<1] = 1
                  upper.ind  = ind + i
                  upper.ind[upper.ind>nrows] = nrows                 
                 
                  band.max.vec = pmax(band.max.vec, column.vec[lower.ind], column.vec[upper.ind])
                  
             }
            
             
             my.ridge.column = rep(0,nrows)
             my.ridge.column[pmax(band.max.vec) == column.vec] = 1
             
             return(my.ridge.column)
 
  }

  Ridge = apply(wavelet.spectrum, 2, ridge.column, band=band)

  Ridge = Ridge * (wavelet.spectrum>min.level)

  return(invisible(Ridge))
}