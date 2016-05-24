sky.colors = function(n){
  nred = floor(n/2)
  nblue = ceiling(n/2) - 1
  sred = seq(0.25, 1, length.out = nred + 1)[1:nred]
  sblue = seq(1, 0.25, length.out = nblue + 1)[1:nblue + 1]

  reds = rgb(1, sred, sred^3)
  white = rgb(1, 1, 1)
  blues = rgb(sblue^3, sblue, 1)
  return(rev(c(reds, white, blues)))
}
