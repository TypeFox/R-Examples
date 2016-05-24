"corrgram" <-
function(x){
  gpairs(x,
         upper.pars = list(conditional = 'barcode', scatter = 'corrgram',
                           mosaic = 'mosaic'),
         lower.pars = list(conditional = 'boxplot', scatter = 'loess', mosaic = 'mosaic'),
         scatter.pars = list(pch = NA),
         outer.labels = 'none',
         gap = 0,
         diag.pars = list(show.hist = FALSE))
}

