output.ci <-
function(ci,lvl){
  paste(switch(ci,"Arithmetic mean","Geometric mean","Median","Median",
    "Arithmetic mean","Geometric mean"), " and its ", lvl*100,"% CI using ",
    switch(ci,"the normal approximation","the (log) normal approximation","the KC Method","bootstrap","2*SD","(log) 2*SD"),sep="")
}
