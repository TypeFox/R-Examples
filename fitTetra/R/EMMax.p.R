EMMax.p <-
function(Z, w, p, ptype) {
  switch(ptype,
    p.free  = apply(Z,2,weighted.mean, w),  # estimate free p's by taking the average per columns over all samples
    p.fixed = p,                            # p's are fixed
    p.HW    = EMMax.p.HW(Z,w),              # HW
    p.part  = apply(Z,2,weighted.mean, w) ) # if ptype=4: partly fixed ===> nog aanpassen, voor geval dat pi=0
}
