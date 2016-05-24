`spearman.brown` <-
function(r.xx,input=2, n.or.r = "n"){

if(n.or.r == "n"){
  r.new <- input*r.xx/(1+(input-1)*r.xx)
  out <- list(r.new=r.new)
}
if(n.or.r == "r"){
  n.new <- input*(1-r.xx)/(r.xx*(1-input))
  out <- list(n.new=n.new)
}
out
}

