"dist.list" <-
function(type=NULL) {
  if(is.null(type)) {
     dist <- c("aep4", "cau", "emu", "exp", "texp", "gam", "gep",
               "gev",  "gld", "glo", "gno", "gov",  "gpa", "gum",
               "kap",  "kmu", "kur", "lap", "lmrq", "ln3", "nor",
               "pe3",  "ray", "revgum",     "rice", "sla", "st3",
               "tri", "wak",  "wei");
     return(dist)
  } else {
    switch(type,
      aep4=4,  cau=2,   emu=2,   exp=2,  gam=2,  gep=3,  gev=3,
       gld=4,  glo=3,   gno=3,   gov=3,  gpa=2,  gum=2,  kap=4,     kmu=2,
       kur=2,  lap=2,  lmrq=2,   ln3=3,  nor=2,  pe3=3,  ray=2,  revgum=2,
      rice=2,  sla=2,   st3=3,  texp=2,  tri=3,  wak=5,  wei=3,
      warning("The given type argument does not match a distribution")
    )
  }
}

