Vardsgr <-
function (ndpts,kvar1,kdv1,rdes) {
       # res holds the results for the fortran call
         res<-.Fortran("vdg",as.integer(ndpts),as.integer(kvar1),as.integer(kdv1),
         as.double(rdes),vdgr=double(84))
         return(res$vdgr)
      }

