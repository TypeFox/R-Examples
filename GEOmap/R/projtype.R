`projtype` <-
function(proj=list())
  {
    if(missing(proj)) { proj = NULL }
    print("Projection Types", quote = FALSE)
    print("0 = None", quote = FALSE)
    print("1 = merc.sphr", quote = FALSE)
    print("2 = utm.sphr", quote = FALSE)
    print("3 = lambert.cc", quote = FALSE)
    print("4 = stereo.sphr", quote = FALSE)
    print("5 = utm.elps", quote = FALSE)
    print("6 = equid.cyl", quote = FALSE)
    
    print("99 = old crosson projection", quote = FALSE)

    if(!is.null(proj))
      {
        print("", quote = FALSE)    
        print(paste(sep=" ", "Current:", proj$type, proj$name), quote = FALSE)
      }
    

  }

