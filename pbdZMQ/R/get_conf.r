### This file is only called by
###   "pbd*/src/Makevars.in" or "pbd*/src/Makevar.win"
### to find the default configurations from
###   "pbd*/etc${R_ARCH}/Makconf".

get.path.lib <- function(arch, fn.in, debug = FALSE){
  ### For the nm outputs to check 32- and 64-bits libraries.
  if(arch == "/i386"){
    n.zero <- 8
  } else if(arch == "/x64"){
    n.zero <- 16 
  } else{
    stop(paste(arch, " is not found.", sep = ""))
  }

  ### Find which path gcc.exe is located.
  path.env <- Sys.getenv("PATH")
  path <- unlist(strsplit(path.env, ";"))
  check.gcc <- FALSE
  for(i.path in 1:length(path)){
    dir.gcc <- gsub("\\\\", "/", path[i.path])
    path.gcc <- paste(dir.gcc, "./gcc.exe", sep = "")
    check.gcc <- file.exists(path.gcc)
    if(check.gcc){
      break
    }
  }

  ### Get root for the Rtools
  if(!check.gcc){
    if(debug){
      print(path)
    }
    stop("gcc is not found.")
  } else{
    if(debug){
      print(path.gcc)
    }
    path.rtools <- gsub("/bin.*", "", path.gcc)
  }

  ### Check libraries.
  # fn.in <- "libiphlpapi.a"
  cmd <- paste("find ", path.rtools, " | grep '", fn.in, "$'", sep = "")
  libs <- shell(cmd, intern = TRUE, ignore.stderr = TRUE)
  if(debug){
    print(libs)
  }
  if(length(libs) == 0){
    stop(paste(fn.in, " is not found in ", path.rtools, sep = ""))
  } else{
    check.bits <- FALSE
    for(i.lib in 1:length(libs)){
      fn.out <- libs[i.lib]
      cmd <- paste("nm ", fn.out, " | grep '^0' | head -1 | sed -e 's/ .*//'",
                   sep = "")
      bits <- shell(cmd, intern = TRUE, ignore.stderr = TRUE)
      zeros <- bits[1]
      if(debug){
        print(zeros)
      } 
      if(nchar(zeros) == n.zero){
        check.bits <- TRUE
        break
      }
    }

    if(check.bits){
      fn.in <- gsub("\\+", "\\\\+", fn.in)
      path.lib <- gsub(paste("/", fn.in, sep = ""), "/", fn.out)
      if(debug){
        print(path.lib)
      }
    } else{
      stop("check.bits is neither 8 nor 16.")
    }
  }

  path.lib
} # End of get.path.lib().

get.path.lib.330 <- function(arch, binpref, fn.in, debug = FALSE){
  ### For the nm outputs to check 32- and 64-bits libraries.
  if(arch == "/i386"){
    arch <- gsub("^.", "", arch)
  } else if(arch == "/x64"){
    arch <- gsub("^.", "", arch)
  } else{
    stop(paste(arch, " is not found.", sep = ""))
  }

  ### Find which path gcc.exe is located.
  if(binpref != ""){
    path.gcc <- paste(binpref, "gcc.exe", sep = "")
  } else{
    ### For Rtools33 or newer version.
    ### This dose not work.
    cmd <- paste("R --arch ", arch, " CMD config CC", sep = "")
    path.gcc <- shell(cmd, intern = TRUE, ignore.stderr = TRUE)
    path.gcc <- paste(path.gcc, ".exe", sep = "")
  }
  check.gcc <- file.exists(path.gcc)

  ### Get root for the Rtools
  if(debug){
    print(path.gcc)
  }
  path.rtools <- gsub("/bin.*", "", path.gcc)

  ### Check libraries.
  # fn.in <- "libiphlpapi.a"
  cmd <- paste("find ", path.rtools, " | grep '", fn.in, "$'", sep = "")
  libs <- shell(cmd, intern = TRUE, ignore.stderr = TRUE)
  fn.out <- libs[1]
  fn.in <- gsub("\\+", "\\\\+", fn.in)
  path.lib <- gsub(paste("/", fn.in, sep = ""), "/", fn.out)
  if(debug){
    print(path.lib)
  }

  path.lib
} # End of get.path.lib().


### For libiphlpapi.a, librpcrt4.a, and libws2_32.a
### C:\Rtools32\gcc-4.6.3\i686-w64-mingw32\lib
### C:\Rtools32\gcc-4.6.3\i686-w64-mingw32\lib64
### C:\Rtools33\mingw_32\i686-w64-mingw32\lib
### C:\Rtools33\mingw_64\x86_64-w64-mingw32\lib
get.mingw.lib <- function(arch = '', binpref = '', debug = FALSE){
  if(getRversion() < '3.3.0'){
    fn.in <- "libiphlpapi.a"
    path.lib <- get.path.lib(arch, fn.in, debug = debug)

    fn.in <- "librpcrt4.a"
    path.fn <- paste(path.lib, fn.in, sep = "")
    check.fn <- file.exists(path.fn)
    if(!check.fn){
      stop(paste(path.fn, " is not found.", sep = ""))
    }

    fn.in <- "libws2_32.a"
    path.fn <- paste(path.lib, fn.in, sep = "")
    check.fn <- file.exists(path.fn)
    if(!check.fn){
      stop(paste(path.fn, " is not found.", sep = ""))
    }

    ### Cat back to "Makefile.win".
    cat(path.lib)
  }

  invisible()
} # End of get.mingw.lib().

### For libstdc++.a
### C:\Rtools32\gcc-4.6.3\lib
### C:\Rtools32\gcc-4.6.3\lib64
### C:\Rtools33\mingw_32\lib\gcc\i686-w64-mingw32\4.9.3
### C:\Rtools33\mingw_64\lib\gcc\x86_64-w64-mingw32\4.9.3
get.stdcxx.lib <- function(arch = '', binpref = '', debug = FALSE){
  if(getRversion() < '3.3.0'){
    fn.in <- "libstdc++.a"
    path.lib <- get.path.lib(arch, fn.in, debug = debug)

    ### Cat back to "Makefile.win".
    cat(path.lib)
  }

  invisible()
} # End of get.stdcxx.lib().
