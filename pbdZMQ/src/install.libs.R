### Modified from Rserve/src/install.libs.R
### For libs
files <- c("pbdZMQ.so", "pbdZMQ.so.dSYM", "pbdZMQ.dylib", "pbdZMQ.dll",
           "symbols.rds",
           "libzmq.so", "libzmq.so.dSYM", "libzmq.4.dylib",
           "libzmq.dll")
files <- files[file.exists(files)]
if(length(files) > 0){
  libsarch <- if (nzchar(R_ARCH)) paste("libs", R_ARCH, sep='') else "libs"
  dest <- file.path(R_PACKAGE_DIR, libsarch)
  dir.create(dest, recursive = TRUE, showWarnings = FALSE)
  file.copy(files, dest, overwrite = TRUE, recursive = TRUE)

  ### For Mac OSX 10.10 Yosemite and when "internal ZMQ" is asked.
  ### Overwrite RPATH from the shared library installed to the destination.
  if(Sys.info()[['sysname']] == "Darwin" && ("libzmq.4.dylib" %in% files)){
    cmd.int <- system("which install_name_tool", intern = TRUE)
    fn.pbdZMQ.so <- file.path(dest, "pbdZMQ.so")
    fn.libzmq.4.dylib <- file.path(dest, "libzmq.4.dylib")

    if(length(grep("install_name_tool", cmd.int)) > 0 &&
       file.exists(fn.pbdZMQ.so) &&
       file.exists(fn.libzmq.4.dylib)){

      cmd.ot <- system("which otool", intern = TRUE) 
      if(length(grep("otool", cmd.ot)) > 0){
        rpath <- system(paste(cmd.ot, " -L ", fn.pbdZMQ.so, sep = ""),
                        intern = TRUE)
        cat("\nBefore install_name_tool:\n")
        print(rpath)
      }

      org <- file.path(getwd(), "zmq/lib/libzmq.4.dylib")
      cmd <- paste(cmd.int, " -change ", org, " ", fn.libzmq.4.dylib, " ",
                   fn.pbdZMQ.so, sep = "")
      cat("\nIn install_name_tool:\n")
      print(cmd) 
      system(cmd)

      if(length(grep("otool", cmd.ot)) > 0){
        rpath <- system(paste(cmd.ot, " -L ", fn.pbdZMQ.so, sep = ""),
                        intern = TRUE)
        cat("\nAfter install_name_tool:\n")
        print(rpath)
      }
    }
  }
}

### For etc
file <- "Makeconf"
if(file.exists(file)){
  etcarch <- if (nzchar(R_ARCH)) paste("etc", R_ARCH, sep='') else "etc"
  dest <- file.path(R_PACKAGE_DIR, etcarch)
  dir.create(dest, recursive = TRUE, showWarnings = FALSE)
  file.copy(file, dest, overwrite = TRUE)
}

### For zmq
# dir.zmq <- "./zmq"
# if(file.exists(dir.zmq)){
#   libarch <- if (nzchar(R_ARCH)) paste("lib", R_ARCH, sep='') else "lib"
#   dest <- file.path(R_PACKAGE_DIR, libarch)
#   dir.create(dest, recursive = TRUE, showWarnings = FALSE)
#   files <- paste(dir.zmq, c("/include", "/lib") , sep = "")
#   file.copy(files, dest, overwrite = TRUE, recursive = TRUE)
# }

