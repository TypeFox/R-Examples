### For work flows only.

get.workflow <- function(dir.name = "workflow/",
    model = c("roc"), pkg = "cubfits"){
  dir.name <- paste(dir.name, model[1], "/", sep = "")
  file.path <- tools::file_path_as_absolute(
                 system.file(dir.name, package = "cubfits"))
  file.path
} # End of get.workflow().

cp.workflow <- function(flow = c("wphi", "wophi", "simu", "wphi_wophi"),
    model = c("roc"), pkg = "cubfits", to = NULL, code = FALSE){
  ### Check OS.
  if(Sys.info()["sysname"] == "Windows"){
    stop("Non of the workflows supports windows system.")
  }

  ### Check flow.
  if(!any(flow[1] %in% c("wphi", "wophi", "simu", "wphi_wophi"))){
    stop("Work flow is not found.")
  }

  ### Set new path.
  if(interactive()){
    cat("This function is supposed to run in batch ...\n")
    if(is.null(to)){
      to <- tempdir()
    }
  } else{
    if(is.null(to)){
      to <- "." 
    }
  }
  path.current <- to
  if(!file.exists(path.current)){
    dir.create(path.current, mode = "0755")
  }

  ### Original work flow path.
  path.workflow <- get.workflow(model = model[1], pkg = pkg)

  ### Make new copies for basic work flows.
  cat("Copy ", flow[1], "_run_0.sh ...\n", sep = "")
  path.file <- paste(path.workflow, "/script/", flow[1], "_run_0.sh",
                     sep = "")
  path.file.new <- paste(path.current, "/run_0.sh", sep = "")
  file.copy(path.file, path.file.new, overwrite = TRUE)

  cat("Copy and rename ", flow[1], "_run_1_.sh ...\n", sep = "")
  path.file <- paste(path.workflow, "/script/", flow[1], "_run_1.sh",
                     sep = "")
  path.file.new <- paste(path.current, "/run_1.sh", sep = "")
  file.copy(path.file, path.file.new, overwrite = TRUE)

  cat("Copy and rename ", flow[1], "_run_2.sh ...\n", sep = "")
  path.file <- paste(path.workflow, "/script/", flow[1], "_run_2.sh",
                     sep = "")
  path.file.new <- paste(path.current, "/run_2.sh", sep = "")
  file.copy(path.file, path.file.new, overwrite = TRUE)

  cat("Copy and rename ", flow[1], "_run_2_ps.sh ...\n", sep = "")
  path.file <- paste(path.workflow, "/script/", flow[1], "_run_2_ps.sh",
                     sep = "")
  path.file.new <- paste(path.current, "/run_2_ps.sh", sep = "")
  file.copy(path.file, path.file.new, overwrite = TRUE)

  cat("Copy and rename 00-set_env_", flow[1], ".r ...\n", sep = "")
  path.file <- paste(path.workflow, "/code/00-set_env_", flow[1], ".r",
                     sep = "")
  path.file.new <- paste(path.current, "/00-set_env.r", sep = "")
  file.copy(path.file, path.file.new, overwrite = TRUE)

  if(flow[1] == "simu"){
    ### simu uses default param.
    cat("Copy param/ ...\n", sep = "")
    path.file <- paste(path.workflow, "/../param", sep = "")
    file.copy(path.file, path.current, overwrite = TRUE, recursive = TRUE)
  } else{
    ### Rest cases use example param.
    cat("Make param/ ...\n", sep = "")
    path.file.new <- paste(path.current, "/param/", sep = "")
    unlink(path.file.new, recursive = TRUE, force = TRUE)
    dir.create(path.file.new, mode = "0755")

    cat("Copy param/genome.fasta ... (fake)\n", sep = "")
    path.file <- get.expath("seq_200.fasta")
    path.file.new <- paste(path.current, "/param/genome.fasta", sep = "")
    file.copy(path.file, path.file.new, overwrite = TRUE)

    if(flow[1] %in% c("wphi", "wphi_wophi")){
      cat("Copy param/genome.phi.tsv ... (fake)\n", sep = "")
      path.file <- get.expath("phi_200.tsv")
      path.file.new <- paste(path.current, "/param/genome.phi.tsv", sep = "")
      file.copy(path.file, path.file.new, overwrite = TRUE)
    }
  }

  ### Copy to local for customizing only.
  if(code){
    cat("Copy code/ ...\n", sep = "")
    path.file <- paste(path.workflow, "/code", sep = "")
    file.copy(path.file, path.current, overwrite = TRUE, recursive = TRUE)
  }

  ### Done.
  cat("Done ...\n\n")
  cat("Please check \"", path.current, "/00-set_env.r\" ...\n",
      sep = "")
  if(flow[1] != "simu"){
    cat("Please update new \"", path.current, "/param/genome.fasta\" ...\n",
        sep = "")
    if(flow[1] %in% c("wphi", "wphi_wophi")){
      cat("Please update new \"", path.current, "/param/genome.phi.tsv\" ...\n",
          sep = "")
    }
  }
  cat("Please adjust all \"run_*.sh\" if needed ...\n")
  cat("\n")

  invisible()
} # End of cp.workflow().

