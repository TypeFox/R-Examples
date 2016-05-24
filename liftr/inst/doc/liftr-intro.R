## ---- eval = FALSE-------------------------------------------------------
#  dir_docker = "~/liftr_docker/"
#  dir.create(dir_docker)
#  file.copy(system.file("docker.Rmd", package = "liftr"), dir_docker)

## ---- eval = FALSE-------------------------------------------------------
#  library("liftr")
#  docker_input = paste0(dir_docker, "docker.Rmd")
#  lift(docker_input)

## ---- eval = FALSE-------------------------------------------------------
#  drender(docker_input)

## ---- eval = FALSE-------------------------------------------------------
#  dir_rabix  = "~/liftr_rabix/"
#  dir.create(dir_rabix)
#  file.copy(system.file("rabix.Rmd", package = "liftr"), dir_rabix)

## ---- eval = FALSE-------------------------------------------------------
#  library("liftr")
#  rabix_input = paste0(dir_rabix, "rabix.Rmd")
#  lift(rabix_input)
#  drender(rabix_input)

