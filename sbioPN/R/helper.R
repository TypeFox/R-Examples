## -*- ess-indent-level: 2; ess-basic-offset: 2; tab-width: 8 -*-
##
## Copyright (C) 2009-2014 Roberto Bertolusso and Marek Kimmel
##
## This file is part of sbioPN.
##
## sbioPN is free software: you can redistribute it and/or modify it
## under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 2 of the License, or
## (at your option) any later version.
##
## sbioPN is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with sbioPN. If not, see <http://www.gnu.org/licenses/>.

init <- function(place) {
  par.frm <- parent.frame()
  model <- list(transitions=0,
                place=place,
                places=length(place),
                h=list(),
                trans.name=list())
  assign("model", model, envir=par.frm)

  assign("L", array(0, dim <- c(1, model$places)), envir=par.frm)
  assign("R", array(0, dim <- c(1, model$places)), envir=par.frm)

  for (n in 1:model$places) {
    assign(place[n], n, envir = par.frm)
  }
}

atr <- function(trans.name=NULL) {
  par.frm <- parent.frame()
  model <- get("model",envir = par.frm)
  model$transitions <- model$transitions + 1
  if (!is.null(trans.name)) {
    assign(trans.name, model$transitions, envir=par.frm)
    model$trans.name[[model$transitions]] <- trans.name
  }
  if (model$transitions == 1) {
    model$pre <- get("L", envir=par.frm)
    model$post <- get("R", envir=par.frm)
  } else {
    model$pre <- rbind(model$pre,get("L", envir=par.frm))
    model$post <- rbind(model$post,get("R", envir=par.frm))
  }
  model$h[[model$transitions]] <- get("h",envir=par.frm)
  assign("model", model, envir = par.frm)

  assign("L", array(0, dim <- c(1, model$places)), envir=par.frm)
  assign("R", array(0, dim <- c(1, model$places)), envir=par.frm)
}

load.cfn <- function(place, code) {
  for (i in 1:length(place)) {
    if (i == 1) {
        concat.places <- place[i]
    } else {
      concat.places <- paste0(concat.places, ", ", place[i])
    }
  }
  this.tempfile <- basename(tempfile(pattern="tmp-cfn-", tmpdir="."))

  c.file <- paste0("typedef enum {", concat.places, "} enum_places;\n", code,"\n")
  cat(c.file, file=paste(this.tempfile, ".c", sep= ""))
  cat(c.file)

  system.str <- paste0("R CMD SHLIB ", this.tempfile, ".c")
  if (.Platform$OS.type == "windows") {
    shell(system.str)
  } else {
    system(system.str)
  }
  dyn.load(paste0(this.tempfile, .Platform$dynlib.ext))
  getNativeSymbolInfo("cfn", PACKAGE=this.tempfile)$address
}

unload.cfns <- function()  {
  files <- list.files(path = ".",
                      pattern = "^tmp-cfn-[[:xdigit:]]*\\.(so|dll|c|o)$",
                      full.names = F, recursive = F)
  for (file in files) {
    split.file <- strsplit(file, "\\.")[[1]]
    file.extension <- paste(".", split.file[2], sep="")
    if (file.extension == .Platform$dynlib.ext &
        is.loaded("cfn", split.file[1])) {
      cat(paste(file, "unloaded\n"))
      dyn.unload(file) 
    }
    cat(paste(file, "deleted\n"))
    unlink(file)
  }
}
