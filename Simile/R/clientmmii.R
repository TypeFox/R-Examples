is.dummy <- function(path.or.handle) {
  any(c("dummy/path", "dummy.dll", "dummy_mh", "dummy_ih", "dummy_ph") ==
      as.character(path.or.handle)[1])
}

use.simile.at <- function(path.to.installation) {
  tcl("set", "::loadedFromR", 1) # lets Tcl client know R is using it
  if (!is.dummy(path.to.installation)) {
    tcl("source", file.path(find.package(package = "Simile"), "exec",
        "client5d.tcl"))
    tcl("UseSimileAt", path.to.installation)
  }
}

load.model <- function(path.to.binary) {
  if (is.dummy(path.to.binary)) {
    "dummy_mh"
  } else {
    tcl("loadmodel", path.to.binary, "R")
  }
}

list.objects <- function(model.handle) {
  if (is.dummy(model.handle)) {
    c("/sector","/sector/output","/sector/flow1","/sector/variation")
  } else {
    as.character(tcl("ListObjPaths", model.handle))
  }
}

get.model.property <- function(model.handle, caption.path, requested.property) {
  if (is.dummy(model.handle)) {
    if (requested.property == "Dims") {
      c(2,5)
    } else {
      "COMPARTMENT"
    }
  } else {
    tcl.result <- tcl("GetModelProperty", model.handle, caption.path,
                      requested.property)
    if (any(c("Dims")==requested.property)) {
# may be more integer cases
      with.trailing.zero <- as.integer(tcl.result)
      with.trailing.zero[-length(with.trailing.zero)] # removes it
    } else if (any(c("MinVal","MaxVal")==requested.property)) {
# may be more integer cases
      as.double(tcl.result)
    } else {
      as.character(tcl.result)
    }
  }
}

create.model <- function(model.handle) {
  if (is.dummy(model.handle)) {
    "dummy_ih"
  } else {
    tcl("CreateModel", model.handle)
  }
}

set.model.step <- function(instance.handle, step.index, step.size) {
  if (!is.dummy(instance.handle)) {
    tcl("c_setstepmodel", instance.handle, step.size, step.index)
  }
}

create.param.array <- function(instance.handle, param.name) {
  if (is.dummy(instance.handle)) {
    "dummy_ph"
  } else {
    tcl("CreateParamArray", instance.handle, param.name)
  }
}

set.model.parameter <- function(param.handle, data, as.enum.types = FALSE) {
  if (!is.dummy(param.handle)) {
  tcl("SetParamArrayFromFlatList", param.handle, data, as.enum.types, 
      dim(data))
  }
}

consult.parameter.metafile <- function(instance.handle, param.file,
                                       target.submodel = "") {
  if (!is.dummy(instance.handle)) {
    tcl("ConsultParameterMetafile", instance.handle, param.file, 
        target.submodel)
  }
}

reset.model <- function(instance.handle, depth, integration.method = "Euler",
                        starting.time = 0) {
  if (!is.dummy(instance.handle)) {
    tcl("DoResetModel", instance.handle, starting.time, integration.method, 
        depth)
  }
}

execute.model <- function(instance.handle, finish.time,
                          integration.method = "Euler", start.time = NA,
                          error.limit = 0, pause.on.event = FALSE) {
  if (!is.dummy(instance.handle)) {
    if (is.na(start.time)) {
      start.time <- get.model.time(instance.handle)
    }
    as.integer(tcl("DoExecuteModel", instance.handle, integration.method,
                   start.time, finish.time, error.limit, pause.on.event))
  }
}

get.model.time <- function(instance.handle) {
  if (is.dummy(instance.handle)) {
    0.0
  } else {
    tcl("GetModelTime", instance.handle)
  }
}

tcl.paired.to.list <- function(paired, as.enum.types) {
  length <- as.integer(tcl("llength", paired))
  if (length==1) {
    if (as.enum.types) {
      as.character(paired)
    } else {
      as.double(paired)
    }
  } else {
    result <- list() # sets none
    for (posn in seq(1,length,by=2)) {
      index <- tcl("lindex", paired, posn-1)
      if (as.enum.types) {
        index <- as.character(index)
      } else {
        index <- as.integer(index)
      }    
      result[[index]] <- tcl.paired.to.list(tcl("lindex", paired, posn),
                                            as.enum.types)
    }
    result
  }
}

tcl.paired.to.array <- function(paired, dims, as.enum.types) {
  # note indices in value from model are ignored, so may be enumerated type
  if (length(dims)) {
    result <- {}
    subDims <- dims[-1] # removes first element
    for (posn in 1:dims[1]) {
      idx <- 2*posn-1
      member <- tcl("lindex", paired, idx)
      result <- c(result, tcl.paired.to.array(member, subDims, as.enum.types))
    }
    array(result,dim=rev(dims))
  } else if (as.enum.types) {
    as.character(paired)
  } else {
    as.double(paired)
  }
}

get.value.list <- function(instance.handle, value.name, as.enum.types = FALSE) {
  if (is.dummy(instance.handle)) {
    paired <- "1 40.76667783660071 2 37.52906643918561 3 33.820213413335914 4 29.694134498874085";
  } else {
    paired <- tcl("GetPairedValues", instance.handle, value.name, as.enum.types)
  }
  tcl.paired.to.list(paired, as.enum.types)
}

get.value.array <- function(instance.handle, value.name, as.enum.types = FALSE) {
  if (is.dummy(instance.handle)) {
    dims <- 10
    paired <- "1 0.8414710 2 0.9092974 3 0.1411200 4 -0.7568025 5 -0.9589243 6 -0.2794155 7 0.6569866 8 0.9893582 9 0.4121185 10 -0.5440211"
  } else {
    i.m.list <- tcl("array", "get", "::modelTypes", instance.handle)
    dims <- get.model.property(tcl("lindex", i.m.list, 1), value.name, "Dims")
    if (any(is.na(dims))) {
      stop("This value is in a variable-membership submodel --
use get.value.list instead")
    }
    paired <- tcl("GetPairedValues", instance.handle, value.name, as.enum.types)
  }
  tcl.paired.to.array(paired, dims, as.enum.types)
}
