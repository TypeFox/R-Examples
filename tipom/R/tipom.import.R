tipom.import <-
  function(imported.data, name, units, description=NULL) {
    imd <- imported.data
    attr(imd, "tipom") <- TRUE
    attr(imd, "name") <- name
    attr(imd, "units") <- units
    attr(imd, "description") <- description
    imd
  }
