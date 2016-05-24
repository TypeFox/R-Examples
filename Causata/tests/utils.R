# Testing utilities
# 
# Author David Barker
###############################################################################
library(stringr)

has.causata.connection.function <- function(){
  result <- try(source('/tmp/causata_environment.R'), silent=TRUE)
  if(class(result) == "try-error"){
    return(FALSE)
  } else {
    return(TRUE)
  }
}

local.connection <- function() {
  Connect(hostname="127.0.0.1", username="superuser", password="changeme")
}

close.connection <- function(conn) {
  if (!is.null(conn)) {
    Close(conn)
  }
}

with.local.connection <- function(function.with.connection) {
  conn <- NULL
  tryCatch({
    conn <- local.connection()
    print(paste("Connected to localhost"))
    function.with.connection(conn)
  }, finally={
    print("cleaning up from with.local.connection(): ")
    close.connection(conn)
    print("Disconnected from localhost")
  })
}

create.primary.variables <- function(conn, config, variables) {
  for (name in names(variables)) {
    result <- Config.CreatePrimaryVariable(config, variable.name=name, variable.expression=variables[[name]])
    stopifnot(result)
    print(paste("Created variable", name))
  }
}

with.primary.variables <- function(conn, config, variables, code) {
  variables.to.delete <- c()
  tryCatch({
    for (name in names(variables)) {
      result <- Config.CreatePrimaryVariable(config, variable.name=name, variable.expression=variables[[name]])
      stopifnot(result)
      print(paste("Created variable", name))
      variables.to.delete <- c(variables.to.delete, name)
    }
    
    # Wait (via loop and sleep) until all the variables appear in GetVariables(query)
    for (name in names(variables)) {
      wait.for.variable(conn, name)
    }
    
    # yield control to the caller.
    code
  }, finally={
    for (variable in variables.to.delete) {
      if (Config.DeleteVariable(config, variable.name=variable)) {
        print(paste("Deleted variable", variable))
      } else {
        print(paste("FAILED to delete variable", variable))
      }
    }
  })
}

delete.variable.at.end <- function(config, variable.name, code) {
  tryCatch({
    code
  }, finally={
    print(paste("cleaning up from delete.variable.at.end(", variable.name, ")"))
    if (Config.DeleteVariable(config, variable.name=variable.name)) {
      print(paste("Deleted variable", variable.name))
    } else {
      print(paste("FAILED to delete variable", variable.name))
    }
  })
}

wait.for.variable <- function(conn, variable.name, timeout.ms=20000) {
  start <- as.integer(1000*as.numeric(format(Sys.time(),"%H%M%OS3")))
  print(paste("Waiting for variable", variable.name))
  while (!variable.exists(conn, variable.name) && as.integer(1000*as.numeric(format(Sys.time(),"%H%M%OS3"))) < (start + timeout.ms)) {
    Sys.sleep(0.1)
  }
  if (!variable.exists(conn, variable.name)) {
    actual <- paste(GetMetadata(conn)[,"system.name"], collapse=",")
    print(paste("FAILURE. Expecting to see variable", variable.name, "only see", actual, collapse=" "))
    Sys.sleep(3600)
  }
  stopifnot(variable.exists(conn, variable.name))
  return(TRUE)
}

variable.exists <- function(conn, variable.name) {
  return(variable.name %in% GetMetadata(conn)$system.name)
}

# Normalizes whitespace in a string to make it one line, with no double, leading or trailing spaces.
# 1) replace newlines with spaces,
# 2) then replace 2 or more consecutive spaces with 1 space,
# 3) then remove leading and trailing space
strip <- function(s) {
  s <- str_replace_all(s, "\\n", " ")
  s <- str_replace_all(s, " +", " ")
  s <- str_replace_all(s, "^ +", "")
  s <- str_replace_all(s, " +$", "")
  return(s)
}

