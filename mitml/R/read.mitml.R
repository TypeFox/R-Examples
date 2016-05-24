read.mitml <- function(filename){
# read mitml objects from file

  env <- new.env(parent=parent.frame())
  load(filename, env)
  obj <- ls(env)
  eval(parse(text=obj), env)

}

