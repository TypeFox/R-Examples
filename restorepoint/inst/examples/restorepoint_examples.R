examples.test.fun = function() {
  library(restorepoint)
  init.restore.point() 
  
  f = function(x=1) {
    restore.point("f")
    new.x = sample(1:4,1)
    if (new.x>2) return(new.x)
    f(x)
  }
 
#  add.restore.point.var.test({
#    if (x<4) {
#      cat("x == ",x)
#    }
#  }) 
  
  add.restore.point.test(x.test = function(env,...) {
    if (exists("x",env)) {
      x = env$x
      if (x<4) {
        cat("\nx == ",env$x,"\n")
      }
    }
  })

  set.restore.point.options(display.restore.point=TRUE)

  f()  
}

examples.env.console = function() {

  # Generate an environment and start a console that can be used to evaluate 
  # that environment
  a = 5
  env = new.env()
  env$b = "Hi"
  env.console(env=env, dots = list("A"="ABC"))
  
  # Try typing the following into the new console
  1+
  1
  a
  b
  list(...)
  
  # Investigate the environment of the function
  f = function(a,b="B",...) {
    env.console(environment(), dots=list(...))
  }
  # Run from the standard R console. Nested calls to env.console don't work well
  f(a=5,c=10,d=20)
  # Try typing the following into the new console
  list(...)
} 
# 
# examples.break.point = function() {
#   library(restorepoint)
#   f = function(x=5) {
#     breakpoint()
#     x*2
#     f(x*3)
#   }
#   f()
# }

# 
# examples.clone.list = function() {
#   # Test if data.table inside is cloned correctly
#   library(data.table)
#   li0 = list(dt = data.table(col=1:2),name="NAME")
#   li1 = clone.list(li)
#   # Change value in clone
#   li1$dt[1,col] = 100
#   # Values should be different
#   li0$dt
#   li1$dt
#   
#   env = clone.environment(as.environment(li))
#   env$dt[1,col:="A"]      
#   env$dt
# } 
# 


examples.break.point = function () {  
  library(restorepoint)
  set.restore.point.options(break.point.to.global = FALSE)
  # A function that shall swap the left and right part of a vector
  swap.in.vector = function(vec,swap.ind) {
    break.point()
    left  = vec[1:(swap.ind-1)]
    right = vec[swap.ind:nrow(vec)]
    c(right,left)
  }
  swap.in.vector(1:10,4)
  
}  
  


examples.restore.point.browser = function () {  
  library(restorepoint)
  
  set.restore.point.options(to.global = FALSE)
  
  # A function that shall swap the left and right part of a vector
  swap.in.vector = function(vec,swap.ind) {
    restore.point("swap.in.vector")
    left  = vec[1:(swap.ind-1)]
    right = vec[swap.ind:nrow(vec)]
    c(right,left)
  }
  swap.in.vector(1:10,4)
  
  # You could call
  restore.point.browser("swap.in.vector")
  
  # But there is no need to. Just running in the R Console (global environment) the command
  restore.point("swap.in.vector")
  # has the same effect (at least once you have set restore.point.options(to.global=FALSE))
  
  # Basically restore.point.browser() is rather an internal function
}  
  

examples.restore.point = function () {  
  
  # See the vignette for a detailed description
  library(restorepoint)
  init.restore.point() 
  set.restore.point.options(to.global = TRUE, display.restore.point=TRUE, trace.calls=TRUE)
  # A function that shall swap the left and right part of a vector
  swap.in.vector = function(vec,swap.ind) {
    restore.point("swap.in.vector")
    left  = vec[1:(swap.ind-1)]
    right = vec[swap.ind:nrow(vec)]
    c(right,left)
  }
  
  swap.in.vector(1:10,4)
  
  # Paste the body of the function
  restore.point("swap.in.vector")
  left  = vec[1:(swap.ind-1)]
  right = vec[swap.ind:nrow(vec)]
  c(right,left)

  
  # Investigate the error line
  swap.ind:nrow(vec)
  nrow(vec)
  
  length(vec)
  
  # Correct the function
  swap.in.vector = function(vec,swap.ind) {
    restore.point("swap.in.vector")
    left  = vec[1:(swap.ind-1)]
    right = vec[swap.ind:length(vec)]
    c(right,left)
  }
  swap.in.vector(1:10,4)

  ###############################################################
  # Nested function calls
  ###############################################################
    
  f = function(vec=1:5) {
    restore.point("f")
    for (i in 1:20) {
      swap.point = sample(1:length(vec),1)
      sw = swap.in.vector(vec,swap.point)
      print(sw)
      if (length(sw)>length(vec)) 
        stop("Error output too long...")
    }
    return(NULL)
  }
  f()  
  # See the tutorial on how to correct this function
  

  ###############################################################
  # Dealing with ...
  ###############################################################
  g = function(...) {
    # restore.points can deal with function calls that use ...
    # only with the manual browser, i.e. one needs to set
    # to.global = FALSE
    restore.point("g", to.global = FALSE) 
    print(paste(...))
    list(...)
  }
  g(a="Hello", b="World!")
  
  ######################################################
  # Restoring also enclosing environment of a function
  #####################################################
  
  library(restorepoint)
  make.fun = function(x="make.fun.x",y="make.fun.y") {
    f = function(y="local.y") {
      restore.point("f") 
      c(x,y,z)
    }
    return(f)
  }
  g = make.fun()
  x = "glob.x"
  z = "glob.z"
  g()
  x = "changed.x"
  z = "changed.z"
  # x and y will be restored, but no variables from .GlobalEnv
  restore.point("f") 
  c(x,y,z)
  
  
  ###############################################################
  # Check when environments are cloned correctly
  ###############################################################
  
  env <- new.env()
  parent.env(env)
  env$x = 1:3
  li <- list(env=env,test="test")
  li$env$x
  
  g = function(env,li) {
    restore.point("g", deep.copy=TRUE)
    as.list(env)
    parent.env(env)
    env$x = c("A","B")
    li$env$x
  }
  g(env,li)
  # Environment has been changed
  env$x
  
  
  # Check if environments are correctly restored
  restore.point("g")
  env$x
  li$env$x
  
  # Are both environments still references to the same object?
  env$x = "new x"
  li$env$x
  
  
  
#   # data.tables are also copied by reference
#   library(data.table)
#   dt <- data.table(col=1:2)
#   init.restore.point()
#   g = function(dt) {
#     restore.point("g")
#     dt[,col:=c("A","B")]
#   }
#   g(dt)
#   dt  # dt has been changed
#   
#   restore.point("g")
#   dt # the original value has been restored
} 

examples.eval.with.error.trace = function() {
  # A function that has an error
  f = function(x) {
    1:3[0]
  }
  g = function(x=4) {
    f(x)
  }
  # Usually no traceback is available for errors that are caught with tryCatch
  h = function() {
    tryCatch(
      g(25),
      error = identity
    )
  }
  h()
  
  # The function eval.with.error.trace adds trace information if an error is thrown
  h = function() {
    tryCatch(
      eval.with.error.trace(g(25)),
      error = identity
    )
  }
  h()
}

test.dots = function() {
  dots = list(b=4)
  f = function(a,b) {
    a*b
  }
  
  as.character(as.call((parse(text="f (a=3,b=3)")))[[1]][[1]])
  as.character(as.call((parse(text="x=5"))))
  
  tryCatch(
    f(a=5,...),
    error = function(e){
      print(names(e))
      print(str(e))
      print(str(e$call))
    }
  ) 
}

examples.calls.to.trace = function() {
  f = function(a=5) calls.to.trace(rev(sys.calls()))
  g = function(x=2) f(x)
  g(10)
}
