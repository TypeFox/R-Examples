# library(devtools);install(pkg = "C:/libraries/restorepoint/restorepoint", reload = TRUE, quick = FALSE,args = NULL, quiet = FALSE)

# Debug with the R console by setting restore points.

rpglob <- new.env()
.onLoad <- function(libname, pkgname) {
  init.restore.point()  
}

init.restore.point = function() {
  rpglob$options = list(storing=TRUE,to.global = TRUE,multi.line.parse.error = get.multi.line.parse.error(), deep.copy=FALSE, break.point.to.global=FALSE, display.restore.point=FALSE, trace.calls=TRUE)
  rpglob$OBJECTS.LIST <- list()
  rpglob$CALLS.LIST <- list()
  rpglob$TESTS.LIST <- list()
}

clear.restore.point.tests = function() {
  rpglob$TESTS.LIST <- list()
}

#' Add one or several test functions
#'
#' A test function is called after a restore point
#' has stored data. It must have an argument env and name.
#' It can check whether certain conditions are satisfied by the variables
#' @param ... a slist of test functions that will be called with the stored arguments
#' @export  
add.restore.point.test = function(...) {
  tests = list(...)
  rpglob$TESTS.LIST[names(tests)] = tests
}



#' Set global options for restore points
#' 
#' @param options a list of options that shall be set. Possible options are listed below
#' @param ... options can also directly be passed. The following options can be set:
#' - storing Default=TRUE enable or disable storing of options, setting storing = FALSE basicially turns off debugging via restore points
#' - deep.copy Default = FALSE. If TRUE then when storing and restoring tries to make a deep copy of R objects that are by default copied by reference, like environments. deep.copy = FALSE substantially speeds up restore.point.
#' - to.global Default=TRUE. If  TRUE then when options are restored, they are simply copied into the global environment and the R console is directly used for debugging. If FALSE a browser mode will be started instead. It is still possible to parse all R commands into the browser and to use copy and paste. To quit the browser press ESC in the R console. The advantage of the browser is that all objects are stored in a newly generated environment that mimics the environemnt of the original function, i.e. global varariables are not overwritten. Furthermore in the browser mode, one can pass the ... object to other functions, while this does not work in the global environment. The drawback is that the browser is still not as convenient as the normal R console, e.g. pressing arrow up does not restore the previous command. Also, one has to press Esc to leave the browser mode.
#' @export 
restore.point.options = set.restore.point.options = function(options=NULL,...) {
  options = c(options,list(...))
  unknown.options = setdiff(names(options),names(get.restore.point.options())) 
  if (length(unknown.options)>0) {
    warning(paste("unknown options", paste(unknown.options, collapse=","),"ignored"))
    options = options[setdiff(names(options),unknown.options)]
  }
  if (!exists("rpglob"))
    init.restore.point()  
  if (length(options)>0)
    rpglob$options[names(options)] = options
  invisible(rpglob$options)
}

#' Get global options for restore points
#' 
#' @export
get.restore.point.options = function() {
  if (!exists("rpglob"))
    init.restore.point()  
  rpglob$options
  
} 
#get.restore.point.options()
#' Retrieves the list of all restore.points with the stored objects
#' 
#' @export
get.stored.object.list = function() {
  rpglob$OBJECTS.LIST
}

#' Set whether objects shall be stored or not
#' 
#' @param storing if FALSE don't store objects if restore.point or store.objects is called. May save time. If TRUE (default) turn on storage again.
#' @export
set.storing <- function(storing=TRUE) {
  set.restore.point.options(storing=storing)
}


#' Check whether objects currently are stored or not
#' 
#' @export
is.storing <- function() {
  return(get.restore.point.options()$storing)
}




#' Sets a restore point
#'
#' The function behaves different when called from a function or when called from the global environemnt. When called from a function, it makes a backup copy of all local objects and stores them internally under a key specified by name. When called from the global environment, it restores the previously stored objects by copying them into the global environment. See the package Vignette for an illustration of how this function can facilitate debugging.
#'
#' @param name key under which the objects are stored. For restore points at the beginning of a function, I would suggest the name of that function.
#' @param to.global if TRUE (default) objects are restored by simply copying them into the global environment. If FALSE a new environment will be created and the restore point browser will be invoked. 
#' @param deep.copy if TRUE try to make deep copies of  objects that are by default copied by reference. Works so far for environments (recursivly). The function will search lists whether they contain reference objects, but for reasons of speed not yet in other containers. E.g. if an evironment is stored in a data.frame, only a shallow copy will be made. Setting deep.copy = FALSE (DEFAULT) may be useful if storing takes very long and variables that are copied by reference are not used or not modified.
#' @param force store even if set.storing(FALSE) has been called
#' @param dots by default a list of the ... argument of the function in whicht restore.point was called
#' @param display.restore.point shall a text be shown in the console if restore.point is called. Can be useful when informative tracebacks are not readily availbale, e.g. when debugging shiny apps.
#' @param indent.level when display.restore.point=TRUE shall level of nestedness be illustrated by identation
#' @param trace.calls when objects are restored, shall a traceback be shown
#' @param max.trace.lines if trace.calls=TRUE how many lines shall be shown at most in the traceback.
#' @export
restore.point = function(name,to.global = get.restore.point.options()$to.global,deep.copy = get.restore.point.options()$deep.copy, force=FALSE,display.restore.point = get.restore.point.options()$display.restore.point, indent.level = TRUE, trace.calls = get.restore.point.options()$trace.calls,max.trace.lines=10, dots = eval(substitute(list(...), env = parent.frame()))) {

  envir = sys.frame(-1)
  if (isTRUE(display.restore.point)) {
    if (!indent.level) {
      cat("\nrestore.point: ", name)
    } else {
      level = max(length(sys.calls())-2,0)
      if (level>9) {
        indent = paste0(level," .......")
      } else {
        indent = paste0(rep(".", level),collapse="")
      }
      cat(paste0("\n",indent,"restore.point: ", name))
    }
  }
  
  # restore objects if called from the global environment
  # when called from a function store objects
  restore = identical(.GlobalEnv,envir)  
  if (restore) {
    if (trace.calls) {
      calls = rpglob$CALLS.LIST[[name]]
      if (length(calls)>0) {
        tr = calls.to.trace(calls,max.lines = max.trace.lines)
        cat(paste0("\n",paste0(tr,collapse="\n")))
      }
    }
    
    
    if (!to.global) {
      restore.point.browser(name,was.forced=force, deep.copy=deep.copy)
    } else {
      restore.objects(name=name,was.forced=force,deep.copy=deep.copy)
    }
  } else {
    if (trace.calls) {
      rpglob$CALLS.LIST[[name]] <- sys.calls()  
    }
    store.objects(name=name,parent.num=-2, deep.copy=deep.copy, force=force,dots=dots)
    
    # call test funs
    for (test.fun in rpglob$TESTS.LIST) {
      test.fun(env = rpglob$OBJECTS.LIST[[name]],name=name)  
    }
  }
}


#' Sets a break point that can be debugged like a restore point
#'
#'  This function can be used as an alternative to browser(). When called inside a function, break.point stores all local objects and then does the following. i) If to global=FALSE (the default for break.point) starts the restore.point.browser for the local objects. ii) if to.global=TRUE copies the local objects to the global environment and stops execution.
#' 
#' An alternative to break points are restore points. In the tutorial on GitHub, I provide some arguments how restore points can facilitate debugging compared to break points.
#'
#' @param name key under which the objects are stored. For restore points at the beginning of a function, I would suggest the name of that function.
#' @param to.global if TRUE (default) objects are restored by simply copying them into the global environment. If FALSE a new environment will be created and the restore point browser will be invoked. 
#' @param deep.copy if TRUE try to make deep copies of  objects that are by default copied by reference. Works so far for environments (recursivly). The function will search lists whether they contain reference objects, but for reasons of speed not yet in other containers. E.g. if an evironment is stored in a data.frame, only a shallow copy will be made. Setting deep.copy = FALSE (DEFAULT) may be useful if storing takes very long and variables that are copied by reference are not used or not modified.
#' @param force store even if set.storing(FALSE) has been called
#' @param dots by default a list of the ... argument of the function in whicht restore.point was called
#' @export

break.point = function(name="BREAK_POINT___",to.global = get.restore.point.options()$break.point.to.global,deep.copy = get.restore.point.options()$deep.copy, force=FALSE, dots = eval(substitute(list(...), env = parent.frame()))) {

  envir = sys.frame(-1);
  
  # break function when called from a function, otherwise ignore
  from.global = identical(.GlobalEnv,envir)  
  # store and restore
  if (!from.global) {
    store.objects(name=name,parent.num=-2, deep.copy=deep.copy, force=force,dots=dots)
    if (!to.global) {
      str = paste0("Break point encountered in function ", as.character(sys.call(-1)[[1]]),". Start restore point browser.")
      message(str)
      restore.point.browser(name,was.forced=force, deep.copy=deep.copy, message.text = NULL)
    } else {
      str = paste0("Break point encountered in function ", as.character(sys.call(-1)[[1]]),". Write objects to global environment")
      restore.objects(name=name,was.forced=force,deep.copy=deep.copy)
      stop(str)
    }    
  }
}


parent.env.to.store = function(penv) {
  if(identical(penv,globalenv()) |
     identical(penv,baseenv())   |
     identical(penv,emptyenv())) {
    return(FALSE)
  }
  if (nchar(environmentName(penv))>0)
    return(FALSE)
  return(TRUE)
}

#' Stores all local objects of the calling environment to be able to restore them later when debugging. Is used by restore.point 
#' 
#' @param name key under which the objects are stored, typical the name of the calling function. If name is NULL by default the name of the calling function is chosen
#' @param deep.copy if TRUE (default) variables that are copied by reference (in the moment environments)  will be stored as deep copy. May take long for large variables but ensures that the value of the stored variable do not change
#' @param force store even if do.store(FALSE) has been called
#' @param store.if.called.from.global if the function is called from the global environment and store.if.called.from.global FALSE (default) does not store objects when called from the global environment but does nothing instead.
#' @param envir the environment from which objects shall be stored. By default the local environemnt of the calling function.
#' @param parent.num can be used to specify envir=sys.frame(parent.num)
#' @param store.parent.env shall objects from enclosing environments of envir also be stored? So far this happens for all enclosing environments except for the global environment or baseenv.
#' @param dots by default a list of the ... argument of the function in whicht restore.point was called
#' @return returns nothing, just called for side effects
#' @export
store.objects = function(name=NULL,parent.num=-1,deep.copy = get.restore.point.options()$deep.copy, force=FALSE, store.if.called.from.global = FALSE, envir = sys.frame(parent.num),store.parent.env = "all.but.global", dots = eval(substitute(list(...), env = parent.frame()))
) {
   
  if (!(is.storing()) & !force) {
    return(NULL)
  }
  
	if (sys.nframe() < 2 & !store.if.called.from.global) {
		warning(paste("store.objects(\"",name,"\") ignored since called from global environment."),sep="")
		return()
	}

  #envir = sys.frame(parent.num);
  
  # Assign name of the calling function
  fun.name = all.names(sys.call(parent.num))[1]
  if (is.null(name)) {
  	name = fun.name
  }
  if (force) {
    #warning(paste("store.objects called by ", fun.name, " with force!"))
  }
  if (deep.copy) {
    rpglob$copied.ref = NULL
    copy.fun = function(envir)
      clone.environment(envir, use.copied.ref = TRUE)
  } else {
    copy.fun = function(envir)
      as.environment(as.list(envir, all.names=TRUE))
  }
  
  copied.env = copy.fun(envir)
  # Generate copies of parent enclosing environments
  if (store.parent.env == "all.but.global") {
    cenv = copied.env
    penv = parent.env(envir)
    while (parent.env.to.store(penv)) {
      copied.penv = copy.fun(penv)
      parent.env(cenv) <- copied.penv
      cenv = copied.penv
      penv = parent.env(penv)
    }
  }
  
  ev.dots <- NULL
  try(ev.dots <- force(dots), silent=TRUE)
  dots <- ev.dots
  if (!is.null(dots) & deep.copy) {
    dots = clone.list(dots)
  }
  attr(copied.env,"dots") <- dots
  
  rpglob$OBJECTS.LIST[[name]] <- copied.env
  return()  
}


#' Restore stored objects by copying them into the specified environment. Is used by restore.point
#' 
#' @param name name under which the variables have been stored
#' @param dest environment into which the stored variables shall be copied. By default the global environment.
#' @param was.forced flag whether storage of objects was forced. If FALSE (default) a warning is shown if restore.objects is called and is.storing()==FALSE, since probably no objects have been stored.
#' @param deep.copy when storing or restoring tries to make a deep copy of R objects that are by default copied by reference, like environments. Setting deep.copy = FALSE can substantially speed up restore.point, however.
#' @return returns nothing but automatically copies the stored variables into the global environment
#' @export
restore.objects = function(name, dest=globalenv(), was.forced=FALSE, deep.copy=get.restore.point.options()$deep.copy) {
  if ((!is.storing()) & (!was.forced)) 
    warning("is.storing() == FALSE\nPossible objects were not correctly stored. Call set.storing(TRUE) to enable storing.")
  
  env =   rpglob$OBJECTS.LIST[[name]]
  
  if (is.null(env)) {
    stop(paste0("No objects stored under name ", name))
  }
  
  if (!deep.copy) {
    # Reference objects will just be taken in their actual state
    copy.into.env(source=env,dest=dest,from.restore.objects=TRUE)
    restored = ls(envir=env, all.names=TRUE)
    # Copy all objects from stored parent environments into dest
    penv = parent.env(env)
    penv.list = list()
    count = 1
    while (parent.env.to.store(penv)) {
      copy.into.env(source=penv,dest=dest,from.restore.objects=TRUE, exclude=restored)
      restored = c(restored,ls(envir=penv, all.names=TRUE))
      penv = parent.env(penv)    
    }
    restored = unique(restored)

  } else if (deep.copy) {
    # Clone stored environment in order to guarantee that the restore point can be used several times even if reference objects are used
    rpglob$copied.ref = NULL
    cenv = clone.environment(env,use.copied.ref = TRUE)
    # Copy the stored objects into the enviornment specified by dest (usually the global environment)
    
    # Restore copies of parent enclosing environments
    penv = parent.env(env)
    penv.list = list()
    count = 1
    while (! (identical(penv,globalenv())|
                identical(penv,baseenv())  |
                identical(penv,emptyenv())  )) {
      penv.list[[count]] = clone.environment(penv,use.copied.ref = TRUE)
      penv = parent.env(penv)    
    }
    restored = NULL
    # Simply copy variables of none-global parent environments into dest
    # The hierachy of enclosing environments is not replicated!!!
    for (i in rev(seq_along(penv.list))) {
      copy.into.env(source=penv.list[[i]],dest=dest,from.restore.objects=TRUE)
      restored = c(restored,ls(envir=penv.list[[i]], all.names=TRUE))
    }  
    copy.into.env(source=cenv,dest=dest,from.restore.objects=TRUE)
    restored = c(restored,ls(envir=cenv, all.names=TRUE))
  }
  message(paste("Restored: ", paste(restored,collapse=",")))
}


clone.list = function(li, use.copied.ref = FALSE) {          
  ret.li = lapply(li,copy.object,use.copied.ref = use.copied.ref)
  return(ret.li)
}

#' Deep copy of an environment
#' @param env the environment to be cloned
#' @param use.copied.ref internal 
#' @param all.names passed to eapply
#' @export
clone.environment = function(env, use.copied.ref = FALSE, all.names=TRUE) {
  #print(as.list(env))
  #browser()
  li = eapply(env,copy.object,use.copied.ref = use.copied.ref, all.names=all.names)
  cloned.env = as.environment(li)
  # Set same enclosing environment as env
  parent.env(cloned.env) <- parent.env(env)
  cloned.env
}


copy.object = function(obj, use.copied.ref = FALSE) {
  #print("copy.object")
  #print(paste("missing: ",missing(obj), "class(obj) ", class(obj)))
  #browser()        
  # Dealing with missing values
  if (is.name(obj)) {
    return(obj)
  }
  oclass =class(obj) 
  
  # If the objects has already been copied, just return the reference of the copied version, don't create an additional copy
  
  if (any(oclass %in% c("environment"))) {
    if (use.copied.ref & !is.null(rpglob$copied.ref)) {
      ind = which(sapply(rpglob$copied.ref[,1],identical,y=obj))
      if (length(ind)>0) {
        return(rpglob$copied.ref[ind[1],2][[1]])
      }
    }
    if ("environment" %in% oclass) {
      copy = clone.environment(obj,use.copied.ref = use.copied.ref)
      
    #} else if ("data.table" %in% oclass) {
    #  copy = data.table::copy(obj)
    }    
    # Store a copy of the reference
    rpglob$copied.ref = rbind(rpglob$copied.ref,c(obj,copy))
  } else {
    if (is.list(obj) & !(is.data.frame(obj))) {
      copy = clone.list(obj,use.copied.ref = use.copied.ref)
    } else {
      copy = obj
    }
  }
  return(copy)
}

#' Checks whether for the installed R version the function env.console is able to correctly parse R expressions that extend over more than a line
#' 
#'  The current implementation of env.console is quite dirty in so far that it parses an error message of the parse() function to check whether a given R expression is assumed to be continued in the next line. That process may not work in R distributions that have error messages that are not in English. The function can.parse.multi.line() tries to check whether that process works or not
#'  @export
can.parse.multi.line = function() {
  is.multi.line("1+")
}

# Tries to get the error message for a multilin.parse.error
get.multi.line.parse.error = function() {
  err= tryCatch(
    parse(text="1+"),
    error = function(e) {
      as.character(e)
    }
  )
  err = strsplit(err,"\n")[[1]][1]
  ret = strsplit(err,":")[[1]]
  ret[length(ret)]
}

# Check whether the supplied code is a multiline input
is.multi.line = function(code, multi.line.parse.error = get.restore.point.options()$multi.line.parse.error) {
  ret = FALSE
  tryCatch(
    parse(text=code),
    error = function(e) {
      str = as.character(e)
      #print(str)
      if (length(grep(multi.line.parse.error,str,fixed=TRUE))>0)
        ret<<-TRUE
      
    }
  )
  return(ret)  
}

#' Examing a previously stored restore point by invoking the browser.
#' 
#' The function is mainly for internal use by restore.point.
#' 
#' @param name name under which the variables have been stored
#' @param was.forced flag whether storage of objects was forced. If FALSE (default) a warning is shown if restore.objects is called and is.storing()==FALSE, since probably no objects have been stored.
#' @param message.text initial shown message
#' @param deep.copy when storing or restoring tries to make a deep copy of R objects that are by default copied by reference, like environments. Setting deep.copy = FALSE can substantially speed up restore.point, however.
#' @return returns nothing
#' @export
restore.point.browser = function(name,was.forced=FALSE, message.text=paste("restore point",name, ", press ESC to return."), deep.copy=get.restore.point.options()$deep.copy) {
  if (!is.null(message.text))
    message(message.text)

  # Generate environment in which the console shall be called
  enclos.env=.GlobalEnv # may store an enclosing environment instead
  env <- new.env(parent=enclos.env)
  # Populate environment with stored variables
  restore.objects(name,dest=env,was.forced=was.forced, deep.copy=deep.copy)
  # Get ... from original function
  dots = get.stored.dots(name, deep.copy=deep.copy)

  local.variables=as.list(env)
  env.console(env=env,dots=dots, startup.message=NULL, local.variables=local.variables)
}


#' Emulates an R console that evaluates expressions in the specified environement env. You return to the standard R console by pressing ESC
#' 
#' @param env The environment in which expressions shall be evaluated. If not specified then a new environment with the given parent.env is created.
#' @param parent.env If env is not specified the parent environemnt in which the new environment shall be created
#' @param dots a list that contains values for the ellipsies ... that will be used if you call other functions like fun(...) from within the console. You can access the values inside the console by typing list(...)
#' @param prompt The prompt that shall be shown in the emulated console. Default = ": "
#' @param startup.message The text that is shown when env.console is started
#' @param multi.line.parse.error A substring used to identify an error by parse that is due to parsing the beginning of a multi-line expression. The substring can depend on the language of R error messages. The packages tries to find a correct substring automatically as default.
#' @param local.variables additional variables that shall be locally available 
#' @return Returns nothing since the function must be stopped by pressing ESC.

#' @export
env.console = function(env = new.env(parent=parent.env), parent.env = parent.frame(), dots=NULL,prompt=": ", startup.message = "Press ESC to return to standard R console", multi.line.parse.error = get.restore.point.options()$multi.line.parse.error, local.variables = NULL) {
  
  
  parse.fun <- function(..., .LOCAL.VARIABLES = NULL) {
    .IS.ENV.CONSOLE.ENVIRONMENT__ <- TRUE
    
    .CONSOLE.INTERNAL$prev.code = ""
    .CONSOLE.INTERNAL$prompt = .CONSOLE.INTERNAL$normal.prompt
  
    if (!is.null(.LOCAL.VARIABLES)) {
      copy.into.env(source=.LOCAL.VARIABLES)
    }
    #message("multi.line.parse.error:")
    #print(.CONSOLE.INTERNAL$multi.line.parse.error)
    while(TRUE) {
      # Read 1 line of code
      .CONSOLE.INTERNAL$code = paste(.CONSOLE.INTERNAL$prev.code,readline(prompt=.CONSOLE.INTERNAL$prompt),sep="")

      # Skip empty lines and comments
      .CONSOLE.INTERNAL$stripped.code = gsub(" ","",.CONSOLE.INTERNAL$code,fixed=TRUE)
      if (nchar(.CONSOLE.INTERNAL$stripped.code)==0)
        next
      if (substring((.CONSOLE.INTERNAL$stripped.code),1,1)=="#")
        next
      
      .CONSOLE.INTERNAL$multi.line = FALSE
      
      # Try to parse the code 
      .CONSOLE.INTERNAL$expr = tryCatch(
        parse(text=.CONSOLE.INTERNAL$code),
        error = function(e) {
          str = as.character(e)
          #print(str)
          if (length(grep(.CONSOLE.INTERNAL$multi.line.parse.error,str,fixed=TRUE))>0) {
            .CONSOLE.INTERNAL$multi.line<<-TRUE
          } else {
            if (length(str)>0)
              message(str)
          }
        }
      )
      
      # A multiline expression need to continue parsing before evaluating
      if (.CONSOLE.INTERNAL$multi.line) {
        .CONSOLE.INTERNAL$prev.code = paste(.CONSOLE.INTERNAL$code, "\n")
        .CONSOLE.INTERNAL$prompt = "+ "
        next
      } else {
        .CONSOLE.INTERNAL$prompt = .CONSOLE.INTERNAL$normal.prompt
        .CONSOLE.INTERNAL$prev.code = ""
      }
      
      # Check if input shall be stopped
      .CONSOLE.INTERNAL$fun.name <- tryCatch({
        .CONSOLE.INTERNAL$fun.name <- as.character(as.call((.CONSOLE.INTERNAL$expr))[[1]][[1]])
        },
        error = function(e) {
          ""
        }
      ) 
      # Check if a function like restore.point or source was called that stops the
      # console
      if (.CONSOLE.INTERNAL$fun.name %in% .CONSOLE.INTERNAL$stop.functions) {
        return(list(fun.name=.CONSOLE.INTERNAL$fun.name, expr=.CONSOLE.INTERNAL$expr))
      }
      # Try to evaluate the expression
      tryCatch(
        eval.with.error.trace({
          if (!is.null(.CONSOLE.INTERNAL$expr)) {
            .CONSOLE.INTERNAL$expr.out <- capture.output(eval(.CONSOLE.INTERNAL$expr))
            if (length(.CONSOLE.INTERNAL$expr.out)>0) {
              cat(paste(.CONSOLE.INTERNAL$expr.out,collapse="\n"),"\n")
            }
          }
        }, remove.early.calls = 7),
        error = function(e) {
          str = conditionMessage(e)
          if (substring(str,1,35)=="Error in eval(expr, envir, enclos):")
            str = paste("Error:",substring(str,36),sep="")
          if (length(str)>0) {
            message(str, appendLF=FALSE)
            #cat("Traceback:\n")
            #cat(paste(e$calls,collapse="\n"))
          }
        }
      )
    }
    return(NULL)
  }
  
  dots = c(dots)
  env$.CONSOLE.INTERNAL = list(prompt=prompt, normal.prompt=prompt, stop.functions = c("restore.point","env.console","source"), multi.line.parse.error = multi.line.parse.error)

  environment(parse.fun) <- env
  
  if (!is.null(startup.message))
    message(startup.message)
  if (!is.null(dots)) {
   # ret = do.call(parse.fun,dots)
    ret = do.call(parse.fun,c(dots,list(.LOCAL.VARIABLES=local.variables))) 
  } else {
    ret = parse.fun(.LOCAL.VARIABLES=local.variables)
  }

  # Call the expression that should be called in the global environment
  if(!is.null(ret)) {
    if (ret$fun.name=="source") {
      message("Stop debugger because new file is sourced")
    }
    eval(ret$expr, envir=.GlobalEnv)
  }
}

#' A default error string function for eval with error trace
#' 
#' @param e the error object
#' @param tb a character vector of the traceback
default.error.string.fun = function(e,tb) {
  if (length(tb)>0) {
    paste0(as.character(e),"\nCall sequence:\n", paste(tb,collapse = "\n"),"\n")
  } else {
    paste0(as.character(e),"\n")
  }
}

#' Evals the expression such that if an error is encountered a traceback is added to the error message.
#' 
#' This function is mostly useful within a tryCatch clause 
#' Adapted from code in tools:::.try_quietly
#' as suggested by Kurt Hornik in the following message
#' https://stat.ethz.ch/pipermail/r-devel/2005-September/034546.html
#' 
#' @param expr the expression to be evaluated
#' @param max.lines as in traceback()
#' @param remove.early.calls an integer specifying a number of calls that won't be shown in the trace.
#' @param error.string.fun a function(e,tb) that takes as arguments an error e and a string vector tb of the stack trace resulting from a call to calls.to.trace() and returns a string with the extended error message
#' @return If no error occurs the value of expr, otherwise an error is thrown with an error message that contains the stack trace of the error.
#' @export
eval.with.error.trace = function(expr, max.lines=4, remove.early.calls =  0, error.string.fun = default.error.string.fun) {
  withRestarts(
    withCallingHandlers(
      eval(expr),
      error = {
        function(e) invokeRestart("grmbl", e, sys.calls())
      }
    ),
    grmbl = function(e, calls) {
      n <- length(sys.calls())
      calls <- calls[-seq.int(length.out = n - 1L)]
      remove = pmax(2,length(calls)-remove.early.calls):length(calls)
      calls <- rev(calls)[-c(1L,2L, remove)]
      if (length(calls)>0) {
        tb <- calls.to.trace(calls)
        stop(error.string.fun(e,tb), call. = FALSE)
      } else {
        stop(error.string.fun(e,tb=NULL), call. = FALSE)
      }
    }
  )    
}

#' Transforms a list returned by sys.calls into a vector of strings that looks like a result of traceback()
#' @param calls a list of calls, e.g. returned by sys.calls
#' @param max.lines as in traceback()
#' @return a character vector with one element for each call formated in a similar fashion as traceback() does
#' @export
calls.to.trace = function(calls=sys.calls(),max.lines=4) {
  x <- lapply(calls, deparse)
  n <- length(x)
  if (n == 0L) 
    str = paste("No traceback", "\n")
  else {
    str = character(n)
    for (i in 1L:n) {
      label <- paste0(n - i + 1L, ": ")
      m <- length(x[[i]])
      if (!is.null(srcref <- attr(x[[i]], "srcref"))) {
        srcfile <- attr(srcref, "srcfile")
        x[[i]][m] <- paste0(x[[i]][m], " at ", basename(srcfile$filename), 
                            "#", srcref[1L])
      }
      if (m > 1) 
        label <- c(label, rep(substr("          ", 1L, 
                                     nchar(label, type = "w")), m - 1L))
      
      if (is.numeric(max.lines) && max.lines > 0L && max.lines < 
            m) {
        str[i] = paste0(paste0(paste0(label[1L:max.lines], x[[i]][1L:max.lines]), 
                               collapse = "\n"),
                        paste("\n",label[max.lines + 1L], " ...\n"))
      }
      else str[i] = paste0(paste0(label, x[[i]]), collapse = "\n")
    }
  }    
  str
}  


#' Returns the ellipsis (...) that has been stored in restore.point name as a list
#' 
#' @param name the name whith which restore.point or store.objects has been called.
#' @param deep.copy shall a deep copy of stored objects be made
#' @export
get.stored.dots = function(name,deep.copy=FALSE) {
  env  = rpglob$OBJECTS.LIST[[name]]
  dots = attributes(env)$dots
  if (length(dots)>0 & deep.copy)
    dots = clone.list(dots)
  dots
}

#check.global.vars()

#' Copies all members of a list or environment into an environment
#' 
#' @param source a list or environment from which objects are copied
#' @param dest the environment into which objects are copied
#' @param names optionally a vector of names that shall be copied. If null all objects are copied
#' @param exclude optionally a vector of names that shall not be copied
#' @param from.restore.objects internal paramater keep FALSE
#' @param overwrite should existing objects in dest with same name be overwritten?
#' @param all.names if TRUE copy all objects if names=NULL, if FALSE omit variables starting with .
#' @export
copy.into.env = function(source=sys.frame(sys.parent(1)),dest=sys.frame(sys.parent(1)),names = NULL, exclude=NULL, from.restore.objects=FALSE, overwrite = TRUE, all.names=TRUE) {

  if (is.null(names)) {
    if (is.environment(source)) {
      names = ls(envir=source, all.names=all.names)
    } else {
      names = names(source)
    }
  }
  if (!overwrite) {
    exclude = c(exclude, ls(envir=dest))
  }
  
  names = setdiff(names,exclude)
  
  if (is.environment(source)) {
    for (na in names) {
      if (!from.restore.objects) {
        assign(na,get(na,envir=source), envir=dest)
      } else {
        tryCatch (
          assign(na,get(na,envir=source), envir=dest),
          error = function(e) {
            message(paste("Variable ", na, " was missing."))
          }
        )
      }
    }
  } else if (is.list(source)) {
    for (na in names) {
      assign(na,source[[na]], envir=dest)
    }
  }
}

#' Checks whether cond holds true if not throws an error
#' 
#' Can be used for checking for errors in functions
#' @param cond a condition that is checked
#' @export
assert = function(cond) {
  if (!all(cond)) {
    label=as.character(match.call()[2])
    calls = sys.calls()
    if (length(calls)>1) {
      trace = paste0(" in " ,as.character(calls[[length(calls)-1]])[1])
    } else {
      trace = ""
    }
    #restore.point("assert")
    stop(paste0("The assertion '",label,"' failed", trace,"."),call.=FALSE)
  }
}
