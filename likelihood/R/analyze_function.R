analyze_function<-function(fun, par, parname, datasets, ...) {

  ##########################################
  ##
  ## Analyzes a function and extracts its arguments
  ## from a parameter list to make a function call.  
  ## This will recursively create a function call for 
  ## every function argument that is itself a function.
  ##
  ## Arguments:
  ## fun = function to analyze and find arguments for
  ## par = list, where item "varname" is the variable name and
  ## "value" is another list in which to look for arguments
  ## parname = name of function in parameter list (NULL if
  ## top-level function)
  ## dataset = array of lists; for each array item (which is a
  ## list), item "varname" is the variable name and "value" is a 
  ## data frame in which to look for data sources
  ## of function arguments
  ## ... = any additional arguments in which to look for
  ## function arguments
  ##
  ## To find the sources of function arguments, this will
  ## look first in par, then in dataset, then finally in ... if
  ## present.  If the function argument has a default, no error
  ## is thrown if a source of that argument is not found. If
  ## there is a default provided for a function argument that
  ## indicates it is of type character, then any character strings
  ## provided in par will not be resolved but will be passed as-is.
  ##
  ## Returns: a nested list representing function
  ## dependencies.  Each function call is packaged
  ## up into a list with its members as follows:
  ## parname = parameter name in list
  ## fun = function call location
  ## call = function call
  ## pre_eval = any function calls to evaluate first
  ##
  ## For example: if the model function is "linear(a,b,N)",
  ## N is listed in the par list as sum(), then the
  ## returned list would be as follows:
  ## results$call[[1]]$parname = NULL
  ## results$call[[1]]$fun = "linear"
  ## results$call[[1]]$call = list(a=par$a,b=par$b,N=eval_results$N)
  ## results$call[[1]]$pre_eval[[1]]$parname = N
  ## results$call[[1]]$pre_eval[[1]]$fun = par$N
  ## results$call[[1]]$pre_eval[[1]]$call = 1:5(or whatever)
  ## results$call[[1]]$pre_eval[[1]]$pre_eval = NULL
  ##
  ##########################################
  
  
  # Get the list of function arguments
  parnames<-names(par$value)
  default_checker <- formals(fun)
  fun_args <- names(default_checker)
  cols<-NULL
  for (i in 1:length(datasets)) {
    cols[[i]]<-names(datasets[[i]]$value)
  }
  extra_args<-list(...)
  extra_arg_names<-names(extra_args)
  
  # Here's where we'll assemble the arguments for calling the function
  results<-list(parname=parname, fun=fun, pre_eval=NULL)
  results$call<-list()
  results$pre_eval<-NULL

  if (is.null(fun_args)) {
    return(results)
  } 

  # Build the call to model by finding each argument in the parameter list
  for (i in 1:length(fun_args)) {
    if (fun_args[[i]] == "...") {
      if (length(extra_args) > 0) {
        results$call[["..."]]<-quote(...)
      }
    }
    else {
      used <- fun_args[[i]]==parnames
      if(any(used)) {
        if(is.function(par$value[used][[1]])) {
          
          # One of the arguments is a function - recursively analyze it
          results$pre_eval[[length(results$pre_eval)+1]]<-analyze_function(par$value[used][[1]], par, fun_args[[i]], datasets, ...)
          results$call[[fun_args[[i]]]]<-parse(text=paste("eval_results$", fun_args[[i]],sep=""))[[1]]
          
        }
        else {
          # This parameter is not a function - see if we can find it amongst
          # the parameters or source data    

          if (is.character(par$value[used][[1]])) {

            # Check to see if this is referencing the predicted value
            if (par$value[used][[1]]=="predicted") {
              results$call[[fun_args[[i]]]]<-quote(predicted)
            } 

            else {

              # If the function argument expects a character value, pass as-is
              if (mode(default_checker[[i]]) == "character") {
                results$call[[fun_args[[i]]]]<-parse(text=paste(par$varname, "$", fun_args[[i]],sep=""))[[1]]
              }
              else {               

                # Look for this in the dataset(s)
                found <- FALSE
                for (j in 1:length(cols)) { 
                  if(any(par$value[used][[1]]==cols[[j]])) {
                    # We found the argument in the data - add it referencing the data column
                    results$call[[fun_args[[i]]]]<-parse(text=paste(datasets[[j]]$varname, "$", par$value[used][[1]], sep=""))[[1]]
                    found<-TRUE
                    break 
                  }
                }
                if (!found) {
                  # This is a character argument and we haven't found it anywhere;
                  # so pass it as-is
                  results$call[[fun_args[[i]]]]<-parse(text=paste(par$varname, "$", fun_args[[i]],sep=""))[[1]]
                }
              }
            } 
          }
          else {
            # Non-character value - presumably numeric and to be passed as-is from par
            results$call[[fun_args[[i]]]]<-parse(text=paste(par$varname, "$", fun_args[[i]],sep=""))[[1]]
          }
        }
      }
      else {
        # Is this in the extra list? Only parse out if the function doesn't
        # take "..." as an argument
        used <- fun_args[[i]]==extra_arg_names
        if(any(used)) {
          if (!any(fun_args=="...")) {
            if(is.function(extra_args[used][[1]])) {
              results$pre_eval[[length(results$pre_eval)+1]]<-analyze_function(extra_args[used][[1]], par, fun_args[[i]], datasets, ...)
              results$call[[fun_args[[i]]]]<-parse(text=paste("eval_results$", fun_args[[i]],sep=""))[[1]]              
            }
            else {
              # This parameter is not a function - see if we can find it amongst
              # the parameters or source data    

              if (is.character(extra_args[used][[1]])) {

                # Check to see if this is referencing the predicted value
                if (extra_args[used][[1]]=="predicted") {
                  results$call[[fun_args[[i]]]]<-quote(predicted)
                } 

                else {

                  # If the function argument expects a character value, pass as-is
                  if (mode(default_checker[[i]]) == "character") {
                    results$call[[fun_args[[i]]]]<-parse(text=paste("list(...)$", fun_args[[i]],sep=""))[[1]]
                  }
                  else {               

                    # Look for this in the dataset

                    # Look for this in the dataset(s)
                    found <- FALSE
                    for (j in 1:length(cols)) { 
                      if(any(extra_args[used][[1]]==cols[[j]])) {
                        # We found the argument in the data - add it referencing the data column
                        results$call[[fun_args[[i]]]]<-parse(text=paste(datasets[[j]]$varname, "$", extra_args[used][[1]], sep=""))[[1]]
                        found<-TRUE
                        break 
                      }
                    }
                    if (!found) {
                      # This is a character argument and we haven't found it anywhere;
                      # so pass it as-is
                      results$call[[fun_args[[i]]]]<-parse(text=paste("list(...)$", fun_args[[i]],sep=""))[[1]]
                    }
                  }
                } 
              }
              else {
                # Non-character value - presumably numeric and to be passed as-is from extra_args
                results$call[[fun_args[[i]]]]<-parse(text=paste("list(...)$", fun_args[[i]],sep=""))[[1]]
              }
            }
          }
        }
        else {
          # Last chance - does this argument have a default?
          # There's gotta be a better way to find this but I can't figure it out
          if (mode(default_checker[[i]]) == "name") {
            # No default, and couldn't find the parameter anywhere
            stop("analyze_function:\nArgument \"",fun_args[[i]],"\" to function \"fun\" not found in par.\n")
          }
        }
      }
    }
  }
  results
}