



###############################################################################
#
# idbg.init - Initialize debugger data structure  
#
# parameters: force - force initialization even if debugger data already exists
#           
# the debugger data is stored in an environment called idbg.data under the environment
# of function idbg. Use function idbg() to get that environment
#
###############################################################################

idbg.init <- function(force=FALSE)
{
#  if (!exists("idbg.data") || force) 
#  {
#    idbg.data <- new.env() 
#    idbg.data[["break_frame"]] <- -1
#    idbg.data[["call_stack"]] <- list()
#    idbg.data[["debug_frame"]] <- -1
#    idbg.data[["ifunc_names"]] <- c()
#    idbg.data[["gui"]] <- FALSE 	#suppressWarnings(library("tcltk", character.only =TRUE, logical.return=TRUE))
#    idbg.data[["gui_toplevel"]] <- NULL
#  }

}

idbg <- function() {
  return(idbg.data)
}

idbg.call_stack_top <- function(frame_id, func_name, pos)
{
  idbg.data$call_stack[[frame_id]] <- list(func_name, pos)
  length(idbg.data$call_stack) <- frame_id
}


idbg.add_ifunc <- function(fname)
{
  if (! (fname %in% idbg.data$ifunc_names))
    idbg.data$ifunc_names <- c(idbg.data$ifunc_names, fname)
}

idbg.gui_mode <- function(is_gui=NA)
{
  return(FALSE)
#  is_gui <- as.logical(is_gui)
#  if (is.na(is_gui))
#    return(idbg.data$gui)
#  idbg.data$gui <- is_gui
#  return(idbg.data$gui)
}


###############################################################################
idbg.bp <- function( func_name, line_number=NA, condition=TRUE)
{
  #func_name <- as.character(substitute(func_name))
  f <- ifunc(func_name)

  return(breakpoint.ifunc(f, line_number, condition))
}
###############################################################################
idbg.get_breakpoints <- function()
{
  
  breakpoints <- NULL
  for (fname in idbg()[["ifunc_names"]])
  {
    f <-  idbg.match.ifunc(fname)
    if (is.null(f) || ! is.ifunc(f))
      next
    fbp <- list_breakpoints.ifunc(f)
    if (nrow(fbp))
      breakpoints <- rbind(breakpoints, data.frame("function_name"=fname,line=fbp$line, condition=fbp$condition, stringsAsFactors=FALSE))
  }
  return(breakpoints)
}
###############################################################################  
idbg.print_breakpoints <- function()
{  
  breakpoints <- idbg.get_breakpoints()
  for (line in capture.output(breakpoints))
    idbg.cat(line,"\n")
  idbg.cat("\n")
  #print(breakpoints)
}
###############################################################################
idbg.source <- function(fname)
{
  breakpoints <- idbg.get_breakpoints()
  source(fname)
  if (! is.null(breakpoints))
  {
    nbreakpoints<- nrow(breakpoints)
    for (i in seq_len(nbreakpoints))
    {
      idbg.bp(breakpoints$function_name[[i]], breakpoints$line[[i]])
    }
  }
}
###############################################################################
idbg.clear_all_breakpoints <- function()
{
  breakpoints <- idbg.get_breakpoints()
  if (! is.null(breakpoints))
  {
    nbreakpoints<- nrow(breakpoints)
    for (i in seq_len(nbreakpoints))
      idbg.bp(breakpoints$function_name[[i]], breakpoints$line[[i]], FALSE)
  }
}
###############################################################################
idbg.interact <- function(pos, func_name)
{
  gui <- idbg.gui_mode(NA)
  if (gui && is.null(idbg()$gui_toplevel))
    idbg.gui()
  
  debug_loop = TRUE

  frame_id <- sys.nframe()
  
  idbg.call_stack_top(frame_id-1, func_name, pos)
  last_debug_frame <- idbg()$debug_frame
  assign("debug_frame", frame_id-1, envir=idbg())
  
  #get("call_stack", envir=idbg())[[frame_id-1]] <- list(func_name, pos)
  #length(idbg()$call_stack) <- frame_id-1


  func <- match.fun(func_name)
  breakpoints <- attr(func, "data")$breakpoints
  
  # eval the breakpoint at the parent - for conditional breakpoint support
  if (! eval.parent(parse(text=breakpoints[[pos]])))
  {
    
    # no breakpoint at this point
    # we may stop in two cases
    # 1. step in command assigned the break_frame and we got to that frame
    # 2. we returned from a debugged function to the caller
    if (! is.na(idbg()$break_frame) && frame_id-1 != idbg()$break_frame) # && last_debug_frame <= frame_id-1)
      return(NULL)
    #browser()
    assign("break_frame", -1, envir=idbg())
  }

# GUI not implemented in current version
#  if (gui)
#    idbg.gui.set_entry_text(func_name,list_source_extended.ifunc(func, pos), TRUE,TRUE)
#  else
    cat(list_source.ifunc(func, pos))


  while (debug_loop)
  {
    # GUI not implemented in current version
	#if (gui)
    #  line <- idb.gui.wait_for_usr_cmd()
    #else
      line <- readline(sprintf("debug %d: ",pos));

    if (line == "")
      line <- "n"

    words = strsplit(line, "\ +")[[1]];
    if (length(words) == 0)
      next
    if (words[1] == "")
      words = words[2:length(words)];

    cmd = words[1];
    words = words[-1];
    if (cmd == "n")
    {
      assign("break_frame", sys.nframe() -1, envir=idbg())
      debug_loop <- FALSE
    }
    else
    if (cmd == "s")
    {
      addr <- attr(func, "data")$key2addr[[pos]]
      e <- body(func)
      has_error <- FALSE
      for (h in addr)
        if (h <= length(e))
          e <- e[[h]]
        else
        {
          has_error <- TRUE
          break
        }

      if (! has_error)
        idbg.prepare_step(e)
      assign("break_frame", NA, envir=idbg())
      debug_loop <- FALSE
    }
    else
    if (cmd == "c")
    {
      debug_loop <- FALSE
    }
    else
    if (cmd == "o")
    {
      assign("break_frame", sys.nframe()-2, envir=idbg())
      debug_loop <- FALSE
    }
    else
    if (cmd == "l")
    {
      q <- idbg()$call_stack[[idbg()$debug_frame]]
      lfunc <- q[[1]]
      # GUI not implemented in current version
	  #if (gui)
      #  idbg.gui.set_entry_text(lfunc,list_source_extended.ifunc(lfunc, NA), TRUE,TRUE)
      #else
      {
        lpos <- q[[2]]
        if (length(words) == 0)
          cat(list_source.ifunc(lfunc, lpos))
        else
        if (length(words) == 1)
          cat(list_source.ifunc(lfunc, lpos, words[[1]], words[[1]]))
        else
        if (length(words) >= 2)
          cat(list_source.ifunc(lfunc, lpos, words[[1]], words[[2]]))
      }
    }
    else
    if (cmd == "b")
    {
      if (length(words) == 0)
        idbg.print_breakpoints()
      else
	  if (words[[1]]=="FALSE")
        idbg.clear_all_breakpoints()
      else
      {
        # if first argument is a number than break at current func otherwise it is a func name string
        bp_line <-  suppressWarnings(as.numeric(words[[1]]))
        bp_cond_arg <- 2
        if (is.na(bp_line))
        {
          # this must be a name of a function, break on the 1st line 
          bp_func_name <- words[[1]]
          if (length(words) > 1)
          {  
            bp_line <-  as.numeric(words[[2]])
            if (! is.na(bp_line))
              bp_cond_arg <- 3
          }
        } 
        else
        {
          # this is a line number to break in current proc
          bp_func_name <- func_name
        }

        if (length(words) >= bp_cond_arg)
		  bp_condition <- do.call("paste", as.list(words[bp_cond_arg:length(words)]))
		else
          bp_condition <- TRUE
        idbg.bp(bp_func_name, bp_line, bp_condition) 
      }
      
      # GUI not implemented in current version
	  #if (gui)
      #  idbg.gui.set_entry_text(func_name,list_source_extended.ifunc(func, pos), TRUE,TRUE)

    }
    else
    if (cmd == "w") 
    {
      # print the stack
      # todo - if the stack include non debugable functions (eg. apply) print data about them
      widx <- 0
      for (q in idbg()$call_stack)
      {
        widx <- widx + 1
        if (is.null(q))
          idbg.cat("???\n")
        else
        {
          if (widx == idbg()$debug_frame)
            wchar <- "*"
          else
            wchar <- " "
          idbg.cat(wchar, widx, q[[1]]," ")
          idbg.cat(list_source.ifunc(q[[1]], q[[2]], 0, 0))
        }
      }
    }
    else
    if (cmd == "u") 
    {
      # up in the stack
      if (idbg()$debug_frame > 1)
        assign("debug_frame", idbg()$debug_frame-1, envir=idbg())
    }
    else
    if (cmd == "d") 
    {
      # down in the stack
      if (idbg()$debug_frame < frame_id)
        assign("debug_frame", idbg()$debug_frame+1, envir=idbg())
    }
    else
    if (cmd == "q")
    {
      invokeRestart(findRestart("abort"))
    }
    else
    if (cmd == "h")
    {
      idbg.cat("h - help. print this message\n")
      idbg.cat("n - next. Empty line is the same as 'n'\n")
      idbg.cat("s - step into\n")
      idbg.cat("o - step out\n")
      idbg.cat("c - continue\n")
      idbg.cat("q - quit\n")
      idbg.cat("b - print breakpoints\n")
	  idbg.cat("b FALSE - clear all breakpoints\n")
      idbg.cat("b <func_name> [FALSE] - set/unset a breakpoint in first line of function\n")
      idbg.cat("b <line_number> [FALSE] - set/unset a breakpoint in current function\n")
      idbg.cat("b <func_name> <line_number> [FALSE] - set/unset a breakpoint in function at at line_number\n")
      idbg.cat("w - print the stack\n")
      idbg.cat("u - go up the stack\n")
      idbg.cat("d - go down the stack\n")
      idbg.cat("l [nlines] - print nlines of source before and after current position\n")
      idbg.cat("l [nlines_back] [nlines_forward] - print source around current position\n")
      idbg.cat("x expr - execute expression. Any expression that doesn't match the above options will also be executed\n")
    }
    else
    {
      if (cmd == "x")
      {
        expr <- ""
        for (w in words)
          expr = paste(expr, w)
      }
      else
        expr <- line
      e<- try(
        #print(eval.parent(parse(text=expr))),
        capture.output(eval(parse(text=expr),envir=sys.frame(idbg()$debug_frame))),
        silent = TRUE
      )
      if (class(e) == "try-error")
        idbg.cat(geterrmessage())
      else
      {
		for (line in e)
          idbg.cat(line,"\n")
        idbg.cat("\n")
      }
      timestamp(expr,prefix="",suffix="",quiet =TRUE)
    }
  }
}
###############################################################################
idbg.match.ifunc <- function(fname)
{
  f <-match.fun(fname)
  if (is.ifunc(f))
    return(f)

  return(NULL)
}
###############################################################################
idbg.prepare_step <- function(expr)
{
  expr_len <- length(expr)

  
  e <- expr[[1]]
  ifunc(as.character(expr[[1]]))
  if (expr_len > 1)
  {
    for (i in 2:expr_len) 
    { 
      if (class(expr[[i]]) == "call")
        idbg.prepare_step(expr[[i]])
    }
  }
}
###############################################################################
idbg.instrument_if <- function(if_expr, func_name, ienv)
{
  l <- list(if_expr[[1]], if_expr[[2]])
  ienv$addr <- c(ienv$addr, NA)

  for (idx in 3:4)
  {
    if (length(if_expr) >= idx)
    {
      ienv$addr[[length(ienv$addr)]] <- idx
      ikey <- ienv$key

      if (class(if_expr[[idx]]) == "{")
        l <- c(l, idbg.instrument_expr_list(if_expr[[idx]], func_name, ienv))
      else
        l <- c(l, idbg.instrument_expr_list(as.call(c(as.name("{"), if_expr[[idx]])), func_name, ienv))
    }
  }
  ienv$addr <- ienv$addr[-length(ienv$addr)]
  return(as.call(l))
}
###############################################################################
idbg.instrument_for <- function(for_expr, func_name, ienv)
{
  key <- ienv$key
  l <- list(for_expr[[1]], for_expr[[2]], for_expr[[3]])
  ienv$addr <- c(ienv$addr, 4)

  if (class(for_expr[[4]]) == "{")
    l <- c(l, idbg.instrument_expr_list(for_expr[[4]], func_name, ienv))
  else
    l <- c(l, idbg.instrument_expr_list(as.call(c(as.name("{"), for_expr[[4]])), func_name, ienv))
  ienv$addr <- ienv$addr[-length(ienv$addr)]
  return(as.call(l))
}
###############################################################################
idbg.instrument_while <- function(while_expr, func_name, ienv)
{
  key <- ienv$key
  l <- list(while_expr[[1]], while_expr[[2]])
  ienv$addr <- c(ienv$addr, 3)

  if (class(while_expr[[3]]) == "{")
    l <- c(l, idbg.instrument_expr_list(while_expr[[3]], func_name, ienv))
  else
    l <- c(l, idbg.instrument_expr_list(as.call(c(as.name("{"), while_expr[[3]])), func_name, ienv))
  return(as.call(l))
}
###############################################################################
idbg.instrument_repeat <- function(repeat_expr, func_name, ienv)
{
  key <- ienv$key
  l <- list(repeat_expr[[1]])
  ienv$addr <- c(ienv$addr, 2)

  if (class(repeat_expr[[2]]) == "{")
    l <- c(l, idbg.instrument_expr_list(repeat_expr[[2]], func_name, ienv))
  else
    l <- c(l, idbg.instrument_expr_list(as.call(c(as.name("{"), repeat_expr[[2]])), func_name, ienv))
  return(as.call(l))
}
###############################################################################
idbg.instrument_expr_list <- function(expr, func_name, ienv)
{
  expr_len <- length(expr)
  l<-list()
  ikey <- ienv$key
  ienv$addr <- c(ienv$addr, NA)

  if (expr_len > 0 && class(expr[[1]]) == "name" && expr[[1]] != "{")
  {
    return(idbg.instrument_expr_list(as.call(c(as.name("{"), expr)), func_name, ienv))
  }
  
  for (i in seq_len(expr_len)) 
  { 
    ienv$addr[[length(ienv$addr)]] <- 2*(i-1)+1
    e <- expr[[i]]
    #ikey<- paste(key,i,sep=".")

    if (class(e) == "{")
      e<- idbg.instrument_expr_list(e, func_name, ienv)
    else
    if (class(e) == "if")
      e<- idbg.instrument_if(e, func_name, ienv)
    else
    if (class(e) == "for")
      e<- idbg.instrument_for(e, func_name, ienv)
    else
    if (class(e) == "while")
      e<- idbg.instrument_while(e, func_name, ienv)
    else    
    if (class(e) == "repeat" || (is.call(e) && e[[1]] == "repeat"))
      e<- idbg.instrument_repeat(e, func_name, ienv)
      
        
    ienv$key <- ienv$key + 1
    
    ikey <- ienv$key
    l <- c(l, e, substitute(idbg.interact(i,name), list(i=ikey,name=func_name))) 
    
    # increase by two in order to get the addess of the line following the bp and not the command address
    ienv$addr[[length(ienv$addr)]] <- ienv$addr[[length(ienv$addr)]] + 2
    ienv$key2addr[[ikey]] <- ienv$addr
    

  }

  ienv$addr <- ienv$addr[-length(ienv$addr)]

  return(as.call(l))
}
###############################################################################
idbg.gen_source <- function(instrumented_body)
{
  s <- format(instrumented_body)

  bp_pos <- regexpr("idbg.interact\\([0-9]+", format(s))
  #bp_len <- attr(bp_pos, "match.length") - 3
  is_bp <- bp_pos != -1
  nbp <- sum(is_bp)
  #bp_key <- as.numeric(substr(s[is_bp], bp_pos[is_bp] + 3, bp_pos[is_bp] + 3+ bp_len[is_bp]-1))
  key2line <- which(c(FALSE, is_bp)) - seq_len(nbp)
  s <- s[! is_bp]
  return(list(src=s, key2line=key2line))
}
###############################################################################
ifunc <- function(fname)
{
  if (! is.character(fname))
    return(NULL)

#  if (fname %in% builtins())
#    return(NULL)

  
  #func <- try (match.fun(fname))
  #if (length(class(func)) == 1 && class(func) == "try-error")
  #  return(NULL)

  # search for the func and the frame that it belongs to
  found <- FALSE
  for (frame_id in seq(sys.nframe(),0))
  {
    envir <- sys.frame(frame_id)
    if (exists(fname, envir = envir, mode = "function", inherits = FALSE))
    {
      func <- get(fname, envir = envir, mode = "function", inherits = FALSE)
      found <- TRUE
      break
    }
  }
  
  # looking for the function in packages doesn't work as there is no way to assign to the instrumented copy
  # cannot change value of locked binding for 'fname'
  #if (! found)
  #{
  #  for (envir in search()[-1])
  #	{
  #	  if (exists(fname, where = envir, mode = "function", inherits = FALSE))
  #	  {
  #	    func <- get(fname, pos = envir, mode = "function", inherits = FALSE)
  #		found <- TRUE
  #		break
  #	  }
  #	}  
  #}

  if (! found)
    return(NULL)

  if (is.primitive(func))
    return(func)
  
  if (is.ifunc(func))
    return(func)

  ienv <- new.env()
  ienv$key <- 0
  ienv$addr <- c()
  ienv$key2addr <- list()

  l <- idbg.instrument_expr_list(body(func),fname,ienv)
  ret <- func
  body(ret) <- as.call(l)  
  attr(ret,"orig") <- func
  q <-idbg.gen_source(l)
  data <- new.env()
  attr(ret,"data") <- data
  data[["src"]] <- q$src
  data[["key2line"]] <- q$key2line
  data[["key2addr"]] <- ienv$key2addr
  data[["breakpoints"]] <- rep(FALSE,length(q$key2line))


  class(ret) <- c("ifunc", class(ret))
  
  
  idbg.add_ifunc(fname)
  if (is.character(envir))
    assign(fname, ret, pos=envir)
  else	
    assign(fname, ret, envir=envir)
}
###############################################################################
is.ifunc <- function(x)
{
  return(inherits(x, "ifunc"))
}
###############################################################################
print.ifunc <- function(x)
{
  src <- attr(x, "data")$src
  cat(format(args(x))[[1]],"\n")
  for (line in src)
    cat(line,"\n")
}
###############################################################################
# to clear a breakpoint set expr to FALSE
# to toggle TRUR/FALSE set exor to NA
breakpoint.ifunc <- function(f, line_number, expr=TRUE)
{
  if (! is.ifunc(f) )
    return(FALSE)

  key2line <- attr(f, "data")$key2line
  if (length(key2line) == 0)
    return(FALSE)
  
  if (is.na(line_number))
    d <- 1
  else
    d <- which.min(abs(line_number - key2line))
  if (length(d) == 0)
    return(FALSE)
  d <- d[[1]]	
	
  if (is.na(expr))
    attr(f,"data")$breakpoints[[d]] <- ifelse( attr(f,"data")$breakpoints[[d]] == FALSE, TRUE, FALSE) 
  else
    attr(f,"data")$breakpoints[[d]] <- expr
  return(TRUE)
}
###############################################################################
line_breakpoint.ifunc <- function(f, line_number)
{
  if (! is.ifunc(f) )
    return(FALSE)

  key2line <- attr(f, "data")$key2line
  if (length(key2line) == 0)
    return(FALSE)
  
  d <- which(line_number - key2line==0)
  if (length(d) != 1)
    return(FALSE)
  return(attr(f,"data")$breakpoints[[d]])
}
###############################################################################
list_breakpoints.ifunc <- function(f)
{
  if (! is.ifunc(f))
    return(NULL)
  keys <- which(attr(f,"data")$breakpoints != FALSE)
  lines <- (attr(f, "data")$key2line)[keys]
  conditions <- attr(f,"data")$breakpoints[keys]
  return(data.frame(line=lines, condition=conditions, stringsAsFactors=FALSE))
}
###############################################################################
list_source.ifunc <- function(func, pos, back=10, forward=10)
{
  back <- suppressWarnings(as.integer(back))
  forward <- suppressWarnings(as.integer(forward))

  if (is.na(back) || back < 0)
    back <- 10
  if (is.na(forward) || forward < 0)
    forward <- 10

  if (is.character(func))
    func <- match.fun(func)

  txt_buffer <- ""
  if (is.ifunc(func))
  {  
    src <- attr(func, "data")$src
    key2line <- attr(func, "data")$key2line
    pos <- key2line[[pos]]
    start <- pos - back
    end <- pos + forward
  
    if (start < 1 )
      start <- 1
    if (end > length(src))
      end <- length(src)

    if (start == 1)
      txt_buffer <- paste(format(args(func))[[1]],"\n",sep="")

    for (i in start:end ) 
    { 
      line <- sprintf("%04d %s",i,src[[i]])
      if (i == pos)
        txt_buffer <-  paste(txt_buffer,"=>",sep="")
      else
        txt_buffer <-  paste(txt_buffer,"  ",sep="")
      txt_buffer <-  paste(txt_buffer,line,"\n",sep="")
    }
  }
  return(txt_buffer)
}
###############################################################################
list_source_extended.ifunc <- function(func, pos)
{
  if (is.character(func))
    func <- match.fun(func)
  
  if (!is.ifunc(func))
    return(NULL)
    
  src <- attr(func, "data")$src
  key2line <- attr(func, "data")$key2line
  
  start <- 1
  end <- length(src)

  
  # TODO:   set the format to be sprintf("%%0%dd %s" , ceiling(log10(end+1)))
  result <- data.frame(SRC=rep("", end+1), IP=FALSE, BP=FALSE, stringsAsFactors=FALSE)
  result$SRC[[1]] <- paste("     ", format(args(func))[[1]], sep="")
  result$SRC[2:(end+1)] <- sprintf("%04d %s",seq_len(end),src)
  if (! is.na(pos))
  {
    pos <- key2line[[pos]]
    result$IP[[pos+1]] <- TRUE
  }
  result$BP[attr(func, "data")$key2line +1] <- attr(func,"data")$breakpoints

  return(result)
}
###############################################################################
idbg.cat <- function(...)
{
  str1 <- do.call("paste", list(...))
  # GUI not implemented in current version
  #if (idbg.gui_mode())
  #  idb.gui.bottom.cat(str1)
  #else
    cat(str1)
}
###############################################################################
idbg.gui <- function()
{
#
#tt <- tktoplevel()
#tkwm.protocol(tt,"WM_DELETE_WINDOW", function()idbg.gui.close())
##tkbind(tt,"<Destroy>", function(W)cat("Bye\n"))
#
#tkbind( tt, "<KeyPress-F5>", function(K)idb.gui.key_press(K) )
#tkbind( tt, "<KeyPress-F6>", function(K)idb.gui.key_press(K) )
#tkbind( tt, "<KeyPress-F8>", function(K)idb.gui.key_press(K) )
#tkbind( tt, "<KeyPress-F9>", function(K)idb.gui.key_press(K) )
#
#idb.gui.send_cmd <<- function(cmd)
#{
#  if (cmd != "")
#    tcl("set", "idb_gui_user_cmd",cmd)
#}
#
#idbg.gui.close <- function()
#{
#  idb.gui.send_cmd("q") 
#  tkdestroy(tt)
#  assign("gui_toplevel", NULL, envir=idbg())
#}
#
#
#tt.menu <- tkmenu(tt, tearoff=0)
#
#tt.menu.file <-tkmenu(tt.menu, tearoff=0)
#tkadd(tt.menu, "cascade", label="File", menu=tt.menu.file, underline=0)
#tkadd(tt.menu.file, "command", label="Exit", command=function()idbg.gui.close())
#tkconfigure(tt,menu=tt.menu)
##$m add separator
#
#
#tt.pane <- ttkpanedwindow(tt,orient="vertical")
#tt.pane.top <- ttkframe(tt.pane)
#tt.pane.bottom <- ttkframe(tt.pane)
#
#
### Make the notebook and set up Ctrl+Tab traversal
#tt.pane.top.note <- ttknotebook(tt.pane.top)
#tkpack(tt.pane.top.note, fill="both", expand=1, padx=2, pady=3)
##ttk::notebook::enableTraversal $w.note
#
#idb.gui.left_click <<- function(note_entry.text,x,y, func_name)
#{
#  addr <- strsplit(as.character(tkindex(note_entry.text,paste("@",x,",",y,sep=""))),"\\.")
#  row <- as.numeric(addr[[1]][[1]])
#  col <- as.numeric(addr[[1]][[2]])
#
#  #cat("left click",x,y, row,col,"\n")
#  if (col ==0)
#    idbg.bp( func_name, row-1, NA)
#  idbg.gui.set_entry_bp(note_entry.text, row, line_breakpoint.ifunc(ifunc(func_name), row-1))
#}
#
#idb.gui.right_click <<- function(note_entry.text,x,y)
#{
#  row.col <- tkindex(w,paste("@",note_entry.text,",",y,sep=""))
#  cat("right click",x,y, as.character(tkindex(w,paste("@",note_entry.text,",",y,sep=""))),"\n")
#
#}
#
#
#idb.gui.wait_for_usr_cmd <<- function()
#{
#  tcl("set", "idb_gui_user_cmd","")
#  tkwait.variable("idb_gui_user_cmd")
#  cmd <- as.character(tclvalue("idb_gui_user_cmd"))
#  #cat("idb.gui.wait_for_usr_cmd: cmd='",cmd,"'\n",sep="")
#  return(cmd)
#}
#
#
#idb.gui.key_press <<- function(K)
#{
#  #cat("K='",K,"'\n",sep="")
#  switch(K,
#    F5={ cmd<- "s" },
#    F6={ cmd<- "n" },
#    F8={ cmd<- "o" },
#    F9={ cmd<- "c" }
#  )
#
#  idb.gui.send_cmd(cmd)
#}
#
#
#
#idb.gui.bottom.key_press <<- function(bottom.text, K)
#{
#  if (K == "space")
#    K <- " "
#  if (K == "Tab")
#    K <- "\t"
#
#  cat("K='",K,"'\n",sep="")  
#
#  pos <-strsplit(as.character(tkindex(bottom.text, "end")), "\\.")
#  last_char <- as.numeric(pos[[1]][[2]])
#  last_line <- as.numeric(pos[[1]][[1]])-1
#
#  if (nchar(K) == 1 || K == "Delete" || K == "BackSpace" || K == "Return" || K == "Control-x" || K == "Control-v" || K == "Shift-Insert")
#  {
#    # if we are not at the last line set the insert position to the end of text
#    pos <-strsplit(as.character(tkindex(bottom.text, "insert")), "\\.")
#    line <- as.numeric(pos[[1]][[1]])
#    cat("line=",line,"\n")
#    cat("last_line=",last_line,"\n")
#    if (line != last_line)
#    {
#      if (K == "Control-x" || K == "Control-v" || K == "Shift-Insert" || K == "BackSpace")
#      {
#        # clear the selection
#        v <-tktag.ranges(bottom.text, "sel")
#        nranges <- as.numeric(as.character(tcl("llength", v)))
#        for (i in (seq_len(nranges/2)-1)*2)
#        {
#          index1<-tcl("lindex", v,i)
#          index2<-tcl("lindex", v,i+1)
#
#          tktag.remove(bottom.text, "sel", index1, index2)
#        }
#      }
#      tkmark.set(bottom.text, "insert", "end")
#    } 
#    else
#    if (K == "BackSpace")
#    {
#      if (last_char == 0)
#        tkinsert(bottom.text, "end", "\n")
#    }
#    else
#    if (K == "Return")
#    {
#      # if enter is pressed while in the middle of the text right to the insert will me moved to next line
#      # avoid that by moving to the end of line anyway
#      tkmark.set(bottom.text, "insert", "end")
#    }
#  }
#  
#  if (K == "Return")
#  {
#    cmd <- as.character(tkget(bottom.text, sprintf("%d.0", last_line), "end"))
#    if (last_char == 0 && last_line > 0)
#    {
#      cat("last_char=",last_char,"\n")
#      cat("last_line=",last_line,"\n")
#      # if we just got an ENTER don't print an empty line -> remove the last enter (==last character)
#      cat("kkkkkkkkkkkkkkkkkkk\n")
#      tkinsert(bottom.text, "end", "n")
#      #tkdelete(bottom.text,sprintf("%d.end",last_line-1) )
#      #tkmark.set(bottom.text, "insert", sprintf("%d.end",last_line-1))
#    }
#
#    tcl("set", "idb_gui_user_cmd",cmd)
#  }
#}
#
#
#idbg.gui.get_tab <<- function(tab_name, b_create=FALSE, b_select=FALSE)
#{
#  ntabs <- as.integer(as.character(tkindex(tt.pane.top.note,"end")))
#  for (tabid in seq_len(ntabs)-1)
#  {
#    tab_text <- as.character(tcl(tt.pane.top.note, "tab",tabid,"-text"))  
#    if (tab_text == tab_name)
#    {
#      if (b_select)
#        tkselect(tt.pane.top.note, tkindex(tt.pane.top.note,tabid))
#
#      return(tabid)
#    }
#  } 
#  if (b_create)
#  {
#    note_entry <- ttkframe(tt.pane.top.note)
#    note_entry.vscrollbar <- tkscrollbar(note_entry,command=function(...)tkyview(note_entry.text,...)) 
#    note_entry.hscrollbar <- tkscrollbar(note_entry,command=function(...)tkxview(note_entry.text,...), orient="horiz") 
#    note_entry.text <- tktext(note_entry, setgrid=1, height=20, undo=1, autosep=1, wrap="none", state="disabled",exportselection=1, yscrollcommand=function(...)tkset(note_entry.vscrollbar,...), xscrollcommand=function(...)tkset(note_entry.hscrollbar,...))
#    tkpack(note_entry.vscrollbar,side="right",fill="y")
#    tkpack(note_entry.hscrollbar,side="bottom",expand="no",fill="both")
#    tkpack(note_entry.text, expand="yes", fill="both") 
#    tktag.configure(note_entry.text, "ip_color", foreground="#00aa00")
#    tktag.configure(note_entry.text, "bp_color", foreground="red")
#    tkadd(tt.pane.top.note, note_entry, text=tab_name, underline=0, padding=2)
#    tkbind( note_entry.text, "<Button-3>", function(x,y)idb.gui.right_click(note_entry.text,x, y, tab_name) )
#    tkbind( note_entry.text, "<Button-1>", function(x,y)idb.gui.left_click(note_entry.text, x, y, tab_name) )
#    tkbind( note_entry.text, "<KeyPress-F5>", function(K)idb.gui.key_press(K) )
#    tkbind( note_entry.text, "<KeyPress-F6>", function(K)idb.gui.key_press(K) )
#    tkbind( note_entry.text, "<KeyPress-F8>", function(K)idb.gui.key_press(K) )
#    tkbind( note_entry.text, "<KeyPress-F9>", function(K)idb.gui.key_press(K) )
#
#    if (b_select)
#      tkselect(tt.pane.top.note, tkindex(tt.pane.top.note,ntabs))
#    return(ntabs)
#  }
#
#  return(-1)
#}
#
#idbg.gui.create_tab <<- function(tab_name)
#{
#  idbg.gui.get_tab(tab_name, TRUE)
#}
#
#idbg.gui.get_tab_obj <<- function(tab_name, b_create=FALSE, b_select= FALSE)
#{
#  tabid <- idbg.gui.get_tab(tab_name, b_create, b_select)  
#  if (tabid == -1)
#    return(NULL)
#
#  tab_obj <- tcl("lindex",tcl(tt.pane.top.note, "tabs"),tabid)
#}
#
#idbg.gui.set_entry_text <<- function(tab_name, text_df, b_create=FALSE, b_select= FALSE, incremental=TRUE)
#{
#  obj <- idbg.gui.get_tab_obj(tab_name,FALSE, b_select)
#  update_source <- incremental && is.null(obj)
#  #cat("update_source=",update_source,"\n")
#  if (is.null(obj))
#    obj <- idbg.gui.get_tab_obj(tab_name,b_create, b_select)
#  if (is.null(obj))
#    return()
#
#  text_entry_obj <- tcl("lindex", tkpack.slaves(obj),2)
#  tkconfigure(text_entry_obj,state="normal")
#  
#  if (update_source)
#  {
#    tkdelete(text_entry_obj,"0.0","end")
#
#    n <- nrow(text_df)
#    for (i in seq_len(n))
#    {
#      if (text_df$BP[[i]])
#        tkinsert(text_entry_obj, "end", "O", "bp_color")
#      else
#        tkinsert(text_entry_obj, "end", " ")
#
#      if (text_df$IP[[i]])
#        tkinsert(text_entry_obj, "end", "->", "ip_color")
#      else
#        tkinsert(text_entry_obj, "end", "  ")
#        tkinsert(text_entry_obj, "end", paste(text_df$SRC[[i]],"\n",sep=""))
#    }
#  }
#  else
#  {
#    # just update the ip
#    v <-tktag.ranges(text_entry_obj, "ip_color")
#    nranges <- as.numeric(as.character(tcl("llength", v)))
#    for (i in (seq_len(nranges/2)-1)*2)
#    {
#      index1<-tcl("lindex", v,i)
#      index2<-tcl("lindex", v,i+1)
#
#      tktag.remove(text_entry_obj, "ip_color", index1, index2)
#      tcl(text_entry_obj, "replace", index1, index2, "  ")
#    }
#    ip_row <- which(text_df$IP)
#    tcl(text_entry_obj, "replace", paste(ip_row,1,sep="."), paste(ip_row,3,sep="."), "->", "ip_color")
#
#  }
#  tkconfigure(text_entry_obj,state="disabled")
#}
#
#
#
#idbg.gui.set_entry_bp <<- function(text_entry_obj, row, value)
#{
#  tkconfigure(text_entry_obj,state="normal")
#  
#  pos <- sprintf("%d.0",row)
#  tkdelete(text_entry_obj,pos)
#  if (is.logical(value) && value == FALSE)
#    tkinsert(text_entry_obj, pos, " ")
#  else
#    tkinsert(text_entry_obj, pos, "O", "bp_color")
#    
#  tkconfigure(text_entry_obj,state="disabled")
#}
#
#
#
#idbg.gui.create_tab("help")
#
##tkselect(tt.pane.top.note, tkindex(tt.pane.top.note,"2"))
##txt<-tkcget(tt.pane.top.note,"text")
##txt <- as.character(tcl(tt.pane.top.note, "tab","2","-text"))
##cat(txt,"\n")
#
#
#
#tt.pane.bottom.vscrollbar <- tkscrollbar(tt.pane.bottom,command=function(...)tkyview(tt.pane.bottom.text,...)) 
#tt.pane.bottom.hscrollbar <- tkscrollbar(tt.pane.bottom,command=function(...)tkxview(tt.pane.bottom.text,...), orient="horiz") 
#tt.pane.bottom.text <- tktext(tt.pane.bottom, setgrid=1, height=7, undo=1, autosep=1, wrap="none", exportselection=1, yscrollcommand=function(...)tkset(tt.pane.bottom.vscrollbar,...), xscrollcommand=function(...)tkset(tt.pane.bottom.hscrollbar,...))
#tktag.configure(tt.pane.bottom.text, "cat_color", foreground="blue")
#tkpack(tt.pane.bottom.vscrollbar,side="right",fill="y")
#tkpack(tt.pane.bottom.hscrollbar,side="bottom",expand="no", fill="both")
#tkpack(tt.pane.bottom.text, expand="yes", fill="both") 
#tkbind(tt.pane.bottom.text, "<Control-x>", function(K)idb.gui.bottom.key_press(tt.pane.bottom.text, "Control-x") )
#tkbind(tt.pane.bottom.text, "<Control-v>", function(K)idb.gui.bottom.key_press(tt.pane.bottom.text, "Control-v") )
#tkbind(tt.pane.bottom.text, "<Shift-Insert>", function(K)idb.gui.bottom.key_press(tt.pane.bottom.text, "Shift-Insert") )
#tkbind(tt.pane.bottom.text, "<KeyPress>", function(K)idb.gui.bottom.key_press(tt.pane.bottom.text, K) )
#
#
#tkadd(tt.pane, tt.pane.top)
#tkadd(tt.pane, tt.pane.bottom)
#
#
#tkpack(tt.pane,side="top",expand="yes",fill="both",pady=2,padx="2m")
#
#
#idb.gui.bottom.cat <<- function(str)
#{
#  tkinsert(tt.pane.bottom.text, "end", str, "cat_color")
#  tkyview(tt.pane.bottom.text, "moveto","1.0")
#}
#
#assign("gui_toplevel", tt, envir=idbg())
#
#return(invisible())
#
}


.First.lib <- function(lib, pkg)
{
  #idbg.init()
}
  
if (!exists("idbg.data") || force) 
{
    idbg.data <- new.env() 
    idbg.data[["break_frame"]] <- -1
    idbg.data[["call_stack"]] <- list()
    idbg.data[["debug_frame"]] <- -1
    idbg.data[["ifunc_names"]] <- c()
    idbg.data[["gui"]] <- FALSE 	#suppressWarnings(library("tcltk", character.only =TRUE, logical.return=TRUE))
    idbg.data[["gui_toplevel"]] <- NULL
}
