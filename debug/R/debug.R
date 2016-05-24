# This is package debug 

".end.incarnation" <-
function() do.in.envir( envir=find.debug.HQ( FALSE), {
#  cat( 'Ending', .frames.[nrow( .frames.), 'function.name'], '\n')
# In R 1.8.0, seem to need to de-link the listboxes before calling 'tkdestroy'
  if( nrow( .frames.) && .frames.$has.window.yet[ nrow( .frames.)])  
    try( evalq( {
        tkconfigure( line.list.win, yscroll=function(...) {}, xscroll=function(...) {})
        tkconfigure( bp.win, yscroll=function(...){})
        tkdestroy( tcl.win) },
      sys.frame( .frames.$debug[ nrow( .frames.)]) ) # evalq
    ) # try

  .frames. <<- .frames.[ -nrow( .frames.),]
  
  # Command recall: clear up temp history file if no more debug windows
  if( !nrow( .frames.) && debug.command.recall)
    unlink( debug.hist.file)
  .nothing.
})


".onAttach" <-
function( libname, pkgname){
  f <- function( val) blah-blah-blah
  for( x in cq( tracees)) {
    body( f) <- substitute( if( missing( val)) x else x <<- val, list( x=as.name( x)))
    environment( f) <- asNamespace( 'debug')
    makeActiveBinding( x, f, as.environment( 'package:debug'))
    dont.lockBindings( x, 'debug', namespace.=FALSE)
  }    
}


".onLoad" <-
function (libname, pkgname) {
  set.presave.hook.mvb(untracer.env)
  
  # Stuff for namespace:
  evalq({
    tracees <- list()
    if (getRversion() >= "2.12") {
      my.index <- function(x, i) if (length(i)) 
        x[[i]]
      else x
      my.index.assign <- function(x, i, value) if (length(i)) 
        `[[<-`(x, i, value)
      else value
    }
    environment(my.index) <- baseenv()
    
    # Idiocy of next 3 lines imposed by CRANkiness
    balloonIsTethered <- mvbutils:::balloonIsTethered
    tetherBalloon <- mvbutils:::tetherBalloon
    untetherBalloon <- mvbutils:::untetherBalloon
    
    step.intos <- c(with = TRUE, within = TRUE, eval = TRUE, 
      evalq = TRUE, try = TRUE, suppressWarnings = TRUE)
  }, asNamespace(pkgname))
  
  dont.lockBindings(c("tracees", "step.intos"), pkgname)
  
  # Set 'try2' in a semi-future-proof way
  try2 <- try
  environment(try2) <- asNamespace(pkgname)
  Rcmd.check.cat <- cat
  if (!identical(body(try)[[2]][[1]], quote(tryCatch))) 
    warning("Can't catch ESC interrupts-- 'try' format has changed-- please notify MVB")
  else {
    # 
    interrupt.fun <- function(e) {
      dog <- cat # anti-CRANal: doesn't like cat() inside onLoad...
      dog("<Interrupted!>\n", file = find.debug.HQ()$debug.catfile())
      invisible(structure("Interrupt", class = "try-error"))
    }
    body(try2)[[2]]$interrupt <- interrupt.fun
  }
  assign("try2", try2, envir = asNamespace(pkgname))
}


".update.debug.window" <-
function( nlocal=sys.parent(), l) mlocal({
  l <- screen.line( lno)
  if( l != old.l) {
    if( textaroo) {
      if( old.l>0) # not first time
        tktag.remove( line.list.win, 'nexect', sprintf( '%i.%i', old.l+1, 0), 
            sprintf( '%i.%i', old.l+1, nchar( line.list[ old.l+1])))
        tktag.add( line.list.win, 'nexect', sprintf( '%i.%i', l+1, 0), 
            sprintf( '%i.%i', l+1, nchar( line.list[ l+1])))
    } else { # listaroo
      if( old.l>0) # not first time
        tkitemconfigure( line.list.win, old.l, background='White', selectforeground='White')
      tkitemconfigure( line.list.win, l, background='Green', selectforeground='Green')
    }
    old.l <- l
  } # l has moved

  tksee( line.list.win, l %&% if( textaroo) '.0') 
  # Also nice to see 2 lines ahead, unless window is too low... dunno how to find that, so use default height
  tksee( line.list.win, min( l+min( 2, getOption( 'debug.height', 10)), length( line.list)) %&% if( textaroo) '.0')

# Make sure breakpoints show OK
  for( l in seq_along( breakpoints)) {
    colour <- if( mark.bp( breakpoints[[l]])) 'red' else 'white'
    tkitemconfigure( bp.win, screen.line( l), foreground=colour, selectforeground=colour)
  }
})


".update.window.with.on.exit" <-
function( nlocal=sys.parent(), i, l) mlocal({
  l <- length( line.list)
  
  if( textaroo) {
  
  } else {
    tkdelete( line.list.win, orig.line.list.length-1, prev.line.list.length-1)
    #  tkinsert( line.list.win, 'end', line.list[ orig.line.list.length:l])
    #  .Tcl( .Tcl.args( line.list.win, 'insert', 'end', line.list[ orig.line.list.length:l]))
    do.call( 'tkinsert', c( list( line.list.win, 'end'), 
        as.vector( line.list[ orig.line.list.length %upto% l])))

    if( prev.line.list.length > l)
      tkdelete( bp.win, l-1, prev.line.list.length-2)
    else if( prev.line.list.length < l)
      do.call( 'tkinsert', c( list( bp.win, 'end'), rep( '*', l-prev.line.list.length)))

    # .Tcl( .Tcl.args( bp.win, 'insert', 'end', rep( '*', l-prev.line.list.length)))
    # tkinsert( bp.win, 'end', rep( '*', l-prev.line.list.length))

    # Blank any old breakpoints that would fall on non-breakpointable on-exit lines
    for( i in orig.line.list.length:l)
      tkitemconfigure( bp.win, i-1, foreground='white', selectforeground='white')
  } # if textaroo
})


"add.numbers" <-
function (expr, width = options()$width, numbering = TRUE, cat.on.exit = FALSE, 
    expr.offset = 0, line.number.offset = 0, preamble = character(0)) 
{
    src <- attr(expr, "source")
    if (is.null(src)) {
        nullf <- expr
        attributes(nullf) <- list()
        environment(nullf) <- .GlobalEnv
        src <- suppressWarnings(deparse(nullf, control = "all", 
            width.cutoff = 120L))
    }
    src.expr <- parse(text = src)
    re.bloody.quire <- require
    using.parser <- getOption("debug.src", FALSE) && suppressWarnings(re.bloody.quire("parser"))
    if (using.parser) {
        src.expr <- parser(text = src)
        ld <- attr(src.expr, "data")
    }
    old.width <- options(width = 1000)$width
    if (cat.on.exit) 
        cat.on.exit <- expression(cat(paste(format(names(line.list)), 
            line.list, sep = ":"), sep = "\n"))[[1]]
    on.exit({
        eval(cat.on.exit)
        options(width = old.width)
    })
    tab.width <- option.or.default("tab.width", 4)
    spaces <- function(n = 1) paste(rep(" ", n), collapse = "")
    tab.sp <- spaces(tab.width)
    prefix <- ""
    tabs <- function(more.or.less = 0, exact = n.tabs + more.or.less) {
        answer <- paste(rep(tab.sp, max(exact, 0)), collapse = "") %&% 
            prefix
        prefix <<- ""
        answer
    }
    assign("[[", my.index)
    assign("[[<-", my.index.assign)
    deparse1 <- function(x) {
        x <- deparse(x, width.cutoff = 500)
        if (length(x) > 3) 
            x <- c(x[1], paste(x[-c(1, length(x))], collapse = "; "), 
                x[length(x)])
        paste(x, collapse = " ")
    }
    ch <- function(...) {
        stuff <- c(...)
        stuff[1] <- stuff[1] + expr.offset
        paste(c(stuff, "."), collapse = ",")
    }
    add.to.last.line <- function(x) if (!is.na(x)) 
        line.list[length(line.list)] <<- paste(line.list[length(line.list)], 
            x, sep = "")
    make.indent <- function(n.in = 1, loop = FALSE) {
        indents[ch(i)] <<- n.in
        n.tabs <<- n.tabs + n.in
        if (loop) 
            last.loop.tabs <<- c(last.loop.tabs, n.tabs)
    }
    if ("debug" %in% loadedNamespaces()) {
        debuggify.system.call <- function(step.check = FALSE, 
            nlocal = sys.parent()) mlocal(expr[[c(i, 1)]] <- if (step.check) 
            substitute({
                if (debug:::stepping(call.type)) debug:::debug.fun.name else orig.fun.name
            }, list(call.type = call.type, debug.fun.name = as.name("debug." %&% 
                as.character(expr[[c(i, 1)]])), orig.fun.name = expr[[c(i, 
                1)]]))
        else call(":::", quote(debug), as.name("debug." %&% as.character(expr[[c(i, 
            1)]]))))
    }
    else {
        debuggify.system.call <- function(step.check = FALSE, 
            nlocal = sys.parent()) mlocal(expr[[c(i, 1)]] <- if (step.check) 
            substitute({
                if (stepping(call.type)) debug.fun.name else orig.fun.name
            }, list(call.type = call.type, debug.fun.name = as.name("debug." %&% 
                as.character(expr[[c(i, 1)]])), orig.fun.name = expr[[c(i, 
                1)]]))
        else as.name("debug." %&% as.character(expr[[c(i, 1)]])))
    }
    default.update.line.list <- function(nlocal = sys.parent()) mlocal(line.list <- c(line.list, 
        tabs() %&% deparse1(expr[[i]])))
    if (is.a.function <- is.function(expr)) {
        line.list <- clip(deparse(do.call("args", list(expr))))
        line.list <- paste(c(line.list, preamble), collapse = " ")
        names(line.list) <- ""
        expr <- body(expr)
    }
    else line.list <- ""
    expr <- do.call("expression", list(expr))
    src.expr <- do.call("expression", list(src.expr[[c(1, length(src.expr[[1]]) - 
        1)]]))
    if (using.parser) {
        id <- ld$id
        numbered.id <- numeric(0)
    }
    breakpoint <- vector("list", 0)
    n <- line.number.offset
    i <- 1
    n.tabs <- 1
    last.loop.tabs <- 0
    suffix <- structure("", names = ch(1))
    indents <- numeric(0)
    while (length(i)) {
        needs.a.number <- TRUE
        next.i <- numeric(0)
        call.above <- as.character(expr[[c(clip(i), 1)]])
        if (i[length(i)] != 1 && call.above == "switch") {
            add.to.last.line(",")
            if (i[length(i)] <= length(expr[[clip(i)]])) 
                line.list <- c(line.list, tabs(-1) %&% "'" %&% 
                  names(expr[[clip(i)]])[i[length(i)]] %&% "' = ")
            if (i[length(i)] == length(expr[[clip(i)]]) && is.name(expr[[i]]) && 
                !nzchar(expr[[i]])) 
                expr[[clip(i)]][tail(i, 1)] <- list(NULL)
        }
        if (mode(expr[[i]]) == "(" || !is.call(expr[[i]])) {
            anything.to.add <- deparse1(expr[[i]])
            if (nzchar(anything.to.add)) 
                line.list <- c(line.list, tabs() %&% anything.to.add)
        }
        else {
            call.type <- expr[[c(i, 1)]]
            if (!is.name(call.type)) 
                call.type <- "default"
            else call.type <- as.character(call.type)
            switch(call.type, `{` = {
                if ((length(i) == 1) || call.above %in% c("if", 
                  "switch", "for", "while", "repeat")) add.to.last.line(" {") else {
                  line.list <- c(line.list, tabs() %&% "{")
                  make.indent()
                }
                needs.a.number <- FALSE
                automove <- FALSE
                suffix[ch(i)] <- " }"
                if (length(expr[[i]]) == 1) expr[[i]] <- call("{", 
                  NULL)
                next.i <- c(i, 2)
            }, `(` = {
            }, `if` = {
                if (call.above != "if" || i[length(i)] != 4) {
                  line.list <- c(line.list, tabs())
                  make.indent()
                }
                add.to.last.line("if( " %&% deparse1(expr[[c(i, 
                  2)]]) %&% ")")
                next.i <- c(i, 3)
            }, switch = {
                line.list <- c(line.list, tabs() %&% "switch( " %&% 
                  deparse1(expr[[c(i, 2)]]))
                make.indent(2)
                suffix[ch(i)] <- ")"
                next.i <- c(i, 3)
            }, `for` = {
                line.list <- c(line.list, paste(tabs(), "for( ", 
                  expr[[c(i, 2)]], " in ", deparse1(expr[[c(i, 
                    3)]]), " )", sep = ""))
                make.indent(loop = TRUE)
                next.i <- c(i, 4)
            }, `while` = {
                line.list <- c(line.list, paste(tabs(), "while( ", 
                  deparse1(expr[[c(i, 2)]]), " )", sep = ""))
                make.indent(loop = TRUE)
                next.i <- c(i, 3)
            }, `repeat` = {
                line.list <- c(line.list, tabs() %&% "repeat")
                make.indent(loop = TRUE)
                needs.a.number <- FALSE
                next.i <- c(i, 2)
            }, `<-` = {
                temp.expr <- expr[[i]]
                temp.expr[[3]] <- 0
                temp.expr <- paste(deparse(temp.expr, width.cutoff = 500), 
                  collapse = " ")
                temp.expr <- substring(temp.expr, 1, nchar(temp.expr) - 
                  1)
                prefix <- prefix %&% temp.expr
                needs.a.number <- FALSE
                next.i <- c(i, 3)
            }, `break` = , `next` = {
                line.list <- c(line.list, tabs(exact = last.loop.tabs[length(last.loop.tabs)] - 
                  1) %&% call.type)
                debuggify.system.call()
            }, return = , invisible = {
                line.list <- c(line.list, deparse1(expr[[i]]))
                debuggify.system.call()
            }, with = , within = , try = , suppressWarnings = , 
                eval = , evalq = {
                  default.update.line.list()
                  debuggify.system.call(step.check = TRUE)
                }, local.on.exit = , on.exit = {
                  default.update.line.list()
                  debuggify.system.call()
                }, default.update.line.list())
        }
        if (needs.a.number) {
            n <- n + 1
            if (using.parser) 
                numbered.id[n] <- attr(src.expr[[i]], "id")
            names(line.list)[length(line.list)] <- format(c(n, 
                1000))[1]
            names(line.list)[is.na(names(line.list))] <- ""
            if (numbering) 
                breakpoint[[ch(i)]] <- is.a.function && (length(breakpoint) == 
                  0)
        }
        if (length(next.i)) 
            i <- next.i
        else addnum.move.to.next.expr()
    }
    if (using.parser) {
        num.ld <- ld[match(numbered.id, id, 0) %such.that% (. > 
            0), ]
        pop <- lapply(1:n, function(ni) {
            popi <- num.ld$id[ni]
            while (popi[1] != 0) popi <- c(ld$parent[ld$id == 
                popi[1]], popi)
            return(popi)
        })
        cpop <- list(numeric(0))
        for (i in 2 %upto% n) cpop[[i]] <- unique(c(cpop[[i - 
            1]], pop[[i - 1]]))
        nonfirst <- index(duplicated(num.ld$line1))
        nonlast <- index(rev(duplicated(rev(num.ld$line1))))
        zap.upto <- rep(0, n)
        zap.from <- nchar(src[num.ld$line1]) + 1
        for (izap in nonfirst) {
            ancid <- (pop[[izap]] %except% cpop[[izap]])[1]
            zap.upto[izap] <- ld$col1[ld$id == ancid]
        }
        zap.from[nonlast] <- zap.upto[nonlast + 1] + 1
        num.srcline <- src[num.ld$line1]
        substring(num.srcline, 1, zap.upto) <- spaces(max(zap.upto))
        substring(num.srcline, zap.from, nchar(num.srcline)) <- spaces(max(nchar(num.srcline)))
        names(num.srcline) <- 1:n
        ll.numbered <- split(num.srcline, num.ld$line1)
        oll <- rep("", length(src))
        rep.spots <- as.numeric(names(ll.numbered))
        line.list <- massrep(src, rep.spots, ll.numbered)
        oll <- massrep(oll, rep.spots, lapply(ll.numbered, names))
        names(line.list) <- oll
    }
    if (is.a.function) {
        if (length(preamble)) 
            line.list[length(line.list)] <- line.list[length(line.list)] %&% 
                ")"
        line.list <- c(line.list, "###### ON EXIT ######")
        n <- n + 1
        line.list <- c(line.list, structure("NULL", names = format(c(n, 
            1000))[1]))
        breakpoint[[ch(2)]] <- FALSE
    }
    expr[[1]] <- debug.mvb.subst(expr[[1]])
    ll <- names(line.list)
    ll[is.na(ll)] <- ""
    ll <- ifelse(nzchar(ll), sprintf("%4s: ", ll), spaces(6))
    ll <- paste(ll, line.list, sep = "")
    names(ll) <- names(line.list)
    class(ll) <- "cat"
    invisible(list(expr = expr[[1]], breakpoint = breakpoint, 
        line.list = ll, n = n))
}


"addnum.move.to.next.expr" <-
function( nlocal=sys.parent(), final) mlocal({
  # cat( 'AMTNE: ', i, '\n')
  final <- i[ length( i)]
  if( length( expr[[ clip( i)]]) <= i[ length( i)]) { # BACK UP 1 LEVEL
    i <- i[ -length( i)]
    if( !length( i))
return( local.return())

    add.to.last.line( suffix[ ch(i)])
    if( !is.na( this.indent <- indents[ ch( i)]))
      n.tabs <- n.tabs - this.indent
    if( is.call( expr[[ i]]) && as.character( expr[[ c( i, 1)]]) %in% c( 'for', 'while', 'repeat') )
      last.loop.tabs <- last.loop.tabs[ -length( last.loop.tabs) ]

    addnum.move.to.next.expr()
  } else { # NORMAL
    if( final==3 & as.character( expr[[ c( clip( i), 1)]])=='if' && length( expr[[ clip( i)]])==4) # "ELSE"
      line.list <- c( line.list, tabs(-1) %&% 'else ')

    i[ length( i)] <- i[ length( i)] + 1
  }
})


"an" <-
function (expr, width = options()$width, numbering = TRUE, cat.on.exit = FALSE, 
    expr.offset = 0, line.number.offset = 0, preamble = character(0)) 
{
    src <- attr(expr, "source")
    if (is.null(src)) {
        nullf <- expr
        attributes(nullf) <- list()
        environment(nullf) <- .GlobalEnv
        src <- deparse(nullf, control = "all", width.cutoff = 120L)
    }
    re.bloody.quire <- require
    using.parser <- getOption("debug.src", TRUE) && suppressWarnings(re.bloody.quire("parser"))
    if (using.parser) {
        src.expr <- parser(text = src)
        ld <- attr(src.expr, "data")
    }
    old.width <- options(width = 1000)$width
    if (cat.on.exit) 
        cat.on.exit <- expression(cat(paste(format(names(line.list)), 
            line.list, sep = ":"), sep = "\n"))[[1]]
    on.exit({
        eval(cat.on.exit)
        options(width = old.width)
    })
    tab.width <- option.or.default("tab.width", 4)
    spaces <- function(n = 1) paste(rep(" ", n), collapse = "")
    tab.sp <- spaces(tab.width)
    prefix <- ""
    tabs <- function(more.or.less = 0, exact = n.tabs + more.or.less) {
        answer <- paste(rep(tab.sp, max(exact, 0)), collapse = "") %&% 
            prefix
        prefix <<- ""
        answer
    }
    assign("[[", my.index)
    assign("[[<-", my.index.assign)
    deparse1 <- function(x) {
        x <- deparse(x, width.cutoff = 500)
        if (length(x) > 3) 
            x <- c(x[1], paste(x[-c(1, length(x))], collapse = "; "), 
                x[length(x)])
        paste(x, collapse = " ")
    }
    ch <- function(...) {
        stuff <- c(...)
        stuff[1] <- stuff[1] + expr.offset
        paste(c(stuff, "."), collapse = ",")
    }
    add.to.last.line <- function(x) if (!is.na(x)) 
        line.list[length(line.list)] <<- paste(line.list[length(line.list)], 
            x, sep = "")
    make.indent <- function(n.in = 1, loop = FALSE) {
        indents[ch(i)] <<- n.in
        n.tabs <<- n.tabs + n.in
        if (loop) 
            last.loop.tabs <<- c(last.loop.tabs, n.tabs)
    }
    if ("debug" %in% loadedNamespaces()) {
        debuggify.system.call <- function(step.check = FALSE, 
            nlocal = sys.parent()) mlocal(expr[[c(i, 1)]] <- if (step.check) 
            substitute({
                if (debug:::stepping(call.type)) debug:::debug.fun.name else orig.fun.name
            }, list(call.type = call.type, debug.fun.name = as.name("debug." %&% 
                as.character(expr[[c(i, 1)]])), orig.fun.name = expr[[c(i, 
                1)]]))
        else call(":::", quote(debug), as.name("debug." %&% as.character(expr[[c(i, 
            1)]]))))
    }
    else {
        debuggify.system.call <- function(step.check = FALSE, 
            nlocal = sys.parent()) mlocal(expr[[c(i, 1)]] <- if (step.check) 
            substitute({
                if (stepping(call.type)) debug.fun.name else orig.fun.name
            }, list(call.type = call.type, debug.fun.name = as.name("debug." %&% 
                as.character(expr[[c(i, 1)]])), orig.fun.name = expr[[c(i, 
                1)]]))
        else as.name("debug." %&% as.character(expr[[c(i, 1)]])))
    }
    default.update.line.list <- function(nlocal = sys.parent()) mlocal(line.list <- c(line.list, 
        tabs() %&% deparse1(expr[[i]])))
    if (is.a.function <- is.function(expr)) {
        line.list <- clip(deparse(do.call("args", list(expr))))
        line.list <- paste(c(line.list, preamble), collapse = " ")
        names(line.list) <- ""
        expr <- body(expr)
    }
    else line.list <- ""
    expr <- do.call("expression", list(expr))
    src.expr <- do.call("expression", list(src.expr[[c(1, length(src.expr[[1]]) - 
        1)]]))
    id <- ld$id
    numbered.id <- numeric(0)
    breakpoint <- vector("list", 0)
    n <- line.number.offset
    i <- 1
    n.tabs <- 1
    last.loop.tabs <- 0
    suffix <- structure("", names = ch(1))
    indents <- numeric(0)
    while (length(i)) {
        needs.a.number <- TRUE
        next.i <- numeric(0)
        call.above <- as.character(expr[[c(clip(i), 1)]])
        if (i[length(i)] != 1 && call.above == "switch") {
            add.to.last.line(",")
            if (i[length(i)] <= length(expr[[clip(i)]])) 
                line.list <- c(line.list, tabs(-1) %&% "'" %&% 
                  names(expr[[clip(i)]])[i[length(i)]] %&% "' = ")
            if (i[length(i)] == length(expr[[clip(i)]]) && is.name(expr[[i]]) && 
                !nzchar(expr[[i]])) 
                expr[[clip(i)]][tail(i, 1)] <- list(NULL)
        }
        if (mode(expr[[i]]) == "(" || !is.call(expr[[i]])) {
            anything.to.add <- deparse1(expr[[i]])
            if (nzchar(anything.to.add)) 
                line.list <- c(line.list, tabs() %&% anything.to.add)
        }
        else {
            call.type <- expr[[c(i, 1)]]
            if (!is.name(call.type)) 
                call.type <- "default"
            else call.type <- as.character(call.type)
            switch(call.type, `{` = {
                if ((length(i) == 1) || call.above %in% c("if", 
                  "switch", "for", "while", "repeat")) add.to.last.line(" {") else {
                  line.list <- c(line.list, tabs() %&% "{")
                  make.indent()
                }
                needs.a.number <- FALSE
                automove <- FALSE
                suffix[ch(i)] <- " }"
                if (length(expr[[i]]) == 1) expr[[i]] <- call("{", 
                  NULL)
                next.i <- c(i, 2)
            }, `(` = {
            }, `if` = {
                if (call.above != "if" || i[length(i)] != 4) {
                  line.list <- c(line.list, tabs())
                  make.indent()
                }
                add.to.last.line("if( " %&% deparse1(expr[[c(i, 
                  2)]]) %&% ")")
                next.i <- c(i, 3)
            }, switch = {
                line.list <- c(line.list, tabs() %&% "switch( " %&% 
                  deparse1(expr[[c(i, 2)]]))
                make.indent(2)
                suffix[ch(i)] <- ")"
                next.i <- c(i, 3)
            }, `for` = {
                line.list <- c(line.list, paste(tabs(), "for( ", 
                  expr[[c(i, 2)]], " in ", deparse1(expr[[c(i, 
                    3)]]), " )", sep = ""))
                make.indent(loop = TRUE)
                next.i <- c(i, 4)
            }, `while` = {
                line.list <- c(line.list, paste(tabs(), "while( ", 
                  deparse1(expr[[c(i, 2)]]), " )", sep = ""))
                make.indent(loop = TRUE)
                next.i <- c(i, 3)
            }, `repeat` = {
                line.list <- c(line.list, tabs() %&% "repeat")
                make.indent(loop = TRUE)
                needs.a.number <- FALSE
                next.i <- c(i, 2)
            }, `<-` = {
                temp.expr <- expr[[i]]
                temp.expr[[3]] <- 0
                temp.expr <- paste(deparse(temp.expr, width.cutoff = 500), 
                  collapse = " ")
                temp.expr <- substring(temp.expr, 1, nchar(temp.expr) - 
                  1)
                prefix <- prefix %&% temp.expr
                needs.a.number <- FALSE
                next.i <- c(i, 3)
            }, `break` = , `next` = {
                line.list <- c(line.list, tabs(exact = last.loop.tabs[length(last.loop.tabs)] - 
                  1) %&% call.type)
                debuggify.system.call()
            }, return = , invisible = {
                line.list <- c(line.list, deparse1(expr[[i]]))
                debuggify.system.call()
            }, with = , within = , try = , suppressWarnings = , 
                eval = , evalq = {
                  default.update.line.list()
                  debuggify.system.call(step.check = TRUE)
                }, local.on.exit = , on.exit = {
                  default.update.line.list()
                  debuggify.system.call()
                }, default.update.line.list())
        }
        if (needs.a.number) {
            n <- n + 1
            if (using.parser) {
                iddo <- attr(src.expr[[i]], "id")
                if (is.null(iddo)) {
                  parid <- attr(src.expr[[clip(i)]], "id")
                  iddo <- ld$id[ld$parent == parid & ld$token.desc == 
                    "NULL_CONST"][1]
                  if (is.na(iddo)) 
                    stop(sprintf("Can't mtrace source; no valid id attr at expr[[c(%s)]]-- try with debug.src=FALSE", 
                      paste(i, collapse = ",")))
                }
                numbered.id[n] <- iddo
            }
            names(line.list)[length(line.list)] <- format(c(n, 
                1000))[1]
            names(line.list)[is.na(names(line.list))] <- ""
            if (numbering) 
                breakpoint[[ch(i)]] <- is.a.function && (length(breakpoint) == 
                  0)
        }
        if (length(next.i)) 
            i <- next.i
        else addnum.move.to.next.expr()
    }
    if (using.parser) {
        num.ld <- ld[match(numbered.id, id, 0) %such.that% (. > 
            0), ]
        pop <- lapply(1:n, function(ni) {
            popi <- num.ld$id[ni]
            while (popi[1] != 0) popi <- c(ld$parent[ld$id == 
                popi[1]], popi)
            return(popi)
        })
        cpop <- list(numeric(0))
        for (i in 2 %upto% n) cpop[[i]] <- unique(c(cpop[[i - 
            1]], pop[[i - 1]]))
        nonfirst <- index(duplicated(num.ld$line1))
        nonlast <- index(rev(duplicated(rev(num.ld$line1))))
        zap.upto <- rep(0, n)
        zap.from <- nchar(src[num.ld$line1]) + 1
        for (izap in nonfirst) {
            ancid <- (pop[[izap]] %except% cpop[[izap]])[1]
            zap.upto[izap] <- ld$col1[ld$id == ancid]
        }
        zap.from[nonlast] <- zap.upto[nonlast + 1] + 1
        num.srcline <- src[num.ld$line1]
        substring(num.srcline, 1, zap.upto) <- spaces(max(zap.upto))
        substring(num.srcline, zap.from, nchar(num.srcline)) <- spaces(max(nchar(num.srcline)))
        names(num.srcline) <- 1:n
        ll.numbered <- split(num.srcline, num.ld$line1)
        oll <- rep("", length(src))
        rep.spots <- as.numeric(names(ll.numbered))
        line.list <- massrep(src, rep.spots, ll.numbered)
        oll <- massrep(oll, rep.spots, lapply(ll.numbered, names))
        names(line.list) <- oll
    }
    if (is.a.function) {
        if (length(preamble)) 
            line.list[length(line.list)] <- line.list[length(line.list)] %&% 
                ")"
        line.list <- c(line.list, "###### ON EXIT ######")
        n <- n + 1
        line.list <- c(line.list, structure("NULL", names = format(c(n, 
            1000))[1]))
        breakpoint[[ch(2)]] <- FALSE
    }
    expr[[1]] <- debug.mvb.subst(expr[[1]])
    ll <- names(line.list)
    ll[is.na(ll)] <- ""
    ll <- ifelse(nzchar(ll), sprintf("%4s: ", ll), spaces(6))
    ll <- paste(ll, line.list, sep = "")
    names(ll) <- names(line.list)
    class(ll) <- "cat"
    invisible(list(expr = expr[[1]], breakpoint = breakpoint, 
        line.list = ll, n = n))
}


"backtrack.to.loop" <-
function( expr, i) {
  i.try <- clip( i)
  while( length( i.try) && (!is.call( expr[[ i.try]]) ||
      !(paste( as.character( expr[[ c( i.try, 1)]]), collapse=' ') %in% c( 'for', 'while', 'repeat')) ) )
    i.try <- clip( i.try)

  i.try
}


"bp" <-
function(line.no, expr = TRUE, fname) do.in.envir( envir=find.debug.HQ( TRUE), { # NB can maybe call "bp" while not debugging
  .system. <<- TRUE
  expr <- do.call( 'substitute' , list( substitute( expr), list( F=FALSE))) # so displayed version will have no star

  repeat { # only to allow break
    # Get info on function being called
    if( missing( fname)) {
      this <- .frames.$actual[ nrow( .frames.)]
      if( !length( this)) {
        cat( 'bp: don\'t know which function to set breakpoint in\n')
  break }

      fname <- .frames.$function.name[ nrow( .frames.)]
      max.line.no <- length( get( 'breakpoints', envir=sys.frame( .frames.$debug[ nrow( .frames.)])))
      general <- line.no <= length( tracees[[ fname]]$breakpoint) # else it's in the 'on.exit' section-- specific to incarnation
    } else {
      general <- TRUE
      max.line.no <- length( tracees[[ fname]]$breakpoint) }

    if( is.null( tracees[[ fname]])) {
      cat( 'bp: no trace info for', fname, '\n')
  break }

    if( !(line.no %in.range% c( 1, max.line.no))) {
      cat( 'bp: out-of-range line number', line.no, 'for', fname, '\n')
  break }

    if( nrow( .frames.)) # ie live
      for( i in .frames.$debug[ .frames.$function.name==fname])
        set.a.breakpoint( expr, line.no, frame.number=i)

    if( general)
      tracees[[ fname]]$breakpoint[[line.no]] <<- expr # used to be substitute( expr)-- surely not??

  break
  }

  .nothing.
})


"check.for.tracees" <-
function( where=1) {
  o <- find.funs( where)
  if( !length( o))
return( character( 0))

  where <- as.environment( where)
  is.tracee <- function( x) {
    x <- get( x, envir=where)
    idx <- c(2,2,1)
    repeat{
      if( !my.index.exists( idx, body( x)))
  return( FALSE)

      # Can run into terrible grief if this is "missing" (0-length name)-- won't assign
      if( is.name( bod <- my.index( body(x), idx)) && !nzchar( as.character( my.index( body(x), idx))))
  return( FALSE)

      if( identical( quote( mlocal), bod) && length( idx)==3)
        idx <- c( 2, 2, 2, 1) # try another
      else
  return( identical( quote( debug:::evaluator), bod))
    }
  }

  o[ sapply( o, is.tracee)]
}


"check.legality" <-
function( thing, call.type) do.in.envir( envir=find.debug.HQ( FALSE), {
# Trap non-logical first arguments to "if", "while", and 
# ...non-subsettable arguments to "for" e.g. for( i in call( 'abc'))
# cat( 'Checking legality in call.type', call.type, 'of', thing, '\n')

# Can issue warning, or cause error, or neither, depending on getOption( 'warn')

  if( call.type %in% c( 'if', 'while') &&
      ( (typeof( thing) %in% dodgy.if.while.types || is.na( as.logical( thing)[ 1])) ||
        ( (length( thing) > 1) && 
        (try( list( eval( substitute( if( thing) TRUE), envir=parent.frame()))) %is.a% 'try-error'))))
    message <- 'illegal if/while test'            
  else if( call.type == 'for' && typeof( thing) %in% dodgy.for.counter.types)
    message <- 'illegal for-loop counter'
  else if( call.type == 'switch' && ( ! (typeof( thing) %in% valid.switch.types) || length( thing) != 1 ) )
    message <- 'illegal switch control argument'
  else
    return( TRUE)

  structure( .Data=FALSE, message=message)
})


"debug.break" <-
function() do.in.envir( envir=sys.frame( find.active.control.frame()), {
  .system. <<- TRUE
  .print.result. <<- FALSE

  i.try <- backtrack.to.loop( expr, i)
  if( !length( i.try))
stop( 'Not in a loop!')

  .evaluated.OK. <<- TRUE
  i <<- i.try
  move.to.next.expression( sorted.out=FALSE, nlocal=find.active.control.frame())

# Slightly weird construction: tell it to skip to itself!
  .skip. <<- TRUE
  .skipto. <<- i
  j # return previous value
})


"debug.C" <-
function( ...){
  mc <- as.list( match.call( expand.dots=TRUE))
  for( i in seq_along( mc)[ nzchar( names( mc))]) {
    scatn( '**** %s ****', names( mc)[i])
    print( str( eval( mc[[i]], parent.frame())))
  }
NULL
}


"debug.do.call" <-
function( what, args) {
  funs.to.replace <- named( cq( break, next, return, on.exit, sys.on.exit)) # should this include "do.call"?!
  funs.to.replace[] <- 'debug.' %&% funs.to.replace
  new.what <- funs.to.replace[ what]
  if( !is.na( new.what))
    what <- new.what
  mc <- match.call()
  mc$what <- what
  mc[[1]] <- do.call
  eval( mc, parent.frame())
}


"debug.eval" <-
function( expr, envir = parent.frame(), enclos = if (is.list(envir) || 
    is.pairlist(envir)) parent.frame() else baseenv()) {
  mc <- mvb.match.call()
  mc[[1]] <- quote( debug:::debug.eval.guts)
  mc$expr <- call( 'quote', expr)
  eval( mc, mvb.parent.frame())
}


"debug.eval.guts" <-
function( expr, envir = parent.frame(), enclos = if (is.list(envir) || 
    is.pairlist(envir)) parent.frame() else baseenv(), fun.name='eval( ..., expr={') {
  f.eval <- function( nlocal=sys.parent()) 9
  body( f.eval) <- call( 'mlocal', expr)
  
  # Used to have find.debug.HQ( FALSE), assuming this is called inside a function being debugged
  # Now changed to TRUE to allow direct command-line use: e.g. debug.eval( myexpr)
  # ... hopefully without side-effects...
  
  eval.name <- tail( find.debug.HQ( TRUE)$.frames.$window.name, 1) %&% '.' %&% 
      sub( '[(].*', '', fun.name)
  assign( eval.name, f.eval)
  mtrace( char.fname=eval.name)
  on.exit( try( mtrace( char.fname=eval.name, tracing=FALSE), silent=TRUE))

  # Is it worth a separate window?
  mc <- mvb.match.call()  
  if( tracees[[ eval.name]]$n > 2) {
    tracees[[ eval.name]]$line.list[ 1] <<- fun.name
    force( envir)
    force( enclos)
  
    mc$expr <- call( eval.name)
    ans <- eval( substitute( get( eval.name, envir=sys.frame( n))(), 
        list( eval.name=eval.name, n=sys.nframe())), envir=envir, enclos=enclos)
  } else {
    mc[[1]] <- quote( eval)
    mc$fun.name <- NULL
    ans <- eval( mc, parent.frame())
  }
  
return( ans)
}


"debug.evalq" <-
function( expr, envir = parent.frame(), enclos = if (is.list(envir) || 
    is.pairlist(envir)) parent.frame() else baseenv()) {
  mc <- mvb.match.call()
  mc[[1]] <- quote( debug:::debug.eval.guts)
  mc$fun.name <- 'evalq( ..., expr={'
  mc$expr <- call( 'quote', substitute( expr))
  eval( mc, mvb.parent.frame())
}


"debug.invisible" <-
function( x) do.in.envir( envir=sys.frame( find.active.control.frame()), {
  .print.result. <<- FALSE
  if( .step.)
    scatn( '<invisible>', file=debug.catfile())
  if( missing( x))
    base::invisible()
  else
    base::invisible( x)
})


"debug.local.on.exit" <-
function( new.expr, add=FALSE) # avoid duplicating code of 'local.on.exit'-- not easy!
  do.call( 'do.in.envir', list( envir=sys.frame( find.active.control.frame()), 
      fbody= match.call( do.in.envir, call=body( debug.on.exit))$fbody))


"debug.mvb.subst" <-
function( expr) {
  sublist <- named( cq( nargs, sys.call, sys.parent, sys.function, sys.nframe, 
      parent.frame, eval.parent, match.call))
  # Next from when debug Depended on mvbutils, rather than Importing it, so these mvb.... funs were on search path
  # sublist[] <- 'mvb.' %&% sublist
  #  sublist <- lapply( sublist, as.name)
  sublist <- lapply( sublist, function( x) 
      call( ':::', quote( mvbutils), as.name( 'mvb.' %&% x)))
      
  # Next added because I moved mvb.sys.on.exit into debug--- it was in mvbutils but had lots of debug stuff.
  # But what about debug.sys.on.exit??! Which is "right"?
  # Well, who really cares-- it works...
  sublist <- c( sublist, list( sys.on.exit=quote( debug:::mvb.sys.on.exit))) 
  sublist <- c( sublist, list( Recall=quote( debug:::debug.Recall)))
  expro <- substitute( substitute( e, sublist), list( e=expr, sublist=sublist))
  eval( expro)
}


"debug.next" <-
function() do.in.envir( envir=sys.frame( find.active.control.frame()), {
  .system. <<- TRUE
  .print.result. <<- FALSE

  i.try <- backtrack.to.loop( expr, i)
  if( !length( i.try))
stop( 'Not in a loop!')

  .evaluated.OK. <<- TRUE
  i.try <- c( i.try, length( expr[[ i.try]]))
  i <<- i.try # from debug.break
  move.to.next.expression( sorted.out=FALSE, nlocal=find.active.control.frame()) # from debug.break

  .skip. <<- TRUE
  .skipto. <<- i
  j # return previous value
})


"debug.on.exit" <-
function( new.expr, add=FALSE) do.in.envir( envir=sys.frame( find.active.control.frame()), {
  if( i[1]>1)
stop( "Can't usefully call 'on.exit' while 'on.exit' is running!") # and return to D()> prompt

  if( !add)
    old <- list()
  else {
    old <- as.list( expr[[ 2]])
    if( !identical( old, list()))
      old <- old[ -c( 1, length( old))] # remove the starting brace and the final NULL
  }

  if( missing( new.expr))
    new.expr <- list()
  else
    new.expr <- substitute( new.expr)

  # Next line clunkily adds new.expr to existing on.exit code, all inside braces
  # Spare NULL at the end is to allow a stop after on.exit code has run
  new.expr <- as.call( c( list( as.name( '{')), old, new.expr, list( NULL)))
#  new.expr <- call( '{', old, new.expr, list( NULL))

  # Could redo the entire function here, but that can be slow
  adn <- add.numbers( new.expr, line.number.offset=orig.breakpoint.length-1, expr.offset=1)
  prev.line.list.length <- length( line.list)
  line.list <<- c( line.list[ 1:(orig.line.list.length-1)], adn$line.list)

  if( stopped.yet)
    .update.window.with.on.exit()

  # Check if any breakpoints were set in previous on.exit code. If so, set bp at start of new on.exit
  set.bp <- !all( sapply( FUN=identical, y=FALSE, breakpoints[ -(1:(orig.breakpoint.length-1))]))
  breakpoints <<- c( breakpoints[1:(orig.breakpoint.length-1)], adn$breakpoint)
  if( set.bp)
    set.a.breakpoint( TRUE, orig.breakpoint.length, find.active.control.frame())
  
  expr[[2]] <<- adn$expr # which has sorted out any breaks, null braces, etc.
  .nothing.
})


"debug.q" <-
function() do.in.envir( envir=find.debug.HQ( FALSE), {
#  .quit.debug. <<- TRUE
  cat( "To quit the debugger, type 'qqq()'\n")
  .nothing.
})


"debug.Recall" <-
function(...){
  mc <- match.call( expand.dots=TRUE)
  fh <- debug:::find.debug.HQ( FALSE)
  w <- which.max( fh$.frames.$actual)
  mc[[1]] <- as.name( fh$.frames.$function.name[ w])
  # print( mc)
  eval( mc, parent.frame())
}


"debug.return" <-
function( ... ) do.in.envir( envir=sys.frame( find.active.control.frame() ), {
  if( i[1]>1)
stop( "Can't \"return\": no function to return from because \"on.exit\" code is executing")
  orig.mc <- mc <- as.list( match.call())[ -1]
  
  if( length( mc)) {
    if( length( mc)==1)
      mc <- eval( mc[[1]], envir=frame)
    else { # multiple arguments, so return as named list
      if( is.null( names( mc))) {
        which <- rep( TRUE, length( mc))
        names( mc) <- rep( '', length( mc))
      } else
        which <- names( mc)==''

      for( i in index( which))
        if( is.symbol( orig.mc[[ i]]))
          names( mc)[ i] <- as.character( orig.mc[[ i]] )
      mc <- lapply( mc, eval, envir=frame)
    }
  } else
    mc <- NULL

# Admin comes last, in case there's a crash in setting "mc"
  .skipto. <<- 2 # start the 'on.exit'
  .skip. <<- TRUE
#  .step. <<- FALSE # no hanging about at the return-value breakpoint
  .evaluated.OK. <<- TRUE

  mc
})


"debug.retval" <-
function() get( 'j', envir=sys.frame( find.active.control.frame()))


"debug.suppressWarnings" <-
function( expr) {
  mc <- mvb.match.call()
  mc[[1]] <- quote( debug:::debug.eval.guts)
  mc$fun.name <- 'suppressWarnings({'
  mc$expr <- call( 'quote', substitute( expr))
  suppressWarnings( eval( mc, mvb.parent.frame()))
}


"debug.sys.on.exit" <-
function() do.in.envir( envir=sys.frame( find.active.control.frame()), {
  ex <- expr[[ 2]]
  if( length( ex)>1) { # it's been set via "debug.on.exit"
    ex <- as.list( ex)
    ex <- ex[ -c( 1, length( ex))] # brace & NULL
    if( length( ex)==1)
      ex <- ex[[ 1]]
    else
      ex <- as.call( c( list( as.name( '{')), ex))
  }
  ex
})


"debug.try" <-
function( expr, silent=FALSE) {
  mc <- mvb.match.call()
  mc[[1]] <- quote( debug:::debug.eval.guts)
  mc$silent <- NULL
  mc$fun.name <- 'try( silent=..., {'
  mc$expr <- call( 'quote', substitute( expr))
  eval( mc, mvb.parent.frame())
}


"debug.with" <-
function( data, expr, ...){
  # Don't bother if non-default
  if( any( ('with.' %&% class( data)) %in% methods( 'with')))
UseMethod( 'with') # auto-returns

  f.with <- function(...) 9
  body( f.with) <- substitute( expr)
  e <- new.env( parent=parent.frame())
  for( i in names( data))
    e[[ i]] <- data[[i]]
  environment( f.with) <- e
  
  with.name <- tail( find.debug.HQ( FALSE)$.frames.$window.name, 1) %&% '.with'
  assign( with.name, f.with)
  mtrace( char.fname=with.name)
  on.exit( try( mtrace( char.fname=with.name, tracing=FALSE), silent=TRUE))
  
  if( tracees[[ with.name]]$n > 2) {
    tracees[[ with.name]]$line.list[ 1] <<- 'with( data, {'
return( do.call( with.name, list( ...)))
  } else {
    mc <- match.call( expand.dots=TRUE)
    mc[[1]] <- quote( with)
return( eval( mc, parent.frame()))
  }
}


"debug.within" <-
function( data, expr, ...){
  if( !is.list( data)) # within has methods for data.frames (which are lists) and lists
UseMethod( within) # crash if no method; if there is, debug doesn't know how to handle it

  # Create a dummy function to execute 'expr
  f.within <- function(...) 9
  body( f.within) <- substitute( expr)
  e <- new.env( parent=parent.frame())
  for( i in names( data))
    e[[ i]] <- data[[i]]
  environment( f.within) <- e
  
  within.name <- tail( find.debug.HQ( FALSE)$.frames.$window.name, 1) %&% '.within'
  assign( within.name, f.within)
  mtrace( char.fname=within.name)
  # Don't stop for single statements, except on error
  if( tracees[[ within.name]]$n <= 2) 
    bp( 1, FALSE, fname=within.name)
  on.exit( try( mtrace( char.fname=within.name, tracing=FALSE), silent=TRUE))
  
  # Now create a dummy 'within' function, that calls 'debug.eval' instead of 'debug'...
  # ... NB *without* checking whether we are stopping for evals, which would be the norm
  
  # Could just hardwire 'within...' code here, but hopefully future-proof it by...
  # ... substituting the 'eval' call. Hopefully only one of them, gulp.
  withindef <- if( is.data.frame( data)) within.data.frame else within.list
  body( withindef) <- do.call( 'substitute', list( body( withindef), 
    list( eval=quote( debug:::debug.eval))))

  within.wrapper.name <- 'withindef.' %&% within.name
  assign( within.wrapper.name, withindef)
  mtrace( char.fname=within.wrapper.name)
  
  # Inside the wrapper, stop only at the end, unless there's an error
  tracees[[ within.wrapper.name]] <<- within( tracees[[ within.wrapper.name]], {
    # Tidy up display
    oline <- grep( 'debug:::debug.eval', line.list, fixed=TRUE)[1]
    line.list[ oline] <- sub( 'debug:::debug.eval', 'eval', line.list[ oline], fixed=TRUE)
    breakpoint[[1]] <- quote( {go(); FALSE}) # put into go-mode at start and...
    breakpoint[[1+sum(nzchar(names(line.list)[1:oline]))]] <- 
        quote( {go(); FALSE}) # ... and after eval 
    breakpoint[[ length( breakpoint)]] <- TRUE #... back to step-mode at the end.
    rm( oline)
  })
  on.exit( try( mtrace( char.fname=within.wrapper.name, tracing=FALSE), silent=TRUE), add=TRUE)  

  mc <- match.call(expand.dots=TRUE)
  mc[[1]] <- get( within.wrapper.name)
return( eval( mc, parent.frame()))
}


"enact.command.r" <-
function( command, frame) do.in.envir( envir=find.debug.HQ( FALSE), {
  .evaluated.OK. <<- FALSE
  .system. <<- FALSE
  .print.result. <<- TRUE

  if( command=='') {
    .evaluated.OK. <<- TRUE
    .system. <<- TRUE
return() }

  .skip. <<- FALSE
  command <- try( list( parse( text=command, srcfile=NULL)))
  if( command %is.not.a% 'try-error') {
    command <- command[[1]] # unwrap list in try
    if( length( command)) {
      command <- command[[1]] # parse returns expression(...)
      command <- do.call( 'substitute', list( command, list.of.command.subs))
      command <- debug.mvb.subst( command)
    #  print( command)
    #  cat( 'Mode=', mode( command), '\n')
      command <- try2( list( eval( command, envir=frame)))

      if( missing( command))
        cat( '<missing>\n', file=debug.catfile())
      else if( command %is.not.a% 'try-error') 
        printIfSmall( command[[1]], ofile=debug.catfile()) # unwrap list from 'try'
      else
        .evaluated.OK. <<- FALSE # paranoid safety net; eOK _might_ have been set to TRUE before a crash!
    } else # length-0, presumably from line starting with a hash
      .system. <<- TRUE # so retval doesn't get set
  } # else try-error will be picked up on return
  
return( command) # to be unwrapped from 'try' list by 'evaluator'
})


"enter.on.exit" <-
function() get( 'orig.breakpoint.length', envir=sys.frame( find.active.control.frame()))


"eval.bp" <-
function(ex, envir) {
  break.time <- try2( eval( ex, envir=envir))
  if( break.time %is.a% 'try.error') {
    cat( '\nInvalid breakpoint expression!\n')
    break.time <- TRUE }

  if( length( break.time) != 1)
    break.time <- TRUE
  else
    break.time <- as.logical( break.time)
    
  if( is.na( break.time)) {
    cat( '\nBreakpoint evaluates to NA!\n')
    break.time <- TRUE
  }

  break.time
}


"eval.catching.errors" <-
function(i, envir) do.in.envir( envir=find.debug.HQ( FALSE), {
#  if( option.or.default( 'erk', FALSE)) {
#    cat( "Warn: ", getOption( 'warn'), '\n')
#    print( i)
#  }

  # Thanks to Luke Tierney for this trick, which avoids "restart"
  j <- try2( list( value=eval( i, envir=envir))) # try2 traps interrupts as well as errors
  .evaluated.OK. <<- j %is.not.a% 'try-error'
  if( .evaluated.OK.) {
    j <- j$value
    if( missing( j)) 
return( formals( evaluator)$fname) # avoid immediate trouble
    # else not missing
return( j) } 
  else
return( NULL)
})


"evaluator" <-
function( fname) do.in.envir( envir=find.debug.HQ( TRUE), {
  on.exit( .end.incarnation())
  next.incarnation()

# Main loop through the expressions
  repeat {
    if( i[1]==2 && in.body.code)
      retval <- j # we have finished the real body of the function, and will return this value
    else if( i[ 1]==3)
break # and then return

    if( i[ 1]!=1)
      in.body.code <- FALSE

    ch.i <- ch( i)
    .skip. <<- FALSE 
    .skipto. <- 0
    
    repeat { # Figure out if user input is possible, and if so keep processing commands until told to continue debugger
      if( .quit.debug.) { # set by user's call to q(), mapped to 'q.debug'; this may be inside child
        cat( '\rNo ', file=debug.catfile())
#        opt <- options(show.error.messages=FALSE)
#        cat( 'SEM=', opt$show.error.messages, '\n', file=debug.catfile())
#       on.exit( options(opt), add=TRUE)

        # If func is being called inside 'try', 'stop' alone won't quit...
        # ... so make up a fake error
        stoppo <- simpleError( "merely quitting mvb's debugger")
        class( stoppo) <- c( 'stoppo.debug', 'condition')
stop( stoppo)      
        #stop( 'merely quitting mvb\'s debugger', call.=FALSE)
      } # if .quit.

      .evaluated.OK. <<- TRUE
      .print.result. <<- FALSE
      find.line <- match( ch.i, names( breakpoints))
      if( !is.na( find.line))
        lno <- find.line
      
      if( !stop.here() || !.step.)
    break

      if( !stopped.yet) { # debug windows aren't launched until we actually stop for input...
        # ... maybe never if no bp's
        # NB must also launch mtrace'd parents that didn't have bp's and so haven't displayed yet
        launch.debug.windows() # win=sdebug.window.name, fun=fname)
        stopped.yet <- TRUE }
      command <- interact()
      try.j <- enact.command.r( command, frame)
      if( try.j %is.not.a% 'try-error' && !.system.)
        j <- try.j[[1]]
      if( .evaluated.OK.) # mostly set to F even with valid command, to force repeat
    break

    } # user commands

    if( !.skip.) { # deal with next debuggee statement
#cat( ch( i), '  ') # these statements can be dehashed if you want to see which statements get looked at
      call.type <- get.call.type( expr[[ i]])
#      afei <- augment.for.eval( i, call.type)
#      if( !identical( afei, i))
#cat( ch( afei))
#cat( '\n')

      if( call.type %in% c( 'normal', 'if', 'for', 'while', 'switch') ) { # alternatives are {} break next
        try.j <- eval.catching.errors( expr[[ augment.for.eval( i, call.type) ]], envir=frame)
        if( .evaluated.OK.) {
          if( missing( try.j)) {
            j <- formals( evaluator)$fname # missing val
            if( .step.)
              cat( '<missing>\n', file=debug.catfile())
          } else {
            j <- try.j
            if(.print.result.)
              printIfSmall(j, ofile=debug.catfile())
          }
          .evaluated.OK. <<- check.legality( j, call.type) # illegal arg to if/while/for/switch
          if( !.evaluated.OK.)
            cat( 'Problem:', attr( .evaluated.OK., 'message'), '\n', file=debug.catfile())
        }

        if( !.evaluated.OK.) { # back to outer loop and user's commands
          .step. <<- TRUE
  next } # jump straight back to getting user input
      } # not normal/if/for/while/switch
    } # not user-driven skip

    if( .skip.) # user OR next/break
      skipto.debug()
    else {
      move.to.next.expression()
      if( !.evaluated.OK.) { # only possible here if assignment failed
        .step. <<- TRUE
  next
      }
    } # if not .skip.

  } # master loop

  if( in.body.code) { # I think this can only happen if the user intervenes heavily with "skip"
    cat( 'Function exited via "skip": return value may be strange\n', file=debug.catfile())
    retval <- j }

  retval
})


"exit.on.exit" <-
function() length( get( 'breakpoints', envir=sys.frame( find.active.control.frame())))


"find.active.control.frame" <-
function() {
  dhq <- find.debug.HQ( FALSE)
  .frames. <- get( '.frames.', envir=dhq)
  .frames.[ nrow( .frames.), 'debug']
}


"find.debug.HQ" <-
function( create.if.not=TRUE ) {
  n.debug.HQ <- index( sapply( sys.frames(), 
      function( x) !is.null( attr( x, 'I.am.the.debug.HQ'))))
  if( length( n.debug.HQ))
    debug.HQ <- sys.frame( n.debug.HQ[ 1])
  else if( create.if.not) {
#    cat( 'Setting up debug HQ: frame=', sys.nframe(), '\nparent=', sys.parent(), 
#        '\nparent.frame=')
#    print( sys.frame( sys.parent()))
#    cat( '\nparent.frame contents=')
#    print( ls( sys.frame( sys.parent())))

    ## Needs changing to make this a separate inheriting frame ##
    
    debug.HQ <- parent.frame()
#cat( "Setting debug HQ:\n")
#cat( 'sys.parent=', sys.parent())
#print( sys.frame( sys.parent()))
#print( debug.HQ)
#cat( lsall( debug.HQ), '\n', sep=' ')
    
    attr( debug.HQ, 'I.am.the.debug.HQ') <- TRUE
#cat( 'after attr:\n')
#print( debug.HQ)
    n.debug.HQ <- sys.parent()
    setup.debug.admin( nlocal=n.debug.HQ) 
#cat( 'after admin:\n')
#print( debug.HQ)
    
  } else # doesn't exist, not supposed to create it
return( FALSE)

  # Create an "active copy" of 'tracees' in debug.HQ
  
  if( exists( 'tracees', envir=debug.HQ)) 
    rm( tracees, envir=debug.HQ)
  f <- function( val) if( missing( val)) tracees else tracees <<- val
  environment( f) <- asNamespace( 'debug')
  makeActiveBinding( 'tracees', f, debug.HQ)

  debug.HQ
}


"find.from" <-
function (char.fname, from = mvb.sys.parent(), look.for.generics = TRUE) 
{
    if (typeof(from) == "closure") 
        from <- environment(from)
    else if (is.numeric(from)) 
        from <- if (from > 0) 
            sys.frames()[[from]]
        else .GlobalEnv
    orig.from <- from
    repeat {
        found <- exists(char.fname, envir = from, inherits = FALSE)
        if (found || is.null(from)) 
            break
        from <- parent.env(from)
    }
    if (!found && look.for.generics) {
        if (length(grep("\\.", char.fname))) {
            parts <- strsplit(char.fname, "\\.")[[1]]
            for (i in 2 %upto% length(parts)) {
                gen <- paste(parts[1:(i - 1)], collapse = ".")
                ff <- find.from(gen, orig.from, look.for.generics = FALSE)
                if (is.logical(ff) || !exists(".__S3MethodsTable__.", 
                  ff, inherits = FALSE)) 
                  next
                S3 <- get(".__S3MethodsTable__.", ff)
                if (!exists(char.fname, envir = S3, inherits = FALSE)) 
                  next
                found <- S3
            }
        }
        from <- found
    }
    from
}


"find.S3.dispatch" <-
function( obj, gen){
  # Work out which method will be called
  poss <- gen %&% '.' %&% c( class( obj), 'default')
  mm <- match( poss, methods( gen), 0)
  poss[ mm>0][1] # NA if missing
}


"fun.locator" <-
function (fname, from = .GlobalEnv, mode = "function") {
######################
# Formatting got lost at some point... this is pretty ancient code
  if( typeof(from) == "closure") {
    from <- environment(from)
  } else if( is.numeric(from)) {
    from <- (if (from > 0) 
        sys.frame(from)
      else 
        .GlobalEnv)
  } else if( isS4(from) && is.environment(from)) {
    try( do.call( "$", list(from, as.symbol(fname))), silent = TRUE)
    return( if( exists( fname, from, inherits = FALSE)) list(from) else list())
  }
  
  is.here <- function(env) 
      exists(fname, envir = env, inherits = FALSE, mode = mode) && 
      (!length(ff) || identical(env[[fname]], ff[[1]][[fname]]))

  orig.from <- from
  ff <- list()
  search.envs <- lapply(search(), as.environment)

  while (!any(sapply(search.envs, identical, y = from))) {
    if (is.here(from)) 
      ff <- c(ff, from)
    if (environmentName(from) == environmentName(emptyenv())) 
  break
    from <- parent.env(from)
  }
  
  for (se in search.envs) {
    if (is.here(se)) 
      ff <- c(ff, list(se))
  }
  
  ln <- lapply(loadedNamespaces(), asNamespace)
  ln <- ln[!sapply(ln, identical, y = orig.from)]
  for (lni in ln) {
    if (is.here(lni)) 
      ff <- c(ff, list(lni))
    if (is.here(parent.env(lni))) 
      ff <- c(ff, list(parent.env(lni)))
  }
  
  S3 <- lapply(ln, function(x) 
      if (exists(".__S3MethodsTable__.", x, inherits = FALSE)) 
        x$.__S3MethodsTable__.
      else 
        0
    )
  S3 <- S3[ !sapply(S3, is.numeric)]
  for (S3i in S3) if (is.here(S3i)) 
    ff <- c(ff, list(S3i))
    
  ff
}


"get.mtraced.callers" <-
function(){
  # Return CONTROL frames of callers that are mtrace'd
  splist <- sys.parent( 1)
  while( tail( splist, 1)>0) {
    splist <- c( splist, sys.parent( length( splist)+1))
    if( tail( splist, 1) %in% clip( splist)) # found already, e.g. if nlminb is calling
  break
}
  frames <- find.debug.HQ( FALSE)$.frames.
return( evalq( debug[ actual %in% splist], frames))
}


"get.retval" <-
function() do.in.envir( envir=sys.frame( find.active.control.frame()), 
 j
)


"go" <-
function(line.no) do.in.envir( envir=sys.frame( find.active.control.frame()), {
  set.global.debug.vars( .system.=TRUE)

  if( missing( line.no))
    temp.bp <<- ''
  else if( !is.numeric( line.no) || line.no[1] < 1) {
    cat("Go how far?\n")
return( .nothing.) }
  else {
    line.no <- trunc(line.no[1])
    l <- length( breakpoints)
    if( line.no>l ) {
      cat("Max. line number=", l, "\n")
return( .nothing.) }

    temp.bp <<- names(breakpoints)[line.no]
#    cat("Temp bp set at", temp.bp, "***\n")
  }

  set.global.debug.vars( .step.=FALSE, .evaluated.OK.=TRUE)
  .nothing.
})


"interact" <-
function( nlocal=sys.parent(), input, i) mlocal({
  for( i in get.mtraced.callers())
    .update.debug.window(nlocal=i)

  # These just don't work-- very-long-standing R tcl loop bug
  # .Tcl( 'update idletasks') # maybe now...
  # .Tcl( 'update') # maybe now...
  
  
  cat( '\nD(' %&% frame.number %&% ')> ', file=debug.catfile())
  input <- try2( readLines( n=1, ok=FALSE))
  if( (!missing( input)) && (input %is.not.a% 'try-error')) {
    if( nzchar( input) && debug.command.recall) {
      cat( input, '\n', sep='', append=TRUE, file=debug.hist.file)
      loadhistory( debug.hist.file) 
    } # if input was adequate
  } else
    input <- 'stop( "Bad input!")' # ctrl-Z in windows, for example
  
  input
})


"is.mtraced" <-
function( f){
  if( is.null( f))
return( FALSE) # mtrace( <<function.that.is.out.of.scope>>, FALSE)

  if( is.function( f))
    f <- list( f)
  else if( is.character( f)) {
    nope <- !sapply( f, exists, where=1)
    flist <- vector( 'list', length( f))
    flist[ !nope] <- lapply( f[!nope], get, pos=1)
    f <- flist
  }

  if( !is.list( f)){
    print( f)
stop()}
stopifnot( is.list( f))

  nope <- sapply( f, function( x) length( body( x))<4)
  if( !all( nope)) {
    bf2 <- sapply( f[ !nope], function( x) paste( deparse( body( x)[[2]]), collapse=' '))

    xret <- '^ *return\\( *'
    xeval <- '(debug:::)?evaluator\\('
    nomatch <- rep( TRUE, length( bf2))
    for( midbodi in c( '', 'mlocal\\( *', 'do.in.envir\\( *envir *=.*, *fbody *= *'))
      nomatch <- nomatch & regexpr( xret %&% midbodi %&% xeval, bf2) < 0

    nope[ !nope][ nomatch] <- TRUE
  }

  !nope
}


"last.try.error" <-
function() structure( geterrmessage(), class='try-error')


"launch.debug.windows" <-
function() do.in.envir( envir=find.debug.HQ( FALSE), {
  for( i.win in index( !.frames.$has.window.yet)) { # i.e. not launched yet
#  cat( 'Launching', .frames.$window.name[ i.win], '\n')
    if( !is.null( debug.first.window.hook <- getOption( 'debug.first.window.hook')))
      debug.first.window.hook()

    setup.tcltk.in.control.frame( title=.frames.$window.name[ i.win], nlocal=.frames.$debug[ i.win]) 
    
    if( getOption( 'debug.window.hack', FALSE) && R.version$os=='mingw32' && exists( 'set.window.state') && 
        is.loaded( 'enumerate_running')) { 
      # Circumvent TCL/TK bug that makes windows appear without text. Hopefully it's Windows-specific
      win.num <- windows.running( .frames.$window.name[ i.win])
      for( ij in c( 2, 1, 5))
        set.window.state( win.num, ij) 
      set.focus.win( r.window.handle)
    }
        
    # Record successful launch
    .frames.$has.window.yet[ i.win] <<- TRUE
    
    if( !is.null( debug.post.window.launch.hook <- getOption( 'debug.post.window.launch.hook')))
      debug.post.window.launch.hook()
  } # for i.win in not-launched-yet

#  cat( 'Launched OK\n')
})


"make.locs" <-
function( namespace.dest=NULL) {
  # returns list of command subs
  list.of.command.subs <- named( cq( break, next, return, invisible, q, on.exit, local.on.exit, 
      sys.on.exit, do.call))
  where.to.look <- 'package:debug'
  if( where.to.look %!in% search()) # 'debug' has been 'cd'ed up
    where.to.look <- match( 'debug', 
        sapply( 1:length( search()), 
            function( i) (names( attr( pos.to.env( i), 'path')) %&% '')[1]))

  if( exists( 'debug.next', where.to.look)) # normal
    subfun <- function( x) as.name( 'debug.' %&% x)
  else { # NAMESPACE 
    subfun <- function( x) call( ':::', quote( debug), as.name( 'debug.' %&% x))

    if( !is.null( namespace.dest)) {
      # Copy unexported functions into debug.HQ
      # This is really a pretty dreadful way of doing things in modern R, and...
      # ... I should change things to make debug.HQ a separate env inheriting asNamespace( debug)
      # ... but still "owned" by the call stack and thus auto-zapped on return to prompt
      
      debug.namespace <- asNamespace( 'debug')
      funs <- lsall( env=debug.namespace) %except% find.funs( 'package:debug')
      funs <- funs %except% lsall( env=namespace.dest)
      funs <- funs[ sapply( funs, function( x) !is.environment( debug.namespace[[x]]))]
      for( ifun in funs)
        assign( x=ifun, value=debug.namespace[[ ifun]], envir=namespace.dest)
    }
  }
  
  lapply( list.of.command.subs, subfun)
}


"move.to.next.expression" <-
function( sorted.out=TRUE, nlocal=sys.parent(), original.i, tryout) mlocal({
# sorted.out=FALSE is used only after user calls break/next/skip
# Start with control decisions

  original.i <- i # in case of failed assignment; this is new for 1.1
  if( sorted.out) { # Normal behaviour: try to move into a loop etc.
    if(call.type == "for") { # NB documented R behaviour when looping over a factor...
      if( j %is.a% 'factor') 
        j <- as.character( j) # as.integer before ~R2.9
      for.counters[[ ch.i]] <- j
      i <- c( i, 4)
      sorted.out <- FALSE } # so we set the counter
    else if( call.type == 'repeat')
      i <- c( i, 2)
    else if( call.type %in% c( 'if', 'while')) {
      if( as.logical( j)[ 1])
        i <- c(i, 3)
      else if( length( expr[[ i]])==4) # there's an ELSE
        i <- c( i, 4)
      else {
        if( call.type=='while')
          j <- NULL # R 2.10; value of while-loop is NULL
        sorted.out <- FALSE 
      }
    } else if( call.type == 'switch') { # NB we have already guaranteed that length(j)==1 && j is char or numeric
      swlen <- length( expr[[ i]])
      if( is.numeric( j))
        jj <- floor( as.double( j)) # avoid complex floor
      else { # character match
        jj <- match( j, names( expr[[ i]]), NA)-2
        # If no match, move to "otherwise" if there is one; else leave it as NA
        # Note that displayed value might be misleading
        if( is.na( jj) && names( expr[[ i]])[ swlen] == '')
          jj <- swlen-2
      }

      if( is.na( jj) || jj<1 || jj>swlen-2) {
        j <- NULL
        sorted.out <- FALSE } # move to expr after switch
      else { # find first non-missing expression at or after matched position
        jj <- jj+2
        # Have to avoid evaluating missing args
        while( jj<swlen && is.name( expr[[ c( i, jj)]]) &&
            !nzchar( expr[[ c(i,jj)]]))
          jj <- jj+1
        i <- c( i, jj)
      } }
    else if( call.type == '{') # }
      i <- c( i, 2)
    else if( call.type == '<-')
      i <- c( i, 3)
    else
      sorted.out <- FALSE }

  if( sorted.out)
return( local.return())

# Now we try to move
  while( !my.all.equal( i, 1)) { # loop terminates when "sorted.out" and at viable statement
    #cat( 'Main move loop: i=', ch( i), '\n' )
    #cat( 'sorted.out', sorted.out, '\n')
    i.parent <- clip( i)
    parent.call.type <- get.call.type( expr[[ i.parent]])
    #cat( 'parent call: ', i, '\n  ')
    #print( expr[[ i.parent]])
    #cat( 'pct', parent.call.type, '\n')

    if( sorted.out && i[ length( i)] <= length( expr[[ i.parent]]) )
return( local.return())

    if( parent.call.type %in% c( '{', 'expression') ) { # }
      if( sorted.out <- i[ length( i)] < length( expr[[ i.parent]]))
        i[ length( i)] <- i[ length( i)]+1
      else
        i <- i.parent } # now sorted.out= FALSE, so force a move
    else if( parent.call.type %in% c( 'while', 'repeat')) {
      i <- i.parent
      sorted.out <- TRUE } # and will definitely exit
    else if( parent.call.type == 'for') {
      ch.ip <- ch( i.parent)
      if( sorted.out <- length( for.counters[[ ch.ip]])) {
        val <- for.counters[[ ch.ip]][[ 1]]
        assign( as.character( expr[[ c( i.parent, 2)]]), val, envir=frame)
        if( .step.) {
          cat( '\nFor-loop counter: ', expr[[ c( i.parent, 2) ]], '\nValue:\n',
              file=debug.catfile())
          printIfSmall( val, ofile=debug.catfile())
        }
        for.counters[[ ch.ip]] <- for.counters[[ ch.ip]][ -1 ] }
      else { # if loop counter exhausted
        j <- NULL # R 2.10; value of for-loop is NULL
        i <- i.parent 
      }
      # sorted.out=FALSE so force a move
    } else if( parent.call.type == '<-') { 
      # This has to be a bit weird to prevent the call to '<-' from over-evaluating
      if( missing( j)){ # but presumably legally, e.g. via formals( fun)$arg.with.no.default
        callo <- call( '<-', expr[[ c( i.parent, 2)]], quote( formals( glm)$subset))
        eval( callo, envir=frame)
        .evaluated.OK. <<- TRUE
      } else {
        callo <- call( '<-', expr[[ c( i.parent, 2)]], call( 'quote', j))
        tryout <- try2( list( eval( callo, envir=frame))) # list wrap new in 10/2008
        .evaluated.OK. <<- (tryout %is.not.a% 'try-error') 
      }
      # was: .evaluated.OK. <<- missing( tryout) || (tryout %is.not.a% 'try-error')
      # not sure what the missing test was for, or how to re-write
      if( !.evaluated.OK.) {
        i <- original.i # back out
return( local.return())
      } # if not eval OK
      
      i <- i.parent
      sorted.out <- FALSE
    } else if( parent.call.type %in% c( 'if', 'switch')) { # can't see what else it might be
      i <- i.parent # but must move on
      sorted.out <- FALSE 
    } # parent call type
  } # while anything left to do

# If we are here, then i==1 and we are meant to move to return-value breakpoint
  i <- 2
})


"mtrace" <-
function( fname=NULL, tracing=TRUE, char.fname=as.character( substitute( fname)),
    from=mvb.sys.parent(), update.tracees=TRUE, return.envs=FALSE) {
# do.in.envir call REMOVED... gulp
# mtrace is "do.in.envir" (of its caller) so it can be called WHILE debugging another function
  assign( '[[', my.index)
  lcf <- length( char.fname)
  if( lcf > 1) {
    dunno <- TRUE
    if( lcf==3) {
      dunno <- FALSE
      from <- switch( char.fname[1],
        ':::' = asNamespace( char.fname[2]),
        '::' = as.environment( 'package:' %&% char.fname[2]),
        '$' = eval( as.name( char.fname[2]), envir=from),
          { 
            dunno <- TRUE;
            NULL
          }
        )
      char.fname <- char.fname[3]
    }
  } else
    dunno <- FALSE
  
  if( dunno)
stop( "Dunno wot to do with " %&% deparse( substitute( fname)))    

  fname <- char.fname

  ff <- fun.locator( fname, from) # will force ref classes
  
  if( !tracing && !length( ff)) # couldn't find
    f <- NULL # not completely useless; zap entry in 'tracees'
  else {
    if( !length( ff))
stop( "Can't find " %&% fname)
    f <- get( char.fname, ff[[1]])
    old.env <- environment( f)
    old.attr <- attributes( f)
    old.isS4 <- isS4( f)
  }

  if( tracing) {
    if( is.mtraced( f)) {
      cat( 'Re-applying trace...\n')
      f <- unmtrace( f)
    }

    preamble <- character(0)
    orig.body <- body( f)
    # normal, mlocal, or do.in.envir?
    if( is.recursive( body( f)) && body( f)[[1]]=='mlocal') {
      cc <- substitute(
          return( mlocal( debug:::evaluator( fname=this.fun.name))),
          list( this.fun.name=fname))
      preamble <- 'mlocal(' #)
      body( f) <- body( f)[[2]] }
    else if( is.recursive( body( f)) && body( f)[[1]]=='do.in.envir') {
      mc <- match.call( definition=do.in.envir, call=body( f))
      if( any( names( mc) == 'envir'))
        cc <- substitute(
            return( do.in.envir( envir=this.envir, fbody=debug:::evaluator( fname=this.fun.name))),
            list( this.fun.name=fname, this.envir=mc$envir))
      else
        cc <- substitute(
            return( do.in.envir( fbody=debug:::evaluator( fname=this.fun.name))),
            list( this.fun.name=fname))
      body( f) <- mc$fbody
      preamble <- 'do.in.envir( envir=' %&% paste( deparse( mc$envir), collapse=' ') %&% ',' } # )
    else # normal
      cc <- substitute(
          return( debug:::evaluator( fname=this.fun.name)), list( this.fun.name=fname) )

    if( getOption( 'debug.src', FALSE))
      add.numbers <- an # with source-code capability

    this.tracee <- add.numbers( f, preamble=preamble)
    orig.args <- args( f)
    body( f) <- call( '{', cc, orig.args, orig.body ) # } to keep matcher happy

    # Now substitute delicate stuff in formal args:
    list.of.command.subs <- make.locs( NULL)
    for( arg.name in names( formals( f))) {
      this.arg <- formals( f)[[ arg.name]]
      if( !missing( this.arg) && ((mode( this.arg) != 'name') || nchar( as.character( this.arg)))) {
        this.arg <- do.call( 'substitute', list( this.arg, list.of.command.subs)) # next etc.
        # Next line has 'list' wrapper so NULL doesn't delete
        formals( f)[ arg.name] <- list( debug.mvb.subst( this.arg)) # sys.nframe etc.
      }
    }

    tracees <<- tracees %without.name% fname
    tracees <<- c( tracees, structure( .Data=list(this.tracee), names=fname))
  } else { # untracing
    if( is.mtraced( f))
      f <- unmtrace( f)
    tracees <<- tracees %without.name% fname
  }

  if( !is.null( f)) { # IE we need to save f
    environment( f) <- old.env
    attributes( f) <- old.attr
    if( old.isS4)
      f <- asS4( f)

    # Now check whether the version being saved is in the temporary frame stack
    if( any( sapply( sys.frames(), identical, y=ff[[1]])))
      ff <- ff[ 1] # don't change any deeper copies or permanent copies

    for( this.ff in ff) {
      locko <- balloonIsTethered( fname, this.ff)
      if( locko)
        untetherBalloon( fname, this.ff)
      assign( fname, f, envir=this.ff)
      if( locko) {
        ow <- getOption( 'warn')
        try( {
          options( warn=-1)
          tetherBalloon( fname, this.ff)
          }
        )
        options( warn=ow)
      }
    }
  }

  if( return.envs)
return( ff)
  else
return( invisible( f))
}


"mtrace.off" <-
function() {
  mtrace.off.quietly <- function( fname) try( mtrace( char.fname=fname, tracing=FALSE), silent=TRUE)
  sapply( names( tracees), mtrace.off.quietly)
  invisible( NULL)
}


"mtrace.S4" <-
function( fname=NULL, signature=NULL, egargs=NULL, tracing=TRUE, char.fname=as.character( substitute( fname)),
    update.tracees=TRUE) {
  fname <- char.fname
  
  if( !is.null( egargs)) {
    mc <- match.call( fname, as.call( c( list( char.fname), egargs)), expand.dots=FALSE)
    mc <- as.list( mc)[-1]
    mc <- sapply( mc, class)
    mc <- mc[ environment( fname)$.SigArgs]
    mc[ is.na( mc)] <- 'missing'
    target <- attr( selectMethod( char.fname, signature=mc), 'target')
    mtrace( char.fname=paste( target, collapse='#'), from=environment( fname)$.AllMTable, tracing=tracing)
  }

  if( is.null( signature)) {
    if( !tracing)
return( mtrace( char.fname=char.fname, tracing=FALSE)) # de-mtrace generic

    # Trace all calls via the generic... gulp!
    # Body of generic needs to be something like this:
    # local( {mc <- match.call( expand.dots=F)
    #   siggo <- named( environment( sys.function())$.SigArgs)
    #   siggo[] <- 'missing'
    #   mc <- mc[ names( mc) %that.are.in% siggo]
    #   siggo[ names( mc)] <- sapply( mc, function( x) class( eval( x, parent.frame())))
    #   target <- attr( selectMethod( MEMEME, signature=siggo), 'target')
    #   mtrace( char.fname=target, from=environment( sys.function())$.AllMTable)
    # })
    # standardGeneric( MEMEME)
  }
  
return( invisible( NULL))  
}


"mvb.sys.on.exit" <-
function() {
  p <- mvb.sys.parent()
  f <- find.debug.HQ()
  f <- get( '.frames.', envir=f)
  dbg <- f$debug[ match( p, f$actual, 0)]
  if( !length( dbg))
stop( "sys.on.exit won't work reliably here when the debugger is being used; must be 'unusual' code!")

  get( 'expr', envir=sys.frame(dbg))[[2]][[2]]
}


"my.simple.func" <-
function() sys.parent()


"next.incarnation" <-
function( nlocal=sys.parent()) mlocal({
# Set up debug admin. Already environmentalized so it'll find .frames. etc.
  trf <- try( get( 'tracees', 'package:debug', inherits=FALSE))
  if( (trf %is.a% 'try-error') || !( fname %in% names( trf)))
stop( "No mtrace info for " %&% fname %&% 
    "; maybe saved before being un-mtraced?")
  trf <- trf[[ fname]] # match must be exact

  frame.number <- mvb.sys.parent( 1)
  frame <- sys.frame( frame.number)

  subframe <- sum( .frames.$actual==frame.number) # this is to do with "mlocal" functions
  if( subframe)
    subframe <- LETTERS[ subframe]
  else
    subframe <- ''
  .frames. <<- rbind( .frames., list( actual=frame.number, debug=mvb.sys.parent(0), function.name=fname, subframe=subframe,
      has.window.yet=FALSE, window.name= fname %&% '(' %&% frame.number %&% subframe %&% ')'))

# Doesn't handle .Last.debug in R version. In the S version, this attempts to shut down the Pascal window and then call the
# top .Last in the global environment, if it exists. This is of session-long duration (doesn't work in R as no frame 0)

  for.counters <- list()
  old.l <- -1 # shows which line is currently highlit
  lno <- 1
  breakpoints <- tracees[[ fname]]$breakpoint
  line.list <- tracees[[ fname]]$line.list
  orig.line.list.length <- length( line.list)
  orig.breakpoint.length <- length( breakpoints)
  temp.bp <- ''
#  expr_ do.call( 'substitute', list( tracees[[ fname]]$expr, list.of.command.subs)) # this should have been done by 'add.numbers'
  expr <- tracees[[ fname]]$expr
  expr <- do.call( 'expression', list( expr, quote( NULL), quote( NULL))) # wrap it up one level; matches 'add.numbers'
  stopped.yet <- FALSE
  in.body.code <- TRUE
  i <- 1
  j <- NULL

# Returns name of the window-- but I think this is obsolete
  .frames.$window.name[ nrow( .frames.)]
})


"printIfSmall" <-
function(x, ..., ofile=stdout()) {
  osx <- try( object.size( x), silent=TRUE)
  if( osx %is.a% 'try-error')
    osx <- NA
  
  if( !is.na( osx) && osx < option.or.default( 'threshold.debug.autoprint.size', 8192)) {
    # Also check for too long
    # Code below is bad because "Error..." is printed even when it's just taking too long
    # Also, would be nice to preserve existing time limits if any... somehow...
    
    tl <- getOption( 'debug.print.time.limit', 0.5)
    if( is.finite( tl))
      setTimeLimit( Inf, tl, TRUE)
    try.to.print <- try2( list( capture.output( print(x, ...), file=ofile)))
    if( is.finite( tl))
      setTimeLimit( Inf, Inf, TRUE)
      
    print.obj.info <- FALSE
    if( try.to.print %is.a% 'try-error') {
      if( grepl( 'reached elapsed time limit', c( try.to.print))) {
        print.obj.info <- TRUE
        cat( "<<Printing is taking too long... truncated>>\n", file=ofile)
      } else 
        cat( "!! Couldn't successfully print result: ", c( try.to.print), "\n", file=ofile)
    }
  } else 
    print.obj.info <- TRUE
    
  if( print.obj.info) {
    if(!is.null(class <- attr(x, "class")))
      cat("Class: ", class, " ", file=ofile)
    cat("Mode: ", mode(x), " ", file=ofile)
    cat("Length: ", length <- length(x), " ", file=ofile)
    if(!is.null(dim <- dim(x)))
      cat("dim: ", dim, " ", file=ofile)
    if(is.numeric(x))
      cat("Storage: ", storage.mode(x), " ", file=ofile)
    cat("Size: ", osx, "\n", file=ofile)
  }
}


"qqq" <-
function() do.in.envir( envir=find.debug.HQ( FALSE), {
  .quit.debug. <<- TRUE
  .nothing.
})


"screen.line" <-
function( l, nlocal=sys.parent()) mlocal(
  index( nzchar( names( line.list)))[ l]-1 # tcl starts at 0, aaargh
)


"set.a.breakpoint" <-
function( bp.expr,line.no, frame.number=sys.parent()) do.in.envir( envir=sys.frame( frame.number), {
  breakpoints[[ line.no]] <<- bp.expr

  if( exists( 'bp.win', envir=sys.frame( frame.number), inherits=FALSE)) { # having cake + eating it
    colour <- if( mark.bp( bp.expr)) 'red' else 'white'
    tkitemconfigure( bp.win, screen.line( line.no), foreground=colour, selectforeground=colour)
  }
})


"set.global.debug.vars" <-
function( ...) {
  env <- find.debug.HQ( FALSE)
  l <- list( ...)
  for( i in names( l))
    assign( i, l[[ i]], envir=env)
}


"set.wordstar.key.bindings" <-
function( nlocal=sys.parent()) mlocal({
  # Wordstar key bindings. 
  bindo <- function( keys, fun) {
    # Allow all combos of caps-nocaps, and ctrl-noctrl on 2nd key
    keyseq <- c( toupper( keys[1]), tolower( keys[1]))
    keyseq <- sprintf( '<Control-KeyPress-%s>', keyseq)

    if( length( keys)==2) {
      keyseq2 <- c( toupper( keys[2]), tolower( keys[2]))
      keyseq2a <- sprintf( '<Control-KeyPress-%s>', keyseq2)
      keyseq2b <- sprintf( '<KeyPress-%s>', keyseq2)
      keyseq2 <- c( keyseq2a, keyseq2b)
      keyseq <- c( outer( keyseq, keyseq2, paste, sep=''))
    }

    lapply( keyseq, function( keys) tkbind( tcl.win, keys, fun))
  }

  bindo( 'z', function(){ tkyview( line.list.win, 'scroll', +1, 'units') })
  bindo( 'w', function(){ tkyview( line.list.win, 'scroll', -1, 'units') })
  bindo( 'd', function(){ tkxview( line.list.win, 'scroll', +1, 'units') })
  bindo( 's', function(){ tkxview( line.list.win, 'scroll', -1, 'units') })    
  bindo( 'r', function(){ tkyview( line.list.win, 'scroll', -1, 'pages') })
  bindo( 'c', function(){ tkyview( line.list.win, 'scroll', +1, 'pages') })    
  
  if( !textaroo) {
    # These are listbox-specific
    bindo( 'x', function(){ 
        cursel <- as.integer( tkcurselection( line.list.win))
        if( cursel+1 < length( line.list)) {
          tkselection.clear( line.list.win, cursel)
          tkselection.set( line.list.win, 1+cursel)
          # Should use tkyview to get window position (2 fractions)
          # ... then scroll if new cursel will be outside limits
          # Could bind same to down-arrow to ensure nicely visible
          tkyview( line.list.win, 'scroll', 0, 'units')
        }
      })
    bindo( 'e', function(){ 
        cursel <- as.integer( tkcurselection( line.list.win))
        if( cursel > 0) {
          tkselection.clear( line.list.win, cursel)
          tkselection.set( line.list.win, cursel-1)
          tkyview( line.list.win, 'scroll', 0, 'units')          
        }
      })
    bindo( c( 'k', 'c'), function(){ 
        # Dunno how to get the value direcly out of the tcltk object 'line.list.win'
        cursel <- as.integer( tkcurselection( line.list.win))
        tkclipboard.clear( )
        tkclipboard.append( line.list[ cursel+1])
      })
    } # else if textaroo...
})


"setup.debug.admin" <-
function( nlocal=sys.parent()) mlocal({
# "sys.parent()" is almost certainly NOT the appropriate value for 'nlocal'!
# Intended to be called only by 'find.debug.HQ'
  assign( '[[', my.index)
  assign( '[[<-', my.index.assign)
  
  .frames. <- empty.data.frame( actual=, debug=0, function.name=, window.name=, subframe='', has.window.yet=FALSE)
  .nothing. <- structure( 0, class='nullprint') # invisible object
  .step. <- .skip. <- .evaluated.OK. <- .system. <- .print.result. <- .in.users.commands. <- .quit.debug. <- FALSE
  .end.debug. <- NULL
  debug.catfile <- if( option.or.default( 'debug.catfile', 'stderr')=='stderr') stderr else stdout
  
  # TCL/TK stuff
  back.colour <- c( ' '='White', '*'='Red')
  select.colour <- c( ' '='Blue', '*'='Red')

  values.of.typeof <- cq( symbol, pairlist, closure, environment, promise, language, special,
      builtin, logical, integer, double, complex, character, '...', any, expression, list,
      externalptr)
  rogue.types <- c( 'for', 'while', 'repeat', 'if', 'switch', 'break', 'next', 'return', '{', '<-') # }
  dodgy.for.counter.types <-  cq( language, symbol, promise, environment, closure, '...', any, externalptr)
  dodgy.if.while.types <- cq( 'NULL', pairlist, closure, environment, promise, language,
      special, builtin, '...', any, expression, list, externalptr)
  valid.switch.types <- cq( character, logical, integer, double, complex)

#  cat( 'Before make.locs:\n')
#  print( sys.frames())
#  print( ls( sys.frame(sys.nframe())))
#  print( sys.frame( sys.nframe()))  
#
#  nsd <- sys.frame( mvb.sys.nframe())
#  print( nsd)

  # Command subs
  list.of.command.subs <- make.locs( namespace.dest=sys.frame( mvb.sys.nframe()))
  
#  cat( 'After make.locs:\n')
#  print( ls( sys.frame(sys.nframe())))
#  print( sys.frames())
#  print( sys.frame( sys.nframe()))  

  # Command recall
  # Could implement several options-- recall debug commands only while in debugger,
  # recall debug and non-debug but only while in debugger,
  # merge all debug commands with command-line commands.
  # For now, only the "merge all" option is implemented
  history.available <- function() {
      df <- tempfile()
      sh <- try( savehistory( df), silent=TRUE)
      unlink( df)
      sh %is.not.a% 'try-error'
    }

  if( debug.command.recall <- option.or.default( 'debug.command.recall',
      history.available())) {
    debug.hist.file <- tempfile()
    savehistory( debug.hist.file)
    # Next line would only be used by "recall debug commands only while debugging"
    # original.hist <- scan( debug.hist.file, what='', quiet=TRUE, sep='\n')
  }

  ch <- function( ...) paste( c( ..., '.'), collapse=',')
  star.or.space <- function( bpex) if( is.logical( bpex) && length( bpex)==1 && !is.na( bpex) && !bpex) ' ' else '*'
  mark.bp <- function( bpex) !( is.logical( bpex) && length( bpex)==1 && !is.na( bpex) && !bpex )
  augment.for.eval <- function( i, call.type) c( i, switch( call.type, 'for'=3, 'if'=, 'switch'=, 'while'=2, numeric( 0)) )
#        get.call.type_ function( ex) if( !is.call( ex) || mode( ex)=='(' || !( i_ match( as.character( ex[[ 1]]), rogue.types, 0))) 'normal' else rogue.types[ i]
  get.call.type <- function( ex) {
    if( !is.call( ex) || mode( ex)=='(' || # ')'
        !( i <- match( paste( as.character( ex[[1]]), collapse=' '), rogue.types, 0)) ) {
      if( !is.call( ex) && is.expression( ex))
        'expression'
      else
        'normal' }
    else
      rogue.types[ i]
  }

  # push.to.foreground( r.window.handle)  # BELIEVED OBSOLETE
})


"setup.tcltk.in.control.frame" <-
function( title, nlocal=sys.parent(), nl, screen.pos, font, height, width, i) mlocal({
    tcl.win <- tktoplevel( )
    tktitle( tcl.win) <- title

    nl <- length( line.list)
    # TODO: define def.screen.pos to be a bit smarter about window placement
    # NB tkwinfo( 'screenwidth', tcl.win)
    screen.pos <- option.or.default( 'debug.screen.pos', '+5-5') # was +5-5
    font <- option.or.default( 'debug.font', 'Courier')
    height <- option.or.default( 'debug.height', 10)
    width <- option.or.default( 'debug.width', 120)
    
    # First create the objects, then link them with 'tkconfigure'
    line.list.win <- tklistbox( tcl.win, font=font, bg='white', 
        fg=option.or.default( 'debug.fg', 'black'),
        height=height, width=width, setgrid=TRUE, borderwidth=0)
        
    textaroo <- FALSE
    if( textaroo) {
      line.list.win <- tktext( tcl.win, font=font, bg='white', 
          fg=option.or.default( 'debug.fg', 'black'),
          height=height, width=width, setgrid=TRUE, wrap='none')    
      tktag.configure( line.list.win, 'nexect', background='Green', foreground='Black') 
      tktag.configure( line.list.win, 'alternate', background='Yellow')
      ### NB NB: need to handle breakpoints differently, too. Use a tag, maybe with a border
      # Not clear what to do about selection. Need State==disabled, but doesn't show cursor. 
      # use tag.raise to give selection priority
    }

    bp.win <- tklistbox( tcl.win, font=font, fg='white', bg='white', selectforeground='blue', 
        height=height, setgrid=TRUE, width=1, borderwidth=0, takefocus=FALSE)
    yscroll.win <- tkscrollbar( tcl.win)
    xscroll.win <- tkscrollbar( tcl.win, orient='horizontal')

    tkconfigure( yscroll.win,
        command=function(...) { 
          tkyview( line.list.win,...); tkyview( bp.win, ...) } )
    tkconfigure( xscroll.win, 
        command=function(...) tkxview( line.list.win, ...))
    tkconfigure( line.list.win, 
        yscroll=function(...) { 
          tkset( yscroll.win,...); tkyview( bp.win, 'moveto', list(...)[[1]]) },
        xscroll=function(...) 
          tkset( xscroll.win, ...)) 
    tkconfigure( bp.win, 
        yscroll=function(...) { 
          tkset( yscroll.win,...); tkyview( line.list.win, 'moveto', list(...)[[1]]) })

    if( textaroo) {
      tkinsert( line.list.win, 'end', paste( line.list, collapse='\n'))
    } else {
      # tkinsert doesn't allow length>1 vectors any more; also need to strip names. Hence:
      do.call( 'tkinsert', c( list( line.list.win, 'end'), as.vector( line.list))) 
    }
    
    do.call( 'tkinsert', c( list( bp.win, 'end'), rep( '*', nl)))

    if( getOption( 'debug.wordstar.keys', FALSE))
      set.wordstar.key.bindings()

    # tkpack(buttons, side="bottom", fill="x")#, pady="2m")
    # tkpack( go.button, goto.button, side='left', expand=TRUE)

    tkpack( xscroll.win, side='bottom', fill='x')
    tkpack( yscroll.win, side='right', fill='y')
    tkpack( bp.win, line.list.win, side='left', fill='both', expand=TRUE, ipadx=0)
    
    lapply( screen.pos, function( x) tkwm.geometry( tcl.win, x)) # to allow several calls
    
    tkfocus( line.list.win)
    if( textaroo) {
      tkmark.set( line.list.win, "insert", "0.0")
      tkconfigure( line.list.win, state='disabled')
    } else # listaroo
      tkselection.set( line.list.win, 0)

    bp.list <- unlist( lapply( breakpoints, mark.bp), use.names=F)
    bps <- rep( FALSE, nl)
#cat( '***Lines:', paste( names( line.list), line.list, sep='#'), '***', sep='\n')
#cat( "nbreaks", length( breakpoints), "nl=", nl, "length( bp.list)=", length( bp.list), 'n.names=', 
#sum( names( line.list) != ''), '\n')
    bps[ names( line.list) != ''] <- bp.list
  
    for( i in index( bps))
      tkitemconfigure( bp.win, i-1, foreground='red', selectforeground='red')
      
    # Wake up the window manager, avoiding redraw problems
    # Side-effect is to switch focus to the code window
    # May not be necessary on your system
    # Seems to be ACTIVELY UNHELPFUL on Windows XP, c. R 2.13
    
    # Use 'debug.window.hack' if set--- presumably means calls to 'handy' winmanip in 'launch.debug.windows'
    def.shakeup <- (.Platform$OS.type=='windows') && !getOption( 'debug.window.hack', FALSE)
    if( getOption( 'shakeup.debug.windows', def.shakeup)) {
      # Done thus to let people set their own patches
      window.shaker <- getOption( 'debug.window.shaker', quote( {
          tkwm.withdraw( tcl.win)
          geo <- tkwm.geometry( tcl.win) # Try making 1 less of everything, then 1 more
          tkwm.geometry( tcl.win, squiggle.geometry( geo))
          .Tcl( 'update') # not idletasks; worth a try    
          tkwm.geometry( tcl.win, geo) # this was commented out
          
          if( .Platform$OS.type=='windows') {
            tkwm.state( tcl.win, 'zoomed') # not Linux-able :/
            tkwm.state( tcl.win, 'normal')
          }
          tkwm.deiconify( tcl.win)
        }))
      eval( window.shaker)
    }

    #.Tcl( 'update', 'idletasks')    
    .Tcl( 'update') # not idletasks; worth a try    
})


"skip" <-
function( line.no) do.in.envir( envir=sys.frame( find.active.control.frame()), {
# This is pretty easy & all we need to check, is that we don't try to move INTO a for-loop
  .system. <<- TRUE
  if( line.no<0)
    line.no <- lno-line.no
  if( line.no<0 || line.no > length( breakpoints) ) {
    cat( "can't go there!")
return( .nothing.) }

  target <- names( breakpoints)[ line.no]
  target <- as.numeric( strsplit( substring( target, 1, nchar( target)-1), ',')[[1]])

  mm <- min( length( i), length( target))
  id <- index( i[ 1:mm] != target[ 1:mm])[ 1]
  if( !is.na( id)) # mismatch at position id
    mm <- id-1
#  cat( 'MM=', mm, '\n')

  i.try <- i[ 1 %upto% mm]
  for( j in (mm+1) %upto% length( target))
    if( !is.call( expr[[ i.try]]) || expr[[ c( i.try, 1)]] != 'for')
      i.try <- c( i.try, target[j])
    else {
      cat( "Can't skip into a new for-loop: stopping at the beginning of the for-loop instead\n")
return( .nothing.) }

# For some reason, I used to try to put .skipto. just before the expression to skip to, 
# then (in skipto.debug) I would force a move. But this doesn't work. Don't know why-- 
# I'm leaving the old code in case I had a good reason.
#  i.try[ length( i.try)]_ i.try[ length( i.try)] - 1 # we will force a move anyhow-- I think it works

  .skipto. <<- i.try
  .skip. <<- TRUE
  .evaluated.OK. <<- TRUE
return( .nothing.)
})


"skipto.debug" <-
function( nlocal=sys.parent()) mlocal({
  i <- .skipto. # See "skip". Old comment was: which had better be the statement BEFORE where you want to go
#  move.to.next.expression( FALSE) # FORCE a move; sets up loop
})


"squiggle.geometry" <-
function( geo) {
  siz <- sub( '([0-9]x[0-9]+).*', '\\1', geo)
  sizes <- as.integer( strsplit( siz, 'x')[[1]])
  # cat( sizes, file=stderr())
return( sprintf( '%ix%i%s', sizes[1]-1, sizes[2]-1, substring( geo, nchar( siz)+1)))
}


"step.into.sysfuns" <-
function( ...){
  tagvals <- list( ...)
  if( !length( tagvals))
return( step.intos)

  if( !all( names( tagvals) %in% names( step.intos)))
stop( "Illegal step-into tag name(s)")

  ovals <- step.intos[ names( tagvals)]
  step.intos[ names( tagvals)] <<- unlist( tagvals)
return( ovals)
}


"stepping" <-
function( call.type=''){
  HQ <- debug:::find.debug.HQ( FALSE)
  (HQ %is.an% 'environment') && HQ$.step. && (call.type=='' || step.intos[ call.type])
}


"stop.here" <-
function( nlocal=sys.parent()) mlocal({
# cat("i=", i, "\n")
# cat("stop at: ", names(tracees[[fname]]$breakpoint), "\n")
# cat("temp.bp", tracees[[fname]]$temp.bp, "*\n")
#  if( .skip.) {
#    .set.( .skip.=F, .step.=T)
#return( T) }

  which.bp <- names(breakpoints) == ch.i #
  if(any(which.bp)) {
#
#   Check for temporary & permanent breakpoints
#   cat( mode( temp.bp))
    if( temp.bp == ch.i)
      .step. <<- TRUE
    else
      .step. <<- eval.bp( breakpoints[ which.bp][[1]], envir=frame) || .step.
    # previous line changed 7/1/04 so that breakpoints are eval'd in step mode too

    if(.step.) # discard the temporary breakpoint
      temp.bp <- ''

#   cat("new .step.=", .step., "\n")
  }

  any(which.bp)
})


"unmtrace" <-
function( f) {
  env.f <- environment( f)
  attr.f <- attributes( f)
  formals( f) <- formals( body( f)[[3]])
  
  # Next line: used to need the list(...) call, but this behaviour is changed in R 2.9+
  body( f) <- if( getRversion() >= '2.9.0') body( f)[[4]] else list( body( f)[[4]]) 
  attributes( f) <- attr.f
  environment( f) <- env.f
  f
}


"untracer" <-
function( env){
  undo <- names( tracees) %that.are.in% lsall( env)
  undo <- undo[ sapply( undo, function( x) is.mtraced( env[[x]]))]
  orig <- lapply( named( undo), get, envir=env)
  for( i in names( orig)) {
    f <- env[[ i]]
    env.f <- environment( f)
    attr.f <- attributes( f)
    formals( f) <- formals( body( f)[[3]])
    body( f) <- list( body( f)[[4]]) # destroys attr & env
    attributes( f) <- attr.f
    environment( f) <- env.f
    env[[ i]] <- f
  }
  
  orig
}


"untracer.env" <-
function( env){
  undo <- names( tracees) %that.are.in% lsall( env)
  undo <- undo[ sapply( undo, function( x) is.mtraced( env[[x]]))]
  orig <- lapply( named( undo), get, envir=env)
  for( i in names( orig)) 
    env[[ i]] <- unmtrace( env[[ i]])
  
  orig
}


"xmtrace" <-
function( fname=NULL, tracing=TRUE, char.fname=as.character( substitute( fname)), fexpr=NULL,
    from=mvb.sys.parent(), update.tracees=TRUE, return.envs=FALSE) {
# do.in.envir call REMOVED... gulp
# mtrace is "do.in.envir" (of its caller) so it can be called WHILE debugging another function
  assign( '[[', my.index)
  
  fexpr <- substitute( fname)
  if( !is.null( fexpr)) {
    ind <- 1
    repeat{
      whatcall <- deparse( fexpr[[ind]], nlines=1, width.cutoff=10)
      if( whatcall %not.in% c( '[[', '$'))
  stop( "Can only mtrace [[- or $- components")
      ind[ length( ind)] <- 2
      mofex <- mode( fexpr[[ind]])
      if( mofex=='name') {
        fname1 <- as.character( fexpr[[ind]])
    break
      } else if( mofex != 'call')
  stop( "Don't know how to mtrace a " %&% mofex)
      ind <- c( ind, 1)
    }
        
    ff <- fun.locator( fname1, from, mode='any')
  } else {
    ff <- fun.locator( char.fname, from)
    fexpr <- as.name( char.fname)
  }
  
  fname <- sub( '[^A-Za-z0-9._].*', '', deparse( fexpr, nlines=1, width.cutoff=60))
  
  if( !tracing && !length( ff)) # couldn't find
    f <- NULL # not completely useless; zap entry in 'tracees'
  else {
    if( !length( ff))
stop( "Can't find " %&% fname)
    f <- eval( fexpr, ff[[1]]) # get( char.fname, ff[[1]])
    old.env <- environment( f)
    old.attr <- attributes( f)
  }

  if( tracing) {
    if( is.mtraced( f)) {
      cat( 'Re-applying trace...\n')
      f <- unmtrace( f)
    }

    preamble <- character(0)
    orig.body <- body( f)
    # normal, mlocal, or do.in.envir?
    if( is.recursive( body( f)) && body( f)[[1]]=='mlocal') {
      cc <- substitute(
          return( mlocal( debug:::evaluator( fname=this.fun.name))),
          list( this.fun.name=fname))
      preamble <- 'mlocal(' #)
      body( f) <- body( f)[[2]] }
    else if( is.recursive( body( f)) && body( f)[[1]]=='do.in.envir') {
      mc <- match.call( definition=do.in.envir, call=body( f))
      if( any( names( mc) == 'envir'))
        cc <- substitute(
            return( do.in.envir( envir=this.envir, fbody=debug:::evaluator( fname=this.fun.name))),
            list( this.fun.name=fname, this.envir=mc$envir))
      else
        cc <- substitute(
            return( do.in.envir( fbody=debug:::evaluator( fname=this.fun.name))),
            list( this.fun.name=fname))
      body( f) <- mc$fbody
      preamble <- 'do.in.envir( envir=' %&% paste( deparse( mc$envir), collapse=' ') %&% ',' } # )
    else # normal
      cc <- substitute(
          return( debug:::evaluator( fname=this.fun.name)), list( this.fun.name=fname) )

    this.tracee <- add.numbers( f, preamble=preamble)
    orig.args <- args( f)
    body( f) <- call( '{', cc, orig.args, orig.body ) # } to keep matcher happy

    # Now substitute delicate stuff in formal args:
    list.of.command.subs <- make.locs( NULL)
    for( arg.name in names( formals( f))) {
      this.arg <- formals( f)[[ arg.name]]
      if( !missing( this.arg) && ((mode( this.arg) != 'name') || nchar( as.character( this.arg)))) {
        this.arg <- do.call( 'substitute', list( this.arg, list.of.command.subs)) # next etc.
        # Next line has 'list' wrapper so NULL doesn't delete
        formals( f)[ arg.name] <- list( debug.mvb.subst( this.arg)) # sys.nframe etc.
      }
    }

    tracees <<- tracees %without.name% fname
    tracees <<- c( tracees, structure( .Data=list(this.tracee), names=fname))
  } else { # untracing
    if( is.mtraced( f))
      f <- unmtrace( f)
    tracees <<- tracees %without.name% fname
  }

  if( !is.null( f)) { # IE we need to save f
    environment( f) <- old.env
    attributes( f) <- old.attr

    # Now check whether the version being saved is in the temporary frame stack
    if( any( sapply( sys.frames(), identical, y=ff[[1]])))
      ff <- ff[ 1] # don't change any deeper copies or permanent copies

    for( this.ff in ff) {
      locko <- balloonIsTethered( fname, this.ff)
      if( locko)
        untetherBalloon( fname, this.ff)
      assign( fname, f, envir=this.ff)
      if( locko) {
        ow <- getOption( 'warn')
        try( {
          options( warn=-1)
          tetherBalloon( fname, this.ff)
          }
        )
        options( warn=ow)
      }
    }
  }

  if( return.envs)
return( ff)
  else
return( invisible( f))
}

