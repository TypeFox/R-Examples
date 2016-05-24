TkListView <- function(list) {

  if( !requireNamespace('tcltk', quietly = TRUE) ) {
    stop('This function is dependent on the tcltk package')
  }
  if( !have.ttk() ) {
    stop('this function depends on having tcl 8.5 or higher')
  }


  tt <- tcltk::tktoplevel()
  tcltk::tkwm.title(tt, deparse(substitute(list)))

  fr1 <- tcltk::tkframe(tt)

  tcltk::tkpack(fr1, '-side', 'left', '-fill', 'both', '-expand', 0)
  Sys.sleep(.1)  ## needed for some strange reason.
  tree <- tcltk::ttktreeview(fr1, '-selectmode','browse','-columns',1,height=21)

  scrtree1 <- tcltk::tkscrollbar(fr1, command=function(...)tcltk::tkyview(tree,...))
  scrtree2 <- tcltk::tkscrollbar(fr1, command=function(...)tcltk::tkxview(tree,...), orient='horizontal')
  tcltk::tkconfigure(tree, yscrollcommand=function(...)tcltk::tkset(scrtree1,...),
              xscrollcommand=function(...)tcltk::tkset(scrtree2,...))

#  tkpack(scrtree1, side='right', fill='y',expand=1)
#  tkpack(tree, side='right',fill='both',expand=1)
  tcltk::tkgrid(tree,scrtree1,sticky='nsew')
  tcltk::tkgrid(scrtree2,sticky='nsew')
  tcltk::tkgrid.columnconfigure(fr1, 0, weight=1)
  tcltk::tkgrid.rowconfigure(fr1,0,weight=1)

  fr2 <- tcltk::tkframe(tt)
  tcltk::tkpack(fr2, '-side','top','-fill','both','-expand',1)

  txt <- tcltk::tktext(fr2, bg="white", font="courier", wrap='none', width=40)
  scrtxt1 <- tcltk::tkscrollbar(fr2, command=function(...)tcltk::tkyview(txt,...))
  scrtxt2 <- tcltk::tkscrollbar(fr2, command=function(...)tcltk::tkxview(txt,...), orient='horizontal')
  tcltk::tkconfigure(txt, yscrollcommand=function(...)tcltk::tkset(scrtxt1,...),
              xscrollcommand=function(...)tcltk::tkset(scrtxt2,...))

  tcltk::tkgrid(txt,scrtxt1, sticky='nsew')
  tcltk::tkgrid(scrtxt2,sticky='nsew')
  tcltk::tkgrid.columnconfigure(fr2, 0, weight=1)
  tcltk::tkgrid.rowconfigure(fr2, 0, weight=1)



  buildtree <- function(list, tree, parent) {
    str.info <- capture.output( str(list, max.level=1, give.attr=FALSE,
                                    no.list=TRUE) )
    str.info <- gsub(' |\t','\\\\ ',str.info)

    n <- length(list)
    nms <- names(list)
    if( is.null(nms) ) nms <- rep('', n)

    if( n < length(str.info) ) {
      str.info <- tail(str.info, n)
    }

    for( i in seq(length.out=n) ){
      id <- paste(parent, '.', i, sep='')
      nm <- nms[i]
      if(is.na(nm) || nm == '') nm <- paste0('[[',i,']]')
      tcltk::tkinsert(tree, parent, 'end', '-id', id, '-text', nm, '-values', str.info[i])
      if( is.list(list[[i]]) ){
        Recall( list[[i]], tree, id )
      } else if( !is.null(attributes(list[[i]])) ) {
        tcltk::tkinsert(tree, id, 'end','-id', paste(id,'.a',sep=''), '-text', '<<Attributes>>')
        Recall( attributes(list[[i]]), tree, paste(id,'.a',sep='') )
      }
    }
    tmp <- as.list(attributes(list))
    tmp$names <- NULL
    if( length(tmp) ) {
      tcltk::tkinsert(tree, parent, 'end', '-id', paste(parent,'.aa',sep=''), '-text', '<<Attributes>>')
      Recall( tmp, tree, paste(parent,'.aa',sep='') )
    }
  }

  tmpvals <- capture.output(str(list,max.level=0))
  tmpvals <- gsub(' ','\\\\ ',tmpvals)
  tcltk::tkinsert(tree,'','end','-id', 0, '-text', deparse(substitute(list)), '-values', tmpvals)


  buildtree(list, tree, '0')


  getx <- function(list){
    tmp <- tcltk::tclvalue(tcltk::tkselect(tree))
    tmp2 <- strsplit(tmp, '\\.')[[1]][-1]
    if(length(tmp2)==0) return(list)
    sb <- function(y, list) {
      if (any( y %in% c('a','aa') ) ) {
        a <- which(y %in% c('a','aa'))[1]
        tmp <- if( a==1 ) {
          as.list(attributes( list ) )
        } else {
          y1 <- y[ seq(length.out=a-1) ]
          as.list(attributes( list[[ as.numeric(y1) ]] ))
        }
        if( a == length(y) ) return(tmp)
        y2 <- y[ seq( from=a+1, length.out = length(y) - a) ]
        if( y[a] == 'aa' ) tmp$names <- NULL
        Recall(y2,tmp)
      } else {
        tmp <- as.numeric(y)
        list[[tmp]]
      }
    }
    sb(tmp2,list)
  }

  pr <- tcltk::tkbutton(tt, text='print', command=function(...) {
    tmp <- capture.output(print(getx(list)))
    tcltk::tkdelete(txt, '1.0','end')
    tcltk::tkinsert(txt, 'end', paste(tmp, collapse='\n'))
  }
  )
  st <- tcltk::tkbutton(tt, text='str',   command=function(...) {
    tmp <- capture.output(print(str(getx(list))))
    tcltk::tkdelete(txt, '1.0','end')
    tcltk::tkinsert(txt, 'end', paste(tmp, collapse='\n'))
  }
  )

  tcltk::tkpack(pr, side='top', anchor='w')
  tcltk::tkpack(st, side='top', anchor='w')

  fr3 <- tcltk::tkframe(tt)
  tcltk::tkpack(fr3, side='top', expand=1, fill='x')

  cmd <- tcltk::tclVar('summary(x)')
  eve <- tcltk::tkentry(fr3, textvariable=cmd)
  ev  <- tcltk::tkbutton(fr3, text='Eval:', command=function(...) {
    tmp <- capture.output( eval(parse(text=tcltk::tclvalue(cmd)), list(x=getx(list))))
    tcltk::tkdelete(txt, '1.0', 'end')
    tcltk::tkinsert(txt, 'end', paste(tmp, collapse='\n'))
  }
  )

  tcltk::tkpack(ev, side='left')
  tcltk::tkpack(eve, side='left')

  tcltk::tkpack(tcltk::tkbutton(tt, text='Quit', command=function() tcltk::tkdestroy(tt)),
         side='bottom', anchor='e')

  invisible(NULL)
}
