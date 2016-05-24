# ==================================================================================
# TestScorer 1.7.2
# Last modification: 2016.02.18
# ==================================================================================

# ==================================================================================
# TestScorerGUI()        ENV: present.test, subject
#   on.nameFile()
#   initialitze.test()
#   on.item,numEntry()
#   on.answer.out()
#   invalid.answer()
#   on.show.record()
#   on.write()
#     score.results    called by on.show.record & on.write
#     show.results     id
#     write.results    id
#   on.clean.all()
#   on.new.test()
#   on.clean.test()
#   on.exit()
#   on.help()
#   on.manager()
#
# set.tests.directory()
#
# add.new.test()         ENV: test.char, scale.char
#   make.test.window()
#     no.trans()
#     other.trans()
#     Tmean.trans()
#     raw.is.sum()
#     prorrate.yes()
#     validate.test()
#     OnOkTest()
#   write.test.begining()
#   make.scale.window()   <----------------------------|
#        scale.characteristics                         |
#        validate.scale.char                           |
#        write.scale.begining                          |
#            trans.winow      <-|                      | loop for each scale
#            validate.trans     | loop for each trans  |
#            write.trans      <-|                      |
#   write.scale.end       <----------------------------|
#   write.test.end()
#   MAIN
#
# delete.test()
#
# ==================================================================================

# variables to be treated as global in cheking the package
if(getRversion() >= '2.15.1') utils::globalVariables(c('scoring.fun', 'testChar'))

TestScorerGUI <- function() {
  # =======================================================
  # Main window
  # =======================================================
  
  # set a directory
  catalog <- set.tests.directory() # choose directory and looks for existing tests
  
  # fonts for main GUI
  courrier10 <- tkfont.create(family="courier", size=10)   # font to show items
  arial8 <- tkfont.create(family="arial", size=8)   # font to show test info
  
  # environments for exchanging variable values easier
  present.test <- new.env(parent = emptyenv()) # environment for: test
  #                                            total.items
  #                                            valid.characters
  #                                            item.string
  #                                            next.item
  subject <- new.env(parent = emptyenv()) # environment for: id
  #                                       age
  #                                       sex
  #                                       date.test
  #                                       com.sbj
  
  # top window with three columns: test, id, items
  # -----------------------------------------------
  top <- tktoplevel()
  tkwm.title(top, 'TestScorer 1.7.2')
  tkwm.resizable(top, FALSE, FALSE)
  
  testFrame <- tkframe(top, relief="groove",borderwidth=2)  # col 1 for test
  idFrame <- tkframe(top, relief="groove",borderwidth=2)    # col 2 for id
  respFrame <- tkframe(top, relief="groove",borderwidth=2)  # col 3 for items
  tkgrid.configure(testFrame, idFrame, respFrame)
  
  # column 1: choose a test
  # --------------------------
  number.of.tests <- length(catalog)
  tests <- character(0)
  for (i in 1:number.of.tests)
    tests <- c(tests, catalog[[i]][[1]]) # test acronyms
  which.testFrame <- tkframe(testFrame, relief="groove", borderwidth=2)
  scr <- tkscrollbar(which.testFrame, repeatinterval=5,
                     command=function(...)tkyview(test.selBox,...))
  test.selBox <- tklistbox(which.testFrame, height=6, width='31',
                           selectmode='single', yscrollcommand=function(...)tkset(scr,...),
                           background='white', exportselection=F)
  for (i in (1:number.of.tests))
    tkinsert(test.selBox,'end',tests[i])  # test acronyms
  tkselection.set(test.selBox,0)  # Default test (rem: indexing starts at zero)
  tkgrid.configure(which.testFrame, sticky='nw')
  tkgrid(tklabel(which.testFrame, text='Which test would you like to score?'), sticky='nw')
  tkgrid(test.selBox, scr, sticky='nw')
  tkgrid.configure(scr, rowspan=3, sticky='nse')
  
  # column 1: where to save scores/items
  # ----------------------------------------------------
  on.nameFile <- function() {    # button function to change file
    choice.top <- tktoplevel()
    tkwm.deiconify(choice.top)
    tkwm.resizable(choice.top, FALSE, FALSE)
    tkgrab.set(choice.top)
    tktitle(choice.top) <- 'Saving scores'
    tkgrid(tklabel(choice.top, text=' '))
    choiceFrame <- tkframe(choice.top, relief="groove", borderwidth=2)
    tkgrid(choiceFrame)
    
    rb.dont <- tkradiobutton(choice.top)
    rb.new <- tkradiobutton(choice.top)
    rb.exist <- tkradiobutton(choice.top)
    rbChoice <- tclVar("dont")
    tkconfigure(rb.dont, variable=rbChoice ,value="dont")
    tkconfigure(rb.new, variable=rbChoice ,value="new")
    tkconfigure(rb.exist, variable=rbChoice, value="exist")
    tkgrid(tklabel(choiceFrame, text="Save to:"), stick='w')
    tkgrid(tklabel(choiceFrame, text="   Don't save"), rb.dont, sticky='w')
    tkgrid(tklabel(choiceFrame, text="   An existing file"), rb.exist, sticky='w')
    tkgrid(tklabel(choiceFrame, text="   A new file"), rb.new, sticky='w')
    
    tkgrid(tklabel(choice.top, text=' '))
    donechoiceVar <- tclVar(0) # to stop continuing
    OkDoneChoice.but <- tkbutton(choice.top, text='   OK   ',
                                 command=function() tclvalue(donechoiceVar) <- 1)
    tkgrid(OkDoneChoice.but)
    tkgrid(tklabel(choice.top, text=paste('Close this window to cancel the action.',
                                          'Option will be set to "Don\'t save".',
                                          sep='\n')))
    tkbind(choice.top, '<Destroy>', function() tclvalue(donechoiceVar) <- 3)
    tkfocus(choice.top)         # Place the focus to this tk window
    tkwait.variable(donechoiceVar)
    doneChoice <- as.integer(tclvalue(donechoiceVar))
    tkgrab.release(choice.top)
    tkdestroy(choice.top)
    
    nameFile=''
    if (doneChoice==1) {
      # choice don't save
      if (tclvalue(rbChoice)=='dont') {
        tclvalue(nameFileVar) <- "Don't save"
        tclvalue(write.itemsVar) <- 'Yes'
      }
      
      else if (tclvalue(rbChoice)=='new') {
        # choice new
        msg <- paste("First choose a directory.",
                     "Then a name for the file will be asked.",
                     "",
                     "Note: It is not save to mix data with scoring scripts",
                     "           in the same directory.",
                     sep="\n")
        tkmessageBox(message=msg)
        name.dir <- tk_choose.dir()
        
        if (!is.na(name.dir)) {
          # give a name for the file
          a.file <- tktoplevel()
          tkwm.deiconify(a.file)
          tkwm.resizable(a.file, FALSE, FALSE)
          tkgrab.set(a.file)
          tktitle(a.file) <- 'File'
          tkgrid(tklabel(a.file, text=''))
          nameFrame <- tkframe(a.file, relief="groove", borderwidth=2)
          tkgrid(nameFrame)
          
          nameVar <- tclVar('')
          nameEntryWidget <- tkentry(nameFrame, width = 20,
                                     textvariable = nameVar)
          tkgrid(tklabel(nameFrame, text = 'Name for the file:'), nameEntryWidget)
          
          tkgrid(tklabel(a.file, text="DON'T MISS THE EXTENSION!"))
          tkgrid(tklabel(a.file, text=''))
          tkfocus(nameEntryWidget)
          donechoiceVar <- tclVar(0) # to stop continuing
          OkDoneChoice.but <- tkbutton(a.file, text='   OK   ',
                                       command=function() tclvalue(donechoiceVar) <- 1)
          tkgrid(OkDoneChoice.but)
          tkgrid(tklabel(a.file, text='Close this window to cancel the action.'))
          tkbind(a.file, '<Destroy>', function() tclvalue(donechoiceVar) <- 3)
          tkfocus(a.file)         # Place the focus to this tk window
          tkwait.variable(donechoiceVar)
          doneChoice <- as.integer(tclvalue(donechoiceVar))
          new.file <- tclvalue(nameVar)
          tkgrab.release(a.file)
          tkdestroy(a.file)
        } else new.file <- ''
        
        if (new.file!='') {
          nameFile <- paste(name.dir, '/', new.file, sep='', collapse='')
          the.file.exists <- file.exists(nameFile)
          if (the.file.exists) {
            msg <- paste('File already exists.',
                         'New data will be appended.',
                         sep='\n')
            tkmessageBox(message=msg, icon='warning', type='ok')
          }
        } else nameFile <- ''
      }
      
      else if (tclvalue(rbChoice)=='exist') {
        # choice existing file
        nameFile <- tk_choose.files(multi=FALSE)
        if (length(nameFile)==0) nameFile <- ''  # windows canceled
      }
    }
    
    if (nameFile=='') {  # window canceled
      tclvalue(nameFileVar) <- "Don't save"
      tclvalue(write.itemsVar) <- 'Yes'
      tkconfigure(showBut, text='      Show results      ')
      tkconfigure(writeBut, foreground='gray50')
      tkconfigure(YesButton, state='disabled')
      tkconfigure(NoButton, state='disabled')
      msg <- "Choosing a file has been canceled. \nOption is set to \"Don't save\"."
      tkmessageBox(message=msg, icon='warning', type='ok')
      } else if (grepl('TST_', nameFile)) {  # names with TST_ are not allowed
          msg <- paste('Names begining with the string "TST_" are not allowed.',
                       'Choose another name.',
                       sep='\n')
          tkmessageBox(message=msg, icon='error', type='ok')
          tclvalue(nameFileVar) <- "Don't save"
          tclvalue(write.itemsVar) <- 'Yes'
          tkconfigure(showBut, text='     Show results     ')
          tkconfigure(writeBut, foreground='gray50')
          tkconfigure(YesButton, state='disabled')
          tkconfigure(NoButton, state='disabled')
          
        } else if (!file.exists(nameFile)) {  # if it is a new file you can continue
            tclvalue(nameFileVar) <- nameFile
            tkconfigure(showBut, text='Show & write results')
            tkconfigure(writeBut, foreground='black')
            tkconfigure(YesButton, state='normal')
            tkconfigure(NoButton, state='normal')
          } else {                  # is an existent valid file for test scorer?
            options(warn = -1)
            disk.file <- read.csv2(nameFile, nrows=1)  # only col.names are needed
            options(warn = 0)
            disk.file.names <- names(disk.file)
            names.needed <- c('id', 'test', 'age', 'sex', 'date', 'obs')
            if(!all(names.needed %in% disk.file.names)) {
              msg <- paste("The file ",
                           nameFile,
                           " has not a valid structure.",
                           "\nChoose or create an other file.",
                           sep="")
              tkmessageBox(message=msg, icon='error', type='ok')
              tclvalue(nameFileVar) <- "Don't save"
              tclvalue(write.itemsVar) <- 'Yes'
              tkconfigure(showBut, text='     Show results     ')
              tkconfigure(writeBut, foreground='gray50')
              tkconfigure(YesButton, state='disabled')
              tkconfigure(NoButton, state='disabled')
            } else {                 # is for this test?
              present.test <- sub('^TST_', '', present.test[['test']])
              if(disk.file[1, 2] != present.test) {
                msg <- paste("The file ",
                             nameFile,
                             " is for an other test (",
                             disk.file[1, 2],
                             ").",
                             "\nChoose or create an other file.",
                             sep="")
                tkmessageBox(message=msg, icon='error', type='ok')
                tclvalue(nameFileVar) <- "Don't save"
                tclvalue(write.itemsVar) <- 'Yes'
                tkconfigure(showBut, text='     Show results     ')
                tkconfigure(writeBut, foreground='gray50')
                tkconfigure(YesButton, state='disabled')
                tkconfigure(NoButton, state='disabled')
              } else {                   # the file is valid
                  tclvalue(nameFileVar) <- nameFile
                  tkconfigure(showBut, text='Show & write results')
                  tkconfigure(writeBut, foreground='black')
                  tkconfigure(YesButton, state='normal')
                  tkconfigure(NoButton, state='normal')
                }
            }
          }
  }  # end on.nameFile
  
  where.writeFrame <- tkframe(testFrame, relief="groove", borderwidth=2)
  fileLbl <- tklabel(where.writeFrame, 
                     text='Where would you like to save scores in tabular format?  ')
  tkgrid(where.writeFrame, sticky='w')
  tkgrid(fileLbl, sticky='w')
  
  nameFileVar <- tclVar("Don't save")     # initial value
  nameFileEntry <- tklabel(where.writeFrame, text=tclvalue(nameFileVar),
                           width=70)
  
  fileNameBut <- tkbutton(where.writeFrame,
                          text='Change\noption',
                          command=on.nameFile)
  tkconfigure(nameFileEntry, textvariable=nameFileVar, bg='lightgrey',
              relief='groove', borderwidth=2)
  tkgrid(nameFileEntry, fileNameBut, sticky='ew')

  write.itemsFrame<- tkframe(where.writeFrame)
  write.itemsVar <- tclVar('Yes')
  YesButton <- tkradiobutton(write.itemsFrame, variable=write.itemsVar,
                             value='Yes')
  NoButton <- tkradiobutton(write.itemsFrame, variable=write.itemsVar,
                            value='No')
  tkconfigure(YesButton, state='disabled')
  tkconfigure(NoButton, state='disabled')
  tkgrid(write.itemsFrame, sticky='w')
  tkgrid.configure(tklabel(write.itemsFrame,
                           text='Would you like to save also the items?'),
                   tklabel(write.itemsFrame, text='Yes'), YesButton,
                   tklabel(write.itemsFrame, text='No'), NoButton,
                   sticky='w')
  
  # column 1: test details (items, valid answers...)
  # ---------------------------------------------------------------------
  testCharLbl <- tklabel(testFrame, text='Test information')
  testCharEntry <- tktext(testFrame, width=71, height=7,
                          bg='lightgrey', font=arial8)
  tkgrid(testCharLbl, sticky='w')
  tkgrid(testCharEntry, sticky='w')
  
  # column 2: id, age, sex, date & comment
  # ----------------------------
  idVar <- tclVar('')
  idEntry <- tkentry(idFrame, width='30', textvariable=idVar)
  ageVar <- tclVar('')
  ageEntry <- tkentry(idFrame, width='2', textvariable=ageVar)
  sexVar <- tclVar('Male')
  sexFrame <- tkframe(idFrame, relief="groove", borderwidth=2)
  maleButton <- tkradiobutton(sexFrame, variable=sexVar, value='Male')
  femaleButton <- tkradiobutton(sexFrame, variable=sexVar, value='Female')
  dateVar <- tclVar('')
  dateEntry <- tkentry(idFrame, width='10', textvariable=dateVar)
  commVar <- tclVar('')
  commEntry <- tkentry(idFrame, width='30', textvariable=commVar)
  tkgrid.configure(idFrame, sticky='nw')
  tkgrid(tklabel(idFrame, text='Identification'), sticky='nw')
  tkgrid(idEntry, sticky='nw')
  tkgrid(tklabel(idFrame, text='Age'), sticky='nw')
  tkgrid(ageEntry, sticky='nw')
  tkgrid(sexFrame, sticky='nw')
  tkgrid.configure(tklabel(sexFrame, text='Male'),
                   maleButton, tklabel(sexFrame, text='Female'),
                   femaleButton)
  tkgrid(tklabel(idFrame, text='Date'), sticky='nw')
  tkgrid(dateEntry, sticky='nw')
  tkgrid(tklabel(idFrame, text='Comment'), sticky='nw')
  tkgrid(commEntry, sticky='nw')
  
  # column 3: box for entry items
  # -------------------------
  entryFrame <- tkframe(respFrame, relief="groove", borderwidth=2)
  tkgrid(tklabel(entryFrame, text="Entry items window"), sticky='nw')
  item.numVar <- tclVar('Item to enter: 1 -->') # shows number, initialized at 1
  item.numEntry <- tklabel(entryFrame, text=tclvalue(item.numVar), width=20)
  itemVar <- tclVar('')                             # to enter items
  itemEntry <- tkentry(entryFrame, width='1',
                       textvariable=idVar, font=courrier10 )
  tkconfigure(item.numEntry, textvariable=item.numVar, fg='blue')
  tkconfigure(itemEntry, textvariable=itemVar)
  tkgrid.configure(respFrame, sticky='nw')
  tkgrid(item.numEntry, itemEntry)
  tkgrid(entryFrame, sticky='ew')
  
  # column 3: window to show items
  # ---------------------------
  itemsFrame <- tkframe(respFrame) # put together items and marginal numbers
  items <- tktext(itemsFrame, bg='lightgrey', font=courrier10,
                  height=10, width='10')
  tkconfigure(items, state='disabled')
  itemsText <- tktext(itemsFrame, bg='lightgrey', font=courrier10,
                      height=10, width='7')                    # row lbl
  tkinsert(itemsText, 'end',
           paste(seq(1, 91, 10), '-', seq(10, 100, 10), '\n', 
                 sep='', collapse='')) # lbl rows
  tkconfigure(itemsText, state='disabled')
  tkgrid(tklabel(itemsFrame, text='1234567890', 
                 font=courrier10)) # columns lbl
  tkgrid(items, itemsText)
  tkgrid(itemsFrame)
  
  initialize.test <- function(num.order) {
    # functions for entering items
    # ------------------------
    # reads tests characteristics, makes a chain with
    # valid answers & shows this information
    test <- paste('TST_', catalog[[num.order]][[1]],
                  sep='')  # TST_ added
    assign('test', test, envir=present.test)
    total.items <- catalog[[num.order]][[4]]
    assign('total.items', total.items, envir=present.test)
    valid.characters <- c(catalog[[num.order]][[5]],
                          catalog[[num.order]][[6]], 'space')
    assign('valid.characters', valid.characters, envir=present.test)
    # items chain <- blanks & '=' to complete the hundred
    items.string <- paste(c(rep(' ', present.test[['total.items']]),
                            rep('=', (present.test[['total.items']] %/% 100 +1)*100 - present.test[['total.items']])),
                          collapse='')
    assign('items.string', items.string, envir=present.test)
    next.item <- 1   # next item to entry
    assign('next.item', next.item, envir=present.test)
    tkconfigure(testCharEntry, stat='normal')   # updates window
    tkdelete(testCharEntry, '0.0', 'end')
    tkinsert(testCharEntry,"end","Test name: ")
    tkinsert(testCharEntry,'end', paste(catalog[[num.order]][[2]]))
    tkinsert(testCharEntry,"end","\nAuthor: ")
    tkinsert(testCharEntry,'end', paste(catalog[[num.order]][[3]]))
    tkinsert(testCharEntry,"end","\nNumber of items: ")
    tkinsert(testCharEntry,'end', paste(catalog[[num.order]][[4]]))
    tkinsert(testCharEntry,"end","\nValid answers: ")
    tkinsert(testCharEntry,'end', paste(catalog[[num.order]][[5]]))
    tkinsert(testCharEntry,'end','\nMissings: ')
    if (length(catalog[[num.order]][[6]])==0)
      tkinsert(testCharEntry,'end', 'space')
    else tkinsert(testCharEntry,'end', paste(catalog[[num.order]][[6]], '& space'))
    tkinsert(testCharEntry,'end','\nComment: ')
    tkinsert(testCharEntry,'end', paste(catalog[[num.order]][[7]]))
    tkconfigure(testCharEntry, state='disabled')
    tkconfigure(items, state='normal')      # updates window with items chain
    tkdelete(items, '0.0', 'end')
    tkinsert(items, '0.0', present.test[['items.string']])
    tkconfigure(items, state='disabled')
  } # end initialize test
  
  on.item.numEntry <- function(K) {  # identifies the pressed key
    # --------- first the item number & the answer are updated
    if (is.element(K, c(present.test[['valid.characters']],
                        'Up', 'Down', 'Left', 'Right',
                        'Num_Lock'))) {  # valid answer
      if (K=='Up')
        if (present.test[['next.item']] < 11) answer.out()
      else present.test[['next.item']] <- present.test[['next.item']] - 10
      else if (K=='Down')
        if (present.test[['next.item']] > present.test[['total.items']]-10) answer.out()
      else present.test[['next.item']] <- present.test[['next.item']] + 10
      else if (K=='Left')
        if (present.test[['next.item']] < 2) answer.out()
      else present.test[['next.item']] <- present.test[['next.item']] - 1
      else if (K=='Right')
        if (present.test[['next.item']] > present.test[['total.items']]) answer.out()
      else present.test[['next.item']] <- present.test[['next.item']] + 1
      else if (K=='Num_Lock') {} # do nothing
      else if (present.test[['next.item']] > present.test[['total.items']]) answer.out()
      else {
        present.test[['next.item']] <- present.test[['next.item']] + 1
        tclvalue(item.numVar) <- present.test[['next.item']]
        if (K=='space') item.value <- ' ' else item.value <- K
        present.test[['items.string']] <- paste(substr(present.test[['items.string']],
                                                       1,
                                                       present.test[['next.item']]-2),
                                                item.value,
                                                substr(present.test[['items.string']],
                                                       present.test[['next.item']],
                                                       nchar(present.test[['items.string']])),
                                                sep='')
      }
    }
    else invalid.answer(K)       # invalid answer
    # screen is updated an asterisk is put where the new item will be located
    items.string.mod <- paste(substr(present.test[['items.string']],
                                     1,
                                     present.test[['next.item']] - 1),
                              '*',
                              substr(present.test[['items.string']],
                                     present.test[['next.item']] + 1,
                                     nchar(present.test[['items.string']])),
                              sep='')
    # computes hundred to select the segment which will be shown
    hundred <- (present.test[['next.item']]-1) %/% 100
    # computes the values of the rows
    tkconfigure(itemsText, state='normal')
    tkdelete(itemsText, '0.0', 'end')
    tkinsert(itemsText, 'end',
             paste(seq(hundred*100+1, hundred*100+91, 10), '-',
                   seq(hundred*100+10, hundred*100+100, 10), '\n',
                   sep='',
                   collapse=''))
    tkconfigure(itemsText, state='disabled')
    # cuts the hundred to be shown
    items.string.mod <- substring(items.string.mod, hundred * 100 + 1,
                                  (hundred + 1) * 100)
    # updates window with items
    tkconfigure(items, state='normal')
    tkdelete(items, '0.0', 'end')
    tkinsert(items, '0.0', items.string.mod)
    tkconfigure(items, state='disabled')
    # updates the entry items box
    tkdelete(itemEntry, '0', 'end')
    tkinsert(itemEntry, '0', substring(present.test[['items.string']],
                                       present.test[['next.item']],
                                       present.test[['next.item']]))
    if (present.test[['next.item']] == present.test[['total.items']]+1) {
      alarm()
      tclvalue(item.numVar) <- 'END OF TEST'
    }
    else tclvalue(item.numVar) <- paste('Item to enter: ',
                                        present.test[['next.item']],
                                        '-->')
  } # end on.item.numEntry
  
  answer.out <- function(K)
  { msg <- paste('You tried to go out of the range of items.',
                 '\n"Test information" shows the number of items of this test.',
                 sep='\n')
    tkmessageBox(message=msg, icon='error')
    tkfocus(itemEntry)
  }  # end answer.out
  
  invalid.answer <- function(K)
  { msg <- paste('Answer "', K, '" is invalid.',
                 '\n\nSee "Test information" window for valid answers.',
                 sep='')
    tkmessageBox(message=msg, icon='error')
    tkfocus(itemEntry)
  }  # end invalid.answer
  
  # buttons: functions and position
  # -------------------------------
  on.show.and.record <- function(){
    results <- score.results()  # results: list(results.lst, results.df, results.scores, answers)
    show.results(results[[1]], results[[2]])  # show on console
    if (tclvalue(nameFileVar) != "Don't save") write.results(results[[3]], results[[4]])  # record results to file
  }  # end on.show.record
  
  on.write <- function(){
    results <- score.results() # # results: list(results.lst, results.df, results.scores, answers)
    write.results(results[[3]], results[[4]]) # save results to file
  }
  
  score.results <- function() {
    # reads identification and answer from tktcl
    assign('id', tclvalue(idVar), envir=subject)
    assign('age', tclvalue(ageVar), envir=subject)
    assign('sex', tclvalue(sexVar), envir=subject)
    assign('date.test', tclvalue(dateVar), envir=subject)
    assign('comm.sbj', tclvalue(commVar), envir=subject)
    answers <- character(present.test[['total.items']]) # string to char
    for (i in 1:present.test[['total.items']])
      answers[i] <- substring(present.test[['items.string']], i, i)
    # now read test script and calls score function
    eval(parse(text=paste(text="source('",
                          present.test[['test']],
                          ".r')",
                          sep=""))) # reads file with scoring code
    results <- scoring.fun(answers,
                           subject[['sex']],
                           subject[['age']],
                           subject[['id']],
                           subject[['date.test']],
                           subject[['comm.sbj']])  # scoring function
    results.lst <- results$results.lst
    results.df <- results$results.df
    results.scores <- results$results.scores
    list(results.lst=results.lst,
         results.df=results.df, 
         results.scores=results.scores,
         answers=answers)
  }  # end score.results 
    
  show.results <- function(results.lst, results.df){
    options(width=100)  # assures enough wide
    cat('\n\n')
    cat('----------------------',
        sub('^TST_', '',
            present.test[['test']]),
        '-----------------------\n')
    cat('Subject:', subject[['id']], '\n')
    cat('Sex:', subject[['sex']], '\n')
    cat('Age:', subject[['age']], '\n')
    cat('Date:', subject[['date.test']], '\n')
    cat('Comment:', subject[['comm.sbj']], '\n\n')
    for (i in seq_along(results.lst))
      cat(results.lst[[i]], '\n')
    cat('\n')
    print(results.df, right=FALSE, row.names=FALSE, na.print='-')
  }  # end show.results
  
  write.results <- function(results.scores, answers){
    if (tclvalue(nameFileVar) != "Don't save") {
      the.file.exists <- file.exists(tclvalue(nameFileVar))
      
      # case = id vars
      options(warn=-1)  # if age is not numeric
      case <- data.frame(id=subject[['id']],
                         test=sub('^TST_', '', present.test[['test']]),
                         age=as.numeric(subject[['age']]),
                         sex=subject[['sex']],
                         date=subject[['date.test']],
                         obs=subject[['comm.sbj']])
      options(warn=0)
      
      # case = id + results
      if (is.data.frame(results.scores))
        case <- cbind(case, results.scores) # scripts v.1.1
      else case <- cbind(case, t(results.scores)) # scripts v.1.2 and later
      
      # case = id + results + items     
      if (tclvalue(write.itemsVar)=='Yes') {
        if (the.file.exists) { # if file exists, cheks if items are admited
          if ((length(results.scores)+6)==length(read.csv2(tclvalue(nameFileVar), header = TRUE, nrows = 1))) {
            #6 for identification vars
            msg <- paste("This file doesn't accept items.",
                         "Items have not been saved, but scales yes.",
                         "To save items use an other (existing or new created) file,",
                         "with the save items option checked.",
                         sep="\n")
            tkmessageBox(message=msg, icon='warning')
            tclvalue(write.itemsVar) <- 'No'
          } else case <- cbind(case, i=t(as.numeric(answers)))
        } else case <- cbind(case, i=t(as.numeric(answers)))
      }
      
      if (!the.file.exists) {  # if new file, write the vars names
        write.table(case, tclvalue(nameFileVar), row.names=F,
                    sep=';', na='', quote=TRUE)
      } else { write.table(case, tclvalue(nameFileVar), append=TRUE,
                           row.names=F, col.names=F, sep=';', na='', quote=TRUE)
      }
      rm(case, inherits=TRUE)   # to avoid dragging
      msg <- 'Test scores has been saved.'
      tkmessageBox(message=msg, icon='info')
    }  # end if... write to file
    else {msg <- paste('First choose a file where to write the resuts.',
                       'Use the "Change option" button.',
                       sep='\n')
          tkmessageBox(message=msg, icon='error')}
  }  # end write.results
  
  on.clean.all <- function() { 	# cleans id & items
    tkdelete(idEntry, 0, 'end')
    tkdelete(ageEntry, 0, 'end')
    tkdelete(dateEntry, 0, 'end')
    tkdelete(commEntry, 0, 'end')
    tkdelete(itemEntry, 0, 'end')
    tclvalue(sexVar) <- 'Male'
    on.clean.test()
  }  # end on.clean.all
                           
  on.new.test <- function() {
    if (tclvalue(nameFileVar)!="Don't save") {
      tclvalue(nameFileVar) <- "Don't save"
      tclvalue(write.itemsVar) <- 'Yes'
      tkconfigure(showBut, text='     Show results     ')
      tkconfigure(writeBut, foreground='gray50')
      tkconfigure(YesButton, state='disabled')
      tkconfigure(NoButton, state='disabled')
      msg <- paste("The file for saving scores is no longer valid for the new test.",
                   "Saving results has been reseted to \"Don't save\".",
                   "Select an other file if appropiate.",
                   sep='\n')
      tkmessageBox(message=msg, icon='warning')
    }
	
    on.clean.test()
  }  # end on.new.test
        
  on.clean.test <- function() { 	# clean test items
    tkdelete(itemEntry, 0, 'end')
    initialize.test(as.numeric(tkcurselection(test.selBox))+1)
    tclvalue(item.numVar) <- 'Item to enter: 1  -->'
    tkconfigure(itemsText, state='normal')
    tkdelete(itemsText, '0.0', 'end')
    tkinsert(itemsText, 'end',
             paste(seq(1, 91, 10), '-', seq(10, 100, 10), '\n', 
                   sep='', collapse='')) # lbl rows
    tkconfigure(itemsText, state='disabled')
    tkgrid(items, itemsText)
  }  # end on.clean.test
        
  on.exit <- function() {
    msg <- paste('   -------------------------------------------------------',
                 '   | To exit the R session close the window or type: q() |',
                 '   |                                                     |',                 
                 '   | TestScorer do not use the workspace image, so       |',
                 '   | you can safely answer "No" to the next question.    |',
                 '   -------------------------------------------------------',
                 sep='\n')
    message(msg)
    tkdestroy(top)
  }  # end on.exit
        
  on.help <- function() {
    # help window
    modalDialog <- function() {
      dlg <- tktoplevel()
      tkwm.title(dlg, 'Help')
      tkwm.resizable(dlg, FALSE, FALSE)
      tkwm.deiconify(dlg)
      tkgrab.set(dlg)
      tkfocus(dlg)
      scrBar <- tkscrollbar(dlg, repeatinterval=5,
                            command=function(...)tkyview(helpText,...))
      helpText <- tktext(dlg, height=30, width='70', font=courrier10,
                         yscrollcommand=function(...)tkset(scrBar,...))
      txt <- readLines('Help.txt', n=-1)
      tkinsert(helpText, 'end', paste(txt, '\n', sep='', collapse=''))
      tkconfigure(helpText, state='disabled')
    
	  onClose <- function() {
        tkgrab.release(dlg)
        tkdestroy(dlg)
        tkfocus(top)
      }
    
	  Close.but <-tkbutton(dlg, text=" Close ",command=onClose)
      tkgrid(helpText, scrBar, sticky='ns')
      tkgrid(Close.but)
      tkfocus(dlg)
      tkbind(dlg, "<Destroy>", function() {tkgrab.release(dlg);tkfocus(top)})
      tkwait.window(dlg)
    }  # end modalDialog
    
	modalDialog()  # calls help window
  }  # end on.help
        
  on.manager <- function() {
    manager.top <- tktoplevel()
    tkwm.deiconify(manager.top)
    tkwm.resizable(manager.top, FALSE, FALSE)
    tkgrab.set(manager.top)
    tktitle(manager.top) <- 'Test Manager'
    tkgrid(tklabel(manager.top, text=' '))
    managerFrame <- tkframe(manager.top, relief="groove",borderwidth=2)
    tkgrid(managerFrame)
    rb.add <- tkradiobutton(managerFrame)
    rb.del <- tkradiobutton(managerFrame)
    rbValue <- tclVar("Add")
    tkconfigure(rb.add, variable=rbValue ,value="Add")
    tkconfigure(rb.del, variable=rbValue, value="Delete")
    tkgrid(tklabel(managerFrame, text="Make a choice:"), stick='w')
    tkgrid(tklabel(managerFrame, text="   Add a new test            "), rb.add)
    tkgrid(tklabel(managerFrame, text="   Delete an existing test"), rb.del)
          
    tkgrid(tklabel(manager.top, text=' '))
    doneManager <- tclVar(0) # to stop continuing
    OkDoneManager.but <- tkbutton(manager.top, text='   OK   ',
                                  command=function() tclvalue(doneManager) <- 1)
    tkgrid(OkDoneManager.but)
    tkgrid(tklabel(manager.top, text='Close this window to cancel the action.'))
    tkbind(manager.top, '<Destroy>', function() tclvalue(doneManager) <- 3)
    tkfocus(manager.top)         # Place the focus to this tk window
    tkwait.variable(doneManager)
    doneVal <- as.integer(tclvalue(doneManager))
    tkgrab.release(manager.top)
    tkdestroy(manager.top)
    if (doneVal==1) {
      tkdestroy(top)
      if (tclvalue(rbValue)=='Add') add.new.test()
      if (tclvalue(rbValue)=='Delete') delete.test()
    }
  }  # end on.manager
        
  # frame for buttons
  buttFrame <- tkframe(top)
  showBut <- tkbutton(buttFrame, text="      Show results      ", command=on.show.and.record)
  writeBut <- tkbutton(buttFrame, text="    Only save    ",
                        foreground='gray50', command=on.write)
  cleanTestBut <- tkbutton(buttFrame, text=" Clean items ", command=on.clean.test)
  cleanAllBut <- tkbutton(buttFrame, text="  Clean all   ", command=on.clean.all)
  exitBut <- tkbutton(buttFrame, text="Exit TestScorer", command=on.exit)
  helpBut <- tkbutton(buttFrame, text="      Help      ", command=on.help)
  managerBut <- tkbutton(buttFrame, text="Test Manager", command=on.manager)
  tkgrid(buttFrame, columnspan=3)
  tkgrid(showBut, writeBut, cleanTestBut, cleanAllBut, managerBut, helpBut, exitBut)
        
  # frame for instructions
  # ----------------------
  notabeneFrame <- tkframe(top, relief="groove", borderwidth=2)
  tkgrid(notabeneFrame, columnspan=3)
  tkgrid(tklabel(notabeneFrame, text='Gray windows are non-editable.'))
  tkgrid(tklabel(notabeneFrame, text='To introduce the items, put the cursor inside the white box in "Entry items window". Then enter the answers using numerical keys.'))
  tkgrid(tklabel(notabeneFrame, text='You can use the arrow keys to skip between items. Do NOT use "backspace".'))
  tkgrid(tklabel(notabeneFrame, text='Use the R console menu option "File" to print or save the results.'))
  tkgrid(tklabel(notabeneFrame, text='Click help for more detailed instructions.'))
        
  # key interception
  # ----------------
  tkbind(itemEntry, "<KeyRelease>", on.item.numEntry)
  tkbind(test.selBox, "<ButtonRelease-1>", on.new.test)
        
  # begining: test 1 as defect
  # ---------------------------
  initialize.test(1) # test 1 as default at the beginning of the session
  on.new.test()    # cleans to avoid workspace loaded with residual info
  tkfocus(top)
        
} # end TestScorer
      
set.tests.directory <- function() {
  # ===================================================
  # Set test directory and look for existing tests
  # 1 Allows user stopping TestScorer
  # 2 Looks for tests in directory
  # 3 If necessary creates one and copy necessary files
  # 4 Builds catalog of tests
  # ===================================================
  content <- dir()
  existing.tests <- content[grep('^TST_.+\\.r$', content)]
  if (length(existing.tests)==0) {
    # actual directory hasn't TST_ files: make window to choose dir
    choice.dir.top <- tktoplevel()
    tkwm.deiconify(choice.dir.top)
    tkwm.resizable(choice.dir.top, FALSE, FALSE)
    tkgrab.set(choice.dir.top)
    tktitle(choice.dir.top) <- 'Choose a directory'
    tkgrid(tklabel(choice.dir.top, text=''))
    tkgrid(tklabel(choice.dir.top, text='You have to choose or create a directory for working with TestScorer.'), sticky='w')
    tkgrid(tklabel(choice.dir.top, text=''))
    tkgrid(tklabel(choice.dir.top, text='If you create a new directory, six files will be copied to it:'), sticky='w')
    tkgrid(tklabel(choice.dir.top, text='     DASS.r: R script for scoring DASS scales'), sticky='w')
    tkgrid(tklabel(choice.dir.top, text='     MHLC.r: R script for scoring MHLC test'), sticky='w')
    tkgrid(tklabel(choice.dir.top, text='     RAAS.r: R script for scoring RAAS scales'), sticky='w')
    tkgrid(tklabel(choice.dir.top, text='     Help.txt: an ascii file with help instructions'), sticky='w')
    tkgrid(tklabel(choice.dir.top, text='     TestScorerHelp.pdf: PDF with further help'), sticky='w')
    tkgrid(tklabel(choice.dir.top, text='     .Profile: to open TestScorer when R is launched from this directory'), sticky='w')
    tkgrid(tklabel(choice.dir.top, text=''))
    
    choice.dirFrame <- tkframe(choice.dir.top, relief='groove',
                               borderwidth=2)
    tkgrid(choice.dirFrame)
    rb.continue <- tkradiobutton(choice.dir.top)
    rb.stop <- tkradiobutton(choice.dir.top)
    rbdirChoiceVar <- tclVar('continue')
    tkconfigure(rb.continue, variable=rbdirChoiceVar, value='continue')
    tkconfigure(rb.stop, variable=rbdirChoiceVar, value='stop')
    tkgrid(tklabel(choice.dirFrame, text='Choose an option:'), sticky='w')
    tkgrid(tklabel(choice.dirFrame, text='    Select an existing directory or create a new one'), rb.continue, sticky='w')
    tkgrid(tklabel(choice.dirFrame, text='    Stop & exit TestScorer'), rb.stop, sticky='w')
    
    tkgrid(tklabel(choice.dir.top, text=' '))
    doneChoiceDirVar <- tclVar(0) # stop continuing
    OkDoneChoiceDir.but <- tkbutton(choice.dir.top, text='   Ok   ',
                                    command=function() tclvalue(doneChoiceDirVar) <- 1)
    tkgrid(OkDoneChoiceDir.but)
    tkbind(choice.dir.top, '<Destroy>', function() tclvalue(doneChoiceDirVar) <- 3)
    tkfocus(choice.dir.top)
    tkwait.variable(doneChoiceDirVar)
    doneChoiceDir <- tclvalue(doneChoiceDirVar)
    choiceDir <- tclvalue(rbdirChoiceVar)
    tkgrab.release(choice.dir.top)
    tkdestroy(choice.dir.top)
    
    if (doneChoiceDir=='1' && choiceDir=='continue') {
      set.dir <- 'stop TestScorer'
      working.dir <- tk_choose.dir()
      if (working.dir!='')  { # a dir has been choosen or created
        setwd(working.dir)
        content <- dir() # looks for test files, if none exists copy examples
        if (length(grep('TST_.+\\.r$', content))==0) {
          file.copy(from=system.file('some.stuff/TST_DASS.r',
                                     package='TestScorer'), to=working.dir)
          file.copy(from=system.file('some.stuff/TST_MHLC.r',
                                     package='TestScorer'), to=working.dir)
          file.copy(from=system.file('some.stuff/TST_RAAS.r',
                                     package='TestScorer'), to=working.dir)
          file.copy(from=system.file('some.stuff/TST_IPIP50.r',
                                     package='TestScorer'), to=working.dir)
          file.copy(from=system.file('some.stuff/Help.txt',
                                     package='TestScorer'), to=working.dir)
          file.copy(from=system.file('doc/TestScorerHelp.pdf',
                                     package='TestScorer'), to=working.dir)
          # write .Rprofile
          cat(paste('# Profile to launch TestScorer',
                    '\nsuppressPackageStartupMessages(require(TestScorer))',
                    '\nTestScorerGUI()'),
              file='.Rprofile')
          
          if (.Platform$OS.type=='windows') {
            # write Test.Scorer.cmd
            if (.Platform$r_arch=='i386') startR <- paste(R.home(), '/bin/i386/Rgui.exe', sep='')
            else startR <- paste(R.home(), '/bin/x64/Rgui.exe', sep='')
            cat(paste('start ""',
                      startR,
                      '--no-restore --quiet'),
                file='TestScorer.cmd')
          }
        }
        set.dir <- 'test has been copied'
      }
    } else
      set.dir <- 'stop TestScorer'  # choosing a dir has been canceled
  } else
    set.dir <- 'tests already exist'
  
  # reads information to build the catalog list
  if (set.dir!='stop TestScorer') {  # a dir has been choosen
    content <- dir()
    existing.tests <- content[grep('^TST_.+\\.r$', content)]
    if (length(existing.tests)>0) {
      catalog <- list()
      for (i in seq_along(existing.tests)) {
        source(existing.tests[i]) # reads files to extract tests characteristics
        catalog <- c(catalog, list(testChar))  # read testChar load by source
        # checks for duplicated tests acronyms
        cat.length <- length(catalog)
        acronyms <- character()
        for (i in 1:cat.length) acronyms <- c(acronyms, catalog[[i]]$acronym)
        dupli <- acronyms[duplicated(acronyms)]
        if (length(dupli>0))
          tkmessageBox(message=paste('Some test in catalog are duplicated:', paste(dupli, collapse=', '),
                                     '\n\nConsiderer revising the test files in the directory.',
                                     '\nProbably an acronym do not match the file name.'),
                       icon='warning', type='ok')
      }
    }
  } else
    stop(paste('No directory has been choosen.',
               'TestScorer has been stopped by the user.',
               sep='\n'))
  catalog
} # end set.tests.directory
      
      
add.new.test <- function() {
  # =======================================================
  # add.new test
  #
  # returns: ret.values=list(continue.test, list(test.characteristics))
  # =======================================================
  
  test.char <- new.env(parent = emptyenv())  # environment for: acronymTestVal
  #                                          nameTestVal
  #                                          refVal
  #                                          commVal
  #                                          n.itemsVal
  #                                          validVal
  #                                          missVal
  #                                          reversedVal
  #                                          scoringVal
  #                                          prorrateVal
  #                                          n.scalesVal
  #                                          transVal
  #                                          equalVal
  #                                          graphVal
  #                                          n.age.groups
  #                                          name.transVal
  #                                          trans.is.numeric
  assign('trans.is.numeric', 'yes', envir=test.char)  # initialize
  scale.char <- new.env(parent = emptyenv())  # environment for: transTable
  assign('trans.is.numeric', 'yes', envir=scale.char)  # initialize
  
  make.test.window <- function() {
    no.trans <- function() {
      assign('graphVal', tclVar('no'), envir=test.char)
      assign('equalVal', tclVar('no'), envir=test.char)
      assign('name.transVal', tclVar('Centil'), envir=test.char)
      tkconfigure(rb.graph.yes, variable=test.char[['graphVal']], state='disabled')
      tkconfigure(rb.graph.no, variable=test.char[['graphVal']], state='disabled')
      tkconfigure(rb.equal.yes, variable=test.char[['equalVal']], state='disabled')
      tkconfigure(rb.equal.no, variable=test.char[['equalVal']], state='disabled')
      tkconfigure(entry.name, textvariable=test.char[['name.transVal']], state='disabled')
      assign('n.age.groupsVal', tclVar('1'), envir=test.char)
      tkconfigure(entry.n.age.groups, textvariable=test.char[['n.age.groupsVal']], state='disabled')
    } # end no.trans
    
    other.trans <- function() {
      if (tclvalue(test.char[['prorrateVal']])=='yes' || tclvalue(test.char[['scoringVal']])=='mean')
        tkmessageBox(message='NB: Raw scores would be rounded before transformation.\n        Considerer if rounding is appropiate.',
                     icon='warning', type='ok')
      if (tclvalue(test.char[['graphVal']])=='yes')
      { assign('graphVal', tclVar('no'), envir=test.char)
        tkmessageBox(message='NB: No graph for non-T transformation, graph option set to "no".',
                     icon='warning', type='ok')
      }
      tkconfigure(rb.graph.yes, variable=test.char[['graphVal']], state='disabled')
      tkconfigure(rb.graph.no, variable=test.char[['graphVal']], state='disabled')
      tkconfigure(rb.equal.yes, state='normal')
      tkconfigure(rb.equal.no, state='normal')
      assign('name.transVal', tclVar('Centil'), envir=test.char)
      tkconfigure(entry.name, textvariable=test.char[['name.transVal']], state='normal')
      assign('n.age.groupsVal', tclVar('1'), envir=test.char)
      tkconfigure(entry.n.age.groups, textvariable=test.char[['n.age.groupsVal']], state='normal')
    }  # end other.trans
    
    Tmean.trans <- function() {
      tkconfigure(rb.equal.yes, state='normal')
      tkconfigure(rb.equal.no, state='normal')
      tkconfigure(rb.graph.yes, state='normal')
      tkconfigure(rb.graph.no, state='normal')
      assign('name.transVal', tclVar('Centil'), envir=test.char)
      tkconfigure(entry.name, textvariable=test.char[['name.transVal']], state='disabled')
      assign('n.age.groupsVal', tclVar('1'), envir=test.char)
      tkconfigure(entry.n.age.groups, textvariable=test.char[['n.age.groupsVal']], state='normal')
    }  # end Tmean.trans
    
    Ttable.trans <- function() {
      if (tclvalue(test.char[['prorrateVal']])=='yes' || tclvalue(test.char[['scoringVal']])=='mean')
        tkmessageBox(message=paste('NB: Raw scores would be rounded before transformation.',
                                   '        Considerer if rounding is appropiate.',
                                   sep='\n'),
                     icon='warning', type='ok')
      tkconfigure(rb.equal.yes, state='normal')
      tkconfigure(rb.equal.no, state='normal')
      tkconfigure(rb.graph.yes, state='normal')
      tkconfigure(rb.graph.no, state='normal')
      assign('name.transVal', tclVar('Centil'), envir=test.char)
      tkconfigure(entry.name, textvariable=test.char[['name.transVal']], state='disabled')
      assign('n.age.groupsVal', tclVar('1'), envir=test.char)
      tkconfigure(entry.n.age.groups, textvariable=test.char[['n.age.groupsVal']], state='normal')
    }  # end Ttable.trans
    
    # Activate/deactivate prorrating of missings in case of sum
    raw.is.sum <- function() {
      tkconfigure(rb.prorat.no, state='normal')
      tkconfigure(rb.prorat.yes, state='normal')
    }  # end raw.is.sum
    
    raw.is.mean <- function() {
      if(tclvalue(test.char[['prorrateVal']])=='yes') {
        assign('prorrateVal', tclVar('no'), envir=test.char)
        tkconfigure(rb.prorat.no, variable=test.char[['prorrateVal']], state='disabled')
        tkconfigure(rb.prorat.yes, variable=test.char[['prorrateVal']], state='disabled')
        tkmessageBox(message='Prorrate has been reseted to "no".', icon='warning', type='ok')
      }
      if(tclvalue(test.char[['transVal']])!='no') {
        assign('transVal', tclVar('no'), envir=test.char)
        tkconfigure(rb.trans.No, variable=test.char[['transVal']], state='normal')
        tkconfigure(rb.trans.Tmean, variable=test.char[['transVal']], state='normal')
        tkconfigure(rb.trans.Ttable, variable=test.char[['transVal']], state='normal')
        tkconfigure(rb.trans.Other, variable=test.char[['transVal']], state='normal')
        assign('graphVal', tclVar('no'), envir=test.char)
        tkconfigure(rb.graph.no, variable=test.char[['graphVal']], state='disabled')
        tkconfigure(rb.graph.yes, variable=test.char[['graphVal']], state='disabled')
        assign('equalVal', tclVar('no'), envir=test.char)
        tkconfigure(rb.equal.no, variable=test.char[['equalVal']], state='disabled')
        tkconfigure(rb.equal.yes, variable=test.char[['equalVal']], state='disabled')
        tkmessageBox(message='Transform has been reseted to "no".', icon='warning', type='ok')
      }
    }  # end raw.is.mean
    
    prorrate.yes <- function() {
      if(tclvalue(test.char[['transVal']])!='no') {
        assign('transVal', tclVar('no'), envir=test.char)
        tkconfigure(rb.trans.No, variable=test.char[['transVal']], state='normal')
        tkconfigure(rb.trans.Tmean, variable=test.char[['transVal']], state='normal')
        tkconfigure(rb.trans.Ttable, variable=test.char[['transVal']], state='normal')
        tkconfigure(rb.trans.Other, variable=test.char[['transVal']], state='normal')
        assign('graphVal', tclVar('no'), envir=test.char)
        tkconfigure(rb.graph.no, variable=test.char[['graphVal']], state='disabled')
        tkconfigure(rb.graph.yes, variable=test.char[['graphVal']], state='disabled')
        assign('equalVal', tclVar('no'), envir=test.char)
        tkconfigure(rb.equal.no, variable=test.char[['equalVal']], state='disabled')
        tkconfigure(rb.equal.yes, variable=test.char[['equalVal']], state='disabled')
        tkmessageBox(message='Transform has been reseted to "no".', icon='warning', type='ok')
      }
    }  # end prorrate.yes
    
    validate.test <- function() {
      is.OkTest <- 'yes' # initialize as 'yes'
      options(warn=-1) # inhibit warnings
      
      # acronym
      acronymTest <- tclvalue(test.char[['acronymTestVal']])
      acronymTest <- sub(' +$', '', acronymTest)  # delte blanks
      acronymTest <- sub('^ +', '', acronymTest)  # delete blanks
      tmp <- grep('([[:alnum:]]|_)+$', acronymTest) # valid only alfanum & '_'
      if (length(tmp)==0) {
        is.OkTest <- 'no'
        msg <- paste('Acronym is not valid.',
                     '\nOnly alphanumeric and the character "_" are allowed.',
                     sep='\n')
        tkmessageBox(message=msg, icon="error", type="ok")
      }
      # already exists?
      content <- dir()
      existing.tests <- content[grep('TST_.+\\.r$', content)]
      existing.tests <- gsub('TST_', '', existing.tests)
      existing.tests <- gsub('.r$', '', existing.tests)
      if(acronymTest %in% existing.tests) {
        is.OkTest <- 'no'
        tkmessageBox(message="This test already exists. Change the acronym or cancel.",
                     icon="error", type="ok")
      }
      
      # number of items
      n.items <- as.numeric(tclvalue(test.char[['n.itemsVal']]))
      if (is.na(n.items) || length(grep('\\.', n.items)>0)) {  # detects missing or decimal
        is.OkTest <- 'no'
        tkmessageBox(message="Number of items is not valid.", icon="error", type="ok")
      }
      
      # valid answers
      valid <- tclvalue(test.char[['validVal']])
      valid <- as.numeric(unlist(strsplit(valid, ','))) # string to a numerical vector
      if (length(valid) == 0 ||        # detect missing
            any(is.na(valid)) ||         # id erroneus character (not a number)
            any(nchar(valid) > 1))  {    # id more than one character
        is.OkTest <- 'no'
        tkmessageBox(message="Valid answers are not valid.", icon="error", type="ok")
      }
      
      # missings
      miss <- tclvalue(test.char[['missVal']])
      miss <- as.numeric(unlist(strsplit(miss, ','))) # string to a numerical vector
      if (any(is.na(miss)) ||         # detect erroneus character (not a number)
            any(nchar(valid) > 1)) {    # id more than one character
        is.OkTest <- 'no'
        tkmessageBox(message="Missings are not valid.", icon="error", type="ok")
      }
      if (any(miss %in% valid)) {
        is.OkTest <- 'no'
        tkmessageBox(message="Don't include missing values into the valid values.", icon="error", type="ok")
      }
      
      # reversed items
      reversed <- tclvalue(test.char[['reversedVal']])
      reversed <- as.numeric(unlist(strsplit(reversed, ','))) # string to a numerical vector
      if (any(is.na(reversed) | any(grep('\\.', reversed)))) {  # detect erroneous characters or decimals
        is.OkTest <- 'no'
        tkmessageBox(message="Reversed items are not valid.", icon="error", type="ok")
      }
      # not greater than number of items
      if (length(reversed)!=0 && max(reversed, na.rm=TRUE)>n.items) {
        is.OkTest <- 'no'
        tkmessageBox(message="Some reversed item is greather than de number of items of the test.",
                     icon="error", type="ok")
      }
      
      # number of scales
      n.scales <- as.numeric(tclvalue(test.char[['n.scalesVal']]))
      if (is.na(n.scales)) {
        is.OkTest <- 'no'
        tkmessageBox(message="Number of scales is not valid.", icon="error", type="ok")
      }
   
      # number of age groups
      n.age.groups <- as.numeric(tclvalue(test.char[['n.age.groupsVal']]))
      if (is.na(n.age.groups) || grepl('\\.', n.age.groups)) {  # should be an integer
        is.OkTest <- 'no'
        tkmessageBox(message="Age-groups is not an integer.", icon="error", type="ok")
      }
      
      else if (n.age.groups == 0) {
        is.OkTest <- 'no'
        tkmessageBox(message="Age-groups should be at least 1.", icon="error", type="ok")
      }
      ######## NOT IMPLEMENTED YET
      else if (n.age.groups > 1) {
        # is.OkTest <- 'no'
        assign('n.age.groupsVal', tclVar('1'), envir=test.char)
        tkconfigure(entry.n.age.groups, textvariable=test.char[['n.age.groupsVal']], state='normal')
        tkmessageBox(message="Sorry, norms only implementet as yet for one group.", icon="error", type="ok")
      }
      
      options(warn=0)
      is.OkTest
    } # end validate.test
    
    OnOkTest <- function() {
      OkTest <- validate.test()
      if (OkTest=='yes')  tclvalue(doneTest) <- 1 # allow continuing
    }  # end OnOkTest
    
    # ---------- make window
    top.test <- tktoplevel()
    tkwm.title(top.test, "Test characteristics")
    tkwm.resizable(top.test, FALSE, FALSE)
    
    # Acronym
    frameAcronym <- tkframe(top.test)
    assign('acronymTestVal', tclVar(""), envir=test.char)
    entry.acronymTest <- tkentry(frameAcronym ,width="10", textvariable=test.char[['acronymTestVal']])
    tkgrid(frameAcronym, sticky='w')
    tkgrid(tklabel(frameAcronym, text="Enter an acronym for the test"), entry.acronymTest, sticky='w')
    tkgrid(tklabel(top.test, text="Note: Use only letters, numbers and '_'. The first character must be a letter."), sticky='w')
    
    # Name of the test
    tkgrid(tklabel(top.test, text=' '))
    frameName <- tkframe(top.test)
    assign('nameTestVal', tclVar(""), envir=test.char)
    entry.nameTest <- tkentry(frameName , width="40", textvariable=test.char[['nameTestVal']])
    tkgrid(frameName, sticky='w')
    tkgrid(tklabel(frameName, text="Enter the name of the test (optional)"), entry.nameTest, sticky='w')
    
    # Reference
    tkgrid(tklabel(top.test, text=' '))
    frameRef <- tkframe(top.test)
    assign('refVal', tclVar(""), envir=test.char)
    entry.ref <- tkentry(frameRef ,width="40", textvariable=test.char[['refVal']])
    tkgrid(frameRef, sticky='w')
    tkgrid(tklabel(frameRef, text="Enter the reference of the test (optional)"), entry.ref, sticky='w')
    
    # Comment
    tkgrid(tklabel(top.test, text=' '))
    frameComm <- tkframe(top.test)
    assign('commVal', tclVar(""), envir=test.char)
    entry.comm <- tkentry(frameComm ,width="40", textvariable=test.char[['commVal']])
    tkgrid(frameComm, sticky='w')
    tkgrid(tklabel(frameComm, text="Comment (optional)"), entry.comm, sticky='w')
    tkgrid(tklabel(top.test, text=paste("Note: Any optional comment (e.g.: Use only for ages greater than 18).")), sticky='w')
    
    # two frames
    frame.rigth <- tkframe(top.test, relief="groove", borderwidth=2)
    frame.left <- tkframe(top.test, relief="groove", borderwidth=2)
    frame.rigth <- tkframe(top.test, relief="groove", borderwidth=2)
    tkgrid(frame.left, frame.rigth, sticky='nw')
    
    # frame left = define items and number of scales
    # -------------------------
    # Number of items
    tkgrid(tklabel(frame.left, text=' '))
    frameN.items <- tkframe(frame.left)
    assign('n.itemsVal', tclVar(""), envir=test.char)
    entry.n.items <- tkentry(frameN.items ,width="10", textvariable=test.char[['n.itemsVal']])
    tkgrid(frameN.items, sticky='w')
    tkgrid(tklabel(frameN.items, text="Enter the number of items"), entry.n.items, sticky='w')
    
    # Valid answers
    tkgrid(tklabel(frame.left, text=' '))
    frameValid <- tkframe(frame.left)
    assign('validVal', tclVar(""), envir=test.char)
    entry.valid <- tkentry(frameValid ,width="20", textvariable=test.char[['validVal']])
    tkgrid(frameValid, sticky='w')
    tkgrid(tklabel(frameValid, text="Valid answers"), entry.valid, sticky='w')
    tkgrid(tklabel(frame.left, text=paste("Note: Can be any number.",
                                          "Do NOT include the missing values here",
                                          "(e.g.: 1,2,3).")), sticky='w')
    tkgrid(tklabel(frame.left, text="Use only letters separated by commas."), sticky='w')
    
    # Missing
    tkgrid(tklabel(frame.left, text=' '))
    frameMiss <- tkframe(frame.left)
    assign('missVal', tclVar(""), envir=test.char)
    entry.miss <- tkentry(frameMiss ,width="10", textvariable=test.char[['missVal']])
    tkgrid(frameMiss, sticky='w')
    tkgrid(tklabel(frameMiss, text="Missing answers"), entry.miss, sticky='w')
    tkgrid(tklabel(frame.left, text=paste("Note: Any numbers which represent a missing.",
                                          "'space' is automatically added as a missing character.")), sticky='w')
    tkgrid(tklabel(frame.left, text="Type only numbers separated by commas."), sticky='w')
    
    # Reversed items
    tkgrid(tklabel(frame.left, text=' '))
    frameReversed <- tkframe(frame.left)
    assign('reversedVal', tclVar(""), envir=test.char)
    entry.reversed <- tkentry(frameReversed, width="60", textvariable=test.char[['reversedVal']])
    tkgrid(frameReversed, sticky='w')
    tkgrid(tklabel(frameReversed, text="Reversed items"), entry.reversed, sticky='w')
    tkgrid(tklabel(frame.left, text=paste("Note: Left blank if none of the items should be reversed.",
                                          "Otherwise, use commas to separate the numbers.")), sticky='w')
    # Form of scoring the scales
    tkgrid(tklabel(frame.left, text=' '))
    frameScoring <- tkframe(frame.left)
    rb.sum <- tkradiobutton(frameScoring)
    rb.mean <- tkradiobutton(frameScoring)
    rb.prorat.no <- tkradiobutton(frameScoring)
    rb.prorat.yes <- tkradiobutton(frameScoring)
    assign('scoringVal', tclVar("sum"), envir=test.char)
    tkconfigure(rb.sum, variable=test.char[['scoringVal']], value="sum", command=raw.is.sum)
    tkconfigure(rb.mean, variable=test.char[['scoringVal']], value="mean", command=raw.is.mean)
    rb.prorat.no <- tkradiobutton(frameScoring)
    rb.prorat.yes <- tkradiobutton(frameScoring)
    assign('prorrateVal', tclVar("no"), envir=test.char)
    tkconfigure(rb.prorat.no, variable=test.char[['prorrateVal']], value="no")
    tkconfigure(rb.prorat.yes, variable=test.char[['prorrateVal']], value="yes", command=prorrate.yes)
    tkgrid(frameScoring, sticky='w')
    tkgrid(tklabel(frameScoring, text="Scoring Procedure?   "),
           tklabel(frameScoring, text="Mean"), rb.mean,
           tklabel(frameScoring, text="Add"),rb.sum,
           tklabel(frameScoring, text="If add, prorrate missings?   "),
           tklabel(frameScoring, text="No"), rb.prorat.no,
           tklabel(frameScoring, text="Yes"),rb.prorat.yes, sticky='w')
    
    # Number of scales
    tkgrid(tklabel(frame.left, text=' '))
    frameN.scales <- tkframe(frame.left)
    assign('n.scalesVal', tclVar(""), envir=test.char)
    entry.n.scales <- tkentry(frameN.scales ,width="2", textvariable=test.char[['n.scalesVal']])
    tkgrid(frameN.scales, sticky='w')
    tkgrid(tklabel(frameN.scales, text="Enter the number of scales"), entry.n.scales, sticky='w')
    
    # frame rigth = define transformation
    # -------------------------
    rb.trans.No <- tkradiobutton(frame.rigth)
    rb.trans.Tmean <- tkradiobutton(frame.rigth)
    rb.trans.Ttable <- tkradiobutton(frame.rigth)
    rb.trans.Other <- tkradiobutton(frame.rigth)
    rb.graph.yes <- tkradiobutton(frame.rigth)
    rb.graph.no <- tkradiobutton(frame.rigth)
    rb.equal.yes <- tkradiobutton(frame.rigth)
    rb.equal.no <- tkradiobutton(frame.rigth)
    assign('transVal', tclVar('no'), envir=test.char)
    assign('graphVal', tclVar('no'), envir=test.char)
    assign('equalVal', tclVar('no'), envir=test.char)
    assign('name.transVal', tclVar('Centil'), envir=test.char)
    assign('n.age.groupsVal', tclVar('1'), envir=test.char)
    entry.name <- tkentry(frame.rigth, width='5', textvariable=test.char[['name.transVal']], state='disabled')
    entry.n.age.groups <- tkentry(frame.rigth, width='1', textvariable=test.char[['n.age.groupsVal']], state='disabled')
    tkconfigure(rb.trans.No, variable=test.char[['transVal']], value="no", command=no.trans)
    tkconfigure(rb.trans.Tmean, variable=test.char[['transVal']], value="Tmean", command=Tmean.trans)
    tkconfigure(rb.trans.Ttable, variable=test.char[['transVal']], value="Ttable", command=Ttable.trans)
    tkconfigure(rb.trans.Other, variable=test.char[['transVal']], value="Other", command=other.trans)
    tkconfigure(rb.graph.yes, variable=test.char[['graphVal']], value="yes", state='disabled')
    tkconfigure(rb.graph.no, variable=test.char[['graphVal']], value="no", state='disabled')
    tkconfigure(rb.equal.yes, variable=test.char[['equalVal']], value="yes", state='disabled')
    tkconfigure(rb.equal.no, variable=test.char[['equalVal']], value="no", state='disabled')
    
    tkgrid(tklabel(frame.rigth, text="Choose one option for transforming raw scores:"), sticky='w')
    tkgrid(tklabel(frame.rigth, text="   1. Do NOT transform"), rb.trans.No, sticky='w')
    tkgrid(tklabel(frame.rigth, text="   2. T'scores using mean & sd"), rb.trans.Tmean, sticky='w')
    tkgrid(tklabel(frame.rigth, text="   3. T'scores using a table"), rb.trans.Ttable, sticky='w')
    tkgrid(tklabel(frame.rigth, text="   4. Other transformation using a table"), rb.trans.Other, sticky='w')
    tkgrid(tklabel(frame.rigth, text="        Choose a name, e.g.: 'Centil'"),
           entry.name, sticky='w')
    tkgrid(tklabel(frame.rigth, text=""), sticky='w')
    tkgrid(tklabel(frame.rigth, text="    Number of age-groups norms"),
           entry.n.age.groups, sticky='w')
    tkgrid(tklabel(frame.rigth, text="    Transformation is equal for both sexes?"),
           tklabel(frame.rigth, text='Yes'), rb.equal.yes,
           tklabel(frame.rigth, text='No'), rb.equal.no, sticky='w')
    tkgrid(tklabel(frame.rigth, text="    Add a pseudo-graph for T'scores?"),
           tklabel(frame.rigth, text='Yes'), rb.graph.yes,
           tklabel(frame.rigth, text='No'), rb.graph.no, sticky='w')
    
    #Test OK?
    doneTest <- tclVar(0) # to stop continuing
    tkgrid(tklabel(top.test, text=''))
    OkTest.but <-tkbutton(top.test, text="   OK   ", command=OnOkTest) # validate
    tkgrid(OkTest.but, columnspan=2)
    tkgrid(tklabel(top.test,
                   text='Close the window to cancel adding a new test.'),
           columnspan=2)
    tkbind(top.test, "<Destroy>", function() tclvalue(doneTest) <- 3)
    tkfocus(top.test)
    tkwait.variable(doneTest) # stops until ok is pressed and data are valid or windows is clossed
    if (tclvalue(doneTest)==1) {
      OkTest <- 'yes'
      acronymTest <- tclvalue(test.char[['acronymTestVal']])
      acronymTest <- sub(' +$', '', acronymTest)  # delte blanks
      acronymTest <- sub('^ +', '', acronymTest)  # delete blanks
      n.items <- as.numeric(tclvalue(test.char[['n.itemsVal']]))
      valid <- tclvalue(test.char[['validVal']])
      valid <- as.numeric(unlist(strsplit(valid, ','))) # string to a numerical vector
      miss <- tclvalue(test.char[['missVal']])
      miss <- as.numeric(unlist(strsplit(miss, ','))) # string to a numerical vector
      reversed <- tclvalue(test.char[['reversedVal']])
      reversed <- as.numeric(unlist(strsplit(reversed, ','))) # string to a numerical vector
      n.scales <- as.numeric(tclvalue(test.char[['n.scalesVal']]))
      test.characteristics <- list(acronymTest <- acronymTest,
                                   nameTest <- tclvalue(test.char[['nameTestVal']]),
                                   ref <- tclvalue(test.char[['refVal']]),
                                   comm <- tclvalue(test.char[['commVal']]),
                                   n.items <- n.items,
                                   valid <- valid,
                                   miss <- miss,
                                   reversed <- reversed,
                                   scoring <- tclvalue(test.char[['scoringVal']]),
                                   prorrate <- tclvalue(test.char[['prorrateVal']]),
                                   n.scales <- n.scales,
                                   trans <- tclvalue(test.char[['transVal']]),
                                   equal <- tclvalue(test.char[['equalVal']]),
                                   graph <- tclvalue(test.char[['graphVal']]),
                                   name.trans <- tclvalue(test.char[['name.transVal']]),
                                   n.age.agroups <- tclvalue(test.char[['n.age.groupsVal']]))
    } else {
      OkTest <- 'no' # windows has been calcelled
      test.characteristics <- as.list(rep('', 15))
    }
    
    tkdestroy(top.test)
    list(OkTest=OkTest, test.characteristics=test.characteristics)
  } # end make.test.window
  
  write.test.begining <- function(acronymTest, nameTest, ref, n.items, valid,
                                  miss, comm, reversed, trans, graph) {
    cat(paste('#',
              acronymTest,
              'scale scoring script'),
        file='tmp')
    cat(paste('\n# Creation date:',
              Sys.Date()),
        file='tmp', append=TRUE)
    cat('\n# --------------\n', file='tmp', append=TRUE)
    cat(paste('\ntestChar <- list(acronym = "',
              acronymTest,
              '",',
              sep=''),
        file='tmp', append=TRUE)
    
    cat(paste('\n                 name = "',
              nameTest,
              '",',
              sep=''),
        file='tmp', append=TRUE)
    
    cat(paste('\n                 ref = "',
              ref,
              '",',
              sep=''),
        file='tmp', append=TRUE)
    
    cat(paste('\n                 n.items = ',
              n.items,
              ',',
              sep=''),
        file='tmp', append=TRUE)
    
    cat('\n                 valid = c(', file='tmp', append=TRUE)
    
    if (length(valid)>1)
      for (i in 1:(length(valid)-1)) cat(paste(valid[i], ', ', sep=''), file='tmp', append=TRUE)
    cat(paste(valid[length(valid)],
              '),',
              sep=''),
        file='tmp', append=TRUE)
    
    cat('\n                 miss = c(', file='tmp', append=TRUE)
    if (length(miss)>1)
      for (i in 1:(length(miss)-1)) cat(paste(miss[i], ', ', sep=''), file='tmp', append=TRUE)
    cat(paste(miss[length(miss)],
              '),',
              sep=''),
        file='tmp', append=TRUE)
    cat(paste('\n                 comm = "',
              comm,
              '")',
              sep=''),
        file='tmp', append=TRUE)
    
    # begining of scoring function
    cat('\n\nscoring.fun <- function(answers, sex, age, id, date.test, comm) {', file='tmp', append=TRUE)
    cat('\n  # "answer" is a *character* vector as introduced through the keyboard.', file='tmp', append=TRUE)
    cat('\n  # "sex", "age", "id", "date.test" & "comm" can be used if appropiate.', file='tmp', append=TRUE)
    cat('\n  answers <- as.numeric(answers) # transform to numeric for easier scoring', file='tmp', append=TRUE)
    
    if (length(miss)!=0) { # some miss have been defined
      cat(paste('\n  answers[answers %in% c(',
                miss,
                ')] <- NA # missing characters to NA',
                sep=''),
          file='tmp', append=TRUE)
    }
    cat('\n  blanks <- sum(is.na(answers)) # compute number of missing items', file='tmp', append=TRUE)
    cat(paste('\n  pcnt.blanks <- round((blanks / ',
              n.items,
              ') * 100) # compute % of missings',
              sep=''),
        file='tmp', append=TRUE)
    
    if (sum(nchar(reversed))!=0) {  # some rev have been defined
      cat('\n  # Items which should be reversed', file='tmp', append=TRUE)
      cat('\n  reversed.items <- c(', file='tmp', append=TRUE)
      if (length(reversed)>1)
        for (i in 1:(length(reversed)-1)) cat(paste(reversed[i], ', ', sep=''), file='tmp', append=TRUE)
      cat(paste(reversed[length(reversed)],
                ')',
                sep=''),
          file='tmp', append=TRUE)
      max.answ <- max(valid)
      cat(paste('\n  answers[reversed.items] <- (',
                max.answ,
                ' + 1) - answers[reversed.items]',
                sep=''),
          file='tmp', append=TRUE)
    }
    
    cat('\n  results <- data.frame(NULL) # Initialize null data frame for results', file='tmp', append=TRUE)
    
    if (trans=='Tmean') {
      cat('\n\n  toT <- function(raw.score, mean, sd) {  # compute T score' , file='tmp', append=TRUE)
      cat('\n    T.score <- round(((raw.score - mean) / sd) * 10 + 50)' , file='tmp', append=TRUE)
      cat('\n    T.score' , file='tmp', append=TRUE)
      cat('\n  } # end toT' , file='tmp', append=TRUE)
    }
    if (graph=='yes') {
      
      cat('\n\n  makeGraph <- function(T.score) {  # make a graph' , file='tmp', append=TRUE)
      cat('\n    template  <- "|    :    :    |    :    |    :    |    :    :    |"' , file='tmp', append=TRUE)
      cat('\n    options(warn=-1)', file='tmp', append=TRUE)
      cat('\n    T.score <- as.integer(T.score)  # NA in case of non-numerical string', file='tmp', append=TRUE)
      cat('\n    options(warn=0)', file='tmp', append=TRUE)
      cat('\n    if (!is.na(T.score)) {', file='tmp', append=TRUE)
      cat('\n      if (T.score < 0) T.score <- 0' , file='tmp', append=TRUE)
      cat('\n        else if (T.score > 100) T.score <- 100' , file='tmp', append=TRUE)
      cat('\n      position <- round((T.score / 2) + 1)' , file='tmp', append=TRUE)
      cat(paste('\n      graph <- paste(substr(template, 1, position-1),',
                '                     substr(template, position + 1, nchar(template)),',
                '                     sep="o")  # "o" marks the position',
                sep='\n'),
          file='tmp', append=TRUE)
      cat('\n    } else {' , file='tmp', append=TRUE)
      cat('\n      graph <- "Not graphicable"',
          file='tmp', append=TRUE)
      cat('\n    }', file='tmp', append=TRUE)
      cat('\n    graph', file='tmp', append=TRUE)
      cat('\n  } # end makeGraph' , file='tmp', append=TRUE)
    }
  } # end write.test.begining
  
  make.scale.window <- function(scale.num, n.items, valid, trans, equal, graph, name.trans) {
    ValidateScaleChar <- function() {
      is.OkScale='yes'
      options(warn=-1)
      # acronym
      acronymScale <- tclvalue(acronymScaleVal)
      acronymScale <- sub(' +$', '', acronymScale)  # delte blanks
      acronymScale <- sub('^ +', '', acronymScale)  # delete blanks
      tmp <- grep('^[[:alpha:]]([[:alnum:]]|_)*$', acronymScale) # valid only alfanum & '_'
      if (length(tmp)==0) {
        is.OkScale <- 'no'
        msg <- paste('Acronym is not valid.',
                     '\nThe first character should be a letter',
                     'Only alphanumeric and the character "_" are allowed.',
                     sep='\n')
        tkmessageBox(message=msg, icon="error", type="ok")
      }
      
      # items
      tmp <- tclvalue(itemsVal)
      items <- as.numeric(unlist(strsplit(tmp, ','))) # string to a numerical vector
      if (length(items) == 0 ||              # detect missing
            any(is.na(items)) ||               # id erroneous character
            any(grep('\\.', items))) {         # detect decimal point
        is.OkScale <- 'no'
        tkconfigure(entry.items, state='normal')
        tkmessageBox(message="Items are not valid.", icon="error", type="ok")
      }
      if (max(items, na.rm=TRUE)>n.items) {
        is.OkScale <- 'no'
        tkconfigure(entry.items, state='normal') 
        tkmessageBox(message="Some value is greather than the number of items of the test.",
                     icon="error", type="ok")
      }
      
      # no transform table
      if (trans=='no' || trans=='Tmean')
        assign('tableDone', 'no table is needed', envir=scale.char)
      
      # transform Tmean
      if (trans=='Tmean') {
        if (equal=='yes') {
          # mean
          mean.1 <- as.numeric(tclvalue(mean.1Val))
          if(is.na(mean.1)) {
            is.OkScale <- 'no'
            tkmessageBox(message="Mean is not valid.", icon="error", type="ok")
          }
          # sd
          sd.1 <- as.numeric(tclvalue(sd.1Val))
          if(is.na(sd.1)) {
            is.OkScale <- 'no'
            tkmessageBox(message="SD is not valid.", icon="error", type="ok")
          }
        } else {
          # mean males
          mean.males <- as.numeric(tclvalue(mean.malesVal))
          if(is.na(mean.malesVal)) {
            is.OkScale <- 'no'
            tkmessageBox(message="Mean for males is not valid.",
                         icon="error", type="ok")
          }
          # sd males
          sd.males <- as.numeric(tclvalue(sd.malesVal))
          if(is.na(sd.malesVal)) {
            is.OkScale <- 'no'
            tkmessageBox(message="SD for males is not valid.",
                         icon="error", type="ok")
          }
          # mean females
          mean.females <- as.numeric(tclvalue(mean.femalesVal))
          if(is.na(mean.femalesVal)) {
            is.OkScale <- 'no'
            tkmessageBox(message="Mean for females is not valid.",
                         icon="error", type="ok")
          }
          # sd females
          sd.females <- as.numeric(tclvalue(sd.femalesVal))
          if(is.na(sd.femalesVal)) {
            is.OkScale <- 'no'
            tkmessageBox(message="SD for females is not valid.",
                         icon="error", type="ok")
          }
        }
      }
      
      # transTable
      if (scale.char[['tableDone']]=='no') {
        is.OkScale <- 'no'
        tkmessageBox(message="You need open and fill the transformation table.",
                     icon="error", type="ok")
      } else
        if (any(grepl('^-$', unlist(scale.char[['transTable']]))) ||
              any(grepl('^$', unlist(scale.char[['transTable']])))) {
          is.OkScale <- 'no'
          msg <- paste('All or some of the transformations values are missing',
                       'Open and refill the table.',
                       sep='\n')
          tkmessageBox(message=msg, icon="error", type="ok")
        }
      
      options(warn=0)
      is.OkScale
    }  # end ValidateScaleChar
    
    # ---- window for scale
    top.scale <- tktoplevel()
    tkwm.title(top.scale, paste("Definition of scale number", scale.num))
    tkwm.resizable(top.scale, FALSE, FALSE)
    
    # Acronym
    frameAcronym <- tkframe(top.scale)
    acronymScaleVal <- tclVar("")
    entry.acronymScale <- tkentry(frameAcronym ,width="10", textvariable=acronymScaleVal)
    tkgrid(frameAcronym, sticky='w')
    tkgrid(tklabel(frameAcronym, text="Enter an acronym for the scale"), entry.acronymScale, sticky='w')
    tkgrid(tklabel(top.scale, text="Note: Use only letters, numbers and '_'. The first character must be a letter."), sticky='w')
    
    # Description
    tkgrid(tklabel(top.scale, text=' '))
    frameName <- tkframe(top.scale)
    nameScaleVal <- tclVar("")
    entry.name <- tkentry(frameName ,width="40", textvariable=nameScaleVal)
    tkgrid(frameName, sticky='w')
    tkgrid(tklabel(frameName, text="Enter the scale's name (optional)"), entry.name, sticky='w')
    
    # Items making up the scale
    tkgrid(tklabel(top.scale, text=' '))
    frameItems <- tkframe(top.scale)
    itemsVal <- tclVar("")
    entry.items <- tkentry(frameItems, width="60", textvariable=itemsVal)
    tkgrid(frameItems, sticky='w')
    tkgrid(tklabel(frameItems, text="Items making up the scale"), entry.items, sticky='w')
    tkgrid(tklabel(top.scale, text="Note: Use commas to separate the numbers."), sticky='w')
    
    # transformation type
    if (trans=='Tmean') {
      if (equal=='yes') {
        tkgrid(tklabel(top.scale, text=' '))
        frameMean.sd <- tkframe(top.scale)
        mean.1Val <- tclVar("")
        sd.1Val <- tclVar("")
        entry.mean <- tkentry(frameMean.sd ,width="8", textvariable=mean.1Val)
        entry.sd <- tkentry(frameMean.sd ,width="8", textvariable=sd.1Val)
        tkgrid(frameMean.sd, sticky='w')
        tkgrid(tklabel(frameMean.sd, text="Mean"),entry.mean,
               tklabel(frameMean.sd, text="   SD"), entry.sd, sticky='w')
      } else {
        tkgrid(tklabel(top.scale, text=' '))
        frameMean.sd.males <- tkframe(top.scale)
        mean.malesVal <- tclVar("")
        sd.malesVal <- tclVar("")
        entry.mean.males <- tkentry(frameMean.sd.males ,width="8", textvariable=mean.malesVal)
        entry.sd.males <- tkentry(frameMean.sd.males ,width="8", textvariable=sd.malesVal)
        tkgrid(frameMean.sd.males, sticky='w')
        tkgrid(tklabel(frameMean.sd.males, text="Males:     Mean"), entry.mean.males,
               tklabel(frameMean.sd.males, text="   SD"), entry.sd.males,
               sticky='w')
        frameMean.sd.females <- tkframe(top.scale)
        mean.femalesVal <- tclVar("")
        sd.femalesVal <- tclVar("")
        entry.mean.females <- tkentry(frameMean.sd.females ,width="8", textvariable=mean.femalesVal)
        entry.sd.females <- tkentry(frameMean.sd.females ,width="8", textvariable=sd.femalesVal)
        tkgrid(frameMean.sd.females, sticky='w')
        tkgrid(tklabel(frameMean.sd.females, text="Females: Mean"), entry.mean.females,
               tklabel(frameMean.sd.females, text="   SD"), entry.sd.females,
               sticky='w')
        
      }
    } else if (trans=='Ttable' || trans=='Other') {
      #---  transformation table
      make.table.win <- function(valid, items, equal) {
        tkconfigure(entry.items, state='disabled') # disable editing items
        
        n.items <- length(unlist(strsplit(items, ',')))
        min.score <- min(valid)*n.items
        max.score <- max(valid)*n.items
        n.scores <- max.score - min.score + 1
        
        if (equal=='yes') {
          matrixTrans <- matrix(c(as.character(min.score:max.score), rep('-', n.scores)),
                                nrow=n.scores,
                                ncol=2)
          matrixTrans <- rbind(c('Raw', 'Transf'), matrixTrans)
          
          tclArrayTrans <- tclArray()
          for (i in (1:(n.scores+1)))
            for (j in (1:2))
              tclArrayTrans[[i-1,j-1]] <- matrixTrans[i,j]
          nrow <- n.scores+1
          ncol <- 2
          
        } else {
          matrixTrans <- matrix(c(as.character(min.score:max.score), rep('-', 2*(n.scores))),
                                nrow=n.scores,
                                ncol=3)
          matrixTrans <- rbind(c('Raw', 'Male', 'Female'), matrixTrans)
          
          tclArrayTrans <- tclArray()
          for (i in (1:(n.scores+1)))
            for (j in (1:3))
              tclArrayTrans[[i-1,j-1]] <- matrixTrans[i,j]
          nrow <- n.scores+1
          ncol <- 3
        }
        
        # transform table windows
        trans.win <- tktoplevel()
        tkwm.resizable(trans.win, FALSE, FALSE)
        tkgrab.set(trans.win)
        tableFrame <- tkframe(trans.win)
        tkgrid(tableFrame)
        tclRequire("Tktable")
        tkwm.title(trans.win, 'Transform')
        table <- tkwidget(tableFrame, "table",
                          rows=nrow, cols=ncol,
                          titlerows=1, titlecols=1,
                          height=16, width=ncol+1,
                          yscrollcommand=function(...) tkset(yscr,...))
        yscr <- tkscrollbar(tableFrame, command=function(...)tkyview(table,...))
        tkgrid(table, yscr)
        tkgrid.configure(yscr, sticky="nsw")
        tkconfigure(table, variable=tclArrayTrans, selectmode="single", background="white")
        tkgrid(tklabel(trans.win, text='Note: "-" sign will be automatically removed.'))
        
        # Ok trans table
        tkgrid(tklabel(trans.win, text=' '))
        doneTransVal <- tclVar(0) # to stop continuing
        OkTrans.but <- tkbutton(trans.win, text="   OK   ", command=function() tclvalue(doneTransVal) <- 1)
        tkgrid(OkTrans.but)
        tkgrid(tklabel(trans.win, text='Close the window to cancel the table.'))
        tkgrid(tklabel(trans.win, text='Canceling the table allows the reedition of items.'))
        tkbind(trans.win, "<Destroy>", function() tclvalue(doneTransVal) <- 3)
        tkfocus(trans.win)
        tkwait.variable(doneTransVal) # stops until ok is pressed
        
        if (equal=='yes') {
          transTable <- character(n.scores)
          for (i in 1:n.scores) transTable[i] <- tclvalue(tclArrayTrans[[i, 1]])
        }
        else {
          transTableMales <- character(n.scores)
          for (i in 1:n.scores) transTableMales[i] <- tclvalue(tclArrayTrans[[i, 1]])
          transTableFemales <- character(n.scores)
          for (i in 1:n.scores) transTableFemales[i] <- tclvalue(tclArrayTrans[[i, 2]])
        }
        
        if (tclvalue(doneTransVal)==3) {
          tkconfigure(entry.items, state='normal')
          msg <- paste('The table has been canceled.',
                       'All the transformations are set to "-".',
                       'Edition of items is enabled.',
                       sep='\n')
          tkmessageBox(message=msg, icon='warning', type='ok')
        }
        
        tkdestroy(trans.win)
        
        if (equal=='yes') {
          transTable <- list(transTable)
        } else
          transTable <- list(transTableMales, transTableFemales)
        assign('transTable', transTable, envir=scale.char)
        assign('tableDone', 'yes', envir=scale.char)
        
      } # end make.table.win
      
      tkgrid(tklabel(top.scale, text=' '))
      defTableFrame <- tkframe(top.scale)
      assign('tableDone', 'no', envir=scale.char)
      openTable.but <- tkbutton(top.scale, text='Open',
                                command=function() make.table.win(valid, tclvalue(itemsVal), equal))
      tkgrid.configure(defTableFrame, sticky='w')
      tkgrid(tklabel(defTableFrame,
                     text = 'Transformation table:   '),
             openTable.but,
             sticky='w')
      tkgrid(tklabel(top.scale,
                     text='Note: To avoid mismatch between number of items and number of raw scores,'),
             sticky='w')
      tkgrid(tklabel(top.scale,
                     text='           when the table is openned, the edition of items is blocked'),
             sticky='w')
      tkgrid(tklabel(top.scale,
                     text='           and the table is reseted accordingly the number of items.'),
             sticky='w')
    }
    
    OnOkScale <- function() {
      OkScale <- ValidateScaleChar()
      if (OkScale=='yes') tclvalue(doneScale) <- 1 # allow continuing
    }  # end OnOkScale
    
    # OK?
    tkgrid(tklabel(top.scale, text=' '))
    doneScale <- tclVar(0) # to stop continuing
    scale.butFrame <- tkframe(top.scale)
    tkgrid(scale.butFrame)
    OkScale.but <- tkbutton(scale.butFrame,
                            text="        OK        ",
                            command=OnOkScale)
    cancelScale.but <- tkbutton(scale.butFrame,
                                text="Cancel scale",
                                command=function() tclvalue(doneScale) <- 2)
    cancelTest.but<- tkbutton(scale.butFrame,
                              text=" Cancel test ",
                              command=function() tclvalue(doneScale) <- 3)
    tkgrid(OkScale.but, cancelScale.but, cancelTest.but)
    tkbind(top.scale, '<Destroy>', function() tclvalue(doneScale) <- 2)
    tkfocus(top.scale)
    tkwait.variable(doneScale) # stop continuing untill ok is pressed an data are valid
    
    if (tclvalue(doneScale)==1) {  # scale ok
      acronymScale <- tclvalue(acronymScaleVal)
      acronymScale <- sub(' +$', '', acronymScale)  # delte blanks
      acronymScale <- sub('^ +', '', acronymScale)  # delete blanks
      nameScale <- tclvalue(nameScaleVal)
      items <- tclvalue(itemsVal)
      items <- as.numeric(unlist(strsplit(items, ','))) # string to a numerical vector
      mean.1 <- NA   # initialitze to NA
      sd.1 <- NA
      mean.males <- NA
      mean.females <- NA
      sd.males <- NA
      sd.females <- NA
      trans.1 <- NA
      trans.males <- NA
      trans.females <- NA
      if(trans=='Tmean' && equal=='yes') {
        mean.1 <- as.numeric(tclvalue(mean.1Val))
        sd.1 <- as.numeric(tclvalue(sd.1Val))
      }
      else if (trans=='Tmean' && equal=='no') {
        mean.males <- as.numeric(tclvalue(mean.malesVal))
        sd.males <- as.numeric(tclvalue(sd.malesVal))
        mean.females <- as.numeric(tclvalue(mean.femalesVal))
        sd.females <- as.numeric(tclvalue(sd.femalesVal))
      }
      else if ((trans=='Ttable' || trans=='Other') && equal=='yes') {
        trans.1 <- scale.char[['transTable']]
        trans.1 <- trans.1[[1]]
        trans.1 <- gsub('(^-)|(-$)', '', trans.1)
        if (sum(!grepl('^[0-9]+$', trans.1))>0) {
          assign('trans.is.numeric', 'no', envir=test.char)
          assign('trans.is.numeric', 'no', envir=scale.char)
        } else 
          assign('trans.is.numeric', 'yes', envir=scale.char)
      }
      else if ((trans=='Ttable' || trans=='Other') && equal=='no') {
        trans <- scale.char[['transTable']]
        trans.males <- trans[[1]]
        trans.females <- trans[[2]]
        trans.males <- gsub('(^-)|(-$)', '', trans.males)  # del -
        trans.females <- gsub('(^-)|(-$)', '', trans.females)
        if (sum(!grepl('^[0-9]+$', trans.males))>0 || sum(!grepl('^[0-9]+$', trans.females))>0) {
          assign('trans.is.numeric', 'no', envir=test.char)
          assign('trans.is.numeric', 'no', envir=scale.char)
        } else
          assign('trans.is.numeric', 'yes', envir=test.char)
      }
      
      OkScale <- 'yes'
      n.items <- length(items)
      min.score <- min(valid)*n.items
      scale.characteristics <- list(acronymScale=acronymScale,
                                    nameScale=nameScale,
                                    items=items,
                                    mean.1=mean.1,
                                    sd.1=sd.1,
                                    mean.males=mean.males,
                                    mean.females=mean.females,
                                    sd.males=sd.males,
                                    sd.females=sd.females,
                                    trans.1=trans.1,
                                    trans.males=trans.males,
                                    trans.females=trans.females,
                                    min.score=min.score)
      
    } else  # cancel scale
      if (tclvalue(doneScale)==2) {
        OkScale <- 'repeat scale'
        scale.characteristics <- list('')
        
      } else {   # cancel test
        OkScale <- 'no'
        scale.characteristics <- list('')
      }
    
    tkdestroy(top.scale)
    list(OkScale=OkScale, scale.characteristics=scale.characteristics)
  } # end make.scale.window
  
  write.scale <- function(acronymScale, nameScale, items,
                          scale.num, scoring, prorrate, trans, equal, name.trans, graph,
                          mean.1, sd.1,
                          mean.males, mean.females, sd.males, sd.females,
                          trans.1,
                          trans.males, trans.females, min.score) {
    cat(paste('\n\n  #',
              acronymScale,
              'scale scoring commands'),
        file='tmp', append=TRUE)
    cat('\n  # --------------', file='tmp', append=TRUE)
    cat(paste('\n  results[',
              scale.num,
              ', "Acronym"] <- "',
              acronymScale,
              '" # acronym',
              sep=''),
        file='tmp', append=TRUE)
    
    cat(paste('\n  results[',
              scale.num,
              ', "Scale"] <- "',
              nameScale,
              '" # name of scale',
              sep=''),
        file='tmp', append=TRUE)
    
    cat('\n  # Items making up the scale', file='tmp', append=TRUE)
    cat('\n  items <- c(', file='tmp', append=TRUE)
    
    if (length(items)>1)
      for (i in 1:(length(items)-1)) cat(paste(items[i], ', ', sep=''), file='tmp', append=TRUE)
    cat(paste(items[length(items)],
              ')',
              sep=''),
        file='tmp', append=TRUE)
    cat(paste('\n  results[',
              scale.num,
              ', "Miss"] <- sum(is.na(answers[items])) # number of missings',
              sep=''),
        file='tmp', append=TRUE)
    
    # if all answ are missing: RAW <- NA
    cat(paste('\n  if (results[',
              scale.num,
              ', "Miss"] == length(items)) {  # all answers are missings',
              sep=''),
        file='tmp', append=TRUE)
    
    cat(paste('\n    results[',
              scale.num,
              ', "Raw"] <- NA',
              sep=''),
        file='tmp', append=TRUE)
    
    if (trans=='Tmean' || trans=='Ttable') {
      cat(paste('\n    results[',
                scale.num,
                ', "T"] <- NA',
                sep=''),
          file='tmp', append=TRUE)
    }
    
    if (trans=='Other') {
      cat(paste('\n    results[',
                scale.num,
                ', "',
                name.trans,
                '"] <- NA',
                sep=''),
          file='tmp', append=TRUE)
    }
    
    if (graph=='yes')
      cat(paste('\n    results[',
                scale.num,
                ', "Graph"] <- "All items missing"',
                sep=''),
          file='tmp', append=TRUE)
    
    cat('\n  } else {',
        file='tmp', append=TRUE)
    
    # three scoring option: sum no prorrate, sum & prorrate, mean
    if (scoring=='sum' && prorrate=='no')
      cat(paste('\n    results[', 
                scale.num,
                ', "Raw"] <- sum(answers[items], na.rm=TRUE) # sum answered items',
                sep=''),
          file='tmp', append=TRUE)
    else if (scoring=='sum' && prorrate=='yes')
      cat(paste('\n    results[',
                scale.num,
                ', "Raw"] <- round(mean(answers[items], na.rm=TRUE)*length(items), 2) # sum pondering missings',
                sep=''),
          file='tmp', append=TRUE)
    else
      cat(paste('\n    results[', 
                scale.num,
                ', "Raw"] <- round(mean(answers[items], na.rm=TRUE), 2) # mean answered items',
                sep=''),
          file='tmp', append=TRUE)
    
    # transformation type
    if(trans=='Tmean' && equal=='no') {
      cat(paste('\n    if (sex=="Male") results[',
                scale.num,
                ', "T"] <- toT(results[',
                scale.num,
                ', "Raw"], ',
                mean.males,
                ', ',
                sd.males,
                ') # compute T score',
                sep=''),
          file='tmp', append=TRUE)
      cat(paste('\n      else results[',
                scale.num,
                ', "T"] <- toT(results[',
                scale.num,
                ', "Raw"], ',
                mean.females,
                ', ',
                sd.females,
                ')',
                sep=''),
          file='tmp', append=TRUE)
      
    } else if(trans=='Tmean' && equal=='yes') {
      cat(paste('\n    results[',
                scale.num,
                ', "T"] <- toT(results[',
                scale.num,
                ', "Raw"], ',
                mean.1,
                ', ',
                sd.1,
                ') # compute T score',
                sep=''),
          file='tmp', append=TRUE)
      
    } else if ((trans=='Ttable' || trans=='Other') && equal=='no') {
      # build vector trans.males.tabel
      cat('\n    # Vector for score transformation', file='tmp', append=TRUE)
      if (scale.char[['trans.is.numeric']]=='yes') {
        cat('\n    trans.males.table <- c(', file='tmp', append=TRUE)
        if (length(trans.males)>1)
          for (i in 1:(length(trans.males)-1)) cat(paste(trans.males[i], ', ', sep=''), file='tmp', append=TRUE)
        cat(paste(trans.males[length(trans.males)],
                  ')',
                  sep=''),
            file='tmp', append=TRUE)
      } else {
        cat('\n    trans.males.table <- c(', file='tmp', append=TRUE)
        if (length(trans.males)>1)
          for (i in 1:(length(trans.males)-1)) cat(paste('"', trans.males[i], '", ', sep=''), file='tmp', append=TRUE)
        cat(paste('"',
                  trans.males[length(trans.males)],
                  '")',
                  sep=''),
            file='tmp', append=TRUE)
      }
      # build vector trans.females.table
      if (scale.char[['trans.is.numeric']]=='yes') {
        cat('\n    trans.females.table <- c(', file='tmp', append=TRUE)
        if (length(trans.females)>1)
          for (i in 1:(length(trans.females)-1)) cat(paste(trans.females[i], ', ', sep=''), file='tmp', append=TRUE)
        cat(paste(trans.females[length(trans.females)],
                  ')',
                  sep=''),
            file='tmp', append=TRUE)
      } else {
        cat('\n    trans.females.table <- c(', file='tmp', append=TRUE)
        if (length(trans.females)>1)
          for (i in 1:(length(trans.females)-1)) cat(paste('"', trans.females[i], '", ', sep=''), file='tmp', append=TRUE)
        cat(paste('"',
                  trans.females[length(trans.females)],
                  '")',
                  sep=''),
            file='tmp', append=TRUE)
      }
      
      # build index
      if (scoring=='sum' && prorrate=='no')
        cat(paste('\n    index <- results[',
                  scale.num,
                  ', "Raw"]',
                  sep=''),
            file='tmp', append=TRUE)
      else
        cat(paste('\n    index <- round(results[',
                  scale.num,
                  ', "Raw"])  # rounds raw score',
                  sep=''),
            file='tmp', append=TRUE)
      cat(paste('\n    index <- ',
                1 - min.score,
                ' + index  # ',
                1 - min.score,
                ', because raw scores begin at ',
                min.score,
                sep=''),
          file='tmp', append=TRUE)
      
      # build transformation
      cat(paste('\n    if (sex=="Male") results[',
                scale.num,
                ', "',
                sep=''),
          file='tmp', append=TRUE)
      if (trans=='Ttable')
        cat('T', file='tmp', append=TRUE)
      else
        cat(name.trans , file='tmp', append=TRUE)
      cat('"] <- trans.males.table[index]', file='tmp', append=TRUE)
      cat(paste('\n      else results[',
                scale.num,
                ', "',
                sep=''),
          file='tmp', append=TRUE)
      if (trans=='Ttable')
        cat('T', file='tmp', append=TRUE)
      else
        cat(name.trans , file='tmp', append=TRUE)
      cat('"] <- trans.females.table[index]', file='tmp', append=TRUE)
      
    } else if ((trans=='Ttable' || trans=='Other') && equal=='yes') {
      # build vector trans.table
      cat('\n  # Vector for score transformation', file='tmp', append=TRUE)
      if (scale.char[['trans.is.numeric']]=='yes') {
        cat('\n    trans.table <- c(', file='tmp', append=TRUE)
        if (length(trans.1)>1)
          for (i in 1:(length(trans.1)-1)) cat(paste(trans.1[i], ', ', sep=''), file='tmp', append=TRUE)
        cat(paste(trans.1[length(trans.1)], ')', sep=''), file='tmp', append=TRUE)
      } else {
        cat('\n    trans.table <- c(', file='tmp', append=TRUE)
        if (length(trans.1)>1)
          for (i in 1:(length(trans.1)-1)) cat(paste('"', trans.1[i], '", ', sep=''), file='tmp', append=TRUE)
        cat(paste('"', trans.1[length(trans.1)], '")', sep=''), file='tmp', append=TRUE)
      }
      # build index
      if (scoring=='sum' && prorrate=='no')
        cat(paste('\n    index <- results[',
                  scale.num,
                  ', "Raw"]',
                  sep=''),
            file='tmp', append=TRUE)
      else
        cat(paste('\n    index <- round(results[',
                  scale.num,
                  ', "Raw"])  # rounds raw score',
                  sep=''),
            file='tmp', append=TRUE)
      cat(paste('\n    index <- ',
                1 - min.score,
                ' + index  # ',
                1 - min.score,
                ', because raw scores begin at ',
                min.score,
                sep=''),
          file='tmp', append=TRUE)
      
      # build transformation
      cat(paste('\n    results[',
                scale.num,
                ', "',
                sep=''),
          file='tmp', append=TRUE)
      if (trans=='Ttable')
        cat('T', file='tmp', append=TRUE)
      else
        cat(name.trans , file='tmp', append=TRUE)
      cat('"] <- trans.table[index]', file='tmp', append=TRUE)
    }
    
    # graph
    if(graph=='yes')
      cat(paste('\n    results[',
                scale.num,
                ', "Graph"] <- makeGraph(results[',
                scale.num,
                ', "T"]) # makes the graph',
                sep=''),
          file='tmp', append=TRUE)
    
    cat('\n  }', file='tmp', append=TRUE)
    
  } # end write.scale
  
  write.test.end <- function(acronymTest, trans, graph, name.trans) {
    cat('\n\n  # Vector for writing scores to a file', file='tmp', append=TRUE)
    cat('\n  # --------------------', file='tmp', append=TRUE)
    if (graph=='no') {
      cat('\n  results.scores <- unlist(results[-c(1, 2)])  # not Acronym & Scale columns',
          file='tmp', append=TRUE)
    } else {
      cat('\n  results.scores <- unlist(results[-c(1, 2, 6)])  # not Acronym, Scale & Graph columns',
          file='tmp', append=TRUE)
    }
    cat('\n  names <- paste(results$Acronym, names(results.scores), sep=".")',
        file='tmp', append=TRUE)
    cat('\n  names(results.scores) <- sub("[0-9]+$", "", names)  # delete ending numbers',
        file='tmp', append=TRUE)
    
    if (graph=='yes') {
      cat('\n\n  # Ruler for graph column', file='tmp', append=TRUE)
      cat('\n  # --------------------', file='tmp', append=TRUE)
      cat('\n  names(results)[6] <- "0    10   20   30   40   50   60   70   80   90  100"',
          file='tmp', append=TRUE)
    }
    
    cat('\n\n  # Output in form of list', file='tmp', append=TRUE)
    cat('\n  # ------------------', file='tmp', append=TRUE)
    cat('\n  results.lst <- list(paste("Total number of missing items: ",',
        file='tmp', append=TRUE)
    cat('\n                            blanks,',
        file='tmp', append=TRUE)
    cat('\n                            " (",',
        file='tmp', append=TRUE)
    cat('\n                            pcnt.blanks,',
        file='tmp', append=TRUE)
    cat('\n                            "%)",',
        file='tmp', append=TRUE)
    cat('\n                            sep=""))',
        file='tmp', append=TRUE)
    
    cat('\n\n  # Return results', file='tmp', append=TRUE)
    cat('\n  # ------------------', file='tmp', append=TRUE)
    cat('\n  list(results.lst = results.lst,', file='tmp', append=TRUE)
    cat('\n       results.df = results,', file='tmp', append=TRUE)
    cat('\n       results.scores = results.scores)', file='tmp', append=TRUE)
    
    cat('\n\n} # end of scoring.fun', file='tmp', append=TRUE)
    
    # change the name of the temporal file
    file.name <- paste('TST_', acronymTest, '.r', sep='')
    file.rename(from='tmp', to=file.name)
    msg <- paste("Test definition and scoring script have been saved as TST_",
                 acronymTest, '.r', sep='')
    tkmessageBox(message=msg)
  } # end write.test.end
  
  # ----- main add.new.test()
  ret.values <- make.test.window()
  continue.test <- ret.values[[1]]
  test.characteristics <- ret.values[[2]]
  acronymTest <- test.characteristics[[1]]
  nameTest    <- test.characteristics[[2]]
  ref         <- test.characteristics[[3]]
  comm        <- test.characteristics[[4]]
  n.items     <- test.characteristics[[5]]
  valid       <- test.characteristics[[6]]
  miss        <- test.characteristics[[7]]
  reversed    <- test.characteristics[[8]]
  scoring     <- test.characteristics[[9]]
  prorrate    <- test.characteristics[[10]]
  n.scales    <- test.characteristics[[11]]
  trans       <- test.characteristics[[12]]
  equal       <- test.characteristics[[13]]
  graph       <- test.characteristics[[14]]
  name.trans  <- test.characteristics[[15]]
  
  if (continue.test=='yes')
    write.test.begining(acronymTest, nameTest, ref, n.items, valid,
                        miss, comm, reversed, trans, graph)
  
  if (continue.test=='yes') {
    scale.num <-0
    while ((scale.num < n.scales) && (continue.test == 'yes')) {
      scale.num <- scale.num + 1
      ret.values <- make.scale.window(scale.num, n.items, valid, trans,
                                      equal, graph, name.trans)
      OkScale <- ret.values[[1]]
      scale.characteristics <- ret.values[[2]]
      acronymScale <- scale.characteristics['acronymScale']
      nameScale <- scale.characteristics['nameScale']
      items <- unlist(scale.characteristics['items'])
      mean.1 <- scale.characteristics['mean.1']
      sd.1 <- scale.characteristics['sd.1']
      mean.males <- scale.characteristics['mean.males']
      mean.females <- scale.characteristics['mean.females']
      sd.males <- scale.characteristics['sd.males']
      sd.females <- scale.characteristics['sd.females']
      trans.1 <- unlist(scale.characteristics['trans.1'])
      trans.males <- unlist(scale.characteristics['trans.males'])
      trans.females <- unlist(scale.characteristics['trans.females'])
      min.score <- unlist(scale.characteristics['min.score'])
      if (OkScale=='yes')
        write.scale(acronymScale, nameScale, items,
                    scale.num, scoring, prorrate, trans,
                    equal, name.trans, graph,
                    mean.1, sd.1,
                    mean.males, mean.females,
                    sd.males, sd.females,
                    trans.1,
                    trans.males, trans.females, min.score)
      else
        if (OkScale=='repeat scale') {
          scale.num <- scale.num - 1
          continue.test <- 'yes'
        }
      else
        continue.test <- 'no' # test has been canceled
    } # end while
    
  } # end if continue.test==yes (build scales)
  
  if (continue.test=='yes')
    write.test.end(acronymTest, trans, graph, name.trans)
  
  if (continue.test=='no') {
    msg <- "Test definition has been canceled."
    tkmessageBox(message=msg, icon='warning', type='ok')
    options(warn=-1) # inhibit warnings
    try(dummy <- file.remove('tmp'))
    options(warn=0)
  }
  TestScorerGUI() # reopenss the main GUI, updating the catalog
  
} # end add.new.test ======================================
      
delete.test <- function() {
  # ==================================================================================
  # delete.test
  #
  # 1 Looks at directory for existing tests
  # 2 Delete one of them
  # ==================================================================================
  content <- dir() # looks for test files
  existing.tests <- content[grep('TST_.+\\.r$', content)]
  existing.tests <- gsub('TST_', '', existing.tests)
  delete.top <- tktoplevel()
  tkwm.title(delete.top, 'Delete a test')
  tkwm.resizable(delete.top, FALSE, FALSE)
  scr <- tkscrollbar(delete.top, repeatinterval=5,
                     command=function(...)tkyview(test.selBox,...))
  test.selBox <- tklistbox(delete.top, height=3, width='31',
                           selectmode='single', yscrollcommand=function(...)tkset(scr,...),
                           background='white', exportselection=F)
  for (i in (seq_along(existing.tests)))
    tkinsert(test.selBox,'end', existing.tests[i])  # test acronyms
  tkselection.set(test.selBox,0)  # Default test (rem: indexing starts at zero)
  tkgrid(tklabel(delete.top, text='Which test would you like to delete?'), sticky='nw')
  tkgrid(test.selBox, scr, sticky='nw')
  tkgrid.configure(scr, rowspan=5, sticky='nse')
  # Ok delete test
  tkgrid(tklabel(delete.top, text=' '))
  doneDeleteTestVal <- tclVar(0) # to stop continuing
  OkDeleteTest.but <- tkbutton(delete.top, text="   OK   ",
                               command=function() tclvalue(doneDeleteTestVal) <- 1)
  tkgrid(OkDeleteTest.but)
  tkgrid(tklabel(delete.top, text='Close the window to cancel deleting a test.'))
  tkbind(delete.top, "<Destroy>", function() tclvalue(doneDeleteTestVal) <- 3)
  tkfocus(delete.top)
  tkwait.variable(doneDeleteTestVal) # stops until ok is pressed and data are valid or windowa is clossed
  if (tclvalue(doneDeleteTestVal)==1) {
    test.to.delete <- existing.tests[as.numeric(tkcurselection(test.selBox))+1]
    test.to.delete <- sub('.r$', '', test.to.delete)
    file.to.rename <- paste('TST_', test.to.delete, '.r', sep='')
    file.new.name <- paste('TST_', test.to.delete, '.bak', sep='')
    file.rename(from=file.to.rename, to=file.new.name)
    tkdestroy(delete.top)
    tkmessageBox(message=paste('Extension type has been changed to ".bak".',
                               'This test will not appear in TestScorer.',
                               sep='\n'),
                 icon='warning', type='ok')
  }
  
  TestScorerGUI() # reopenss the main GUI, updating the catalog
} # end delete test
      