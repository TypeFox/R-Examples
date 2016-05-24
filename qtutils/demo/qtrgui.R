

## QMainWindow

require("qtutils")

qtrgui <- function()
{
    w <- Qt$QMainWindow()

    mdi <- Qt$QMdiArea()
    w$setCentralWidget(mdi)

    repl <- qrepl()
    grdev <- QT()
    objbrowser <- qbrowser()
    misctab <- Qt$QTabWidget()
    misctab$setTabsClosable(TRUE)
    qconnect(misctab, "tabCloseRequested",
             function(i) {
                 misctab$removeTab(i)
                 gc()
             })

    misctab$addTab(objbrowser, "Objects")

    mdi$addSubWindow(grdev)
    mdi$addSubWindow(misctab)
    mdi$addSubWindow(repl)

    ## FIXME: reset par(oldopt) when QUI is closed (Ctrl+Q shortcut?)
    oldopt <- options(pager =
                      function(...) {
                          page <- qpager(...)
                          ipage <- misctab$addTab(page, title = page$windowTitle)
                          misctab$setCurrentIndex(ipage)
                      },
                      browser =
                      function(...) {
                          web <- qwebbrowser(...)
                          iweb <- misctab$addTab(web, web$title)
                          misctab$setCurrentIndex(iweb)
                          qconnect(web, "titleChanged",
                                   function(s) {
                                       i <- misctab$indexOf(web)
                                       misctab$setTabText(i, s)
                                   })
                      },
                      error = qrecover,
                      help_type = "html")
    help.start()
    mdi$tileSubWindows()
    w$resize(600, 800)
    w
}

qtrgui()$show()
