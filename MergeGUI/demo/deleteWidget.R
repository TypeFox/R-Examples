library(gWidgetsRGtk2)

testEnv = new.env()

testEnv$msg = 'text'

testEnv$g0 = ggroup(horizontal = FALSE, container = gwindow('demo: delete widgets'))

testEnv$g1 = gbutton('delete', container = testEnv$g0, handler = function(h, ...) {
    if (!is.null(testEnv$g2)) {
        delete(testEnv$g0, testEnv$g2)
        testEnv$msg = paste(testEnv$msg, 'text')
        testEnv$g2 = gtext(testEnv$msg, container = testEnv$g0)
    }
})

testEnv$g2 = gtext('delete me', container = testEnv$g0)
