

###################################################
### code chunk number 170: extend-date-dialog
###################################################
qsetClass("DateDialog", Qt$QDialog, 
          function(parent = NULL) {
            super(parent=parent)
            setWindowTitle("Choose a date")
            this$calendar <- Qt$QCalendarWidget()
            #
            btn_box <- 
              Qt$QDialogButtonBox(Qt$QMessageBox$Cancel | 
                                  Qt$QMessageBox$Ok)
            qconnect(btn_box, "accepted", function() {
              this$close()
              this$setResult(Qt$QMessageBox$Ok)    
            })
            qconnect(btn_box, "rejected", 
                     function() this$close())
            #
            layout <- Qt$QVBoxLayout()
            sapply(list(calendar, btn_box), layout$addWidget)
            setLayout(layout)
          })


###################################################
### code chunk number 171: extend-date-get
###################################################
qsetMethod("selectedDate", DateDialog, 
           function(x) calendar$selectedDate$toString())


###################################################
### code chunk number 172: extend-date-exec
###################################################
date_dialog <- DateDialog()
if (date_dialog$exec())
  message(date_dialog$selectedDate())
