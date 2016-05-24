############################################
### code chunk number 173: Wizard
###################################################
wizard <- Qt$QWizard()
wizard$setWindowTitle("A wizard")


###################################################
### code chunk number 174: Dialogs.Rnw:401-407
###################################################
get_age_page <- Qt$QWizardPage(wizard)
get_age_page$setTitle("What is your age?")
layout <- Qt$QFormLayout()
get_age_page$setLayout(layout)
layout$addRow("Age", (age <- Qt$QLineEdit()))
wizard$addPage(get_age_page)


###################################################
### code chunk number 175: Dialogs.Rnw:411-424
###################################################
get_toys_page <- Qt$QWizardPage(wizard)
get_toys_page$setTitle("What toys do you like?")
layout <- Qt$QFormLayout()
get_toys_page$setLayout(layout)
layout$addRow("Toys", (toys <- Qt$QLineEdit()))
wizard$addPage(get_toys_page)
##
get_games_page <- Qt$QWizardPage(wizard)
get_games_page$setTitle("What games do you like?")
layout <- Qt$QFormLayout()
get_games_page$setLayout(layout)
layout$addRow("Games", (games <- Qt$QLineEdit()))
wizard$addPage(get_games_page)


###################################################
### code chunk number 176: Dialogs.Rnw:429-432
###################################################
response <- wizard$exec()
if(response)
  message(toys$text)
