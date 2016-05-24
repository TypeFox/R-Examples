## Adapted from Qt's AddressBook tutorial
library(qtbase)

qsetClass("AddressBook", Qt$QWidget, function(parent = NULL) {
  super(parent)
  
  nameLabel <- Qt$QLabel("Name:")
  this$nameLine <- Qt$QLineEdit()
  nameLine$readOnly <- TRUE

  addressLabel <- Qt$QLabel("Address:")
  this$addressText <- Qt$QTextEdit()
  addressText$readOnly <- TRUE

  this$addButton <- Qt$QPushButton("&Add")
  addButton$show()
  this$submitButton <- Qt$QPushButton("&Submit")
  submitButton$hide()
  this$cancelButton <- Qt$QPushButton("&Cancel")
  cancelButton$hide()
  this$editButton <- Qt$QPushButton("&Edit")
  editButton$enabled <- FALSE
  this$removeButton <- Qt$QPushButton("&Remove")
  removeButton$enabled <- FALSE
  this$findButton <- Qt$QPushButton("&Find")
  findButton$enabled <- FALSE
  this$loadButton <- Qt$QPushButton("&Load...");
  loadButton$toolTip <- "Load contacts from a file"
  this$saveButton <- Qt$QPushButton("Sa&ve...")
  saveButton$toolTip <- "Save contacts to a file"
  saveButton$enabled <- FALSE

  qconnect(addButton, "clicked", addContact)
  qconnect(submitButton, "clicked", submitContact)
  qconnect(cancelButton, "clicked", cancel)
  qconnect(editButton, "clicked", editContact)
  qconnect(removeButton, "clicked", removeContact)
  qconnect(findButton, "clicked", findContact)
  qconnect(loadButton, "clicked", loadFromFile)
  qconnect(saveButton, "clicked", saveToFile)

  buttonLayout1 <- Qt$QVBoxLayout()
  buttonLayout1$addWidget(addButton, Qt$Qt$AlignTop)
  buttonLayout1$addWidget(editButton)
  buttonLayout1$addWidget(removeButton)
  buttonLayout1$addWidget(findButton)
  buttonLayout1$addWidget(submitButton)
  buttonLayout1$addWidget(cancelButton)
  buttonLayout1$addWidget(loadButton)
  buttonLayout1$addWidget(saveButton)
  buttonLayout1$addStretch()
  
  this$nextButton <- Qt$QPushButton("&Next")
  nextButton$enabled <- FALSE
  this$previousButton <- Qt$QPushButton("&Previous")
  previousButton$enabled <- FALSE

  qconnect(nextButton, "clicked", nextContact)
  qconnect(previousButton, "clicked", prevContact)

  buttonLayout2 <- Qt$QHBoxLayout()
  buttonLayout2$addWidget(previousButton)
  buttonLayout2$addWidget(nextButton)
  
  mainLayout <- Qt$QGridLayout()
  mainLayout$addWidget(nameLabel, 0, 0)
  mainLayout$addWidget(nameLine, 0, 1)
  mainLayout$addWidget(addressLabel, 1, 0, Qt$Qt$AlignTop)
  mainLayout$addWidget(addressText, 1, 1)
  mainLayout$addLayout(buttonLayout1, 1, 2)
  mainLayout$addLayout(buttonLayout2, 3, 1)
  
  setLayout(mainLayout)
  setWindowTitle("Simple Address Book")

  this$contacts <- character()
  this$currentMode <-
    factor("navigation", levels = c("navigation", "adding", "editing"))
  this$dialog <- FindDialog()
})

qsetMethod("addContact", AddressBook, function() {
  this$oldName <- nameLine$text
  this$oldAddress <- addressText$plainText

  nameLine$clear()
  addressText$clear()

  updateInterface("adding")
})

qsetMethod("submitContact", AddressBook, function() {
  name <- nameLine$text
  address <- addressText$plainText

  if (name == "" || address == "") {
    Qt$QMessageBox$information(this, "Empty Field",
                               "Please enter a name and address.");
    return;
  }

  if (currentMode == "adding") {
    if (is.na(contacts[name])) {
      contacts[name] <<- address
      msg <- sprintf("\"%s\" has been added to your address book.", name)
      Qt$QMessageBox$information(this, "Add Successful", msg)
    } else {
      msg <- sprintf("Sorry, \"%s\" is already in your address book.", name)
      Qt$QMessageBox$information(this, "Add Unsuccessful", msg)
      return;
    }
  } else if (currentMode == "editing") {
    if (oldName != name) {
      if (is.na(contacts[name])) {
        contacts <<- contacts[names(contacts) != oldName]
        contacts[name] <<- address
        msg <- sprintf("\"%s\" has been added to your address book.", name)
        Qt$QMessageBox$information(this, "Edit Successful", msg)
      } else {
        msg <- sprintf("Sorry, \"%s\" is already in your address book.", name)
        Qt$QMessageBox$information(this, "Edit Unsuccessful", msg)
        return;
      }
    } else if (oldAddress != address) {
      msg <- sprintf("\"%s\" has been edited in your address book.", name)
      Qt$QMessageBox$information(this, "Edit Successful", msg)
      contacts[name] <<- address
    }
  }

  updateInterface("navigation")
})

qsetMethod("cancel", AddressBook, function() {
  nameLine$text <- oldName
  addressText$plainText <- oldAddress

  updateInterface("navigation")
})

qsetMethod("nextContact", AddressBook, function() {
  i <- match(nameLine$text, names(contacts), length(contacts))
  if (i == length(contacts))
    i <- 1
  else i <- i + 1

  nameLine$text <- names(contacts)[i]
  addressText$plainText <- contacts[i]
})

qsetMethod("prevContact", AddressBook, function() {
  i <- match(nameLine$text, names(contacts))

  if (is.na(i)) {
    nameLine$clear()
    addressText$clear()
    return()
  }
    
  if (i == 1)
    i <- length(contacts)
  else i <- i - 1
  
  nameLine$text <- names(contacts)[i]
  addressText$plainText <- contacts[i]
})

qsetMethod("editContact", AddressBook, function() {
  oldName <<- nameLine$text
  oldAddress <<- addressText$plainText

  updateInterface("editing")
})

qsetMethod("removeContact", AddressBook, function() {
  name <- nameLine$text
  address <- addressText$plainText

  if (!is.na(contacts[name])) {
    msg <- sprintf("Are you sure you want to remove \"%s\"?", name)
    button <- Qt$QMessageBox$question(this, "Confirm Remove", msg,
                                      Qt$QMessageBox$Yes | Qt$QMessageBox$No)
    if (button == Qt$QMessageBox$Yes) {
      prevContact()
      contacts <<- contacts[names(contacts) != name]
      msg <- sprintf("\"%s\" has been removed from your address book.", name)
      Qt$QMessageBox$information(this, "Remove Successful", msg)
    }
  }

  updateInterface("navigation")
})

qsetMethod("findContact", AddressBook, function() {
  dialog$show()

  if (dialog$exec() == Qt$QDialog$Accepted) {
    contactName <- dialog$getFindText()

    if (!is.na(contacts[contactName])) {
      nameLine$text <- contactName
      addressText$plainText <- contacts[contactName]
    } else {
      msg <- sprintf("Sorry, \"%s\" is not in your address book.", contactName)
      Qt$QMessageBox$information(this, "Contact Not Found", msg)
      return()
    }
  }

  updateInterface("navigation")
})

qsetMethod("saveToFile", AddressBook, function() {
  fileName <-
    Qt$QFileDialog$getSaveFileName(this, "Save Address Book", "",
                                   "Address Book (*.abk);;All Files (*)")
  if (!nchar(fileName))
    return()
  else save(contacts, file = fileName)
})

qsetMethod("loadFromFile", AddressBook, function() {
  fileName <-
    Qt$QFileDialog$getOpenFileName(this, "Open Address Book", "",
                                   "Address Book (*.abk);;All Files (*)")
  if (!nchar(fileName))
    return()
  contacts <<- get(load(fileName))

  if (!length(contacts)) {
    msg <- "The file you are attempting to open contains no contacts."
    Qt$QMessageBox$information(this, "No contacts in file", msg)
  } else {
    nameLine$text <- names(contacts)[1]
    addressText$plainText <- contacts[1]
  }
  updateInterface("navigation")
})

qsetMethod("updateInterface", AddressBook, function(mode)
{
  currentMode[] <<- mode
  switch (as.character(currentMode),
          adding =,
          editing = {            
            nameLine$readOnly <- FALSE
            nameLine$setFocus()
            addressText$readOnly <- FALSE

            addButton$enabled <- FALSE
            editButton$enabled <- FALSE
            removeButton$enabled <- FALSE

            nextButton$enabled <- FALSE
            previousButton$enabled <- FALSE

            submitButton$show()
            cancelButton$show()

            loadButton$enabled <- FALSE
            saveButton$enabled <- FALSE
          },
          navigation = {
            if (!length(contacts)) {
              nameLine$clear()
              addressText$clear()
            }
            nameLine$readOnly <- TRUE
            addressText$readOnly <- TRUE
            addButton$enabled <- TRUE

            number <- length(contacts)
            editButton$enabled <- number >= 1
            removeButton$enabled <- number >= 1
            findButton$enabled <- number > 2
            nextButton$enabled <- number > 1
            previousButton$enabled <- number > 1

            submitButton$hide()
            cancelButton$hide()

            loadButton$enabled <- TRUE
            saveButton$enabled <- number >= 1
          })
}, "private")

qsetClass("FindDialog", Qt$QDialog, function(parent = NULL) {
  super(parent)
  
  findLabel <- Qt$QLabel("Enter the name of a contact:")
  this$lineEdit <- Qt$QLineEdit()

  this$findButton <- Qt$QPushButton("&Find")
  this$findText <- ""

  layout <- Qt$QHBoxLayout()
  layout$addWidget(findLabel)
  layout$addWidget(lineEdit)
  layout$addWidget(findButton)

  setLayout(layout)
  setWindowTitle("Find a Contact")
  
  qconnect(findButton, "clicked", findClicked)
  qconnect(findButton, "clicked", accept)
})

qsetMethod("findClicked", FindDialog, function() {
  text <- lineEdit$text

  if (!nchar(text)) {
    Qt$QMessageBox$information(this, "Empty Field", "Please enter a name.")
    return;
  } else {
    findText <<- text
    lineEdit$clear()
    hide()
  }
})

qsetMethod("getFindText", FindDialog, function() findText)

ab <- AddressBook()
ab$show()

