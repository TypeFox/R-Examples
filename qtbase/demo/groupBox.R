## Group boxes with different types of buttons

qsetClass("Window", Qt$QWidget, function(parent = NULL) {
  super(parent)
  
  grid <- Qt$QGridLayout()
  ## NOTE: the layout does not take ownership of the widgets, so we
  ## need to assign the layout to our widget up-front. Our widget then
  ## takes ownership of the widgets in the layout.
  setLayout(grid) 
  grid$addWidget(createFirstExclusiveGroup(), 0, 0)
  grid$addWidget(createSecondExclusiveGroup(), 1, 0)
  grid$addWidget(createNonExclusiveGroup(), 0, 1)
  grid$addWidget(createPushButtonGroup(), 1, 1)
  
  setWindowTitle("Group Boxes")
  resize(480, 320)
})

qsetMethod("createFirstExclusiveGroup", Window, function() {
  groupBox <- Qt$QGroupBox("Exclusive Radio Buttons")

  radio1 <- Qt$QRadioButton("&Radio button 1")
  radio2 <- Qt$QRadioButton("R&adio button 2")
  radio3 <- Qt$QRadioButton("Ra&dio button 3")

  radio1$setChecked(TRUE)

  vbox <- Qt$QVBoxLayout()
  vbox$addWidget(radio1)
  vbox$addWidget(radio2)
  vbox$addWidget(radio3)
  vbox$addStretch(1)
  groupBox$setLayout(vbox)

  groupBox
}, "private")

qsetMethod("createSecondExclusiveGroup", Window, function() {
  groupBox <- Qt$QGroupBox("E&xclusive Radio Buttons")
  groupBox$setCheckable(TRUE)
  groupBox$setChecked(FALSE)

  radio1 <- Qt$QRadioButton("Rad&io button 1")
  radio2 <- Qt$QRadioButton("Radi&o button 2")
  radio3 <- Qt$QRadioButton("Radio &button 3")
  radio1$setChecked(TRUE)
  checkBox <- Qt$QCheckBox("Ind&ependent checkbox")
  checkBox$setChecked(TRUE)

  vbox <- Qt$QVBoxLayout()
  vbox$addWidget(radio1)
  vbox$addWidget(radio2)
  vbox$addWidget(radio3)
  vbox$addWidget(checkBox)
  vbox$addStretch(1)
  groupBox$setLayout(vbox)

  groupBox
}, "private")

qsetMethod("createNonExclusiveGroup", Window, function() {
  groupBox <- Qt$QGroupBox("Non-Exclusive Checkboxes")
  groupBox$setFlat(TRUE)

  checkBox1 <- Qt$QCheckBox("&Checkbox 1")
  checkBox2 <- Qt$QCheckBox("C&heckbox 2")
  checkBox2$setChecked(TRUE)
  tristateBox <- Qt$QCheckBox("Tri-&state button")
  tristateBox$setTristate(TRUE)
  tristateBox$setCheckState(Qt$Qt$PartiallyChecked)

  vbox <- Qt$QVBoxLayout()
  vbox$addWidget(checkBox1)
  vbox$addWidget(checkBox2)
  vbox$addWidget(tristateBox)
  vbox$addStretch(1)
  groupBox$setLayout(vbox)

  groupBox
}, "private")

qsetMethod("createPushButtonGroup", Window, function() {
  groupBox <- Qt$QGroupBox("&Push Buttons")
  groupBox$setCheckable(TRUE)
  groupBox$setChecked(TRUE)

  pushButton <- Qt$QPushButton("&Normal Button")
  toggleButton <- Qt$QPushButton("&Toggle Button")
  toggleButton$setCheckable(TRUE)
  toggleButton$setChecked(TRUE)
  flatButton <- Qt$QPushButton("&Flat Button")
  flatButton$setFlat(TRUE)

  popupButton <- Qt$QPushButton("Pop&up Button")
  menu <- Qt$QMenu(this)
  menu$addAction("&First Item")
  menu$addAction("&Second Item")
  menu$addAction("&Third Item")
  menu$addAction("F&ourth Item")
  popupButton$setMenu(menu)

  newAction <- menu$addAction("Submenu")
  subMenu <- Qt$QMenu("Popup Submenu")
  subMenu$addAction("Item 1")
  subMenu$addAction("Item 2")
  subMenu$addAction("Item 3")
  newAction$setMenu(subMenu)

  vbox <- Qt$QVBoxLayout()
  vbox$addWidget(pushButton)
  vbox$addWidget(toggleButton)
  vbox$addWidget(flatButton)
  vbox$addWidget(popupButton)
  vbox$addStretch(1)
  groupBox$setLayout(vbox)

  groupBox
}, "private")

Window()$show()
