## Features of line edits
qsetClass("Window", Qt$QWidget, function() {
  echoGroup <- Qt$QGroupBox("Echo")

  echoLabel <- Qt$QLabel("Mode:")
  echoComboBox <- Qt$QComboBox()
  echoComboBox$addItem("Normal")
  echoComboBox$addItem("Password")
  echoComboBox$addItem("PasswordEchoOnEdit")
  echoComboBox$addItem("No Echo")

  this$echoLineEdit <- Qt$QLineEdit()
  echoLineEdit$setFocus()

  validatorGroup <- Qt$QGroupBox("Validator")

  validatorLabel <- Qt$QLabel("Type:")
  validatorComboBox <- Qt$QComboBox()
  validatorComboBox$addItem("No validator")
  validatorComboBox$addItem("Integer validator")
  validatorComboBox$addItem("Double validator")

  this$validatorLineEdit <- Qt$QLineEdit()

  alignmentGroup <- Qt$QGroupBox("Alignment")

  alignmentLabel <- Qt$QLabel("Type:")
  alignmentComboBox <- Qt$QComboBox()
  alignmentComboBox$addItem("Left")
  alignmentComboBox$addItem("Centered")
  alignmentComboBox$addItem("Right")

  this$alignmentLineEdit <- Qt$QLineEdit()

  inputMaskGroup <- Qt$QGroupBox("Input mask")

  inputMaskLabel <- Qt$QLabel("Type:")
  inputMaskComboBox <- Qt$QComboBox()
  inputMaskComboBox$addItem("No mask")
  inputMaskComboBox$addItem("Phone number")
  inputMaskComboBox$addItem("ISO date")
  inputMaskComboBox$addItem("License key")

  this$inputMaskLineEdit <- Qt$QLineEdit()

  accessGroup <- Qt$QGroupBox("Access")

  accessLabel <- Qt$QLabel("Read-only:")
  accessComboBox <- Qt$QComboBox()
  accessComboBox$addItem("False")
  accessComboBox$addItem("True")

  this$accessLineEdit <- Qt$QLineEdit()

  qconnect(echoComboBox, "activated(int)", echoChanged)
  qconnect(validatorComboBox, "activated(int)", validatorChanged)
  qconnect(alignmentComboBox, "activated(int)", alignmentChanged)
  qconnect(inputMaskComboBox, "activated(int)", inputMaskChanged)
  qconnect(accessComboBox, "activated(int)", accessChanged)

  echoLayout <- Qt$QGridLayout()
  echoLayout$addWidget(echoLabel, 0, 0)
  echoLayout$addWidget(echoComboBox, 0, 1)
  echoLayout$addWidget(echoLineEdit, 1, 0, 1, 2)
  echoGroup$setLayout(echoLayout)

  validatorLayout <- Qt$QGridLayout()
  validatorLayout$addWidget(validatorLabel, 0, 0)
  validatorLayout$addWidget(validatorComboBox, 0, 1)
  validatorLayout$addWidget(validatorLineEdit, 1, 0, 1, 2)
  validatorGroup$setLayout(validatorLayout)

  alignmentLayout <- Qt$QGridLayout()
  alignmentLayout$addWidget(alignmentLabel, 0, 0)
  alignmentLayout$addWidget(alignmentComboBox, 0, 1)
  alignmentLayout$addWidget(alignmentLineEdit, 1, 0, 1, 2)
  alignmentGroup$setLayout(alignmentLayout)

  inputMaskLayout <- Qt$QGridLayout()
  inputMaskLayout$addWidget(inputMaskLabel, 0, 0)
  inputMaskLayout$addWidget(inputMaskComboBox, 0, 1)
  inputMaskLayout$addWidget(inputMaskLineEdit, 1, 0, 1, 2)
  inputMaskGroup$setLayout(inputMaskLayout)

  accessLayout <- Qt$QGridLayout()
  accessLayout$addWidget(accessLabel, 0, 0)
  accessLayout$addWidget(accessComboBox, 0, 1)
  accessLayout$addWidget(accessLineEdit, 1, 0, 1, 2)
  accessGroup$setLayout(accessLayout)

  layout <- Qt$QGridLayout()
  layout$addWidget(echoGroup, 0, 0)
  layout$addWidget(validatorGroup, 1, 0)
  layout$addWidget(alignmentGroup, 2, 0)
  layout$addWidget(inputMaskGroup, 0, 1)
  layout$addWidget(accessGroup, 1, 1)
  setLayout(layout)

  setWindowTitle("Line Edits")
})

qsetMethod("echoChanged", Window, function(index) {
  echoMode <- switch(index + 1L, Qt$QLineEdit$Normal, Qt$QLineEdit$Password,
                     Qt$QLineEdit$PasswordEchoOnEdit, Qt$QLineEdit$NoEcho)
  echoLineEdit$setEchoMode(echoMode)
}, "private")

qsetMethod("validatorChanged", Window, function(index) {
  validator <- switch(index + 1L, 0, Qt$QIntValidator(validatorLineEdit),
                      Qt$QDoubleValidator(-999.0, 999.0, 2, validatorLineEdit))
  validatorLineEdit$setValidator(validator)
  validatorLineEdit$clear()
}, "private")

qsetMethod("alignmentChanged", Window, function(index) {
  alignment <- switch(index + 1L, Qt$Qt$AlignLeft, Qt$Qt$AlignCenter,
                      Qt$Qt$AlignRight)
  alignmentLineEdit$setAlignment(alignment)
}, "private")

qsetMethod("inputMaskChanged", Window, function(index) {
  inputMask <- switch(index + 1L, "", "+99 99 99 99 99;_", "0000-00-00",
                      ">AAAAA-AAAAA-AAAAA-AAAAA-AAAAA;#")
  inputMaskLineEdit$setInputMask(inputMask)
  if (index == 2L) {
    inputMaskLineEdit$setText("00000000")
    inputMaskLineEdit$setCursorPosition(0)
  }
}, "private")

qsetMethod("accessChanged", Window, function(index) {
  accessLineEdit$setReadOnly(index)
}, "private")

Window()$show()
