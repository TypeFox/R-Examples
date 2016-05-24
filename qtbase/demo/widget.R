### From the Qt Widget tutorial

## construct a widget
window <- Qt$QWidget()

## resize and show
window$resize(320, 240)
window$show()

## now a widget that does something
button <- Qt$QPushButton("Press me", window)
button$move(100, 100)
button$show()

qconnect(button, "pressed", function() print("hello world"))

## layout
window <- Qt$QWidget()
label <- Qt$QLabel("Name:")
lineEdit <- Qt$QLineEdit()

layout <- Qt$QHBoxLayout()
layout$addWidget(label)
layout$addWidget(lineEdit)
window$setLayout(layout)
window$show()

## a bit more complex

window <- Qt$QWidget()
queryLabel <- Qt$QLabel("Query:")
queryEdit <- Qt$QLineEdit()
resultView <- Qt$QTableView()

queryLayout <- Qt$QHBoxLayout()
queryLayout$addWidget(queryLabel)
queryLayout$addWidget(queryEdit)

mainLayout <- Qt$QVBoxLayout()
mainLayout$addLayout(queryLayout)
mainLayout$addWidget(resultView)
window$setLayout(mainLayout)
window$show()

