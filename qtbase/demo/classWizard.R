## From the Qt Class Wizard example

imagefile <- function(x) system.file(file.path("images", x), package = "qtbase")

qsetClass("ClassWizard", Qt$QWizard, function(parent = NULL) {
  super(parent)

  addPage(IntroPage())
  addPage(ClassInfoPage())
  addPage(OutputFilesPage())
  addPage(ConclusionPage())

  setPixmap(Qt$QWizard$BannerPixmap, Qt$QPixmap(imagefile("banner.png")))
  setPixmap(Qt$QWizard$BackgroundPixmap,
            Qt$QPixmap(imagefile("background.png")))

  setWindowTitle("Class Wizard")
})

qsetMethod("accept", ClassWizard, function() {
  className <- field("className")
  baseClass <- field("baseClass")
  path <- field("path")

  block <- character()
  if (!nchar(baseClass))
    baseClass <- "Qt$QObject"
  block <- paste(block, "qsetClass(\"", className, "\", ", baseClass, sep = "")
  if (field("qobjectCtor"))
    block <- paste(block, ", function(parent = NULL) {\n\n}", sep = "")
  block <- paste(block, ")\n\n", sep = "")
  
  out <- try(cat(block, file=path))
  if (is(out, "try-error"))
    Qt$QMessageBox$warning(NULL, "Class Wizard",
                           sprintf("Cannot write file %s:\n%s", path, out))
  else super("accept")
})

qsetClass("IntroPage", Qt$QWizardPage, function(parent = NULL) {
  super(parent)

  setTitle("Introduction")
  setPixmap(Qt$QWizard$WatermarkPixmap, Qt$QPixmap(imagefile("watermark1.png")))

  label <- Qt$QLabel(paste("This wizard will generate a skeleton R/Qt class",
                           "definition, including a few functions. You simply",
                           "need to specify the class name and set a few",
                           "options to produce a source file",
                           "for your new R/Qt class."))
  label$wordWrap <- TRUE
  
  layout <- Qt$QVBoxLayout()
  layout$addWidget(label)
  setLayout(layout)
})

qsetClass("ClassInfoPage", Qt$QWizardPage, function(parent = NULL) {
  setTitle("Class Information")
  setSubTitle(paste("Specify basic information about the class for which you",
                    "want to generate skeleton source code files."))
  setPixmap(Qt$QWizard$LogoPixmap, Qt$QPixmap(imagefile("logo1.png")))

  classNameLabel <- Qt$QLabel("&Class name:")
  classNameLineEdit <- Qt$QLineEdit()
  classNameLabel$setBuddy(classNameLineEdit)

  baseClassLabel <- Qt$QLabel("&Base class:")
  baseClassLineEdit <- Qt$QLineEdit()
  baseClassLabel$setBuddy(baseClassLineEdit)
  
  groupBox <- Qt$QGroupBox("C&onstructor")
  
  qobjectCtorRadioButton <- Qt$QRadioButton("&QObject-style constructor")
  defaultCtorRadioButton <- Qt$QRadioButton("&Default constructor")
  
  defaultCtorRadioButton$checked <- TRUE

  registerField("className*", classNameLineEdit)
  registerField("baseClass", baseClassLineEdit)
  registerField("qobjectCtor", qobjectCtorRadioButton)
  registerField("defaultCtor", defaultCtorRadioButton)
  
  groupBoxLayout <- Qt$QVBoxLayout()
  groupBoxLayout$addWidget(qobjectCtorRadioButton)
  groupBoxLayout$addWidget(defaultCtorRadioButton)
  groupBox$setLayout(groupBoxLayout)

  layout <- Qt$QGridLayout()
  layout$addWidget(classNameLabel, 0, 0)
  layout$addWidget(classNameLineEdit, 0, 1)
  layout$addWidget(baseClassLabel, 1, 0)
  layout$addWidget(baseClassLineEdit, 1, 1)
  layout$addWidget(groupBox, 3, 0, 1, 2)
  setLayout(layout)
})

qsetClass("OutputFilesPage", Qt$QWizardPage, function(parent = NULL) {
  super(parent)

  setTitle("Output Files")
  setSubTitle(paste("Specify where you want the wizard to put the generated",
                    "skeleton code."))
  setPixmap(Qt$QWizard$LogoPixmap, Qt$QPixmap(imagefile("logo3.png")))

  pathLabel <- Qt$QLabel("&Output file:")
  this$pathLineEdit <- Qt$QLineEdit()
  pathLabel$setBuddy(pathLineEdit)
  
  registerField("path*", pathLineEdit)
  
  layout <- Qt$QGridLayout()
  layout$addWidget(pathLabel, 0, 0)
  layout$addWidget(pathLineEdit, 0, 1)
  setLayout(layout)
})

qsetMethod("initializePage", OutputFilesPage, function() {
  className <- field("className")
  pathLineEdit$text <- file.path(getwd(), paste(className, "R", sep = "."))
})

qsetClass("ConclusionPage", Qt$QWizardPage, function(parent = NULL) {
  super(parent)

  setTitle("Conclusion")
  setPixmap(Qt$QWizard$WatermarkPixmap, Qt$QPixmap(imagefile("watermark2.png")))

  this$label <- Qt$QLabel()
  label$wordWrap <- TRUE

  layout <- Qt$QVBoxLayout()
  layout$addWidget(label)
  setLayout(layout)
})

qsetMethod("initializePage", ConclusionPage, function() {
  finishText <- sub("&", "", wizard()$buttonText(Qt$QWizard$FinishButton))
  label$text <- sprintf("Click %s to generate the class skeleton.", finishText)
})

ClassWizard()$show()
