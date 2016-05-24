# # Last modified Feb 14, 2008
# # simulations optimized by Tyler Drombosky 2007
# 
# `betaSimulate.ipsur` <-
# function () 
# {
#     initializeDialog(title = gettextRcmdr("Simulate Beta Variates"))
#     parameterFrame <- tkframe(top)
#     locationFrame <- tkframe(top)
#     if (!is.character(ActiveDataSet())) {
#         locVariable <- tclVar("new")
#     }
#     else {
#         locVariable <- tclVar("add")
#     }
#     addtoactiveButton <- tkradiobutton(locationFrame, variable = locVariable, 
#         value = "add")
#     newDataButton <- tkradiobutton(locationFrame, variable = locVariable, 
#         value = "new")
#     samplesVar <- tclVar("1")
#     samplesEntry <- tkentry(top, width = "6", textvariable = samplesVar)
#     shape1Var <- tclVar("1")
#     shape1Entry <- tkentry(top, width = "6", textvariable = shape1Var)
#     shape2Var <- tclVar("1")
#     shape2Entry <- tkentry(top, width = "6", textvariable = shape2Var)
#     ncpVar <- tclVar("0")
#     ncpEntry <- tkentry(top, width = "6", textvariable = ncpVar)
#     onOK <- function() {
#         nsamples <- round(as.numeric(tclvalue(samplesVar)))
#         shape1 <- tclvalue(shape1Var)
#         shape2 <- tclvalue(shape2Var)
#         ncp <- tclvalue(ncpVar)
#         if (is.na(nsamples)) {
#             errorCondition(recall = betaSimulate.ipsur, message = gettextRcmdr("Number of samples must be a positive integer."))
#             return()
#         }
#         if (is.na(shape1)) {
#             errorCondition(recall = betaSimulate.ipsur, message = gettextRcmdr("The shape1 parameter was not specified."))
#             return()
#         }
#         if (is.na(shape2)) {
#             errorCondition(recall = betaSimulate.ipsur, message = gettextRcmdr("The shape2 parameter was not specified."))
#             return()
#         }
#         if (is.na(ncp)) {
#             errorCondition(recall = betaSimulate.ipsur, message = gettextRcmdr("The noncentrality parameter was not specified."))
#             return()
#         }
#         closeDialog()
#         store <- tclvalue(locVariable)
#         if (store == "new") {
#             initializeDialog(title = gettextRcmdr("Simulation Dataset"))
#             dsname <- tclVar("Simset")
#             entryDsname <- tkentry(top, width = "20", textvariable = dsname)
#             newDataSS <- tclVar("100")
#             entryNewDataSS <- tkentry(top, width = "6", textvariable = newDataSS)
#             onOK <- function() {
#                 dsnameValue <- trim.blanks(tclvalue(dsname))
#                 newSS <- round(as.numeric(tclvalue(newDataSS)))
#                 closeDialog()
#                 if (dsnameValue == "") {
#                   errorCondition(recall = betaSimulate.ipsur, 
#                     message = gettextRcmdr("You must enter the name of a data set."))
#                   return()
#                 }
#                 if (!is.valid.name(dsnameValue)) {
#                   errorCondition(recall = betaSimulate.ipsur, 
#                     message = paste("\"", dsnameValue, "\" ", 
#                       gettextRcmdr("is not a valid name."), sep = ""))
#                   return()
#                 }
#                 if (is.element(dsnameValue, listDataSets())) {
#                   if ("no" == tclvalue(checkReplace(dsnameValue, 
#                     gettextRcmdr("Data set")))) {
#                     betaSimulate.ipsur()
#                     return()
#                   }
#                 }
#                 if (is.na(newSS)) {
#                   errorCondition(recall = betaSimulate.ipsur, 
#                     message = gettextRcmdr("Sample Size must be a positive integer."))
#                   return()
#                 }
#                 UpdatebetasimNumber()
#                 justDoIt(paste(dsnameValue, " = data.frame(beta.sim", 
#                   getRcmdr("betasimNumber"), "=1:", newSS, ")", 
#                   sep = ""))
#                 logger(paste(dsnameValue, "has been initialized."))
#                 for (k in getRcmdr("betasimNumber"):(nsamples + 
#                   getRcmdr("betasimNumber") - 1)) {
#                   justDoIt(paste(dsnameValue, "$beta.sim", k, 
#                     " <- rbeta(", newSS, ", shape1=", shape1, 
#                     ", shape2=", shape2, ", ncp=", ncp, ")", 
#                     sep = ""))
#                 }
#                 activeDataSet(dsnameValue)
#                 putRcmdr("betasimNumber", k)
#                 if (nsamples == 1) {
#                   logger(paste("There was 1 beta variate sample stored in ", 
#                     dsnameValue, ".", sep = ""))
#                 }
#                 else {
#                   logger(paste("There were ", nsamples, " beta variate samples stored in ", 
#                     dsnameValue, ".", sep = ""))
#                 }
#             }
#             OKCancelHelp(helpSubject = "rbeta")
#             tkgrid(tklabel(top, text = gettextRcmdr("Enter name for data set:")), 
#                 entryDsname, sticky = "e")
#             tkgrid(tklabel(top, text = gettextRcmdr("Sample Size (rows):")), 
#                 entryNewDataSS, sticky = "e")
#             tkgrid(buttonsFrame, columnspan = "2", sticky = "w")
#             tkgrid.configure(entryDsname, sticky = "w")
#             tkgrid.configure(entryNewDataSS, sticky = "w")
#             tkfocus(CommanderWindow())
#             dialogSuffix(rows = 2, columns = 2, focus = entryDsname)
#         }
#         else {
#             if (!is.character(ActiveDataSet())) {
#                 errorCondition(recall = betaSimulate.ipsur, message = gettextRcmdr("There is no active data set."))
#                 return()
#             }
#             .activeDataSet <- ActiveDataSet()
#             justDoIt(paste("samplesn <- dim(", .activeDataSet, 
#                 ")[1]", sep = ""))
#             UpdatebetasimNumber()
#             for (k in getRcmdr("betasimNumber"):(nsamples + getRcmdr("betasimNumber") - 
#                 1)) {
#                 justDoIt(paste(.activeDataSet, "$beta.sim", k, 
#                   " <- rbeta(", samplesn, ", shape1=", shape1, 
#                   ", shape2=", shape2, ", ncp=", ncp, ")", sep = ""))
#             }
#             activeDataSet(.activeDataSet)
#             putRcmdr("betasimNumber", k)
#             if (nsamples == 1) {
#                 logger(paste("There was 1 beta variate sample stored in ", 
#                   .activeDataSet, ".", sep = ""))
#             }
#             else {
#                 logger(paste("There were ", nsamples, " beta variate samples stored in ", 
#                   .activeDataSet, ".", sep = ""))
#             }
#         }
#         tkfocus(CommanderWindow())
#     }
#     OKCancelHelp(helpSubject = "rbeta")
#     tkgrid(tklabel(top, text = gettextRcmdr("Number of samples (columns):")), 
#         samplesEntry, sticky = "w")
#     tkgrid(tklabel(top, text = gettextRcmdr("Parameters:"), fg = "blue"), 
#         columnspan = 4, sticky = "w")
#     tkgrid(tklabel(top, text = gettextRcmdr("shape1")), shape1Entry, 
#         sticky = "w")
#     tkgrid(tklabel(top, text = gettextRcmdr("shape2")), shape2Entry, 
#         sticky = "w")
#     tkgrid(tklabel(top, text = gettextRcmdr("ncp (noncentrality parameter)")), 
#         ncpEntry, sticky = "w")
#     tkgrid(tklabel(locationFrame, text = gettextRcmdr("Store values in:"), 
#         fg = "blue"), columnspan = 4, sticky = "w")
#     tkgrid(tklabel(locationFrame, text = gettextRcmdr("Active Dataset")), 
#         addtoactiveButton, sticky = "w")
#     tkgrid(tklabel(locationFrame, text = "New Dataset"), newDataButton, 
#         sticky = "w")
#     tkgrid.configure(samplesEntry, sticky = "w")
#     tkgrid.configure(shape1Entry, sticky = "w")
#     tkgrid.configure(shape2Entry, sticky = "w")
#     tkgrid.configure(ncpEntry, sticky = "w")
#     tkgrid(locationFrame, sticky = "w")
#     tkgrid(buttonsFrame, sticky = "w", columnspan = 2)
#     dialogSuffix(rows = 6, columns = 1, focus = samplesEntry)
# }
# 
# 
# `binomialSimulate.ipsur` <-
# function () 
# {
#     initializeDialog(title = gettextRcmdr("Simulate Binomial Variates"))
#     parameterFrame <- tkframe(top)
#     locationFrame <- tkframe(top)
#     if (!is.character(ActiveDataSet())) {
#         locVariable <- tclVar("new")
#     }
#     else {
#         locVariable <- tclVar("add")
#     }
#     addtoactiveButton <- tkradiobutton(locationFrame, variable = locVariable, 
#         value = "add")
#     newDataButton <- tkradiobutton(locationFrame, variable = locVariable, 
#         value = "new")
#     samplesVar <- tclVar("1")
#     samplesEntry <- tkentry(top, width = "6", textvariable = samplesVar)
#     sizeVar <- tclVar("1")
#     sizeEntry <- tkentry(top, width = "6", textvariable = sizeVar)
#     probVar <- tclVar("0.5")
#     probEntry <- tkentry(top, width = "6", textvariable = probVar)
#     onOK <- function() {
#         nsamples <- round(as.numeric(tclvalue(samplesVar)))
#         size <- tclvalue(sizeVar)
#         prob <- tclvalue(probVar)
#         store <- tclvalue(locVariable)
#         if (is.na(nsamples)) {
#             errorCondition(recall = binomialSimulate.ipsur, message = gettextRcmdr("Number of samples must be a positive integer."))
#             return()
#         }
#         if (is.na(size)) {
#             errorCondition(recall = binomialSimulate.ipsur, message = gettextRcmdr("Number of trials was not specified."))
#             return()
#         }
#         if (is.na(prob)) {
#             errorCondition(recall = binomialSimulate.ipsur, message = gettextRcmdr("The success probability was not specified."))
#             return()
#         }
#         closeDialog()
#         if (store == "new") {
#             initializeDialog(title = gettextRcmdr("Simulation Dataset"))
#             dsname <- tclVar("Simset")
#             entryDsname <- tkentry(top, width = "20", textvariable = dsname)
#             newDataSS <- tclVar("100")
#             entryNewDataSS <- tkentry(top, width = "6", textvariable = newDataSS)
#             onOK <- function() {
#                 dsnameValue <- trim.blanks(tclvalue(dsname))
#                 newSS <- round(as.numeric(tclvalue(newDataSS)))
#                 closeDialog()
#                 if (dsnameValue == "") {
#                   errorCondition(recall = binomialSimulate.ipsur, 
#                     message = gettextRcmdr("You must enter the name of a data set."))
#                   return()
#                 }
#                 if (!is.valid.name(dsnameValue)) {
#                   errorCondition(recall = binomialSimulate.ipsur, 
#                     message = paste("\"", dsnameValue, "\" ", 
#                       gettextRcmdr("is not a valid name."), sep = ""))
#                   return()
#                 }
#                 if (is.element(dsnameValue, listDataSets())) {
#                   if ("no" == tclvalue(checkReplace(dsnameValue, 
#                     gettextRcmdr("Data set")))) {
#                     binomialSimulate.ipsur()
#                     return()
#                   }
#                 }
#                 if (is.na(newSS)) {
#                   errorCondition(recall = binomialSimulate.ipsur, 
#                     message = gettextRcmdr("Sample Size must be a positive integer."))
#                   return()
#                 }
#                 UpdatebinomsimNumber()
#                 justDoIt(paste(dsnameValue, " = data.frame(binom.sim", 
#                   getRcmdr("binomsimNumber"), "=1:", newSS, ")", 
#                   sep = ""))
#                 logger(paste(dsnameValue, "has been initialized."))
#                 for (k in getRcmdr("binomsimNumber"):(nsamples + 
#                   getRcmdr("binomsimNumber") - 1)) {
#                   justDoIt(paste(dsnameValue, "$binom.sim", k, 
#                     " <- rbinom(", newSS, ", size=", size, ", prob=", 
#                     prob, ")", sep = ""))
#                 }
#                 activeDataSet(dsnameValue)
#                 putRcmdr("binomsimNumber", k)
#                 if (nsamples == 1) {
#                   logger(paste("There was 1 binomial variate sample stored in ", 
#                     dsnameValue, ".", sep = ""))
#                 }
#                 else {
#                   logger(paste("There were ", nsamples, " binomial variate samples stored in ", 
#                     dsnameValue, ".", sep = ""))
#                 }
#             }
#             OKCancelHelp(helpSubject = "rbinom")
#             tkgrid(tklabel(top, text = gettextRcmdr("Enter name for data set:")), 
#                 entryDsname, sticky = "e")
#             tkgrid(tklabel(top, text = gettextRcmdr("Sample Size (rows):")), 
#                 entryNewDataSS, sticky = "e")
#             tkgrid(buttonsFrame, columnspan = "2", sticky = "w")
#             tkgrid.configure(entryDsname, sticky = "w")
#             tkgrid.configure(entryNewDataSS, sticky = "w")
#             tkfocus(CommanderWindow())
#             dialogSuffix(rows = 2, columns = 2, focus = entryDsname)
#         }
#         else {
#             if (!is.character(ActiveDataSet())) {
#                 errorCondition(recall = binomialSimulate.ipsur, 
#                   message = gettextRcmdr("There is no active data set."))
#                 return()
#             }
#             .activeDataSet <- ActiveDataSet()
#             justDoIt(paste("samplesn <- dim(", .activeDataSet, 
#                 ")[1]", sep = ""))
#             UpdatebinomsimNumber()
#             for (k in getRcmdr("binomsimNumber"):(nsamples + 
#                 getRcmdr("binomsimNumber") - 1)) {
#                 justDoIt(paste(.activeDataSet, "$binom.sim", 
#                   k, " <- rbinom(", samplesn, ", size=", size, 
#                   ", prob=", prob, ")", sep = ""))
#             }
#             activeDataSet(.activeDataSet)
#             putRcmdr("binomsimNumber", k)
#             if (nsamples == 1) {
#                 logger(paste("There was 1 binomial variate sample stored in ", 
#                   .activeDataSet, ".", sep = ""))
#             }
#             else {
#                 logger(paste("There were ", nsamples, " binomial variate samples stored in ", 
#                   .activeDataSet, ".", sep = ""))
#             }
#         }
#         tkfocus(CommanderWindow())
#     }
#     OKCancelHelp(helpSubject = "rbinom")
#     tkgrid(tklabel(top, text = gettextRcmdr("Number of samples (columns):")), 
#         samplesEntry, sticky = "w")
#     tkgrid(tklabel(top, text = gettextRcmdr("Parameters:"), fg = "blue"), 
#         columnspan = 4, sticky = "w")
#     tkgrid(tklabel(top, text = gettextRcmdr("size (number of trials)")), 
#         sizeEntry, sticky = "w")
#     tkgrid(tklabel(top, text = gettextRcmdr("prob (of success)")), 
#         probEntry, sticky = "w")
#     tkgrid(tklabel(locationFrame, text = gettextRcmdr("Store values in:"), 
#         fg = "blue"), columnspan = 4, sticky = "w")
#     tkgrid(tklabel(locationFrame, text = gettextRcmdr("Active Dataset")), 
#         addtoactiveButton, sticky = "w")
#     tkgrid(tklabel(locationFrame, text = "New Dataset"), newDataButton, 
#         sticky = "w")
#     tkgrid.configure(samplesEntry, sticky = "w")
#     tkgrid.configure(sizeEntry, sticky = "w")
#     tkgrid.configure(probEntry, sticky = "w")
#     tkgrid(locationFrame, sticky = "w")
#     tkgrid(buttonsFrame, sticky = "w", columnspan = 2)
#     dialogSuffix(rows = 6, columns = 1, focus = samplesEntry)
# }
# 
# 
# `cauchySimulate.ipsur` <-
# function () 
# {
#     initializeDialog(title = gettextRcmdr("Simulate Cauchy Variates"))
#     parameterFrame <- tkframe(top)
#     locationFrame <- tkframe(top)
#     if (!is.character(ActiveDataSet())) {
#         locVariable <- tclVar("new")
#     }
#     else {
#         locVariable <- tclVar("add")
#     }
#     addtoactiveButton <- tkradiobutton(locationFrame, variable = locVariable, 
#         value = "add")
#     newDataButton <- tkradiobutton(locationFrame, variable = locVariable, 
#         value = "new")
#     samplesVar <- tclVar("1")
#     samplesEntry <- tkentry(top, width = "6", textvariable = samplesVar)
#     locationVar <- tclVar("0")
#     locationEntry <- tkentry(top, width = "6", textvariable = locationVar)
#     scale1Var <- tclVar("1")
#     scale1Entry <- tkentry(top, width = "6", textvariable = scale1Var)
#     onOK <- function() {
#         nsamples <- round(as.numeric(tclvalue(samplesVar)))
#         location <- tclvalue(locationVar)
#         scale1 <- tclvalue(scale1Var)
#         if (is.na(nsamples)) {
#             errorCondition(recall = cauchySimulate.ipsur, message = gettextRcmdr("Number of samples must be a positive integer."))
#             return()
#         }
#         if (is.na(location)) {
#             errorCondition(recall = cauchySimulate.ipsur, message = gettextRcmdr("The location parameter was not specified."))
#             return()
#         }
#         if (is.na(scale1)) {
#             errorCondition(recall = cauchySimulate.ipsur, message = gettextRcmdr("The scale parameter was not specified."))
#             return()
#         }
#         closeDialog()
#         store <- tclvalue(locVariable)
#         if (store == "new") {
#             initializeDialog(title = gettextRcmdr("Simulation Dataset"))
#             dsname <- tclVar("Simset")
#             entryDsname <- tkentry(top, width = "20", textvariable = dsname)
#             newDataSS <- tclVar("100")
#             entryNewDataSS <- tkentry(top, width = "6", textvariable = newDataSS)
#             onOK <- function() {
#                 dsnameValue <- trim.blanks(tclvalue(dsname))
#                 newSS <- round(as.numeric(tclvalue(newDataSS)))
#                 closeDialog()
#                 if (dsnameValue == "") {
#                   errorCondition(recall = cauchySimulate.ipsur, 
#                     message = gettextRcmdr("You must enter the name of a data set."))
#                   return()
#                 }
#                 if (!is.valid.name(dsnameValue)) {
#                   errorCondition(recall = cauchySimulate.ipsur, 
#                     message = paste("\"", dsnameValue, "\" ", 
#                       gettextRcmdr("is not a valid name."), sep = ""))
#                   return()
#                 }
#                 if (is.element(dsnameValue, listDataSets())) {
#                   if ("no" == tclvalue(checkReplace(dsnameValue, 
#                     gettextRcmdr("Data set")))) {
#                     cauchySimulate.ipsur()
#                     return()
#                   }
#                 }
#                 if (is.na(newSS)) {
#                   errorCondition(recall = cauchySimulate.ipsur, 
#                     message = gettextRcmdr("Sample Size must be a positive integer."))
#                   return()
#                 }
#                 UpdatecauchysimNumber()
#                 justDoIt(paste(dsnameValue, " = data.frame(cauchy.sim", 
#                   getRcmdr("cauchysimNumber"), "=1:", newSS, 
#                   ")", sep = ""))
#                 logger(paste(dsnameValue, "has been initialized."))
#                 for (k in getRcmdr("cauchysimNumber"):(nsamples + 
#                   getRcmdr("cauchysimNumber") - 1)) {
#                   justDoIt(paste(dsnameValue, "$cauchy.sim", 
#                     k, " <- rcauchy(", newSS, ", location=", 
#                     location, ", scale=", scale1, ")", sep = ""))
#                 }
#                 activeDataSet(dsnameValue)
#                 putRcmdr("cauchysimNumber", k)
#                 if (nsamples == 1) {
#                   logger(paste("There was 1 Cauchy variate sample stored in ", 
#                     dsnameValue, ".", sep = ""))
#                 }
#                 else {
#                   logger(paste("There were ", nsamples, " Cauchy variate samples stored in ", 
#                     dsnameValue, ".", sep = ""))
#                 }
#             }
#             OKCancelHelp(helpSubject = "rcauchy")
#             tkgrid(tklabel(top, text = gettextRcmdr("Enter name for data set:")), 
#                 entryDsname, sticky = "e")
#             tkgrid(tklabel(top, text = gettextRcmdr("Sample Size (rows):")), 
#                 entryNewDataSS, sticky = "e")
#             tkgrid(buttonsFrame, columnspan = "2", sticky = "w")
#             tkgrid.configure(entryDsname, sticky = "w")
#             tkgrid.configure(entryNewDataSS, sticky = "w")
#             tkfocus(CommanderWindow())
#             dialogSuffix(rows = 2, columns = 2, focus = entryDsname)
#         }
#         else {
#             if (!is.character(ActiveDataSet())) {
#                 errorCondition(recall = cauchySimulate.ipsur, 
#                   message = gettextRcmdr("There is no active data set."))
#                 return()
#             }
#             .activeDataSet <- ActiveDataSet()
#             justDoIt(paste("samplesn <- dim(", .activeDataSet, 
#                 ")[1]", sep = ""))
#             UpdatecauchysimNumber()
#             for (k in getRcmdr("cauchysimNumber"):(nsamples + 
#                 getRcmdr("cauchysimNumber") - 1)) {
#                 justDoIt(paste(.activeDataSet, "$cauchy.sim", 
#                   k, " <- rcauchy(", samplesn, ", location=", 
#                   location, ", scale=", scale1, ")", sep = ""))
#             }
#             activeDataSet(.activeDataSet)
#             putRcmdr("cauchysimNumber", k)
#             if (nsamples == 1) {
#                 logger(paste("There was 1 Cauchy variate sample stored in ", 
#                   .activeDataSet, ".", sep = ""))
#             }
#             else {
#                 logger(paste("There were ", nsamples, " Cauchy variate samples stored in ", 
#                   .activeDataSet, ".", sep = ""))
#             }
#         }
#         tkfocus(CommanderWindow())
#     }
#     OKCancelHelp(helpSubject = "rcauchy")
#     tkgrid(tklabel(top, text = gettextRcmdr("Number of samples (columns):")), 
#         samplesEntry, sticky = "w")
#     tkgrid(tklabel(top, text = gettextRcmdr("Parameters:"), fg = "blue"), 
#         columnspan = 4, sticky = "w")
#     tkgrid(tklabel(top, text = gettextRcmdr("location")), locationEntry, 
#         sticky = "w")
#     tkgrid(tklabel(top, text = gettextRcmdr("scale")), scale1Entry, 
#         sticky = "w")
#     tkgrid(tklabel(locationFrame, text = gettextRcmdr("Store values in:"), 
#         fg = "blue"), columnspan = 4, sticky = "w")
#     tkgrid(tklabel(locationFrame, text = gettextRcmdr("Active Dataset")), 
#         addtoactiveButton, sticky = "w")
#     tkgrid(tklabel(locationFrame, text = "New Dataset"), newDataButton, 
#         sticky = "w")
#     tkgrid.configure(samplesEntry, sticky = "w")
#     tkgrid.configure(locationEntry, sticky = "w")
#     tkgrid.configure(scale1Entry, sticky = "w")
#     tkgrid(locationFrame, sticky = "w")
#     tkgrid(buttonsFrame, sticky = "w", columnspan = 2)
#     dialogSuffix(rows = 6, columns = 1, focus = samplesEntry)
# }
# 
# 
# `chisqSimulate.ipsur` <-
# function () 
# {
#     initializeDialog(title = gettextRcmdr("Simulate Chi-Squared Variates"))
#     parameterFrame <- tkframe(top)
#     locationFrame <- tkframe(top)
#     if (!is.character(ActiveDataSet())) {
#         locVariable <- tclVar("new")
#     }
#     else {
#         locVariable <- tclVar("add")
#     }
#     addtoactiveButton <- tkradiobutton(locationFrame, variable = locVariable, 
#         value = "add")
#     newDataButton <- tkradiobutton(locationFrame, variable = locVariable, 
#         value = "new")
#     samplesVar <- tclVar("1")
#     samplesEntry <- tkentry(top, width = "6", textvariable = samplesVar)
#     dfVar <- tclVar("1")
#     dfEntry <- tkentry(top, width = "6", textvariable = dfVar)
#     ncpVar <- tclVar("0")
#     ncpEntry <- tkentry(top, width = "6", textvariable = ncpVar)
#     onOK <- function() {
#         nsamples <- round(as.numeric(tclvalue(samplesVar)))
#         df <- tclvalue(dfVar)
#         ncp <- tclvalue(ncpVar)
#         if (is.na(nsamples)) {
#             errorCondition(recall = chisqSimulate.ipsur, message = gettextRcmdr("Number of samples must be a positive integer."))
#             return()
#         }
#         if (is.na(df)) {
#             errorCondition(recall = chisqSimulate.ipsur, message = gettextRcmdr("The degrees of freedom were not specified."))
#             return()
#         }
#         if (is.na(ncp)) {
#             errorCondition(recall = chisqSimulate.ipsur, message = gettextRcmdr("The noncentrality parameter was not specified."))
#             return()
#         }
#         closeDialog()
#         store <- tclvalue(locVariable)
#         if (store == "new") {
#             initializeDialog(title = gettextRcmdr("Simulation Dataset"))
#             dsname <- tclVar("Simset")
#             entryDsname <- tkentry(top, width = "20", textvariable = dsname)
#             newDataSS <- tclVar("100")
#             entryNewDataSS <- tkentry(top, width = "6", textvariable = newDataSS)
#             onOK <- function() {
#                 dsnameValue <- trim.blanks(tclvalue(dsname))
#                 newSS <- round(as.numeric(tclvalue(newDataSS)))
#                 closeDialog()
#                 if (dsnameValue == "") {
#                   errorCondition(recall = chisqSimulate.ipsur, 
#                     message = gettextRcmdr("You must enter the name of a data set."))
#                   return()
#                 }
#                 if (!is.valid.name(dsnameValue)) {
#                   errorCondition(recall = chisqSimulate.ipsur, 
#                     message = paste("\"", dsnameValue, "\" ", 
#                       gettextRcmdr("is not a valid name."), sep = ""))
#                   return()
#                 }
#                 if (is.element(dsnameValue, listDataSets())) {
#                   if ("no" == tclvalue(checkReplace(dsnameValue, 
#                     gettextRcmdr("Data set")))) {
#                     chisqSimulate.ipsur()
#                     return()
#                   }
#                 }
#                 if (is.na(newSS)) {
#                   errorCondition(recall = chisqSimulate.ipsur, 
#                     message = gettextRcmdr("Sample Size must be a positive integer."))
#                   return()
#                 }
#                 UpdatechisqsimNumber()
#                 justDoIt(paste(dsnameValue, " = data.frame(chisq.sim", 
#                   getRcmdr("chisqsimNumber"), "=1:", newSS, ")", 
#                   sep = ""))
#                 logger(paste(dsnameValue, "has been initialized."))
#                 for (k in getRcmdr("chisqsimNumber"):(nsamples + 
#                   getRcmdr("chisqsimNumber") - 1)) {
#                   justDoIt(paste(dsnameValue, "$chisq.sim", k, 
#                     " <- rchisq(", newSS, ", df=", df, ", ncp=", 
#                     ncp, ")", sep = ""))
#                 }
#                 activeDataSet(dsnameValue)
#                 putRcmdr("chisqsimNumber", k)
#                 if (nsamples == 1) {
#                   logger(paste("There was 1 chi-squared variate sample stored in ", 
#                     dsnameValue, ".", sep = ""))
#                 }
#                 else {
#                   logger(paste("There were ", nsamples, " chi-squared variate samples stored in ", 
#                     dsnameValue, ".", sep = ""))
#                 }
#             }
#             OKCancelHelp(helpSubject = "rchisq")
#             tkgrid(tklabel(top, text = gettextRcmdr("Enter name for data set:")), 
#                 entryDsname, sticky = "e")
#             tkgrid(tklabel(top, text = gettextRcmdr("Sample Size (rows):")), 
#                 entryNewDataSS, sticky = "e")
#             tkgrid(buttonsFrame, columnspan = "2", sticky = "w")
#             tkgrid.configure(entryDsname, sticky = "w")
#             tkgrid.configure(entryNewDataSS, sticky = "w")
#             tkfocus(CommanderWindow())
#             dialogSuffix(rows = 2, columns = 2, focus = entryDsname)
#         }
#         else {
#             if (!is.character(ActiveDataSet())) {
#                 errorCondition(recall = chisqSimulate.ipsur, 
#                   message = gettextRcmdr("There is no active data set."))
#                 return()
#             }
#             .activeDataSet <- ActiveDataSet()
#             justDoIt(paste("samplesn <- dim(", .activeDataSet, 
#                 ")[1]", sep = ""))
#             UpdatechisqsimNumber()
#             for (k in getRcmdr("chisqsimNumber"):(nsamples + 
#                 getRcmdr("chisqsimNumber") - 1)) {
#                 justDoIt(paste(.activeDataSet, "$chisq.sim", 
#                   k, " <- rchisq(", samplesn, ", df=", df, ", ncp=", 
#                   ncp, ")", sep = ""))
#             }
#             activeDataSet(.activeDataSet)
#             putRcmdr("chisqsimNumber", k)
#             if (nsamples == 1) {
#                 logger(paste("There was 1 chi-squared variate sample stored in ", 
#                   .activeDataSet, ".", sep = ""))
#             }
#             else {
#                 logger(paste("There were ", nsamples, "    chi-squared variate samples stored in ", 
#                   .activeDataSet, ".", sep = ""))
#             }
#         }
#         tkfocus(CommanderWindow())
#     }
#     OKCancelHelp(helpSubject = "rchisq")
#     tkgrid(tklabel(top, text = gettextRcmdr("Number of samples (columns):")), 
#         samplesEntry, sticky = "w")
#     tkgrid(tklabel(top, text = gettextRcmdr("Parameters:"), fg = "blue"), 
#         columnspan = 4, sticky = "w")
#     tkgrid(tklabel(top, text = gettextRcmdr("df (degrees of freedom)")), 
#         dfEntry, sticky = "w")
#     tkgrid(tklabel(top, text = gettextRcmdr("ncp (noncentrality parameter)")), 
#         ncpEntry, sticky = "w")
#     tkgrid(tklabel(locationFrame, text = gettextRcmdr("Store values in:"), 
#         fg = "blue"), columnspan = 4, sticky = "w")
#     tkgrid(tklabel(locationFrame, text = gettextRcmdr("Active Dataset")), 
#         addtoactiveButton, sticky = "w")
#     tkgrid(tklabel(locationFrame, text = "New Dataset"), newDataButton, 
#         sticky = "w")
#     tkgrid.configure(samplesEntry, sticky = "w")
#     tkgrid.configure(dfEntry, sticky = "w")
#     tkgrid.configure(ncpEntry, sticky = "w")
#     tkgrid(locationFrame, sticky = "w")
#     tkgrid(buttonsFrame, sticky = "w", columnspan = 2)
#     dialogSuffix(rows = 6, columns = 1, focus = samplesEntry)
# }
# 
# 
# `disunifSimulate.ipsur` <-
# function () 
# {
#     initializeDialog(title = gettextRcmdr("Simulate Discrete Uniform Variates"))
#     parameterFrame <- tkframe(top)
#     locationFrame <- tkframe(top)
#     if (!is.character(ActiveDataSet())) {
#         locVariable <- tclVar("new")
#     }
#     else {
#         locVariable <- tclVar("add")
#     }
#     addtoactiveButton <- tkradiobutton(locationFrame, variable = locVariable, 
#         value = "add")
#     newDataButton <- tkradiobutton(locationFrame, variable = locVariable, 
#         value = "new")
#     samplesVar <- tclVar("1")
#     samplesEntry <- tkentry(top, width = "6", textvariable = samplesVar)
#     from1Var <- tclVar("1")
#     from1Entry <- tkentry(top, width = "6", textvariable = from1Var)
#     to1Var <- tclVar("10")
#     to1Entry <- tkentry(top, width = "6", textvariable = to1Var)
#     by1Var <- tclVar("1")
#     by1Entry <- tkentry(top, width = "6", textvariable = by1Var)
#     userdefEntry <- tkentry(top, width = "30", textvariable = "")
#     onOK <- function() {
#         nsamples <- round(as.numeric(tclvalue(samplesVar)))
#         from1 <- tclvalue(from1Var)
#         to1 <- tclvalue(to1Var)
#         by1 <- tclvalue(by1Var)
#         if (is.na(nsamples)) {
#             errorCondition(recall = disunifSimulate.ipsur, message = gettextRcmdr("Number of samples must be a positive integer."))
#             return()
#         }
#         if (is.na(from1)) {
#             errorCondition(recall = disunifSimulate.ipsur, message = gettextRcmdr("The from parameter was not specified."))
#             return()
#         }
#         if (is.na(to1)) {
#             errorCondition(recall = disunifSimulate.ipsur, message = gettextRcmdr("The to parameter was not specified."))
#             return()
#         }
#         if (is.na(by1)) {
#             errorCondition(recall = disunifSimulate.ipsur, message = gettextRcmdr("The by parameter was not specified."))
#             return()
#         }
#         closeDialog()
#         command <- paste("support <- seq(", from1, ", ", to1, ", by=", by1, 
#             ")", sep = "")
#         justDoIt(command)
#         store <- tclvalue(locVariable)
#         if (store == "new") {
#             initializeDialog(title = gettextRcmdr("Simulation Dataset"))
#             dsname <- tclVar("Simset")
#             entryDsname <- tkentry(top, width = "20", textvariable = dsname)
#             newDataSS <- tclVar("100")
#             entryNewDataSS <- tkentry(top, width = "6", textvariable = newDataSS)
#             onOK <- function() {
#                 dsnameValue <- trim.blanks(tclvalue(dsname))
#                 newSS <- round(as.numeric(tclvalue(newDataSS)))
#                 closeDialog()
#                 if (dsnameValue == "") {
#                   errorCondition(recall = disunifSimulate.ipsur, 
#                     message = gettextRcmdr("You must enter the name of a data set."))
#                   return()
#                 }
#                 if (!is.valid.name(dsnameValue)) {
#                   errorCondition(recall = disunifSimulate.ipsur, 
#                     message = paste("\"", dsnameValue, "\" ", 
#                       gettextRcmdr("is not a valid name."), sep = ""))
#                   return()
#                 }
#                 if (is.element(dsnameValue, listDataSets())) {
#                   if ("no" == tclvalue(checkReplace(dsnameValue, 
#                     gettextRcmdr("Data set")))) {
#                     disunifSimulate.ipsur()
#                     return()
#                   }
#                 }
#                 if (is.na(newSS)) {
#                   errorCondition(recall = disunifSimulate.ipsur, 
#                     message = gettextRcmdr("Sample Size must be a positive integer."))
#                   return()
#                 }
#                 UpdatedisunifsimNumber()
#                 justDoIt(paste(dsnameValue, " = data.frame(disunif.sim", 
#                   getRcmdr("disunifsimNumber"), "=1:", newSS, 
#                   ")", sep = ""))
#                 logger(paste(dsnameValue, "has been initialized."))
#                 for (k in getRcmdr("disunifsimNumber"):(nsamples + 
#                   getRcmdr("disunifsimNumber") - 1)) {
#                   justDoIt(paste(dsnameValue, "$disunif.sim", 
#                     k, " <- sample(support, size=", newSS, ", replace = TRUE)", 
#                     sep = ""))
#                 }
#                 activeDataSet(dsnameValue)
#                 putRcmdr("disunifsimNumber", k)
#                 if (nsamples == 1) {
#                   logger(paste("There was 1 discrete uniform variate sample stored in ", 
#                     dsnameValue, ".", sep = ""))
#                 }
#                 else {
#                   logger(paste("There were ", nsamples, " discrete uniform variate samples stored in ", 
#                     dsnameValue, ".", sep = ""))
#                 }
#             }
#             OKCancelHelp(helpSubject = "rdisunif")
#             tkgrid(tklabel(top, text = gettextRcmdr("Enter name for data set:")), 
#                 entryDsname, sticky = "e")
#             tkgrid(tklabel(top, text = gettextRcmdr("Sample Size (rows):")), 
#                 entryNewDataSS, sticky = "e")
#             tkgrid(buttonsFrame, columnspan = "2", sticky = "w")
#             tkgrid.configure(entryDsname, sticky = "w")
#             tkgrid.configure(entryNewDataSS, sticky = "w")
#             tkfocus(CommanderWindow())
#             dialogSuffix(rows = 2, columns = 2, focus = entryDsname)
#         }
#         else {
#             if (!is.character(ActiveDataSet())) {
#                 errorCondition(recall = disunifSimulate.ipsur, 
#                   message = gettextRcmdr("There is no active data set."))
#                 return()
#             }
#             .activeDataSet <- ActiveDataSet()
#             justDoIt(paste("samplesn <- dim(", .activeDataSet, 
#                 ")[1]", sep = ""))
#             UpdatedisunifsimNumber()
#             for (k in getRcmdr("disunifsimNumber"):(nsamples + 
#                 getRcmdr("disunifsimNumber") - 1)) {
#                 justDoIt(paste(.activeDataSet, "$disunif.sim", 
#                   k, " <- sample(support, size=", samplesn, ", replace = TRUE)", 
#                   sep = ""))
#             }
#             activeDataSet(.activeDataSet)
#             putRcmdr("disunifsimNumber", k)
#             if (nsamples == 1) {
#                 logger(paste("There was 1 discrete uniform variate sample stored in ", 
#                   .activeDataSet, ".", sep = ""))
#             }
#             else {
#                 logger(paste("There were ", nsamples, " discrete uniform variate samples stored in ", 
#                   .activeDataSet, ".", sep = ""))
#             }
#         }
#         remove(support, envir = .GlobalEnv)
#         tkfocus(CommanderWindow())
#     }
#     OKCancelHelp(helpSubject = "rdisunif")
#     tkgrid(tklabel(top, text = gettextRcmdr("Number of samples (columns):")), 
#         samplesEntry, sticky = "w")
#     tkgrid(tklabel(top, text = gettextRcmdr("Parameters:"), fg = "blue"), 
#         columnspan = 4, sticky = "w")
#     tkgrid(tklabel(top, text = gettextRcmdr("from (lower limit)")), 
#         from1Entry, sticky = "w")
#     tkgrid(tklabel(top, text = gettextRcmdr("to (upper limit)")), 
#         to1Entry, sticky = "w")
#     tkgrid(tklabel(top, text = gettextRcmdr("by (step size)")), 
#         by1Entry, sticky = "w")
#     tkgrid(tklabel(locationFrame, text = gettextRcmdr("Store values in:"), 
#         fg = "blue"), columnspan = 4, sticky = "w")
#     tkgrid(tklabel(locationFrame, text = gettextRcmdr("Active Dataset")), 
#         addtoactiveButton, sticky = "w")
#     tkgrid(tklabel(locationFrame, text = "New Dataset"), newDataButton, 
#         sticky = "w")
#     tkgrid.configure(samplesEntry, sticky = "w")
#     tkgrid.configure(from1Entry, sticky = "w")
#     tkgrid.configure(to1Entry, sticky = "w")
#     tkgrid.configure(by1Entry, sticky = "w")
#     tkgrid(locationFrame, sticky = "w")
#     tkgrid(buttonsFrame, sticky = "w", columnspan = 2)
#     dialogSuffix(rows = 6, columns = 1, focus = samplesEntry)
# }
# 
# 
# `expSimulate.ipsur` <-
# function () 
# {
#     initializeDialog(title = gettextRcmdr("Simulate Exponential Variates"))
#     parameterFrame <- tkframe(top)
#     locationFrame <- tkframe(top)
#     if (!is.character(ActiveDataSet())) {
#         locVariable <- tclVar("new")
#     }
#     else {
#         locVariable <- tclVar("add")
#     }
#     addtoactiveButton <- tkradiobutton(locationFrame, variable = locVariable, 
#         value = "add")
#     newDataButton <- tkradiobutton(locationFrame, variable = locVariable, 
#         value = "new")
#     samplesVar <- tclVar("1")
#     samplesEntry <- tkentry(top, width = "6", textvariable = samplesVar)
#     rateVar <- tclVar("1")
#     rateEntry <- tkentry(top, width = "6", textvariable = rateVar)
#     onOK <- function() {
#         nsamples <- round(as.numeric(tclvalue(samplesVar)))
#         rate <- tclvalue(rateVar)
#         if (is.na(nsamples)) {
#             errorCondition(recall = expSimulate.ipsur, message = gettextRcmdr("Number of samples must be a positive integer."))
#             return()
#         }
#         if (is.na(rate)) {
#             errorCondition(recall = expSimulate.ipsur, message = gettextRcmdr("The rate parameter was not specified."))
#             return()
#         }
#         closeDialog()
#         store <- tclvalue(locVariable)
#         if (store == "new") {
#             initializeDialog(title = gettextRcmdr("Simulation Dataset"))
#             dsname <- tclVar("Simset")
#             entryDsname <- tkentry(top, width = "20", textvariable = dsname)
#             newDataSS <- tclVar("100")
#             entryNewDataSS <- tkentry(top, width = "6", textvariable = newDataSS)
#             onOK <- function() {
#                 dsnameValue <- trim.blanks(tclvalue(dsname))
#                 newSS <- round(as.numeric(tclvalue(newDataSS)))
#                 closeDialog()
#                 if (dsnameValue == "") {
#                   errorCondition(recall = expSimulate.ipsur, 
#                     message = gettextRcmdr("You must enter the name of a data set."))
#                   return()
#                 }
#                 if (!is.valid.name(dsnameValue)) {
#                   errorCondition(recall = expSimulate.ipsur, 
#                     message = paste("\"", dsnameValue, "\" ", 
#                       gettextRcmdr("is not a valid name."), sep = ""))
#                   return()
#                 }
#                 if (is.element(dsnameValue, listDataSets())) {
#                   if ("no" == tclvalue(checkReplace(dsnameValue, 
#                     gettextRcmdr("Data set")))) {
#                     expSimulate.ipsur()
#                     return()
#                   }
#                 }
#                 if (is.na(newSS)) {
#                   errorCondition(recall = expSimulate.ipsur, 
#                     message = gettextRcmdr("Sample Size must be a positive integer."))
#                   return()
#                 }
#                 UpdateexpsimNumber()
#                 justDoIt(paste(dsnameValue, " = data.frame(exp.sim", 
#                   getRcmdr("expsimNumber"), "=1:", newSS, ")", 
#                   sep = ""))
#                 logger(paste(dsnameValue, "has been initialized."))
#                 for (k in getRcmdr("expsimNumber"):(nsamples + 
#                   getRcmdr("expsimNumber") - 1)) {
#                   justDoIt(paste(dsnameValue, "$exp.sim", k, 
#                     " <- rexp(", newSS, ", rate=", rate, ")", 
#                     sep = ""))
#                 }
#                 activeDataSet(dsnameValue)
#                 putRcmdr("expsimNumber", k)
#                 if (nsamples == 1) {
#                   logger(paste("There was 1 exponential variate sample stored in ", 
#                     dsnameValue, ".", sep = ""))
#                 }
#                 else {
#                   logger(paste("There were ", nsamples, "  exponential variate samples stored in ", 
#                     dsnameValue, ".", sep = ""))
#                 }
#             }
#             OKCancelHelp(helpSubject = "rexp")
#             tkgrid(tklabel(top, text = gettextRcmdr("Enter name for data set:")), 
#                 entryDsname, sticky = "e")
#             tkgrid(tklabel(top, text = gettextRcmdr("Sample Size (rows):")), 
#                 entryNewDataSS, sticky = "e")
#             tkgrid(buttonsFrame, columnspan = "2", sticky = "w")
#             tkgrid.configure(entryDsname, sticky = "w")
#             tkgrid.configure(entryNewDataSS, sticky = "w")
#             tkfocus(CommanderWindow())
#             dialogSuffix(rows = 2, columns = 2, focus = entryDsname)
#         }
#         else {
#             if (!is.character(ActiveDataSet())) {
#                 errorCondition(recall = expSimulate.ipsur, message = gettextRcmdr("There is no active data set."))
#                 return()
#             }
#             .activeDataSet <- ActiveDataSet()
#             justDoIt(paste("samplesn <- dim(", .activeDataSet, 
#                 ")[1]", sep = ""))
#             UpdateexpsimNumber()
#             for (k in getRcmdr("expsimNumber"):(nsamples + getRcmdr("expsimNumber") - 
#                 1)) {
#                 justDoIt(paste(.activeDataSet, "$exp.sim", k, 
#                   " <- rexp(", samplesn, ", rate=", rate, ")", 
#                   sep = ""))
#             }
#             activeDataSet(.activeDataSet)
#             putRcmdr("expsimNumber", k)
#             if (nsamples == 1) {
#                 logger(paste("There was 1 exponential variate sample stored in ", 
#                   .activeDataSet, ".", sep = ""))
#             }
#             else {
#                 logger(paste("There were ", nsamples, " exponential variate samples stored in ", 
#                   .activeDataSet, ".", sep = ""))
#             }
#         }
#         tkfocus(CommanderWindow())
#     }
#     OKCancelHelp(helpSubject = "rexp")
#     tkgrid(tklabel(top, text = gettextRcmdr("Number of samples (columns):")), 
#         samplesEntry, sticky = "w")
#     tkgrid(tklabel(top, text = gettextRcmdr("Parameters:"), fg = "blue"), 
#         columnspan = 4, sticky = "w")
#     tkgrid(tklabel(top, text = gettextRcmdr("rate (of arrivals in unit time)")), 
#         rateEntry, sticky = "w")
#     tkgrid(tklabel(locationFrame, text = gettextRcmdr("Store values in:"), 
#         fg = "blue"), columnspan = 4, sticky = "w")
#     tkgrid(tklabel(locationFrame, text = gettextRcmdr("Active Dataset")), 
#         addtoactiveButton, sticky = "w")
#     tkgrid(tklabel(locationFrame, text = "New Dataset"), newDataButton, 
#         sticky = "w")
#     tkgrid.configure(samplesEntry, sticky = "w")
#     tkgrid.configure(rateEntry, sticky = "w")
#     tkgrid(locationFrame, sticky = "w")
#     tkgrid(buttonsFrame, sticky = "w", columnspan = 2)
#     dialogSuffix(rows = 6, columns = 1, focus = samplesEntry)
# }
# 
# 
# `fSimulate.ipsur` <-
# function () 
# {
#     initializeDialog(title = gettextRcmdr("Simulate F Variates"))
#     parameterFrame <- tkframe(top)
#     locationFrame <- tkframe(top)
#     if (!is.character(ActiveDataSet())) {
#         locVariable <- tclVar("new")
#     }
#     else {
#         locVariable <- tclVar("add")
#     }
#     addtoactiveButton <- tkradiobutton(locationFrame, variable = locVariable, 
#         value = "add")
#     newDataButton <- tkradiobutton(locationFrame, variable = locVariable, 
#         value = "new")
#     samplesVar <- tclVar("1")
#     samplesEntry <- tkentry(top, width = "6", textvariable = samplesVar)
#     df1Var <- tclVar("1")
#     df1Entry <- tkentry(top, width = "6", textvariable = df1Var)
#     df2Var <- tclVar("1")
#     df2Entry <- tkentry(top, width = "6", textvariable = df2Var)
#     ncpVar <- tclVar("0")
#     ncpEntry <- tkentry(top, width = "6", textvariable = ncpVar)
#     onOK <- function() {
#         nsamples <- round(as.numeric(tclvalue(samplesVar)))
#         df1 <- tclvalue(df1Var)
#         df2 <- tclvalue(df2Var)
#         ncp <- tclvalue(ncpVar)
#         if (is.na(nsamples)) {
#             errorCondition(recall = fSimulate.ipsur, message = gettextRcmdr("Number of samples must be a positive integer."))
#             return()
#         }
#         if (is.na(df1) || is.na(df2)) {
#             errorCondition(recall = fSimulate.ipsur, message = gettextRcmdr("Degrees of freedom were not specified."))
#             return()
#         }
#         if (is.na(ncp)) {
#             errorCondition(recall = fSimulate.ipsur, message = gettextRcmdr("The noncentrality parameter was not specified."))
#             return()
#         }
#         closeDialog()
#         store <- tclvalue(locVariable)
#         if (store == "new") {
#             initializeDialog(title = gettextRcmdr("Simulation Dataset"))
#             dsname <- tclVar("Simset")
#             entryDsname <- tkentry(top, width = "20", textvariable = dsname)
#             newDataSS <- tclVar("100")
#             entryNewDataSS <- tkentry(top, width = "6", textvariable = newDataSS)
#             onOK <- function() {
#                 dsnameValue <- trim.blanks(tclvalue(dsname))
#                 newSS <- round(as.numeric(tclvalue(newDataSS)))
#                 closeDialog()
#                 if (dsnameValue == "") {
#                   errorCondition(recall = fSimulate.ipsur, message = gettextRcmdr("You must enter the name of a data set."))
#                   return()
#                 }
#                 if (!is.valid.name(dsnameValue)) {
#                   errorCondition(recall = fSimulate.ipsur, message = paste("\"", 
#                     dsnameValue, "\" ", gettextRcmdr("is not a valid name."), 
#                     sep = ""))
#                   return()
#                 }
#                 if (is.element(dsnameValue, listDataSets())) {
#                   if ("no" == tclvalue(checkReplace(dsnameValue, 
#                     gettextRcmdr("Data set")))) {
#                     fSimulate.ipsur()
#                     return()
#                   }
#                 }
#                 if (is.na(newSS)) {
#                   errorCondition(recall = fSimulate.ipsur, message = gettextRcmdr("Sample Size must be a positive integer."))
#                   return()
#                 }
#                 UpdatefsimNumber()
#                 justDoIt(paste(dsnameValue, " = data.frame(f.sim", 
#                   getRcmdr("fsimNumber"), "=1:", newSS, ")", 
#                   sep = ""))
#                 logger(paste(dsnameValue, " has been initialized."))
#                 for (k in getRcmdr("fsimNumber"):(nsamples + 
#                   getRcmdr("fsimNumber") - 1)) {
#                   justDoIt(paste(dsnameValue, "$f.sim", k, " <- rf(", 
#                     newSS, ", df1=", df1, ", df2=", df2, ", ncp=", 
#                     ncp, ")", sep = ""))
#                 }
#                 activeDataSet(dsnameValue)
#                 putRcmdr("fsimNumber", k)
#                 if (nsamples == 1) {
#                   logger(paste("There was 1 F variate sample stored in ", 
#                     dsnameValue, ".", sep = ""))
#                 }
#                 else {
#                   logger(paste("There were ", nsamples, " F variate samples stored in ", 
#                     dsnameValue, ".", sep = ""))
#                 }
#             }
#             OKCancelHelp(helpSubject = "rf")
#             tkgrid(tklabel(top, text = gettextRcmdr("Enter name for data set:")), 
#                 entryDsname, sticky = "e")
#             tkgrid(tklabel(top, text = gettextRcmdr("Sample Size (rows):")), 
#                 entryNewDataSS, sticky = "e")
#             tkgrid(buttonsFrame, columnspan = "2", sticky = "w")
#             tkgrid.configure(entryDsname, sticky = "w")
#             tkgrid.configure(entryNewDataSS, sticky = "w")
#             tkfocus(CommanderWindow())
#             dialogSuffix(rows = 2, columns = 2, focus = entryDsname)
#         }
#         else {
#             if (!is.character(ActiveDataSet())) {
#                 errorCondition(recall = fSimulate.ipsur, message = gettextRcmdr("There is no active data set."))
#                 return()
#             }
#             .activeDataSet <- ActiveDataSet()
#             justDoIt(paste("samplesn <- dim(", .activeDataSet, 
#                 ")[1]", sep = ""))
#             UpdatefsimNumber()
#             for (k in getRcmdr("fsimNumber"):(nsamples + getRcmdr("fsimNumber") - 
#                 1)) {
#                 justDoIt(paste(.activeDataSet, "$f.sim", k, " <- rf(", 
#                   samplesn, ", df1=", df1, ", df2=", df2, ", ncp=", 
#                   ncp, ")", sep = ""))
#             }
#             activeDataSet(.activeDataSet)
#             putRcmdr("fsimNumber", k)
#             if (nsamples == 1) {
#                 logger(paste("There was 1 F variate sample stored in ", 
#                   .activeDataSet, ".", sep = ""))
#             }
#             else {
#                 logger(paste("There were ", nsamples, " F variate samples stored in ", 
#                   .activeDataSet, ".", sep = ""))
#             }
#         }
#         tkfocus(CommanderWindow())
#     }
#     OKCancelHelp(helpSubject = "rf")
#     tkgrid(tklabel(top, text = gettextRcmdr("Number of samples (columns):")), 
#         samplesEntry, sticky = "w")
#     tkgrid(tklabel(top, text = gettextRcmdr("Parameters:"), fg = "blue"), 
#         columnspan = 4, sticky = "w")
#     tkgrid(tklabel(top, text = gettextRcmdr("df1 (num degrees of freedom)")), 
#         df1Entry, sticky = "w")
#     tkgrid(tklabel(top, text = gettextRcmdr("df2 (denom degrees of freedom)")), 
#         df2Entry, sticky = "w")
#     tkgrid(tklabel(top, text = gettextRcmdr("ncp (noncentrality parameter)")), 
#         ncpEntry, sticky = "w")
#     tkgrid(tklabel(locationFrame, text = gettextRcmdr("Store values in:"), 
#         fg = "blue"), columnspan = 4, sticky = "w")
#     tkgrid(tklabel(locationFrame, text = gettextRcmdr("Active Dataset")), 
#         addtoactiveButton, sticky = "w")
#     tkgrid(tklabel(locationFrame, text = "New Dataset"), newDataButton, 
#         sticky = "w")
#     tkgrid.configure(samplesEntry, sticky = "w")
#     tkgrid.configure(df1Entry, sticky = "w")
#     tkgrid.configure(df2Entry, sticky = "w")
#     tkgrid.configure(ncpEntry, sticky = "w")
#     tkgrid(locationFrame, sticky = "w")
#     tkgrid(buttonsFrame, sticky = "w", columnspan = 2)
#     dialogSuffix(rows = 6, columns = 1, focus = samplesEntry)
# }
# 
# 
# `gammaSimulate.ipsur` <-
# function () 
# {
#     initializeDialog(title = gettextRcmdr("Simulate Gamma Variates"))
#     parameterFrame <- tkframe(top)
#     locationFrame <- tkframe(top)
#     if (!is.character(ActiveDataSet())) {
#         locVariable <- tclVar("new")
#     }
#     else {
#         locVariable <- tclVar("add")
#     }
#     addtoactiveButton <- tkradiobutton(locationFrame, variable = locVariable, 
#         value = "add")
#     newDataButton <- tkradiobutton(locationFrame, variable = locVariable, 
#         value = "new")
#     samplesVar <- tclVar("1")
#     samplesEntry <- tkentry(top, width = "6", textvariable = samplesVar)
#     shapeVar <- tclVar("1")
#     shapeEntry <- tkentry(top, width = "6", textvariable = shapeVar)
#     scale1Var <- tclVar("1")
#     scale1Entry <- tkentry(top, width = "6", textvariable = scale1Var)
#     onOK <- function() {
#         nsamples <- round(as.numeric(tclvalue(samplesVar)))
#         shape <- tclvalue(shapeVar)
#         scale1 <- tclvalue(scale1Var)
#         if (is.na(nsamples)) {
#             errorCondition(recall = gammaSimulate.ipsur, message = gettextRcmdr("Number of samples must be a positive integer."))
#             return()
#         }
#         if (is.na(shape)) {
#             errorCondition(recall = gammaSimulate.ipsur, message = gettextRcmdr("The shape parameter was not specified."))
#             return()
#         }
#         if (is.na(scale1)) {
#             errorCondition(recall = gammaSimulate.ipsur, message = gettextRcmdr("The rate parameter was not specified."))
#             return()
#         }
#         closeDialog()
#         store <- tclvalue(locVariable)
#         if (store == "new") {
#             initializeDialog(title = gettextRcmdr("Simulation Dataset"))
#             dsname <- tclVar("Simset")
#             entryDsname <- tkentry(top, width = "20", textvariable = dsname)
#             newDataSS <- tclVar("100")
#             entryNewDataSS <- tkentry(top, width = "6", textvariable = newDataSS)
#             onOK <- function() {
#                 dsnameValue <- trim.blanks(tclvalue(dsname))
#                 newSS <- round(as.numeric(tclvalue(newDataSS)))
#                 closeDialog()
#                 if (dsnameValue == "") {
#                   errorCondition(recall = gammaSimulate.ipsur, 
#                     message = gettextRcmdr("You must enter the name of a data set."))
#                   return()
#                 }
#                 if (!is.valid.name(dsnameValue)) {
#                   errorCondition(recall = gammaSimulate.ipsur, 
#                     message = paste("\"", dsnameValue, "\" ", 
#                       gettextRcmdr("is not a valid name."), sep = ""))
#                   return()
#                 }
#                 if (is.element(dsnameValue, listDataSets())) {
#                   if ("no" == tclvalue(checkReplace(dsnameValue, 
#                     gettextRcmdr("Data set")))) {
#                     gammaSimulate.ipsur()
#                     return()
#                   }
#                 }
#                 if (is.na(newSS)) {
#                   errorCondition(recall = gammaSimulate.ipsur, 
#                     message = gettextRcmdr("Sample Size must be a positive integer."))
#                   return()
#                 }
#                 UpdategammasimNumber()
#                 justDoIt(paste(dsnameValue, " = data.frame(gamma.sim", 
#                   getRcmdr("gammasimNumber"), "=1:", newSS, ")", 
#                   sep = ""))
#                 logger(paste(dsnameValue, "has been initialized."))
#                 for (k in getRcmdr("gammasimNumber"):(nsamples + 
#                   getRcmdr("gammasimNumber") - 1)) {
#                   justDoIt(paste(dsnameValue, "$gamma.sim", k, 
#                     " <- rgamma(", newSS, ", shape=", shape, 
#                     ", scale=", scale1, ")", sep = ""))
#                 }
#                 activeDataSet(dsnameValue)
#                 putRcmdr("gammasimNumber", k)
#                 if (nsamples == 1) {
#                   logger(paste("There was 1 gamma variate sample stored in ", 
#                     dsnameValue, ".", sep = ""))
#                 }
#                 else {
#                   logger(paste("There were ", nsamples, " gamma variate samples stored in ", 
#                     dsnameValue, ".", sep = ""))
#                 }
#             }
#             OKCancelHelp(helpSubject = "rgamma")
#             tkgrid(tklabel(top, text = gettextRcmdr("Enter name for data set:")), 
#                 entryDsname, sticky = "e")
#             tkgrid(tklabel(top, text = gettextRcmdr("Sample Size (rows):")), 
#                 entryNewDataSS, sticky = "e")
#             tkgrid(buttonsFrame, columnspan = "2", sticky = "w")
#             tkgrid.configure(entryDsname, sticky = "w")
#             tkgrid.configure(entryNewDataSS, sticky = "w")
#             tkfocus(CommanderWindow())
#             dialogSuffix(rows = 2, columns = 2, focus = entryDsname)
#         }
#         else {
#             if (!is.character(ActiveDataSet())) {
#                 errorCondition(recall = gammaSimulate.ipsur, 
#                   message = gettextRcmdr("There is no active data set."))
#                 return()
#             }
#             .activeDataSet <- ActiveDataSet()
#             justDoIt(paste("samplesn <- dim(", .activeDataSet, 
#                 ")[1]", sep = ""))
#             UpdategammasimNumber()
#             for (k in getRcmdr("gammasimNumber"):(nsamples + 
#                 getRcmdr("gammasimNumber") - 1)) {
#                 justDoIt(paste(.activeDataSet, "$gamma.sim", 
#                   k, " <- rgamma(", samplesn, ", shape=", shape, 
#                   ", rate=", scale1, ")", sep = ""))
#             }
#             activeDataSet(.activeDataSet)
#             putRcmdr("gammasimNumber", k)
#             if (nsamples == 1) {
#                 logger(paste("There was 1 gamma variate sample stored in ", 
#                   .activeDataSet, ".", sep = ""))
#             }
#             else {
#                 logger(paste("There were ", nsamples, " gamma variate samples stored in ", 
#                   .activeDataSet, ".", sep = ""))
#             }
#         }
#         tkfocus(CommanderWindow())
#     }
#     OKCancelHelp(helpSubject = "rgamma")
#     tkgrid(tklabel(top, text = gettextRcmdr("Number of samples (columns):")), 
#         samplesEntry, sticky = "w")
#     tkgrid(tklabel(top, text = gettextRcmdr("Parameters:"), fg = "blue"), 
#         columnspan = 4, sticky = "w")
#     tkgrid(tklabel(top, text = gettextRcmdr("shape")), shapeEntry, 
#         sticky = "w")
#     tkgrid(tklabel(top, text = gettextRcmdr("rate (= 1/scale)")), 
#         scale1Entry, sticky = "w")
#     tkgrid(tklabel(locationFrame, text = gettextRcmdr("Store values in:"), 
#         fg = "blue"), columnspan = 4, sticky = "w")
#     tkgrid(tklabel(locationFrame, text = gettextRcmdr("Active Dataset")), 
#         addtoactiveButton, sticky = "w")
#     tkgrid(tklabel(locationFrame, text = "New Dataset"), newDataButton, 
#         sticky = "w")
#     tkgrid.configure(samplesEntry, sticky = "w")
#     tkgrid.configure(shapeEntry, sticky = "w")
#     tkgrid.configure(scale1Entry, sticky = "w")
#     tkgrid(locationFrame, sticky = "w")
#     tkgrid(buttonsFrame, sticky = "w", columnspan = 2)
#     dialogSuffix(rows = 6, columns = 1, focus = samplesEntry)
# }
# 
# 
# `geomSimulate.ipsur` <-
# function () 
# {
#     initializeDialog(title = gettextRcmdr("Simulate Geometric Variates"))
#     parameterFrame <- tkframe(top)
#     locationFrame <- tkframe(top)
#     if (!is.character(ActiveDataSet())) {
#         locVariable <- tclVar("new")
#     }
#     else {
#         locVariable <- tclVar("add")
#     }
#     addtoactiveButton <- tkradiobutton(locationFrame, variable = locVariable, 
#         value = "add")
#     newDataButton <- tkradiobutton(locationFrame, variable = locVariable, 
#         value = "new")
#     samplesVar <- tclVar("1")
#     samplesEntry <- tkentry(top, width = "6", textvariable = samplesVar)
#     probVar <- tclVar("0.5")
#     probEntry <- tkentry(top, width = "6", textvariable = probVar)
#     onOK <- function() {
#         nsamples <- round(as.numeric(tclvalue(samplesVar)))
#         prob <- tclvalue(probVar)
#         if (is.na(nsamples)) {
#             errorCondition(recall = geomSimulate.ipsur, message = gettextRcmdr("Number of samples must be a positive integer."))
#             return()
#         }
#         if (is.na(prob)) {
#             errorCondition(recall = geomSimulate.ipsur, message = gettextRcmdr("The probability of success was not specified."))
#             return()
#         }
#         closeDialog()
#         store <- tclvalue(locVariable)
#         if (store == "new") {
#             initializeDialog(title = gettextRcmdr("Simulation Dataset"))
#             dsname <- tclVar("Simset")
#             entryDsname <- tkentry(top, width = "20", textvariable = dsname)
#             newDataSS <- tclVar("100")
#             entryNewDataSS <- tkentry(top, width = "6", textvariable = newDataSS)
#             onOK <- function() {
#                 dsnameValue <- trim.blanks(tclvalue(dsname))
#                 newSS <- round(as.numeric(tclvalue(newDataSS)))
#                 closeDialog()
#                 if (dsnameValue == "") {
#                   errorCondition(recall = geomSimulate.ipsur, 
#                     message = gettextRcmdr("You must enter the name of a data set."))
#                   return()
#                 }
#                 if (!is.valid.name(dsnameValue)) {
#                   errorCondition(recall = geomSimulate.ipsur, 
#                     message = paste("\"", dsnameValue, "\" ", 
#                       gettextRcmdr("is not a valid name."), sep = ""))
#                   return()
#                 }
#                 if (is.element(dsnameValue, listDataSets())) {
#                   if ("no" == tclvalue(checkReplace(dsnameValue, 
#                     gettextRcmdr("Data set")))) {
#                     geomSimulate.ipsur()
#                     return()
#                   }
#                 }
#                 if (is.na(newSS)) {
#                   errorCondition(recall = geomSimulate.ipsur, 
#                     message = gettextRcmdr("Sample Size must be a positive integerr."))
#                   return()
#                 }
#                 UpdategeomsimNumber()
#                 justDoIt(paste(dsnameValue, " = data.frame(geom.sim", 
#                   getRcmdr("geomsimNumber"), "=1:", newSS, ")", 
#                   sep = ""))
#                 logger(paste(dsnameValue, "has been initialized."))
#                 for (k in getRcmdr("geomsimNumber"):(nsamples + 
#                   getRcmdr("geomsimNumber") - 1)) {
#                   justDoIt(paste(dsnameValue, "$geom.sim", k, 
#                     " <- rgeom(", newSS, ", prob=", prob, ")", 
#                     sep = ""))
#                 }
#                 activeDataSet(dsnameValue)
#                 putRcmdr("geomsimNumber", k)
#                 if (nsamples == 1) {
#                   logger(paste("There was 1 geometric variate sample stored in ", 
#                     dsnameValue, ".", sep = ""))
#                 }
#                 else {
#                   logger(paste("There were ", nsamples, " geometric variate samples stored in ", 
#                     dsnameValue, ".", sep = ""))
#                 }
#             }
#             OKCancelHelp(helpSubject = "rgeom")
#             tkgrid(tklabel(top, text = gettextRcmdr("Enter name for data set:")), 
#                 entryDsname, sticky = "e")
#             tkgrid(tklabel(top, text = gettextRcmdr("Sample Size (rows):")), 
#                 entryNewDataSS, sticky = "e")
#             tkgrid(buttonsFrame, columnspan = "2", sticky = "w")
#             tkgrid.configure(entryDsname, sticky = "w")
#             tkgrid.configure(entryNewDataSS, sticky = "w")
#             tkfocus(CommanderWindow())
#             dialogSuffix(rows = 2, columns = 2, focus = entryDsname)
#         }
#         else {
#             if (!is.character(ActiveDataSet())) {
#                 errorCondition(recall = geomSimulate.ipsur, message = gettextRcmdr("There is no active data set."))
#                 return()
#             }
#             .activeDataSet <- ActiveDataSet()
#             justDoIt(paste("samplesn <- dim(", .activeDataSet, 
#                 ")[1]", sep = ""))
#             UpdategeomsimNumber()
#             for (k in getRcmdr("geomsimNumber"):(nsamples + getRcmdr("geomsimNumber") - 
#                 1)) {
#                 justDoIt(paste(.activeDataSet, "$geom.sim", k, 
#                   " <- rgeom(", samplesn, ", prob=", prob, ")", 
#                   sep = ""))
#             }
#             activeDataSet(.activeDataSet)
#             putRcmdr("geomsimNumber", k)
#             if (nsamples == 1) {
#                 logger(paste("There was 1 geometric variate sample stored in ", 
#                   .activeDataSet, ".", sep = ""))
#             }
#             else {
#                 logger(paste("There were ", nsamples, " geometric variate samples stored in ", 
#                   .activeDataSet, ".", sep = ""))
#             }
#         }
#         tkfocus(CommanderWindow())
#     }
#     OKCancelHelp(helpSubject = "rgeom")
#     tkgrid(tklabel(top, text = gettextRcmdr("Number of samples (columns):")), 
#         samplesEntry, sticky = "w")
#     tkgrid(tklabel(top, text = gettextRcmdr("Parameters:"), fg = "blue"), 
#         columnspan = 4, sticky = "w")
#     tkgrid(tklabel(top, text = gettextRcmdr("prob (of success in each trial)")), 
#         probEntry, sticky = "w")
#     tkgrid(tklabel(locationFrame, text = gettextRcmdr("Store values in:"), 
#         fg = "blue"), columnspan = 4, sticky = "w")
#     tkgrid(tklabel(locationFrame, text = gettextRcmdr("Active Dataset")), 
#         addtoactiveButton, sticky = "w")
#     tkgrid(tklabel(locationFrame, text = "New Dataset"), newDataButton, 
#         sticky = "w")
#     tkgrid.configure(samplesEntry, sticky = "w")
#     tkgrid.configure(probEntry, sticky = "w")
#     tkgrid(locationFrame, sticky = "w")
#     tkgrid(buttonsFrame, sticky = "w", columnspan = 2)
#     dialogSuffix(rows = 6, columns = 1, focus = samplesEntry)
# }
# 
# 
# `hyperSimulate.ipsur` <-
# function () 
# {
#     initializeDialog(title = gettextRcmdr("Simulate Hypergeometric Variates"))
#     parameterFrame <- tkframe(top)
#     locationFrame <- tkframe(top)
#     if (!is.character(ActiveDataSet())) {
#         locVariable <- tclVar("new")
#     }
#     else {
#         locVariable <- tclVar("add")
#     }
#     addtoactiveButton <- tkradiobutton(locationFrame, variable = locVariable, 
#         value = "add")
#     newDataButton <- tkradiobutton(locationFrame, variable = locVariable, 
#         value = "new")
#     samplesVar <- tclVar("1")
#     samplesEntry <- tkentry(top, width = "6", textvariable = samplesVar)
#     mVar <- tclVar("1")
#     mEntry <- tkentry(top, width = "6", textvariable = mVar)
#     nVar <- tclVar("1")
#     nEntry <- tkentry(top, width = "6", textvariable = nVar)
#     k1Var <- tclVar("1")
#     k1Entry <- tkentry(top, width = "6", textvariable = k1Var)
#     onOK <- function() {
#         nsamples <- round(as.numeric(tclvalue(samplesVar)))
#         m <- tclvalue(mVar)
#         n <- tclvalue(nVar)
#         k1 <- tclvalue(k1Var)
#         if (is.na(nsamples) || nsamples < 1) {
#             errorCondition(recall = hyperSimulate.ipsur, message = gettextRcmdr("Number of samples must be a positive integer."))
#             return()
#         }
#         if (is.na(m)) {
#             errorCondition(recall = hyperSimulate.ipsur, message = gettextRcmdr("The m parameter was not specified."))
#             return()
#         }
#         if (is.na(n)) {
#             errorCondition(recall = hyperSimulate.ipsur, message = gettextRcmdr("The n parameter was not specified."))
#             return()
#         }
#         if (is.na(k1)) {
#             errorCondition(recall = hyperSimulate.ipsur, message = gettextRcmdr("The k parameter was not specified."))
#             return()
#         }
#         closeDialog()
#         store <- tclvalue(locVariable)
#         if (store == "new") {
#             initializeDialog(title = gettextRcmdr("Simulation Dataset"))
#             dsname <- tclVar("Simset")
#             entryDsname <- tkentry(top, width = "20", textvariable = dsname)
#             newDataSS <- tclVar("100")
#             entryNewDataSS <- tkentry(top, width = "6", textvariable = newDataSS)
#             onOK <- function() {
#                 dsnameValue <- trim.blanks(tclvalue(dsname))
#                 newSS <- round(as.numeric(tclvalue(newDataSS)))
#                 closeDialog()
#                 if (dsnameValue == "") {
#                   errorCondition(recall = hyperSimulate.ipsur, 
#                     message = gettextRcmdr("You must enter the name of a data set."))
#                   return()
#                 }
#                 if (!is.valid.name(dsnameValue)) {
#                   errorCondition(recall = hyperSimulate.ipsur, 
#                     message = paste("\"", dsnameValue, "\" ", 
#                       gettextRcmdr("is not a valid name."), sep = ""))
#                   return()
#                 }
#                 if (is.element(dsnameValue, listDataSets())) {
#                   if ("no" == tclvalue(checkReplace(dsnameValue, 
#                     gettextRcmdr("Data set")))) {
#                     hyperSimulate.ipsur()
#                     return()
#                   }
#                 }
#                 if (is.na(newSS)) {
#                   errorCondition(recall = hyperSimulate.ipsur, 
#                     message = gettextRcmdr("Sample Size must be a positive integer."))
#                   return()
#                 }
#                 UpdatehypersimNumber()
#                 justDoIt(paste(dsnameValue, " = data.frame(hyper.sim", 
#                   getRcmdr("hypersimNumber"), "=1:", newSS, ")", 
#                   sep = ""))
#                 logger(paste(dsnameValue, "has been initialized."))
#                 for (k in getRcmdr("hypersimNumber"):(nsamples + 
#                   getRcmdr("hypersimNumber") - 1)) {
#                   justDoIt(paste(dsnameValue, "$hyper.sim", k, 
#                     " <- rhyper(", newSS, ", m=", m, ", n=", 
#                     n, ", k=", k1, ")", sep = ""))
#                 }
#                 activeDataSet(dsnameValue)
#                 putRcmdr("hypersimNumber", k)
#                 if (nsamples == 1) {
#                   logger(paste("There was 1 hypergeometric variate sample stored in ", 
#                     dsnameValue, ".", sep = ""))
#                 }
#                 else {
#                   logger(paste("There were ", nsamples, " hyergeometric variate samples stored in ", 
#                     dsnameValue, ".", sep = ""))
#                 }
#             }
#             OKCancelHelp(helpSubject = "rhyper")
#             tkgrid(tklabel(top, text = gettextRcmdr("Enter name for data set:")), 
#                 entryDsname, sticky = "e")
#             tkgrid(tklabel(top, text = gettextRcmdr("Sample Size (rows):")), 
#                 entryNewDataSS, sticky = "e")
#             tkgrid(buttonsFrame, columnspan = "2", sticky = "w")
#             tkgrid.configure(entryDsname, sticky = "w")
#             tkgrid.configure(entryNewDataSS, sticky = "w")
#             tkfocus(CommanderWindow())
#             dialogSuffix(rows = 2, columns = 2, focus = entryDsname)
#         }
#         else {
#             if (!is.character(ActiveDataSet())) {
#                 errorCondition(recall = hyperSimulate.ipsur, 
#                   message = gettextRcmdr("There is no active data set."))
#                 return()
#             }
#             .activeDataSet <- ActiveDataSet()
#             justDoIt(paste("samplesn <- dim(", .activeDataSet, 
#                 ")[1]", sep = ""))
#             UpdatehypersimNumber()
#             for (k in getRcmdr("hypersimNumber"):(nsamples + 
#                 getRcmdr("hypersimNumber") - 1)) {
#                 justDoIt(paste(.activeDataSet, "$hyper.sim", 
#                   k, " <- rhyper(", samplesn, ", m=", m, ", n=", 
#                   n, ", k=", k1, ")", sep = ""))
#             }
#             activeDataSet(.activeDataSet)
#             putRcmdr("hypersimNumber", k)
#             if (nsamples == 1) {
#                 logger(paste("There was 1 hypergeometric variate sample stored in ", 
#                   .activeDataSet, ".", sep = ""))
#             }
#             else {
#                 logger(paste("There were ", nsamples, " hypergeometric variate samples stored in ", 
#                   .activeDataSet, ".", sep = ""))
#             }
#         }
#         tkfocus(CommanderWindow())
#     }
#     OKCancelHelp(helpSubject = "rhyper")
#     tkgrid(tklabel(top, text = gettextRcmdr("Number of samples (columns):")), 
#         samplesEntry, sticky = "w")
#     tkgrid(tklabel(top, text = gettextRcmdr("Parameters:"), fg = "blue"), 
#         columnspan = 4, sticky = "w")
#     tkgrid(tklabel(top, text = gettextRcmdr("m (num of white balls in the urn)")), 
#         mEntry, sticky = "w")
#     tkgrid(tklabel(top, text = gettextRcmdr("n (num of black balls in the urn)")), 
#         nEntry, sticky = "w")
#     tkgrid(tklabel(top, text = gettextRcmdr("k (num of balls drawn from the urn)")), 
#         k1Entry, sticky = "w")
#     tkgrid(tklabel(locationFrame, text = gettextRcmdr("Store values in:"), 
#         fg = "blue"), columnspan = 4, sticky = "w")
#     tkgrid(tklabel(locationFrame, text = gettextRcmdr("Active Dataset")), 
#         addtoactiveButton, sticky = "w")
#     tkgrid(tklabel(locationFrame, text = "New Dataset"), newDataButton, 
#         sticky = "w")
#     tkgrid.configure(samplesEntry, sticky = "w")
#     tkgrid.configure(mEntry, sticky = "w")
#     tkgrid.configure(nEntry, sticky = "w")
#     tkgrid.configure(k1Entry, sticky = "w")
#     tkgrid(locationFrame, sticky = "w")
#     tkgrid(buttonsFrame, sticky = "w", columnspan = 2)
#     dialogSuffix(rows = 6, columns = 1, focus = samplesEntry)
# }
# 
# 
# `lnormalSimulate.ipsur` <-
# function () 
# {
#     initializeDialog(title = gettextRcmdr("Simulate Log Normal Variates"))
#     parameterFrame <- tkframe(top)
#     locationFrame <- tkframe(top)
#     if (!is.character(ActiveDataSet())) {
#         locVariable <- tclVar("new")
#     }
#     else {
#         locVariable <- tclVar("add")
#     }
#     addtoactiveButton <- tkradiobutton(locationFrame, variable = locVariable, 
#         value = "add")
#     newDataButton <- tkradiobutton(locationFrame, variable = locVariable, 
#         value = "new")
#     samplesVar <- tclVar("1")
#     samplesEntry <- tkentry(top, width = "6", textvariable = samplesVar)
#     mulogVar <- tclVar("0")
#     mulogEntry <- tkentry(top, width = "6", textvariable = mulogVar)
#     sigmalogVar <- tclVar("1")
#     sigmalogEntry <- tkentry(top, width = "6", textvariable = sigmalogVar)
#     onOK <- function() {
#         nsamples <- round(as.numeric(tclvalue(samplesVar)))
#         mulog <- tclvalue(mulogVar)
#         sigmalog <- tclvalue(sigmalogVar)
#         if (is.na(nsamples)) {
#             errorCondition(recall = lnormalSimulate.ipsur, message = gettextRcmdr("Number of samples must be a positive integer."))
#             return()
#         }
#         if (is.na(mulog)) {
#             errorCondition(recall = lnormalSimulate.ipsur, message = gettextRcmdr("The mean was not specified."))
#             return()
#         }
#         if (is.na(sigmalog)) {
#             errorCondition(recall = lnormalSimulate.ipsur, message = gettextRcmdr("The standard deviation was not specified."))
#             return()
#         }
#         closeDialog()
#         store <- tclvalue(locVariable)
#         if (store == "new") {
#             initializeDialog(title = gettextRcmdr("Simulation Dataset"))
#             dsname <- tclVar("Simset")
#             entryDsname <- tkentry(top, width = "20", textvariable = dsname)
#             newDataSS <- tclVar("100")
#             entryNewDataSS <- tkentry(top, width = "6", textvariable = newDataSS)
#             onOK <- function() {
#                 dsnameValue <- trim.blanks(tclvalue(dsname))
#                 newSS <- round(as.numeric(tclvalue(newDataSS)))
#                 closeDialog()
#                 if (dsnameValue == "") {
#                   errorCondition(recall = lnormalSimulate.ipsur, 
#                     message = gettextRcmdr("You must enter the name of a data set."))
#                   return()
#                 }
#                 if (!is.valid.name(dsnameValue)) {
#                   errorCondition(recall = lnormalSimulate.ipsur, 
#                     message = paste("\"", dsnameValue, "\" ", 
#                       gettextRcmdr("is not a valid name."), sep = ""))
#                   return()
#                 }
#                 if (is.element(dsnameValue, listDataSets())) {
#                   if ("no" == tclvalue(checkReplace(dsnameValue, 
#                     gettextRcmdr("Data set")))) {
#                     lnormalSimulate.ipsur()
#                     return()
#                   }
#                 }
#                 if (is.na(newSS)) {
#                   errorCondition(recall = lnormalSimulate.ipsur, 
#                     message = gettextRcmdr("Sample Size must be a positive integer."))
#                   return()
#                 }
#                 UpdatelnormsimNumber()
#                 justDoIt(paste(dsnameValue, " = data.frame(lnorm.sim", 
#                   getRcmdr("lnormsimNumber"), "=1:", newSS, ")", 
#                   sep = ""))
#                 logger(paste(dsnameValue, "has been initialized."))
#                 for (k in getRcmdr("lnormsimNumber"):(nsamples + 
#                   getRcmdr("lnormsimNumber") - 1)) {
#                   justDoIt(paste(dsnameValue, "$lnorm.sim", k, 
#                     " <- rlnorm(", newSS, ", meanlog=", mulog, 
#                     ", sdlog=", sigmalog, ")", sep = ""))
#                 }
#                 activeDataSet(dsnameValue)
#                 putRcmdr("lnormsimNumber", k)
#                 if (nsamples == 1) {
#                   logger(paste("There was 1 log normal variate sample stored in ", 
#                     dsnameValue, ".", sep = ""))
#                 }
#                 else {
#                   logger(paste("There were ", nsamples, " log normal variate samples stored in ", 
#                     dsnameValue, ".", sep = ""))
#                 }
#             }
#             OKCancelHelp(helpSubject = "rlnorm")
#             tkgrid(tklabel(top, text = gettextRcmdr("Enter name for data set:")), 
#                 entryDsname, sticky = "e")
#             tkgrid(tklabel(top, text = gettextRcmdr("Sample Size (rows):")), 
#                 entryNewDataSS, sticky = "e")
#             tkgrid(buttonsFrame, columnspan = "2", sticky = "w")
#             tkgrid.configure(entryDsname, sticky = "w")
#             tkgrid.configure(entryNewDataSS, sticky = "w")
#             tkfocus(CommanderWindow())
#             dialogSuffix(rows = 2, columns = 2, focus = entryDsname)
#         }
#         else {
#             if (!is.character(ActiveDataSet())) {
#                 errorCondition(recall = lnormalSimulate.ipsur, 
#                   message = gettextRcmdr("There is no active data set."))
#                 return()
#             }
#             .activeDataSet <- ActiveDataSet()
#             justDoIt(paste("samplesn <- dim(", .activeDataSet, 
#                 ")[1]", sep = ""))
#             UpdatelnormsimNumber()
#             for (k in getRcmdr("lnormsimNumber"):(nsamples + 
#                 getRcmdr("lnormsimNumber") - 1)) {
#                 justDoIt(paste(.activeDataSet, "$lnorm.sim", 
#                   k, " <- rlnorm(", samplesn, ", meanlog=", mulog, 
#                   ", sdlog=", sigmalog, ")", sep = ""))
#             }
#             activeDataSet(.activeDataSet)
#             putRcmdr("lnormsimNumber", k)
#             if (nsamples == 1) {
#                 logger(paste("There was 1 log normal variate sample stored in ", 
#                   .activeDataSet, ".", sep = ""))
#             }
#             else {
#                 logger(paste("There were ", nsamples, " log normal variate samples stored in ", 
#                   .activeDataSet, ".", sep = ""))
#             }
#         }
#         tkfocus(CommanderWindow())
#     }
#     OKCancelHelp(helpSubject = "rlnorm")
#     tkgrid(tklabel(top, text = gettextRcmdr("Number of samples (columns):")), 
#         samplesEntry, sticky = "w")
#     tkgrid(tklabel(top, text = gettextRcmdr("Parameters:"), fg = "blue"), 
#         columnspan = 4, sticky = "w")
#     tkgrid(tklabel(top, text = gettextRcmdr("meanlog (mean of dist'n on log scale)")), 
#         mulogEntry, sticky = "w")
#     tkgrid(tklabel(top, text = gettextRcmdr("sdlog (std dev of dist'n on log scale)")), 
#         sigmalogEntry, sticky = "w")
#     tkgrid(tklabel(locationFrame, text = gettextRcmdr("Store values in:"), 
#         fg = "blue"), columnspan = 4, sticky = "w")
#     tkgrid(tklabel(locationFrame, text = gettextRcmdr("Active Dataset")), 
#         addtoactiveButton, sticky = "w")
#     tkgrid(tklabel(locationFrame, text = "New Dataset"), newDataButton, 
#         sticky = "w")
#     tkgrid.configure(samplesEntry, sticky = "w")
#     tkgrid.configure(mulogEntry, sticky = "w")
#     tkgrid.configure(sigmalogEntry, sticky = "w")
#     tkgrid(locationFrame, sticky = "w")
#     tkgrid(buttonsFrame, sticky = "w", columnspan = 2)
#     dialogSuffix(rows = 6, columns = 1, focus = samplesEntry)
# }
# 
# 
# 
# `logisSimulate.ipsur` <-
# function () 
# {
#     initializeDialog(title = gettextRcmdr("Simulate Logistic Variates"))
#     parameterFrame <- tkframe(top)
#     locationFrame <- tkframe(top)
#     if (!is.character(ActiveDataSet())) {
#         locVariable <- tclVar("new")
#     }
#     else {
#         locVariable <- tclVar("add")
#     }
#     addtoactiveButton <- tkradiobutton(locationFrame, variable = locVariable, 
#         value = "add")
#     newDataButton <- tkradiobutton(locationFrame, variable = locVariable, 
#         value = "new")
#     samplesVar <- tclVar("1")
#     samplesEntry <- tkentry(top, width = "6", textvariable = samplesVar)
#     locationVar <- tclVar("0")
#     locationEntry <- tkentry(top, width = "6", textvariable = locationVar)
#     scale1Var <- tclVar("1")
#     scale1Entry <- tkentry(top, width = "6", textvariable = scale1Var)
#     onOK <- function() {
#         nsamples <- round(as.numeric(tclvalue(samplesVar)))
#         location <- tclvalue(locationVar)
#         scale1 <- tclvalue(scale1Var)
#         if (is.na(nsamples) || nsamples < 1) {
#             errorCondition(recall = logisSimulate.ipsur, message = gettextRcmdr("Number of samples must be a positive integer."))
#             return()
#         }
#         if (is.na(location)) {
#             errorCondition(recall = logisSimulate.ipsur, message = gettextRcmdr("The location was not specified."))
#             return()
#         }
#         if (is.na(scale1)) {
#             errorCondition(recall = logisSimulate.ipsur, message = gettextRcmdr("The scale parameter was not specified."))
#             return()
#         }
#         closeDialog()
#         store <- tclvalue(locVariable)
#         if (store == "new") {
#             initializeDialog(title = gettextRcmdr("Simulation Dataset"))
#             dsname <- tclVar("Simset")
#             entryDsname <- tkentry(top, width = "20", textvariable = dsname)
#             newDataSS <- tclVar("100")
#             entryNewDataSS <- tkentry(top, width = "6", textvariable = newDataSS)
#             onOK <- function() {
#                 dsnameValue <- trim.blanks(tclvalue(dsname))
#                 newSS <- round(as.numeric(tclvalue(newDataSS)))
#                 closeDialog()
#                 if (dsnameValue == "") {
#                   errorCondition(recall = logisSimulate.ipsur, 
#                     message = gettextRcmdr("You must enter the name of a data set."))
#                   return()
#                 }
#                 if (!is.valid.name(dsnameValue)) {
#                   errorCondition(recall = logisSimulate.ipsur, 
#                     message = paste("\"", dsnameValue, "\" ", 
#                       gettextRcmdr("is not a valid name."), sep = ""))
#                   return()
#                 }
#                 if (is.element(dsnameValue, listDataSets())) {
#                   if ("no" == tclvalue(checkReplace(dsnameValue, 
#                     gettextRcmdr("Data set")))) {
#                     logisSimulate.ipsur()
#                     return()
#                   }
#                 }
#                 if (is.na(newSS)) {
#                   errorCondition(recall = logisSimulate.ipsur, 
#                     message = gettextRcmdr("Sample Size must be a positive integer."))
#                   return()
#                 }
#                 UpdatelogissimNumber()
#                 justDoIt(paste(dsnameValue, " = data.frame(logis.sim", 
#                   getRcmdr("logissimNumber"), "=1:", newSS, ")", 
#                   sep = ""))
#                 logger(paste(dsnameValue, "has been initialized."))
#                 for (k in getRcmdr("logissimNumber"):(nsamples + 
#                   getRcmdr("logissimNumber") - 1)) {
#                   justDoIt(paste(dsnameValue, "$logis.sim", k, 
#                     " <- rlogis(", newSS, ", location=", location, 
#                     ", scale=", scale1, ")", sep = ""))
#                 }
#                 activeDataSet(dsnameValue)
#                 putRcmdr("logissimNumber", k)
#                 if (nsamples == 1) {
#                   logger(paste("There was 1 logistic variate sample stored in ", 
#                     dsnameValue, ".", sep = ""))
#                 }
#                 else {
#                   logger(paste("There were ", nsamples, " logistic variate samples stored in ", 
#                     dsnameValue, ".", sep = ""))
#                 }
#             }
#             OKCancelHelp(helpSubject = "rlogis")
#             tkgrid(tklabel(top, text = gettextRcmdr("Enter name for data set:")), 
#                 entryDsname, sticky = "e")
#             tkgrid(tklabel(top, text = gettextRcmdr("Sample Size (rows):")), 
#                 entryNewDataSS, sticky = "e")
#             tkgrid(buttonsFrame, columnspan = "2", sticky = "w")
#             tkgrid.configure(entryDsname, sticky = "w")
#             tkgrid.configure(entryNewDataSS, sticky = "w")
#             tkfocus(CommanderWindow())
#             dialogSuffix(rows = 2, columns = 2, focus = entryDsname)
#         }
#         else {
#             if (!is.character(ActiveDataSet())) {
#                 errorCondition(recall = logisSimulate.ipsur, 
#                   message = gettextRcmdr("There is no active data set."))
#                 return()
#             }
#             .activeDataSet <- ActiveDataSet()
#             justDoIt(paste("samplesn <- dim(", .activeDataSet, 
#                 ")[1]", sep = ""))
#             UpdatelogissimNumber()
#             for (k in getRcmdr("logissimNumber"):(nsamples + 
#                 getRcmdr("logissimNumber") - 1)) {
#                 justDoIt(paste(.activeDataSet, "$logis.sim", 
#                   k, " <- rlogis(", samplesn, ", location=", 
#                   location, ", scale=", scale1, ")", sep = ""))
#             }
#             activeDataSet(.activeDataSet)
#             putRcmdr("logissimNumber", k)
#             if (nsamples == 1) {
#                 logger(paste("There was 1 logistic variate sample stored in ", 
#                   .activeDataSet, ".", sep = ""))
#             }
#             else {
#                 logger(paste("There were ", nsamples, " logistic variate samples stored in ", 
#                   .activeDataSet, ".", sep = ""))
#             }
#         }
#         tkfocus(CommanderWindow())
#     }
#     OKCancelHelp(helpSubject = "rlogis")
#     tkgrid(tklabel(top, text = gettextRcmdr("Number of samples (columns):")), 
#         samplesEntry, sticky = "w")
#     tkgrid(tklabel(top, text = gettextRcmdr("Parameters:"), fg = "blue"), 
#         columnspan = 4, sticky = "w")
#     tkgrid(tklabel(top, text = gettextRcmdr("location")), locationEntry, 
#         sticky = "w")
#     tkgrid(tklabel(top, text = gettextRcmdr("scale")), scale1Entry, 
#         sticky = "w")
#     tkgrid(tklabel(locationFrame, text = gettextRcmdr("Store values in:"), 
#         fg = "blue"), columnspan = 4, sticky = "w")
#     tkgrid(tklabel(locationFrame, text = gettextRcmdr("Active Dataset")), 
#         addtoactiveButton, sticky = "w")
#     tkgrid(tklabel(locationFrame, text = "New Dataset"), newDataButton, 
#         sticky = "w")
#     tkgrid.configure(samplesEntry, sticky = "w")
#     tkgrid.configure(locationEntry, sticky = "w")
#     tkgrid.configure(scale1Entry, sticky = "w")
#     tkgrid(locationFrame, sticky = "w")
#     tkgrid(buttonsFrame, sticky = "w", columnspan = 2)
#     dialogSuffix(rows = 6, columns = 1, focus = samplesEntry)
# }
# 
# 
# `nbinomSimulate.ipsur` <-
# function () 
# {
#     initializeDialog(title = gettextRcmdr("Simulate Negative Binomial Variates"))
#     parameterFrame <- tkframe(top)
#     locationFrame <- tkframe(top)
#     if (!is.character(ActiveDataSet())) {
#         locVariable <- tclVar("new")
#     }
#     else {
#         locVariable <- tclVar("add")
#     }
#     addtoactiveButton <- tkradiobutton(locationFrame, variable = locVariable, 
#         value = "add")
#     newDataButton <- tkradiobutton(locationFrame, variable = locVariable, 
#         value = "new")
#     samplesVar <- tclVar("1")
#     samplesEntry <- tkentry(top, width = "6", textvariable = samplesVar)
#     sizeVar <- tclVar("1")
#     sizeEntry <- tkentry(top, width = "6", textvariable = sizeVar)
#     probVar <- tclVar("0.5")
#     probEntry <- tkentry(top, width = "6", textvariable = probVar)
#     onOK <- function() {
#         nsamples <- round(as.numeric(tclvalue(samplesVar)))
#         size <- tclvalue(sizeVar)
#         prob <- tclvalue(probVar)
#         if (is.na(nsamples) || nsamples < 1) {
#             errorCondition(recall = nbinomSimulate.ipsur, message = gettextRcmdr("Number of samples must be a positive integer."))
#             return()
#         }
#         if (is.na(size)) {
#             errorCondition(recall = nbinomSimulate.ipsur, message = gettextRcmdr("The size was not specified."))
#             return()
#         }
#         if (is.na(prob)) {
#             errorCondition(recall = nbinomSimulate.ipsur, message = gettextRcmdr("The probability of success was not specified."))
#             return()
#         }
#         closeDialog()
#         store <- tclvalue(locVariable)
#         if (store == "new") {
#             initializeDialog(title = gettextRcmdr("Simulation Dataset"))
#             dsname <- tclVar("Simset")
#             entryDsname <- tkentry(top, width = "20", textvariable = dsname)
#             newDataSS <- tclVar("100")
#             entryNewDataSS <- tkentry(top, width = "6", textvariable = newDataSS)
#             onOK <- function() {
#                 dsnameValue <- trim.blanks(tclvalue(dsname))
#                 newSS <- round(as.numeric(tclvalue(newDataSS)))
#                 closeDialog()
#                 if (dsnameValue == "") {
#                   errorCondition(recall = nbinomSimulate.ipsur, 
#                     message = gettextRcmdr("You must enter the name of a data set."))
#                   return()
#                 }
#                 if (!is.valid.name(dsnameValue)) {
#                   errorCondition(recall = nbinomSimulate.ipsur, 
#                     message = paste("\"", dsnameValue, "\" ", 
#                       gettextRcmdr("is not a valid name."), sep = ""))
#                   return()
#                 }
#                 if (is.element(dsnameValue, listDataSets())) {
#                   if ("no" == tclvalue(checkReplace(dsnameValue, 
#                     gettextRcmdr("Data set")))) {
#                     nbinomSimulate.ipsur()
#                     return()
#                   }
#                 }
#                 if (is.na(newSS)) {
#                   errorCondition(recall = nbinomSimulate.ipsur, 
#                     message = gettextRcmdr("Sample Size must be a positive integer."))
#                   return()
#                 }
#                 UpdatenbinomsimNumber()
#                 justDoIt(paste(dsnameValue, " = data.frame(nbinom.sim", 
#                   getRcmdr("nbinomsimNumber"), "=1:", newSS, 
#                   ")", sep = ""))
#                 logger(paste(dsnameValue, "has been initialized."))
#                 for (k in getRcmdr("nbinomsimNumber"):(nsamples + 
#                   getRcmdr("nbinomsimNumber") - 1)) {
#                   justDoIt(paste(dsnameValue, "$nbinom.sim", 
#                     k, " <- rnbinom(", newSS, ", size=", size, 
#                     ", prob=", prob, ")", sep = ""))
#                 }
#                 activeDataSet(dsnameValue)
#                 putRcmdr("nbinomsimNumber", k)
#                 if (nsamples == 1) {
#                   logger(paste("There was 1 negative binomial variate sample stored in ", 
#                     dsnameValue, ".", sep = ""))
#                 }
#                 else {
#                   logger(paste("There were ", nsamples, " negative binomial variate samples stored in ", 
#                     dsnameValue, ".", sep = ""))
#                 }
#             }
#             OKCancelHelp(helpSubject = "rnbinom")
#             tkgrid(tklabel(top, text = gettextRcmdr("Enter name for data set:")), 
#                 entryDsname, sticky = "e")
#             tkgrid(tklabel(top, text = gettextRcmdr("Sample Size (rows):")), 
#                 entryNewDataSS, sticky = "e")
#             tkgrid(buttonsFrame, columnspan = "2", sticky = "w")
#             tkgrid.configure(entryDsname, sticky = "w")
#             tkgrid.configure(entryNewDataSS, sticky = "w")
#             tkfocus(CommanderWindow())
#             dialogSuffix(rows = 2, columns = 2, focus = entryDsname)
#         }
#         else {
#             if (!is.character(ActiveDataSet())) {
#                 errorCondition(recall = nbinomSimulate.ipsur, 
#                   message = gettextRcmdr("There is no active data set."))
#                 return()
#             }
#             .activeDataSet <- ActiveDataSet()
#             justDoIt(paste("samplesn <- dim(", .activeDataSet, 
#                 ")[1]", sep = ""))
#             UpdatenbinomsimNumber()
#             for (k in getRcmdr("nbinomsimNumber"):(nsamples + 
#                 getRcmdr("nbinomsimNumber") - 1)) {
#                 justDoIt(paste(.activeDataSet, "$nbinom.sim", 
#                   k, " <- rnbinom(", samplesn, ", size=", size, 
#                   ", prob=", prob, ")", sep = ""))
#             }
#             activeDataSet(.activeDataSet)
#             putRcmdr("nbinomsimNumber", k)
#             if (nsamples == 1) {
#                 logger(paste("There was 1 negative binomial variate sample stored in ", 
#                   .activeDataSet, ".", sep = ""))
#             }
#             else {
#                 logger(paste("There were ", nsamples, " negative binomial variate samples stored in ", 
#                   .activeDataSet, ".", sep = ""))
#             }
#         }
#         tkfocus(CommanderWindow())
#     }
#     OKCancelHelp(helpSubject = "rnbinom")
#     tkgrid(tklabel(top, text = gettextRcmdr("Number of samples (columns):")), 
#         samplesEntry, sticky = "w")
#     tkgrid(tklabel(top, text = gettextRcmdr("Parameters:"), fg = "blue"), 
#         columnspan = 4, sticky = "w")
#     tkgrid(tklabel(top, text = gettextRcmdr("size (target number of successes)")), 
#         sizeEntry, sticky = "w")
#     tkgrid(tklabel(top, text = gettextRcmdr("prob (of success in each trial)")), 
#         probEntry, sticky = "w")
#     tkgrid(tklabel(locationFrame, text = gettextRcmdr("Store values in:"), 
#         fg = "blue"), columnspan = 4, sticky = "w")
#     tkgrid(tklabel(locationFrame, text = gettextRcmdr("Active Dataset")), 
#         addtoactiveButton, sticky = "w")
#     tkgrid(tklabel(locationFrame, text = "New Dataset"), newDataButton, 
#         sticky = "w")
#     tkgrid.configure(samplesEntry, sticky = "w")
#     tkgrid.configure(sizeEntry, sticky = "w")
#     tkgrid.configure(probEntry, sticky = "w")
#     tkgrid(locationFrame, sticky = "w")
#     tkgrid(buttonsFrame, sticky = "w", columnspan = 2)
#     dialogSuffix(rows = 6, columns = 1, focus = samplesEntry)
# }
# 
# 
# `normalSimulate.ipsur` <-
# function () 
# {
#     initializeDialog(title = gettextRcmdr("Simulate Normal Variates"))
#     parameterFrame <- tkframe(top)
#     locationFrame <- tkframe(top)
#     if (!is.character(ActiveDataSet())) {
#         locVariable <- tclVar("new")
#     }
#     else {
#         locVariable <- tclVar("add")
#     }
#     addtoactiveButton <- tkradiobutton(locationFrame, variable = locVariable, 
#         value = "add")
#     newDataButton <- tkradiobutton(locationFrame, variable = locVariable, 
#         value = "new")
#     samplesVar <- tclVar("1")
#     samplesEntry <- tkentry(top, width = "6", textvariable = samplesVar)
#     muVar <- tclVar("0")
#     muEntry <- tkentry(top, width = "6", textvariable = muVar)
#     sigmaVar <- tclVar("1")
#     sigmaEntry <- tkentry(top, width = "6", textvariable = sigmaVar)
#     onOK <- function() {
#         nsamples <- round(as.numeric(tclvalue(samplesVar)))
#         mu <- tclvalue(muVar)
#         sigma <- tclvalue(sigmaVar)
#         if (is.na(nsamples)) {
#             errorCondition(recall = normalSimulate.ipsur, message = gettextRcmdr("Number of samples must be a positive integer."))
#             return()
#         }
#         if (is.na(mu)) {
#             errorCondition(recall = normalSimulate.ipsur, message = gettextRcmdr("The mean was not specified."))
#             return()
#         }
#         if (is.na(sigma)) {
#             errorCondition(recall = normalSimulate.ipsur, message = gettextRcmdr("The standard deviation was not specified."))
#             return()
#         }
#         closeDialog()
#         store <- tclvalue(locVariable)
#         if (store == "new") {
#             initializeDialog(title = gettextRcmdr("Simulation Dataset"))
#             dsname <- tclVar("Simset")
#             entryDsname <- tkentry(top, width = "20", textvariable = dsname)
#             newDataSS <- tclVar("100")
#             entryNewDataSS <- tkentry(top, width = "6", textvariable = newDataSS)
#             onOK <- function() {
#                 dsnameValue <- trim.blanks(tclvalue(dsname))
#                 newSS <- round(as.numeric(tclvalue(newDataSS)))
#                 closeDialog()
#                 if (dsnameValue == "") {
#                   errorCondition(recall = normalSimulate.ipsur, 
#                     message = gettextRcmdr("You must enter the name of a data set."))
#                   return()
#                 }
#                 if (!is.valid.name(dsnameValue)) {
#                   errorCondition(recall = normalSimulate.ipsur, 
#                     message = paste("\"", dsnameValue, "\" ", 
#                       gettextRcmdr("is not a valid name."), sep = ""))
#                   return()
#                 }
#                 if (is.element(dsnameValue, listDataSets())) {
#                   if ("no" == tclvalue(checkReplace(dsnameValue, 
#                     gettextRcmdr("Data set")))) {
#                     normalSimulate.ipsur()
#                     return()
#                   }
#                 }
#                 if (is.na(newSS)) {
#                   errorCondition(recall = normalSimulate.ipsur, 
#                     message = gettextRcmdr("Sample Size must be a positive integer."))
#                   return()
#                 }
#                 UpdatenormsimNumber()
#                 justDoIt(paste(dsnameValue, " = data.frame(norm.sim", 
#                   getRcmdr("normsimNumber"), "=1:", newSS, ")", 
#                   sep = ""))
#                 logger(paste(dsnameValue, "has been initialized."))
#                 for (k in getRcmdr("normsimNumber"):(nsamples + 
#                   getRcmdr("normsimNumber") - 1)) {
#                   justDoIt(paste(dsnameValue, "$norm.sim", k, 
#                     " <- rnorm(", newSS, ", mean=", mu, ", sd=", 
#                     sigma, ")", sep = ""))
#                 }
#                 activeDataSet(dsnameValue)
#                 putRcmdr("normsimNumber", k)
#                 if (nsamples == 1) {
#                   logger(paste("There was 1 normal variate sample stored in ", 
#                     dsnameValue, ".", sep = ""))
#                 }
#                 else {
#                   logger(paste("There were ", nsamples, " normal variate samples stored in ", 
#                     dsnameValue, ".", sep = ""))
#                 }
#             }
#             OKCancelHelp(helpSubject = "rnorm")
#             tkgrid(tklabel(top, text = gettextRcmdr("Enter name for data set:")), 
#                 entryDsname, sticky = "e")
#             tkgrid(tklabel(top, text = gettextRcmdr("Sample Size (rows):")), 
#                 entryNewDataSS, sticky = "e")
#             tkgrid(buttonsFrame, columnspan = "2", sticky = "w")
#             tkgrid.configure(entryDsname, sticky = "w")
#             tkgrid.configure(entryNewDataSS, sticky = "w")
#             tkfocus(CommanderWindow())
#             dialogSuffix(rows = 2, columns = 2, focus = entryDsname)
#         }
#         else {
#             if (!is.character(ActiveDataSet())) {
#                 errorCondition(recall = normalSimulate.ipsur, 
#                   message = gettextRcmdr("There is no active data set."))
#                 return()
#             }
#             .activeDataSet <- ActiveDataSet()
#             justDoIt(paste("samplesn <- dim(", .activeDataSet, 
#                 ")[1]", sep = ""))
#             UpdatenormsimNumber()
#             for (k in getRcmdr("normsimNumber"):(nsamples + getRcmdr("normsimNumber") - 
#                 1)) {
#                 justDoIt(paste(.activeDataSet, "$norm.sim", k, 
#                   " <- rnorm(", samplesn, ", mean=", mu, ", sd=", 
#                   sigma, ")", sep = ""))
#             }
#             activeDataSet(.activeDataSet)
#             putRcmdr("normsimNumber", k)
#             if (nsamples == 1) {
#                 logger(paste("There was 1 normal variate sample stored in ", 
#                   .activeDataSet, ".", sep = ""))
#             }
#             else {
#                 logger(paste("There were ", nsamples, " normal variate samples stored in ", 
#                   .activeDataSet, ".", sep = ""))
#             }
#         }
#         tkfocus(CommanderWindow())
#     }
#     OKCancelHelp(helpSubject = "rnorm")
#     tkgrid(tklabel(top, text = gettextRcmdr("Number of samples (columns):")), 
#         samplesEntry, sticky = "w")
#     tkgrid(tklabel(top, text = gettextRcmdr("Parameters:"), fg = "blue"), 
#         columnspan = 4, sticky = "w")
#     tkgrid(tklabel(top, text = gettextRcmdr("mean (mu)")), muEntry, 
#         sticky = "w")
#     tkgrid(tklabel(top, text = gettextRcmdr("sd (sigma)")), sigmaEntry, 
#         sticky = "w")
#     tkgrid(tklabel(locationFrame, text = gettextRcmdr("Store values in:"), 
#         fg = "blue"), columnspan = 4, sticky = "w")
#     tkgrid(tklabel(locationFrame, text = gettextRcmdr("Active Dataset")), 
#         addtoactiveButton, sticky = "w")
#     tkgrid(tklabel(locationFrame, text = "New Dataset"), newDataButton, 
#         sticky = "w")
#     tkgrid.configure(samplesEntry, sticky = "w")
#     tkgrid.configure(muEntry, sticky = "w")
#     tkgrid.configure(sigmaEntry, sticky = "w")
#     tkgrid(locationFrame, sticky = "w")
#     tkgrid(buttonsFrame, sticky = "w", columnspan = 2)
#     dialogSuffix(rows = 6, columns = 1, focus = samplesEntry)
# }
# 
# 
# `poisSimulate.ipsur` <-
# function () 
# {
#     initializeDialog(title = gettextRcmdr("Simulate Poisson Variates"))
#     parameterFrame <- tkframe(top)
#     locationFrame <- tkframe(top)
#     if (!is.character(ActiveDataSet())) {
#         locVariable <- tclVar("new")
#     }
#     else {
#         locVariable <- tclVar("add")
#     }
#     addtoactiveButton <- tkradiobutton(locationFrame, variable = locVariable, 
#         value = "add")
#     newDataButton <- tkradiobutton(locationFrame, variable = locVariable, 
#         value = "new")
#     samplesVar <- tclVar("1")
#     samplesEntry <- tkentry(top, width = "6", textvariable = samplesVar)
#     lambdaVar <- tclVar("1")
#     lambdaEntry <- tkentry(top, width = "6", textvariable = lambdaVar)
#     onOK <- function() {
#         nsamples <- round(as.numeric(tclvalue(samplesVar)))
#         lambda <- tclvalue(lambdaVar)
#         if (is.na(nsamples)) {
#             errorCondition(recall = poisSimulate.ipsur, message = gettextRcmdr("Number of samples must be a positive integer."))
#             return()
#         }
#         if (is.na(lambda)) {
#             errorCondition(recall = poisSimulate.ipsur, message = gettextRcmdr("The mean parameter was not specified."))
#             return()
#         }
#         closeDialog()
#         store <- tclvalue(locVariable)
#         if (store == "new") {
#             initializeDialog(title = gettextRcmdr("Simulation Dataset"))
#             dsname <- tclVar("Simset")
#             entryDsname <- tkentry(top, width = "20", textvariable = dsname)
#             newDataSS <- tclVar("100")
#             entryNewDataSS <- tkentry(top, width = "6", textvariable = newDataSS)
#             onOK <- function() {
#                 dsnameValue <- trim.blanks(tclvalue(dsname))
#                 newSS <- round(as.numeric(tclvalue(newDataSS)))
#                 closeDialog()
#                 if (dsnameValue == "") {
#                   errorCondition(recall = poisSimulate.ipsur, 
#                     message = gettextRcmdr("You must enter the name of a data set."))
#                   return()
#                 }
#                 if (!is.valid.name(dsnameValue)) {
#                   errorCondition(recall = poisSimulate.ipsur, 
#                     message = paste("\"", dsnameValue, "\" ", 
#                       gettextRcmdr("is not a valid name."), sep = ""))
#                   return()
#                 }
#                 if (is.element(dsnameValue, listDataSets())) {
#                   if ("no" == tclvalue(checkReplace(dsnameValue, 
#                     gettextRcmdr("Data set")))) {
#                     poisSimulate.ipsur()
#                     return()
#                   }
#                 }
#                 if (is.na(newSS)) {
#                   errorCondition(recall = poisSimulate.ipsur, 
#                     message = gettextRcmdr("Sample Size must be a positive integer."))
#                   return()
#                 }
#                 UpdatepoissimNumber()
#                 justDoIt(paste(dsnameValue, " = data.frame(pois.sim", 
#                   getRcmdr("poissimNumber"), "=1:", newSS, ")", 
#                   sep = ""))
#                 logger(paste(dsnameValue, "has been initialized."))
#                 for (k in getRcmdr("poissimNumber"):(nsamples + 
#                   getRcmdr("poissimNumber") - 1)) {
#                   justDoIt(paste(dsnameValue, "$pois.sim", k, 
#                     " <- rpois(", newSS, ", lambda=", lambda, 
#                     ")", sep = ""))
#                 }
#                 activeDataSet(dsnameValue)
#                 putRcmdr("poissimNumber", k)
#                 if (nsamples == 1) {
#                   logger(paste("There was 1 Poisson variate sample stored in ", 
#                     dsnameValue, ".", sep = ""))
#                 }
#                 else {
#                   logger(paste("There were ", nsamples, " Poisson variate samples stored in ", 
#                     dsnameValue, ".", sep = ""))
#                 }
#             }
#             OKCancelHelp(helpSubject = "rpois")
#             tkgrid(tklabel(top, text = gettextRcmdr("Enter name for data set:")), 
#                 entryDsname, sticky = "e")
#             tkgrid(tklabel(top, text = gettextRcmdr("Sample Size (rows):")), 
#                 entryNewDataSS, sticky = "e")
#             tkgrid(buttonsFrame, columnspan = "2", sticky = "w")
#             tkgrid.configure(entryDsname, sticky = "w")
#             tkgrid.configure(entryNewDataSS, sticky = "w")
#             tkfocus(CommanderWindow())
#             dialogSuffix(rows = 2, columns = 2, focus = entryDsname)
#         }
#         else {
#             if (!is.character(ActiveDataSet())) {
#                 errorCondition(recall = poisSimulate.ipsur, message = gettextRcmdr("There is no active data set."))
#                 return()
#             }
#             .activeDataSet <- ActiveDataSet()
#             justDoIt(paste("samplesn <- dim(", .activeDataSet, 
#                 ")[1]", sep = ""))
#             UpdatepoissimNumber()
#             for (k in getRcmdr("poissimNumber"):(nsamples + getRcmdr("poissimNumber") - 
#                 1)) {
#                 justDoIt(paste(.activeDataSet, "$pois.sim", k, 
#                   " <- rpois(", samplesn, ", lambda=", lambda, 
#                   ")", sep = ""))
#             }
#             activeDataSet(.activeDataSet)
#             putRcmdr("poissimNumber", k)
#             if (nsamples == 1) {
#                 logger(paste("There was 1 Poisson variate sample stored in ", 
#                   .activeDataSet, ".", sep = ""))
#             }
#             else {
#                 logger(paste("There were ", nsamples, " Poisson variate samples stored in ", 
#                   .activeDataSet, ".", sep = ""))
#             }
#         }
#         tkfocus(CommanderWindow())
#     }
#     OKCancelHelp(helpSubject = "rpois")
#     tkgrid(tklabel(top, text = gettextRcmdr("Number of samples (columns):")), 
#         samplesEntry, sticky = "w")
#     tkgrid(tklabel(top, text = gettextRcmdr("Parameters:"), fg = "blue"), 
#         columnspan = 4, sticky = "w")
#     tkgrid(tklabel(top, text = gettextRcmdr("lambda (mean)")), 
#         lambdaEntry, sticky = "w")
#     tkgrid(tklabel(locationFrame, text = gettextRcmdr("Store values in:"), 
#         fg = "blue"), columnspan = 4, sticky = "w")
#     tkgrid(tklabel(locationFrame, text = gettextRcmdr("Active Dataset")), 
#         addtoactiveButton, sticky = "w")
#     tkgrid(tklabel(locationFrame, text = "New Dataset"), newDataButton, 
#         sticky = "w")
#     tkgrid.configure(samplesEntry, sticky = "w")
#     tkgrid.configure(lambdaEntry, sticky = "w")
#     tkgrid(locationFrame, sticky = "w")
#     tkgrid(buttonsFrame, sticky = "w", columnspan = 2)
#     dialogSuffix(rows = 6, columns = 1, focus = samplesEntry)
# }
# 
# 
# `RcmdrEnv` <-
# function () 
# {
#     pos <- match("RcmdrEnv", search())
#     if (is.na(pos)) {
#         RcmdrEnv <- list()
#         rm(RcmdrEnv)
#         pos <- match("RcmdrEnv", search())
#     }
#     return(pos.to.env(pos))
# }
# 
# 
# `tSimulate.ipsur` <-
# function () 
# {
#     initializeDialog(title = gettextRcmdr("Simulate t Variates"))
#     parameterFrame <- tkframe(top)
#     locationFrame <- tkframe(top)
#     if (!is.character(ActiveDataSet())) {
#         locVariable <- tclVar("new")
#     }
#     else {
#         locVariable <- tclVar("add")
#     }
#     addtoactiveButton <- tkradiobutton(locationFrame, variable = locVariable, 
#         value = "add")
#     newDataButton <- tkradiobutton(locationFrame, variable = locVariable, 
#         value = "new")
#     samplesVar <- tclVar("1")
#     samplesEntry <- tkentry(top, width = "6", textvariable = samplesVar)
#     dfVar <- tclVar("1")
#     dfEntry <- tkentry(top, width = "6", textvariable = dfVar)
#     ncpVar <- tclVar("0")
#     ncpEntry <- tkentry(top, width = "6", textvariable = ncpVar)
#     onOK <- function() {
#         nsamples <- round(as.numeric(tclvalue(samplesVar)))
#         df <- tclvalue(dfVar)
#         ncp <- tclvalue(ncpVar)
#         if (is.na(nsamples)) {
#             errorCondition(recall = tSimulate.ipsur, message = gettextRcmdr("Number of samples must be a positive integer."))
#             return()
#         }
#         if (is.na(df)) {
#             errorCondition(recall = tSimulate.ipsur, message = gettextRcmdr("The degrees of freedom were not specified."))
#             return()
#         }
#         if (is.na(ncp)) {
#             errorCondition(recall = tSimulate.ipsur, message = gettextRcmdr("The noncentrality parameter was not specified."))
#             return()
#         }
#         closeDialog()
#         store <- tclvalue(locVariable)
#         if (store == "new") {
#             initializeDialog(title = gettextRcmdr("Simulation Dataset"))
#             dsname <- tclVar("Simset")
#             entryDsname <- tkentry(top, width = "20", textvariable = dsname)
#             newDataSS <- tclVar("100")
#             entryNewDataSS <- tkentry(top, width = "6", textvariable = newDataSS)
#             onOK <- function() {
#                 dsnameValue <- trim.blanks(tclvalue(dsname))
#                 newSS <- round(as.numeric(tclvalue(newDataSS)))
#                 closeDialog()
#                 if (dsnameValue == "") {
#                   errorCondition(recall = tSimulate.ipsur, message = gettextRcmdr("You must enter the name of a data set."))
#                   return()
#                 }
#                 if (!is.valid.name(dsnameValue)) {
#                   errorCondition(recall = tSimulate.ipsur, message = paste("\"", 
#                     dsnameValue, "\" ", gettextRcmdr("is not a valid name."), 
#                     sep = ""))
#                   return()
#                 }
#                 if (is.element(dsnameValue, listDataSets())) {
#                   if ("no" == tclvalue(checkReplace(dsnameValue, 
#                     gettextRcmdr("Data set")))) {
#                     tSimulate.ipsur()
#                     return()
#                   }
#                 }
#                 if (is.na(newSS)) {
#                   errorCondition(recall = tSimulate.ipsur, message = gettextRcmdr("Sample Size must be a positive integer."))
#                   return()
#                 }
#                 UpdatetsimNumber()
#                 justDoIt(paste(dsnameValue, " = data.frame(t.sim", 
#                   getRcmdr("tsimNumber"), "=1:", newSS, ")", 
#                   sep = ""))
#                 logger(paste(dsnameValue, "has been initialized."))
#                 for (k in getRcmdr("tsimNumber"):(nsamples + 
#                   getRcmdr("tsimNumber") - 1)) {
#                   justDoIt(paste(dsnameValue, "$t.sim", k, " <- rt(", 
#                     newSS, ", df=", df, ", ncp=", ncp, ")", sep = ""))
#                 }
#                 activeDataSet(dsnameValue)
#                 putRcmdr("tsimNumber", k)
#                 if (nsamples == 1) {
#                   logger(paste("There was 1 Student's t variate sample stored in ", 
#                     dsnameValue, ".", sep = ""))
#                 }
#                 else {
#                   logger(paste("There were ", nsamples, " Student's t variate samples stored in ", 
#                     dsnameValue, ".", sep = ""))
#                 }
#             }
#             OKCancelHelp(helpSubject = "rt")
#             tkgrid(tklabel(top, text = gettextRcmdr("Enter name for data set:")), 
#                 entryDsname, sticky = "e")
#             tkgrid(tklabel(top, text = gettextRcmdr("Sample Size (rows):")), 
#                 entryNewDataSS, sticky = "e")
#             tkgrid(buttonsFrame, columnspan = "2", sticky = "w")
#             tkgrid.configure(entryDsname, sticky = "w")
#             tkgrid.configure(entryNewDataSS, sticky = "w")
#             tkfocus(CommanderWindow())
#             dialogSuffix(rows = 2, columns = 2, focus = entryDsname)
#         }
#         else {
#             if (!is.character(ActiveDataSet())) {
#                 errorCondition(recall = tSimulate.ipsur, message = gettextRcmdr("There is no active data set."))
#                 return()
#             }
#             .activeDataSet <- ActiveDataSet()
#             justDoIt(paste("samplesn <- dim(", .activeDataSet, 
#                 ")[1]", sep = ""))
#             UpdatetsimNumber()
#             for (k in getRcmdr("tsimNumber"):(nsamples + getRcmdr("tsimNumber") - 
#                 1)) {
#                 justDoIt(paste(.activeDataSet, "$t.sim", k, " <- rt(", 
#                   samplesn, ", df=", df, ", ncp=", ncp, ")", 
#                   sep = ""))
#             }
#             activeDataSet(.activeDataSet)
#             putRcmdr("tsimNumber", k)
#             if (nsamples == 1) {
#                 logger(paste("There was 1 Student's t variate sample stored in ", 
#                   .activeDataSet, ".", sep = ""))
#             }
#             else {
#                 logger(paste("There were ", nsamples, " Student's t variate samples stored in ", 
#                   .activeDataSet, ".", sep = ""))
#             }
#         }
#         tkfocus(CommanderWindow())
#     }
#     OKCancelHelp(helpSubject = "rt")
#     tkgrid(tklabel(top, text = gettextRcmdr("Number of samples (columns):")), 
#         samplesEntry, sticky = "w")
#     tkgrid(tklabel(top, text = gettextRcmdr("Parameters:"), fg = "blue"), 
#         columnspan = 4, sticky = "w")
#     tkgrid(tklabel(top, text = gettextRcmdr("df (degrees of freedom)")), 
#         dfEntry, sticky = "w")
#     tkgrid(tklabel(top, text = gettextRcmdr("ncp (noncentrality parameter) ")), 
#         ncpEntry, sticky = "w")
#     tkgrid(tklabel(locationFrame, text = gettextRcmdr("Store values in:"), 
#         fg = "blue"), columnspan = 4, sticky = "w")
#     tkgrid(tklabel(locationFrame, text = gettextRcmdr("Active Dataset")), 
#         addtoactiveButton, sticky = "w")
#     tkgrid(tklabel(locationFrame, text = "New Dataset"), newDataButton, 
#         sticky = "w")
#     tkgrid.configure(samplesEntry, sticky = "w")
#     tkgrid.configure(dfEntry, sticky = "w")
#     tkgrid.configure(ncpEntry, sticky = "w")
#     tkgrid(locationFrame, sticky = "w")
#     tkgrid(buttonsFrame, sticky = "w", columnspan = 2)
#     dialogSuffix(rows = 6, columns = 1, focus = samplesEntry)
# }
# 
# 
# `unifSimulate.ipsur` <-
# function () 
# {
#     initializeDialog(title = gettextRcmdr("Simulate Uniform Variates"))
#     parameterFrame <- tkframe(top)
#     locationFrame <- tkframe(top)
#     if (!is.character(ActiveDataSet())) {
#         locVariable <- tclVar("new")
#     }
#     else {
#         locVariable <- tclVar("add")
#     }
#     addtoactiveButton <- tkradiobutton(locationFrame, variable = locVariable, 
#         value = "add")
#     newDataButton <- tkradiobutton(locationFrame, variable = locVariable, 
#         value = "new")
#     samplesVar <- tclVar("1")
#     samplesEntry <- tkentry(top, width = "6", textvariable = samplesVar)
#     min1Var <- tclVar("0")
#     min1Entry <- tkentry(top, width = "6", textvariable = min1Var)
#     max1Var <- tclVar("1")
#     max1Entry <- tkentry(top, width = "6", textvariable = max1Var)
#     onOK <- function() {
#         nsamples <- round(as.numeric(tclvalue(samplesVar)))
#         min1 <- tclvalue(min1Var)
#         max1 <- tclvalue(max1Var)
#         if (is.na(nsamples)) {
#             errorCondition(recall = unifSimulate.ipsur, message = gettextRcmdr("Number of samples must be a positive integer."))
#             return()
#         }
#         if (is.na(min1)) {
#             errorCondition(recall = unifSimulate.ipsur, message = gettextRcmdr("The lower limit(min) was not specified."))
#             return()
#         }
#         if (is.na(max1)) {
#             errorCondition(recall = unifSimulate.ipsur, message = gettextRcmdr("The upper limit(max) was not specified."))
#             return()
#         }
#         closeDialog()
#         store <- tclvalue(locVariable)
#         if (store == "new") {
#             initializeDialog(title = gettextRcmdr("Simulation Dataset"))
#             dsname <- tclVar("Simset")
#             entryDsname <- tkentry(top, width = "20", textvariable = dsname)
#             newDataSS <- tclVar("100")
#             entryNewDataSS <- tkentry(top, width = "6", textvariable = newDataSS)
#             onOK <- function() {
#                 dsnameValue <- trim.blanks(tclvalue(dsname))
#                 newSS <- round(as.numeric(tclvalue(newDataSS)))
#                 closeDialog()
#                 if (dsnameValue == "") {
#                   errorCondition(recall = unifSimulate.ipsur, 
#                     message = gettextRcmdr("You must enter the name of a data set."))
#                   return()
#                 }
#                 if (!is.valid.name(dsnameValue)) {
#                   errorCondition(recall = unifSimulate.ipsur, 
#                     message = paste("\"", dsnameValue, "\" ", 
#                       gettextRcmdr("is not a valid name."), sep = ""))
#                   return()
#                 }
#                 if (is.element(dsnameValue, listDataSets())) {
#                   if ("no" == tclvalue(checkReplace(dsnameValue, 
#                     gettextRcmdr("Data set")))) {
#                     unifSimulate.ipsur()
#                     return()
#                   }
#                 }
#                 if (is.na(newSS)) {
#                   errorCondition(recall = unifSimulate.ipsur, 
#                     message = gettextRcmdr("Sample Size must be a positive integer."))
#                   return()
#                 }
#                 UpdateunifsimNumber()
#                 justDoIt(paste(dsnameValue, " = data.frame(unif.sim", 
#                   getRcmdr("unifsimNumber"), "=1:", newSS, ")", 
#                   sep = ""))
#                 logger(paste(dsnameValue, "has been initialized."))
#                 for (k in getRcmdr("unifsimNumber"):(nsamples + 
#                   getRcmdr("unifsimNumber") - 1)) {
#                   justDoIt(paste(dsnameValue, "$unif.sim", k, 
#                     " <- runif(", newSS, ", min=", min1, ", max=", 
#                     max1, ")", sep = ""))
#                 }
#                 activeDataSet(dsnameValue)
#                 putRcmdr("unifsimNumber", k)
#                 if (nsamples == 1) {
#                   logger(paste("There was 1 uniform variate sample stored in ", 
#                     dsnameValue, ".", sep = ""))
#                 }
#                 else {
#                   logger(paste("There were ", nsamples, " uniform variate samples stored in ", 
#                     dsnameValue, ".", sep = ""))
#                 }
#             }
#             OKCancelHelp(helpSubject = "runif")
#             tkgrid(tklabel(top, text = gettextRcmdr("Enter name for data set:")), 
#                 entryDsname, sticky = "e")
#             tkgrid(tklabel(top, text = gettextRcmdr("Sample Size (rows):")), 
#                 entryNewDataSS, sticky = "e")
#             tkgrid(buttonsFrame, columnspan = "2", sticky = "w")
#             tkgrid.configure(entryDsname, sticky = "w")
#             tkgrid.configure(entryNewDataSS, sticky = "w")
#             tkfocus(CommanderWindow())
#             dialogSuffix(rows = 2, columns = 2, focus = entryDsname)
#         }
#         else {
#             if (!is.character(ActiveDataSet())) {
#                 errorCondition(recall = unifSimulate.ipsur, message = gettextRcmdr("There is no active data set."))
#                 return()
#             }
#             .activeDataSet <- ActiveDataSet()
#             justDoIt(paste("samplesn <- dim(", .activeDataSet, 
#                 ")[1]", sep = ""))
#             UpdateunifsimNumber()
#             for (k in getRcmdr("unifsimNumber"):(nsamples + getRcmdr("unifsimNumber") - 
#                 1)) {
#                 justDoIt(paste(.activeDataSet, "$unif.sim", k, 
#                   " <- runif(", samplesn, ", min=", min1, ", max=", 
#                   max1, ")", sep = ""))
#             }
#             activeDataSet(.activeDataSet)
#             putRcmdr("unifsimNumber", k)
#             if (nsamples == 1) {
#                 logger(paste("There was 1 uniform variate sample stored in ", 
#                   .activeDataSet, ".", sep = ""))
#             }
#             else {
#                 logger(paste("There were ", nsamples, " uniform variate samples stored in ", 
#                   .activeDataSet, ".", sep = ""))
#             }
#         }
#         tkfocus(CommanderWindow())
#     }
#     OKCancelHelp(helpSubject = "runif")
#     tkgrid(tklabel(top, text = gettextRcmdr("Number of samples (columns):")), 
#         samplesEntry, sticky = "w")
#     tkgrid(tklabel(top, text = gettextRcmdr("Parameters:"), fg = "blue"), 
#         columnspan = 4, sticky = "w")
#     tkgrid(tklabel(top, text = gettextRcmdr("min (lower limit of the distribution)")), 
#         min1Entry, sticky = "w")
#     tkgrid(tklabel(top, text = gettextRcmdr("max (upper limit of the distribution)")), 
#         max1Entry, sticky = "w")
#     tkgrid(tklabel(locationFrame, text = gettextRcmdr("Store values in:"), 
#         fg = "blue"), columnspan = 4, sticky = "w")
#     tkgrid(tklabel(locationFrame, text = gettextRcmdr("Active Dataset")), 
#         addtoactiveButton, sticky = "w")
#     tkgrid(tklabel(locationFrame, text = "New Dataset"), newDataButton, 
#         sticky = "w")
#     tkgrid.configure(samplesEntry, sticky = "w")
#     tkgrid.configure(min1Entry, sticky = "w")
#     tkgrid.configure(max1Entry, sticky = "w")
#     tkgrid(locationFrame, sticky = "w")
#     tkgrid(buttonsFrame, sticky = "w", columnspan = 2)
#     dialogSuffix(rows = 6, columns = 1, focus = samplesEntry)
# }
# 
# 
# `UpdatebetasimNumber` <-
# function (increment = 1) 
# {
#     betasimNumber <- getRcmdr("betasimNumber")
#     putRcmdr("betasimNumber", betasimNumber + increment)
# }
# `UpdatebinomsimNumber` <-
# function (increment = 1) 
# {
#     fsimNumber <- getRcmdr("binomsimNumber")
#     putRcmdr("binomsimNumber", binomsimNumber + increment)
# }
# `UpdatecauchysimNumber` <-
# function (increment = 1) 
# {
#     cauchysimNumber <- getRcmdr("cauchysimNumber")
#     putRcmdr("cauchysimNumber", cauchysimNumber + increment)
# }
# `UpdatechisqsimNumber` <-
# function (increment = 1) 
# {
#     chisqsimNumber <- getRcmdr("chisqsimNumber")
#     putRcmdr("chisqsimNumber", chisqsimNumber + increment)
# }
# `UpdatedisunifsimNumber` <-
# function (increment = 1) 
# {
#     disunifsimNumber <- getRcmdr("disunifsimNumber")
#     putRcmdr("disunifsimNumber", disunifsimNumber + increment)
# }
# `UpdateexpsimNumber` <-
# function (increment = 1) 
# {
#     expsimNumber <- getRcmdr("expsimNumber")
#     putRcmdr("expsimNumber", expsimNumber + increment)
# }
# `UpdatefsimNumber` <-
# function (increment = 1) 
# {
#     fsimNumber <- getRcmdr("fsimNumber")
#     putRcmdr("fsimNumber", fsimNumber + increment)
# }
# `UpdategammasimNumber` <-
# function (increment = 1) 
# {
#     gammasimNumber <- getRcmdr("gammasimNumber")
#     putRcmdr("gammasimNumber", gammasimNumber + increment)
# }
# `UpdategeomsimNumber` <-
# function (increment = 1) 
# {
#     geomsimNumber <- getRcmdr("geomsimNumber")
#     putRcmdr("geomsimNumber", geomsimNumber + increment)
# }
# `UpdatehypersimNumber` <-
# function (increment = 1) 
# {
#     hypersimNumber <- getRcmdr("hypersimNumber")
#     putRcmdr("hypersimNumber", expsimNumber + increment)
# }
# `UpdatelnormsimNumber` <-
# function (increment = 1) 
# {
#     lnormsimNumber <- getRcmdr("lnormsimNumber")
#     putRcmdr("lnormsimNumber", lnormsimNumber + increment)
# }
# `UpdatelogissimNumber` <-
# function (increment = 1) 
# {
#     logissimNumber <- getRcmdr("logissimNumber")
#     putRcmdr("logissimNumber", logissimNumber + increment)
# }
# `UpdatenbinomsimNumber` <-
# function (increment = 1) 
# {
#     nbinomsimNumber <- getRcmdr("nbinomsimNumber")
#     putRcmdr("nbinomsimNumber", nbinomsimNumber + increment)
# }
# `UpdatenormsimNumber` <-
# function (increment = 1) 
# {
#     normsimNumber <- getRcmdr("normsimNumber")
#     putRcmdr("normsimNumber", normsimNumber + increment)
# }
# `UpdatepoissimNumber` <-
# function (increment = 1) 
# {
#     poissimNumber <- getRcmdr("poissimNumber")
#     putRcmdr("poissimNumber", poissimNumber + increment)
# }
# `UpdatetsimNumber` <-
# function (increment = 1) 
# {
#     tsimNumber <- getRcmdr("tsimNumber")
#     putRcmdr("tsimNumber", tsimNumber + increment)
# }
# `UpdateunifsimNumber` <-
# function (increment = 1) 
# {
#     unifsimNumber <- getRcmdr("unifsimNumber")
#     putRcmdr("unifsimNumber", unifsimNumber + increment)
# }
# `UpdateweibullsimNumber` <-
# function (increment = 1) 
# {
#     weibullsimNumber <- getRcmdr("weibullsimNumber")
#     putRcmdr("weibullsimNumber", weibullsimNumber + increment)
# }
# 
# 
# `weibullSimulate.ipsur` <-
# function () 
# {
#     initializeDialog(title = gettextRcmdr("Simulate Weibull Variates"))
#     parameterFrame <- tkframe(top)
#     locationFrame <- tkframe(top)
#     if (!is.character(ActiveDataSet())) {
#         locVariable <- tclVar("new")
#     }
#     else {
#         locVariable <- tclVar("add")
#     }
#     addtoactiveButton <- tkradiobutton(locationFrame, variable = locVariable, 
#         value = "add")
#     newDataButton <- tkradiobutton(locationFrame, variable = locVariable, 
#         value = "new")
#     samplesVar <- tclVar("1")
#     samplesEntry <- tkentry(top, width = "6", textvariable = samplesVar)
#     shapeVar <- tclVar("1")
#     shapeEntry <- tkentry(top, width = "6", textvariable = shapeVar)
#     scale1Var <- tclVar("1")
#     scale1Entry <- tkentry(top, width = "6", textvariable = scale1Var)
#     onOK <- function() {
#         nsamples <- round(as.numeric(tclvalue(samplesVar)))
#         shape <- tclvalue(shapeVar)
#         scale1 <- tclvalue(scale1Var)
#         if (is.na(nsamples)) {
#             errorCondition(recall = weibullSimulate.ipsur, message = gettextRcmdr("Number of samples must be a positive integer."))
#             return()
#         }
#         if (is.na(shape)) {
#             errorCondition(recall = weibullSimulate.ipsur, message = gettextRcmdr("The shape parameter was not specified."))
#             return()
#         }
#         if (is.na(scale1)) {
#             errorCondition(recall = weibullSimulate.ipsur, message = gettextRcmdr("The scale parameter was not specified."))
#             return()
#         }
#         closeDialog()
#         store <- tclvalue(locVariable)
#         if (store == "new") {
#             initializeDialog(title = gettextRcmdr("Simulation Dataset"))
#             dsname <- tclVar("Simset")
#             entryDsname <- tkentry(top, width = "20", textvariable = dsname)
#             newDataSS <- tclVar("100")
#             entryNewDataSS <- tkentry(top, width = "6", textvariable = newDataSS)
#             onOK <- function() {
#                 dsnameValue <- trim.blanks(tclvalue(dsname))
#                 newSS <- round(as.numeric(tclvalue(newDataSS)))
#                 closeDialog()
#                 if (dsnameValue == "") {
#                   errorCondition(recall = weibullSimulate.ipsur, 
#                     message = gettextRcmdr("You must enter the name of a data set."))
#                   return()
#                 }
#                 if (!is.valid.name(dsnameValue)) {
#                   errorCondition(recall = weibullSimulate.ipsur, 
#                     message = paste("\"", dsnameValue, "\" ", 
#                       gettextRcmdr("is not a valid name."), sep = ""))
#                   return()
#                 }
#                 if (is.element(dsnameValue, listDataSets())) {
#                   if ("no" == tclvalue(checkReplace(dsnameValue, 
#                     gettextRcmdr("Data set")))) {
#                     weibullSimulate.ipsur()
#                     return()
#                   }
#                 }
#                 if (is.na(newSS)) {
#                   errorCondition(recall = weibullSimulate.ipsur, 
#                     message = gettextRcmdr("Sample Size must be a positive integer."))
#                   return()
#                 }
#                 UpdateweibullsimNumber()
#                 justDoIt(paste(dsnameValue, " = data.frame(weibull.sim", 
#                   getRcmdr("weibullsimNumber"), "=1:", newSS, 
#                   ")", sep = ""))
#                 logger(paste(dsnameValue, "has been initialized."))
#                 for (k in getRcmdr("weibullsimNumber"):(nsamples + 
#                   getRcmdr("weibullsimNumber") - 1)) {
#                   justDoIt(paste(dsnameValue, "$weibull.sim", 
#                     k, " <- rweibull(", newSS, ", shape=", shape, 
#                     ", scale=", scale1, ")", sep = ""))
#                 }
#                 activeDataSet(dsnameValue)
#                 putRcmdr("weibullsimNumber", k)
#                 if (nsamples == 1) {
#                   logger(paste("There was 1 weibull variate sample stored in ", 
#                     dsnameValue, ".", sep = ""))
#                 }
#                 else {
#                   logger(paste("There were ", nsamples, " weibull variate samples stored in ", 
#                     dsnameValue, ".", sep = ""))
#                 }
#             }
#             OKCancelHelp(helpSubject = "rweibull")
#             tkgrid(tklabel(top, text = gettextRcmdr("Enter name for data set:")), 
#                 entryDsname, sticky = "e")
#             tkgrid(tklabel(top, text = gettextRcmdr("Sample Size (rows):")), 
#                 entryNewDataSS, sticky = "e")
#             tkgrid(buttonsFrame, columnspan = "2", sticky = "w")
#             tkgrid.configure(entryDsname, sticky = "w")
#             tkgrid.configure(entryNewDataSS, sticky = "w")
#             tkfocus(CommanderWindow())
#             dialogSuffix(rows = 2, columns = 2, focus = entryDsname)
#         }
#         else {
#             if (!is.character(ActiveDataSet())) {
#                 errorCondition(recall = weibullSimulate.ipsur, 
#                   message = gettextRcmdr("There is no active data set."))
#                 return()
#             }
#             .activeDataSet <- ActiveDataSet()
#             justDoIt(paste("samplesn <- dim(", .activeDataSet, 
#                 ")[1]", sep = ""))
#             UpdateweibullsimNumber()
#             for (k in getRcmdr("weibullsimNumber"):(nsamples + 
#                 getRcmdr("weibullsimNumber") - 1)) {
#                 justDoIt(paste(.activeDataSet, "$weibull.sim", 
#                   k, " <- rweibull(", samplesn, ", shape=", shape, 
#                   ", scale=", scale1, ")", sep = ""))
#             }
#             activeDataSet(.activeDataSet)
#             putRcmdr("weibullsimNumber", k)
#             if (nsamples == 1) {
#                 logger(paste("There was 1 weibull variate sample stored in ", 
#                   .activeDataSet, ".", sep = ""))
#             }
#             else {
#                 logger(paste("There were ", nsamples, " weibull variate samples stored in ", 
#                   .activeDataSet, ".", sep = ""))
#             }
#         }
#         tkfocus(CommanderWindow())
#     }
#     OKCancelHelp(helpSubject = "rweibull")
#     tkgrid(tklabel(top, text = gettextRcmdr("Number of samples (columns):")), 
#         samplesEntry, sticky = "w")
#     tkgrid(tklabel(top, text = gettextRcmdr("Parameters:"), fg = "blue"), 
#         columnspan = 4, sticky = "w")
#     tkgrid(tklabel(top, text = gettextRcmdr("shape")), shapeEntry, 
#         sticky = "w")
#     tkgrid(tklabel(top, text = gettextRcmdr("scale")), scale1Entry, 
#         sticky = "w")
#     tkgrid(tklabel(locationFrame, text = gettextRcmdr("Store values in:"), 
#         fg = "blue"), columnspan = 4, sticky = "w")
#     tkgrid(tklabel(locationFrame, text = gettextRcmdr("Active Dataset")), 
#         addtoactiveButton, sticky = "w")
#     tkgrid(tklabel(locationFrame, text = "New Dataset"), newDataButton, 
#         sticky = "w")
#     tkgrid.configure(samplesEntry, sticky = "w")
#     tkgrid.configure(shapeEntry, sticky = "w")
#     tkgrid.configure(scale1Entry, sticky = "w")
#     tkgrid(locationFrame, sticky = "w")
#     tkgrid(buttonsFrame, sticky = "w", columnspan = 2)
#     dialogSuffix(rows = 6, columns = 1, focus = samplesEntry)
# }
