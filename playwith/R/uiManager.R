## playwith: interactive plots in R using GTK+
##
## Copyright (c) 2007 Felix Andrews <felix@nfrac.org>
## GPL version 2 or newer

constructUIManager <- function(playState, tools)
{
    menuEntries <-
        list(
             list("FileMenu", NULL, "_File"),
             list("ViewMenu", NULL, "_View"),
             list("StyleMenu", NULL, "_Style"),
             list("ThemeMenu", NULL, "Th_eme"),
             list("LabelsMenu", NULL, "_Labels"),
             list("ToolsMenu", NULL, "_Tools"),
             list("DataMenu", NULL, "_Data"),
             list("OptionsMenu", NULL, "_Options"),
             list("HelpMenu", NULL, "_Help")
             )
    menuGroup <- gtkActionGroupNew("Menus")
    menuGroup$addActions(menuEntries)
#    for (nm in c("ShortcutsMenu", "ThemesMenu")) {
#        tmp <- menuGroup$getAction(nm)
#        if (!is.null(tmp)) tmp["hide-if-empty"] <- FALSE
#    }
    ## ui manager
    manager <- gtkUIManagerNew()
    window <- playState$win
    window$setData("ui-manager", manager)
    manager$insertActionGroup(plotActionGroup(playState), 0)
    manager$insertActionGroup(globalActionGroup(playState), 0)
    ## user-defined actions:
    ## (from playwith.options and the tools argument)
    customGroup <- gtkActionGroupNew("CustomActions")
    if (length(tools) > 0) {
        ## toggle actions must have 7 elements or a name "is_active"
        toggles <- ((sapply(tools, length) == 7) |
                    sapply(tools, function(x) !is.null(x$is_active)))
        if (any(toggles))
            customGroup$addToggleActions(tools[toggles], playState)
        if (any(!toggles))
            customGroup$addActions(tools[!toggles], playState)
        manager$insertActionGroup(customGroup, 1)
    }
    ## the menus themselves
    manager$insertActionGroup(menuGroup, 0)
    window$addAccelGroup(manager$getAccelGroup())
    ## read in structure of menus and toolbars specified in XML
    for (opt in c("ui.menus.xml", "ui.toolbars.xml", "ui.custom.xml")) {
        uifile <- playwith.getOption(opt)
        if (is.character(uifile) && (nchar(uifile) > 0))
            manager$addUiFromFile(uifile)
    }
    manager$ensureUpdate()
    createStyleActions(playState, manager)
    ## construct UI for user-defined actions
    ## (unless ui.custom.xml was specified)
    if (is.null(playwith.getOption("ui.custom.xml"))) {
        customTb <- playwith.getOption("custom.toolbar")
        for (i in seq_along(tools)) {
            actionName <- tools[[i]][[1]]
            action <- customGroup$getAction(actionName)
            tm <- manager$getWidget("/MenuBar/ToolsMenu")
            tm$getSubmenu()$append(action$createMenuItem())
            if (!is.null(customTb)) {
                tb <- manager$getWidget(paste("/", customTb, sep = ""))
                tb$insert(action$createToolItem(), -1)
            }
            ## TODO: (e.g. need for accelerators)
#            manager$addUi(manager$newMergeId(),
#                      path = "/ui/MenuBar/ToolsMenu",
#                      name = actionName,
#                      action = actionName,
#                      type = GtkUIManagerItemType["auto"],
#                      top = FALSE)
#            manager$addUi(manager$newMergeId(),
#                      path = paste("/ui/", customTb, sep = ""),
#                      name = actionName,
#                      action = actionName,
#                      type = GtkUIManagerItemType["auto"],
#                      top = FALSE)
        }
    }
    manager
}

initActions <- function(playState)
{
    playDevSet(playState)
    initClickActions(playState)
    initIdentifyActions(playState)
    initOptionsActions(playState)
    ## custom init actions
    customAct <- c(playwith.getOption("init.actions"),
                   playState$init.actions)
    for (x in customAct) {
        playDevSet(playState)
        if (is.character(x)) x <- get(x)
        if (is.function(x)) x(playState)
        if (is.language(x)) eval(x, playState$env)
    }
}

preplotActions <- function(playState)
{
    customAct <- c(playwith.getOption("preplot.actions"),
                   playState$preplot.actions)
    for (x in customAct) {
        playDevSet(playState)
        if (is.character(x)) x <- get(x)
        if (is.function(x)) x(playState)
        if (is.language(x)) eval(x, playState$env)
    }
}

updateActions <- function(playState)
{
    playDevSet(playState)
    updateGlobalActions(playState)
    updateClickActions(playState)
    updatePlotActions(playState)
    updateGrobActions(playState)
    updateOptionsActions(playState)
    ## custom update actions
    customAct <- c(playwith.getOption("update.actions"),
                   playState$update.actions)
    for (x in customAct) {
        playDevSet(playState)
        if (is.character(x)) x <- get(x)
        if (is.function(x)) x(playState)
        if (is.language(x)) eval(x, playState$env)
    }
    ## these are user-defined, should draw on top:
    updateIdentifyActions(playState)
    updateAnnotationActions(playState)
}
