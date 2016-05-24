GtkAnchorType<-c("center" = 0,
	"north" = 1,
	"north-west" = 2,
	"north-east" = 3,
	"south" = 4,
	"south-west" = 5,
	"south-east" = 6,
	"west" = 7,
	"east" = 8,
	"n" = 9,
	"nw" = 10,
	"ne" = 11,
	"s" = 12,
	"sw" = 13,
	"se" = 14,
	"w" = 15,
	"e" = 16)
storage.mode(GtkAnchorType) <- 'integer'
class(GtkAnchorType) <- 'enums' 

GtkArrowType<-c("up" = 0,
	"down" = 1,
	"left" = 2,
	"right" = 3)
storage.mode(GtkArrowType) <- 'integer'
class(GtkArrowType) <- 'enums' 

GtkButtonBoxStyle<-c("default-style" = 0,
	"spread" = 1,
	"edge" = 2,
	"start" = 3,
	"end" = 4)
storage.mode(GtkButtonBoxStyle) <- 'integer'
class(GtkButtonBoxStyle) <- 'enums' 

GtkButtonsType<-c("none" = 0,
	"ok" = 1,
	"close" = 2,
	"cancel" = 3,
	"yes-no" = 4,
	"ok-cancel" = 5)
storage.mode(GtkButtonsType) <- 'integer'
class(GtkButtonsType) <- 'enums' 

GtkCellRendererMode<-c("inert" = 0,
	"activatable" = 1,
	"editable" = 2)
storage.mode(GtkCellRendererMode) <- 'integer'
class(GtkCellRendererMode) <- 'enums' 

GtkCellType<-c("empty" = 0,
	"text" = 1,
	"pixmap" = 2,
	"pixtext" = 3,
	"widget" = 4)
storage.mode(GtkCellType) <- 'integer'
class(GtkCellType) <- 'enums' 

GtkCListDragPos<-c("none" = 0,
	"before" = 1,
	"into" = 2,
	"after" = 3)
storage.mode(GtkCListDragPos) <- 'integer'
class(GtkCListDragPos) <- 'enums' 

GtkCornerType<-c("top-left" = 0,
	"bottom-left" = 1,
	"top-right" = 2,
	"bottom-right" = 3)
storage.mode(GtkCornerType) <- 'integer'
class(GtkCornerType) <- 'enums' 

GtkCTreeExpanderStyle<-c("none" = 0,
	"square" = 1,
	"triangle" = 2,
	"circular" = 3)
storage.mode(GtkCTreeExpanderStyle) <- 'integer'
class(GtkCTreeExpanderStyle) <- 'enums' 

GtkCTreeExpansionType<-c("expand" = 0,
	"expand-recursive" = 1,
	"collapse" = 2,
	"collapse-recursive" = 3,
	"toggle" = 4,
	"toggle-recursive" = 5)
storage.mode(GtkCTreeExpansionType) <- 'integer'
class(GtkCTreeExpansionType) <- 'enums' 

GtkCTreeLineStyle<-c("none" = 0,
	"solid" = 1,
	"dotted" = 2,
	"tabbed" = 3)
storage.mode(GtkCTreeLineStyle) <- 'integer'
class(GtkCTreeLineStyle) <- 'enums' 

GtkCTreePos<-c("before" = 0,
	"as-child" = 1,
	"after" = 2)
storage.mode(GtkCTreePos) <- 'integer'
class(GtkCTreePos) <- 'enums' 

GtkCurveType<-c("linear" = 0,
	"spline" = 1,
	"free" = 2)
storage.mode(GtkCurveType) <- 'integer'
class(GtkCurveType) <- 'enums' 

GtkDeleteType<-c("chars" = 0,
	"word-ends" = 1,
	"words" = 2,
	"display-lines" = 3,
	"display-line-ends" = 4,
	"paragraph-ends" = 5,
	"paragraphs" = 6,
	"whitespace" = 7)
storage.mode(GtkDeleteType) <- 'integer'
class(GtkDeleteType) <- 'enums' 

GtkDirectionType<-c("tab-forward" = 0,
	"tab-backward" = 1,
	"up" = 2,
	"down" = 3,
	"left" = 4,
	"right" = 5)
storage.mode(GtkDirectionType) <- 'integer'
class(GtkDirectionType) <- 'enums' 

GtkExpanderStyle<-c("collapsed" = 0,
	"semi-collapsed" = 1,
	"semi-expanded" = 2,
	"expanded" = 3)
storage.mode(GtkExpanderStyle) <- 'integer'
class(GtkExpanderStyle) <- 'enums' 

GtkFileChooserAction<-c("open" = 0,
	"save" = 1,
	"select-folder" = 2,
	"create-folder" = 3)
storage.mode(GtkFileChooserAction) <- 'integer'
class(GtkFileChooserAction) <- 'enums' 

GtkFileChooserError<-c("nonexistent" = 0,
	"bad-filename" = 1)
storage.mode(GtkFileChooserError) <- 'integer'
class(GtkFileChooserError) <- 'enums' 

GtkFileChooserConfirmation<-c("confirm" = 0,
	"accept-filename" = 1,
	"select-again" = 2)
storage.mode(GtkFileChooserConfirmation) <- 'integer'
class(GtkFileChooserConfirmation) <- 'enums' 

GtkIconSize<-c("invalid" = 0,
	"menu" = 1,
	"small-toolbar" = 2,
	"large-toolbar" = 3,
	"button" = 4,
	"dnd" = 5,
	"dialog" = 6)
storage.mode(GtkIconSize) <- 'integer'
class(GtkIconSize) <- 'enums' 

GtkIconThemeError<-c("not-found" = 0,
	"failed" = 1)
storage.mode(GtkIconThemeError) <- 'integer'
class(GtkIconThemeError) <- 'enums' 

GtkIconViewDropPosition<-c("no-drop" = 0,
	"drop-into" = 1,
	"drop-left" = 2,
	"drop-right" = 3,
	"drop-above" = 4,
	"drop-below" = 5)
storage.mode(GtkIconViewDropPosition) <- 'integer'
class(GtkIconViewDropPosition) <- 'enums' 

GtkImageType<-c("empty" = 0,
	"pixmap" = 1,
	"image" = 2,
	"pixbuf" = 3,
	"stock" = 4,
	"icon-set" = 5,
	"animation" = 6)
storage.mode(GtkImageType) <- 'integer'
class(GtkImageType) <- 'enums' 

GtkIMPreeditStyle<-c("nothing" = 0,
	"callback" = 1,
	"none" = 2)
storage.mode(GtkIMPreeditStyle) <- 'integer'
class(GtkIMPreeditStyle) <- 'enums' 

GtkIMStatusStyle<-c("nothing" = 0,
	"callback" = 1)
storage.mode(GtkIMStatusStyle) <- 'integer'
class(GtkIMStatusStyle) <- 'enums' 

GtkJustification<-c("left" = 0,
	"right" = 1,
	"center" = 2,
	"fill" = 3)
storage.mode(GtkJustification) <- 'integer'
class(GtkJustification) <- 'enums' 

GtkMatchType<-c("all" = 0,
	"all-tail" = 1,
	"head" = 2,
	"tail" = 3,
	"exact" = 4,
	"last" = 5)
storage.mode(GtkMatchType) <- 'integer'
class(GtkMatchType) <- 'enums' 

GtkMenuDirectionType<-c("parent" = 0,
	"child" = 1,
	"next" = 2,
	"prev" = 3)
storage.mode(GtkMenuDirectionType) <- 'integer'
class(GtkMenuDirectionType) <- 'enums' 

GtkMessageType<-c("info" = 0,
	"warning" = 1,
	"question" = 2,
	"error" = 3)
storage.mode(GtkMessageType) <- 'integer'
class(GtkMessageType) <- 'enums' 

GtkMetricType<-c("pixels" = 0,
	"inches" = 1,
	"centimeters" = 2)
storage.mode(GtkMetricType) <- 'integer'
class(GtkMetricType) <- 'enums' 

GtkMovementStep<-c("logical-positions" = 0,
	"visual-positions" = 1,
	"words" = 2,
	"display-lines" = 3,
	"display-line-ends" = 4,
	"paragraphs" = 5,
	"paragraph-ends" = 6,
	"pages" = 7,
	"buffer-ends" = 8,
	"horizontal-pages" = 9)
storage.mode(GtkMovementStep) <- 'integer'
class(GtkMovementStep) <- 'enums' 

GtkNotebookTab<-c("first" = 0,
	"last" = 1)
storage.mode(GtkNotebookTab) <- 'integer'
class(GtkNotebookTab) <- 'enums' 

GtkOrientation<-c("horizontal" = 0,
	"vertical" = 1)
storage.mode(GtkOrientation) <- 'integer'
class(GtkOrientation) <- 'enums' 

GtkPackDirection<-c("ltr" = 0,
	"rtl" = 1,
	"ttb" = 2,
	"btt" = 3)
storage.mode(GtkPackDirection) <- 'integer'
class(GtkPackDirection) <- 'enums' 

GtkPackType<-c("start" = 0,
	"end" = 1)
storage.mode(GtkPackType) <- 'integer'
class(GtkPackType) <- 'enums' 

GtkPathPriorityType<-c("lowest" = 0,
	"gtk" = 1,
	"application" = 2,
	"theme" = 3,
	"rc" = 4,
	"highest" = 5)
storage.mode(GtkPathPriorityType) <- 'integer'
class(GtkPathPriorityType) <- 'enums' 

GtkPathType<-c("widget" = 0,
	"widget-class" = 1,
	"class" = 2)
storage.mode(GtkPathType) <- 'integer'
class(GtkPathType) <- 'enums' 

GtkPolicyType<-c("always" = 0,
	"automatic" = 1,
	"never" = 2)
storage.mode(GtkPolicyType) <- 'integer'
class(GtkPolicyType) <- 'enums' 

GtkPositionType<-c("left" = 0,
	"right" = 1,
	"top" = 2,
	"bottom" = 3)
storage.mode(GtkPositionType) <- 'integer'
class(GtkPositionType) <- 'enums' 

GtkPreviewType<-c("color" = 0,
	"grayscale" = 1)
storage.mode(GtkPreviewType) <- 'integer'
class(GtkPreviewType) <- 'enums' 

GtkProgressBarOrientation<-c("left-to-right" = 0,
	"right-to-left" = 1,
	"bottom-to-top" = 2,
	"top-to-bottom" = 3)
storage.mode(GtkProgressBarOrientation) <- 'integer'
class(GtkProgressBarOrientation) <- 'enums' 

GtkProgressBarStyle<-c("continuous" = 0,
	"discrete" = 1)
storage.mode(GtkProgressBarStyle) <- 'integer'
class(GtkProgressBarStyle) <- 'enums' 

GtkRcTokenType<-c("invalid" = 0,
	"include" = 1,
	"normal" = 2,
	"active" = 3,
	"prelight" = 4,
	"selected" = 5,
	"insensitive" = 6,
	"fg" = 7,
	"bg" = 8,
	"text" = 9,
	"base" = 10,
	"xthickness" = 11,
	"ythickness" = 12,
	"font" = 13,
	"fontset" = 14,
	"font-name" = 15,
	"bg-pixmap" = 16,
	"pixmap-path" = 17,
	"style" = 18,
	"binding" = 19,
	"bind" = 20,
	"widget" = 21,
	"widget-class" = 22,
	"class" = 23,
	"lowest" = 24,
	"gtk" = 25,
	"application" = 26,
	"theme" = 27,
	"rc" = 28,
	"highest" = 29,
	"engine" = 30,
	"module-path" = 31,
	"im-module-path" = 32,
	"im-module-file" = 33,
	"stock" = 34,
	"ltr" = 35,
	"rtl" = 36,
	"last" = 37)
storage.mode(GtkRcTokenType) <- 'integer'
class(GtkRcTokenType) <- 'enums' 

GtkReliefStyle<-c("normal" = 0,
	"half" = 1,
	"none" = 2)
storage.mode(GtkReliefStyle) <- 'integer'
class(GtkReliefStyle) <- 'enums' 

GtkResizeMode<-c("parent" = 0,
	"queue" = 1,
	"immediate" = 2)
storage.mode(GtkResizeMode) <- 'integer'
class(GtkResizeMode) <- 'enums' 

GtkResponseType<-c("none" = -1,
	"reject" = -2,
	"accept" = -3,
	"delete-event" = -4,
	"ok" = -5,
	"cancel" = -6,
	"close" = -7,
	"yes" = -8,
	"no" = -9,
	"apply" = -10,
	"help" = -11)
storage.mode(GtkResponseType) <- 'integer'
class(GtkResponseType) <- 'enums' 

GtkScrollStep<-c("steps" = 0,
	"pages" = 1,
	"ends" = 2,
	"horizontal-steps" = 3,
	"horizontal-pages" = 4,
	"horizontal-ends" = 5)
storage.mode(GtkScrollStep) <- 'integer'
class(GtkScrollStep) <- 'enums' 

GtkScrollType<-c("none" = 0,
	"jump" = 1,
	"step-backward" = 2,
	"step-forward" = 3,
	"page-backward" = 4,
	"page-forward" = 5,
	"step-up" = 6,
	"step-down" = 7,
	"page-up" = 8,
	"page-down" = 9,
	"step-left" = 10,
	"step-right" = 11,
	"page-left" = 12,
	"page-right" = 13,
	"start" = 14,
	"end" = 15)
storage.mode(GtkScrollType) <- 'integer'
class(GtkScrollType) <- 'enums' 

GtkSelectionMode<-c("none" = 0,
	"single" = 1,
	"browse" = 2,
	"multiple" = 3,
	"extended" = 4)
storage.mode(GtkSelectionMode) <- 'integer'
class(GtkSelectionMode) <- 'enums' 

GtkShadowType<-c("none" = 0,
	"in" = 1,
	"out" = 2,
	"etched-in" = 3,
	"etched-out" = 4)
storage.mode(GtkShadowType) <- 'integer'
class(GtkShadowType) <- 'enums' 

GtkSideType<-c("top" = 0,
	"bottom" = 1,
	"left" = 2,
	"right" = 3)
storage.mode(GtkSideType) <- 'integer'
class(GtkSideType) <- 'enums' 

GtkSizeGroupMode<-c("none" = 0,
	"horizontal" = 1,
	"vertical" = 2,
	"both" = 3)
storage.mode(GtkSizeGroupMode) <- 'integer'
class(GtkSizeGroupMode) <- 'enums' 

GtkSortType<-c("ascending" = 0,
	"descending" = 1)
storage.mode(GtkSortType) <- 'integer'
class(GtkSortType) <- 'enums' 

GtkSpinButtonUpdatePolicy<-c("always" = 0,
	"if-valid" = 1)
storage.mode(GtkSpinButtonUpdatePolicy) <- 'integer'
class(GtkSpinButtonUpdatePolicy) <- 'enums' 

GtkSpinType<-c("step-forward" = 0,
	"step-backward" = 1,
	"page-forward" = 2,
	"page-backward" = 3,
	"home" = 4,
	"end" = 5,
	"user-defined" = 6)
storage.mode(GtkSpinType) <- 'integer'
class(GtkSpinType) <- 'enums' 

GtkStateType<-c("normal" = 0,
	"active" = 1,
	"prelight" = 2,
	"selected" = 3,
	"insensitive" = 4)
storage.mode(GtkStateType) <- 'integer'
class(GtkStateType) <- 'enums' 

GtkSubmenuDirection<-c("left" = 0,
	"right" = 1)
storage.mode(GtkSubmenuDirection) <- 'integer'
class(GtkSubmenuDirection) <- 'enums' 

GtkSubmenuPlacement<-c("top-bottom" = 0,
	"left-right" = 1)
storage.mode(GtkSubmenuPlacement) <- 'integer'
class(GtkSubmenuPlacement) <- 'enums' 

GtkTextDirection<-c("none" = 0,
	"ltr" = 1,
	"rtl" = 2)
storage.mode(GtkTextDirection) <- 'integer'
class(GtkTextDirection) <- 'enums' 

GtkTextWindowType<-c("private" = 0,
	"widget" = 1,
	"text" = 2,
	"left" = 3,
	"right" = 4,
	"top" = 5,
	"bottom" = 6)
storage.mode(GtkTextWindowType) <- 'integer'
class(GtkTextWindowType) <- 'enums' 

GtkToolbarChildType<-c("space" = 0,
	"button" = 1,
	"togglebutton" = 2,
	"radiobutton" = 3,
	"widget" = 4)
storage.mode(GtkToolbarChildType) <- 'integer'
class(GtkToolbarChildType) <- 'enums' 

GtkToolbarSpaceStyle<-c("empty" = 0,
	"line" = 1)
storage.mode(GtkToolbarSpaceStyle) <- 'integer'
class(GtkToolbarSpaceStyle) <- 'enums' 

GtkToolbarStyle<-c("icons" = 0,
	"text" = 1,
	"both" = 2,
	"both-horiz" = 3)
storage.mode(GtkToolbarStyle) <- 'integer'
class(GtkToolbarStyle) <- 'enums' 

GtkTreeViewColumnSizing<-c("grow-only" = 0,
	"autosize" = 1,
	"fixed" = 2)
storage.mode(GtkTreeViewColumnSizing) <- 'integer'
class(GtkTreeViewColumnSizing) <- 'enums' 

GtkTreeViewDropPosition<-c("before" = 0,
	"after" = 1,
	"into-or-before" = 2,
	"into-or-after" = 3)
storage.mode(GtkTreeViewDropPosition) <- 'integer'
class(GtkTreeViewDropPosition) <- 'enums' 

GtkUpdateType<-c("continuous" = 0,
	"discontinuous" = 1,
	"delayed" = 2)
storage.mode(GtkUpdateType) <- 'integer'
class(GtkUpdateType) <- 'enums' 

GtkVisibility<-c("none" = 0,
	"partial" = 1,
	"full" = 2)
storage.mode(GtkVisibility) <- 'integer'
class(GtkVisibility) <- 'enums' 

GtkWidgetHelpType<-c("tooltip" = 0,
	"whats-this" = 1)
storage.mode(GtkWidgetHelpType) <- 'integer'
class(GtkWidgetHelpType) <- 'enums' 

GtkWindowPosition<-c("none" = 0,
	"center" = 1,
	"mouse" = 2,
	"center-always" = 3,
	"center-on-parent" = 4)
storage.mode(GtkWindowPosition) <- 'integer'
class(GtkWindowPosition) <- 'enums' 

GtkWindowType<-c("toplevel" = 0,
	"popup" = 1)
storage.mode(GtkWindowType) <- 'integer'
class(GtkWindowType) <- 'enums' 

GtkWrapMode<-c("none" = 0,
	"char" = 1,
	"word" = 2,
	"word_char" = 3)
storage.mode(GtkWrapMode) <- 'integer'
class(GtkWrapMode) <- 'enums' 

GtkAssistantPageType<-c("content" = 0,
	"intro" = 1,
	"confirm" = 2,
	"summary" = 3,
	"progress" = 4)
storage.mode(GtkAssistantPageType) <- 'integer'
class(GtkAssistantPageType) <- 'enums' 

GtkCellRendererAccelMode<-c("gtk" = 0,
	"other" = 1)
storage.mode(GtkCellRendererAccelMode) <- 'integer'
class(GtkCellRendererAccelMode) <- 'enums' 

GtkSensitivityType<-c("auto" = 0,
	"on" = 1,
	"off" = 2)
storage.mode(GtkSensitivityType) <- 'integer'
class(GtkSensitivityType) <- 'enums' 

GtkPrintPages<-c("all" = 0,
	"current" = 1,
	"ranges" = 2)
storage.mode(GtkPrintPages) <- 'integer'
class(GtkPrintPages) <- 'enums' 

GtkPageSet<-c("all" = 0,
	"even" = 1,
	"odd" = 2)
storage.mode(GtkPageSet) <- 'integer'
class(GtkPageSet) <- 'enums' 

GtkPageOrientation<-c("portrait" = 0,
	"landscape" = 1,
	"reverse-portrait" = 2,
	"reverse-landscape" = 3)
storage.mode(GtkPageOrientation) <- 'integer'
class(GtkPageOrientation) <- 'enums' 

GtkPrintQuality<-c("low" = 0,
	"normal" = 1,
	"high" = 2,
	"draft" = 3)
storage.mode(GtkPrintQuality) <- 'integer'
class(GtkPrintQuality) <- 'enums' 

GtkPrintDuplex<-c("simplex" = 0,
	"horizontal" = 1,
	"vertical" = 2)
storage.mode(GtkPrintDuplex) <- 'integer'
class(GtkPrintDuplex) <- 'enums' 

GtkUnit<-c("pixel" = 0,
	"points" = 1,
	"inch" = 2,
	"mm" = 3)
storage.mode(GtkUnit) <- 'integer'
class(GtkUnit) <- 'enums' 

GtkTreeViewGridLines<-c("none" = 0,
	"horizontal" = 1,
	"vertical" = 2,
	"both" = 3)
storage.mode(GtkTreeViewGridLines) <- 'integer'
class(GtkTreeViewGridLines) <- 'enums' 

GtkPrintStatus<-c("initial" = 0,
	"preparing" = 1,
	"generating-data" = 2,
	"sending-data" = 3,
	"pending" = 4,
	"pending-issue" = 5,
	"printing" = 6,
	"finished" = 7,
	"finished-aborted" = 8)
storage.mode(GtkPrintStatus) <- 'integer'
class(GtkPrintStatus) <- 'enums' 

GtkPrintOperationResult<-c("error" = 0,
	"apply" = 1,
	"cancel" = 2,
	"in-progress" = 3)
storage.mode(GtkPrintOperationResult) <- 'integer'
class(GtkPrintOperationResult) <- 'enums' 

GtkPrintOperationAction<-c("print-dialog" = 0,
	"print" = 1,
	"preview" = 2,
	"export" = 3)
storage.mode(GtkPrintOperationAction) <- 'integer'
class(GtkPrintOperationAction) <- 'enums' 

GtkPrintError<-c("general" = 0,
	"internal-error" = 1,
	"nomem" = 2)
storage.mode(GtkPrintError) <- 'integer'
class(GtkPrintError) <- 'enums' 

GtkRecentSortType<-c("none" = 0,
	"mru" = 1,
	"lru" = 2,
	"custom" = 3)
storage.mode(GtkRecentSortType) <- 'integer'
class(GtkRecentSortType) <- 'enums' 

GtkRecentChooserError<-c("not-found" = 0,
	"invalid-uri" = 1)
storage.mode(GtkRecentChooserError) <- 'integer'
class(GtkRecentChooserError) <- 'enums' 

GtkRecentManagerError<-c("not-found" = 0,
	"invalid-uri" = 1,
	"invalid-encoding" = 2,
	"not-registered" = 3,
	"read" = 4,
	"write" = 5,
	"unknown" = 6)
storage.mode(GtkRecentManagerError) <- 'integer'
class(GtkRecentManagerError) <- 'enums' 

GtkTextBufferTargetInfo<-c("buffer-contents" = 0,
	"rich-text" = 1,
	"text" = 2)
storage.mode(GtkTextBufferTargetInfo) <- 'integer'
class(GtkTextBufferTargetInfo) <- 'enums' 

GtkBuilderError<-c("invalid-type-function" = 0,
	"unhandled-tag" = 1,
	"missing-attribute" = 2,
	"invalid-attribute" = 3,
	"invalid-tag" = 4,
	"missing-property-value" = 5,
	"invalid-value" = 6)
storage.mode(GtkBuilderError) <- 'integer'
class(GtkBuilderError) <- 'enums' 

GtkNumberUpLayout<-c("left-to-right-top-to-bottom" = 0,
	"left-to-right-bottom-to-top" = 1,
	"right-to-left-top-to-bottom" = 2,
	"right-to-left-bottom-to-top" = 3,
	"top-to-bottom-left-to-right" = 4,
	"top-to-bottom-right-to-left" = 5,
	"bottom-to-top-left-to-right" = 6,
	"bottom-to-top-right-to-left" = 7)
storage.mode(GtkNumberUpLayout) <- 'integer'
class(GtkNumberUpLayout) <- 'enums' 

GtkEntryIconPosition<-c("primary" = 0,
	"secondary" = 1)
storage.mode(GtkEntryIconPosition) <- 'integer'
class(GtkEntryIconPosition) <- 'enums' 

GtkAccelFlags<-c("visible" = 1,
	"locked" = 2,
	"mask" = 4)
storage.mode(GtkAccelFlags) <- 'numeric'
class(GtkAccelFlags) <- 'flags' 

GtkAttachOptions<-c("expand" = 1,
	"shrink" = 2,
	"fill" = 4)
storage.mode(GtkAttachOptions) <- 'numeric'
class(GtkAttachOptions) <- 'flags' 

GtkButtonAction<-c("ignored" = 1,
	"selects" = 2,
	"drags" = 4,
	"expands" = 8)
storage.mode(GtkButtonAction) <- 'numeric'
class(GtkButtonAction) <- 'flags' 

GtkCalendarDisplayOptions<-c("show-heading" = 1,
	"show-day-names" = 2,
	"no-month-change" = 4,
	"show-week-numbers" = 8,
	"week-start-monday" = 16)
storage.mode(GtkCalendarDisplayOptions) <- 'numeric'
class(GtkCalendarDisplayOptions) <- 'flags' 

GtkCellRendererState<-c("selected" = 1,
	"prelit" = 2,
	"insensitive" = 4,
	"sorted" = 8,
	"focused" = 16)
storage.mode(GtkCellRendererState) <- 'numeric'
class(GtkCellRendererState) <- 'flags' 

GtkDestDefaults<-c("motion" = 1,
	"highlight" = 2,
	"drop" = 4,
	"all" = 8)
storage.mode(GtkDestDefaults) <- 'numeric'
class(GtkDestDefaults) <- 'flags' 

GtkDialogFlags<-c("modal" = 1,
	"destroy-with-parent" = 2,
	"no-separator" = 4)
storage.mode(GtkDialogFlags) <- 'numeric'
class(GtkDialogFlags) <- 'flags' 

GtkFileFilterFlags<-c("filename" = 1,
	"uri" = 2,
	"display-name" = 4,
	"mime-type" = 8)
storage.mode(GtkFileFilterFlags) <- 'numeric'
class(GtkFileFilterFlags) <- 'flags' 

GtkIconLookupFlags<-c("no-svg" = 1,
	"force-svg" = 2,
	"use-builtin" = 4)
storage.mode(GtkIconLookupFlags) <- 'numeric'
class(GtkIconLookupFlags) <- 'flags' 

GtkRcFlags<-c("fg" = 1,
	"bg" = 2,
	"text" = 4,
	"base" = 8)
storage.mode(GtkRcFlags) <- 'numeric'
class(GtkRcFlags) <- 'flags' 

GtkTargetFlags<-c("app" = 1,
	"widget" = 2)
storage.mode(GtkTargetFlags) <- 'numeric'
class(GtkTargetFlags) <- 'flags' 

GtkTextSearchFlags<-c("visible-only" = 1,
	"text-only" = 2)
storage.mode(GtkTextSearchFlags) <- 'numeric'
class(GtkTextSearchFlags) <- 'flags' 

GtkTreeModelFlags<-c("iters-persist" = 1,
	"list-only" = 2)
storage.mode(GtkTreeModelFlags) <- 'numeric'
class(GtkTreeModelFlags) <- 'flags' 

GtkUIManagerItemType<-c("auto" = 1,
	"menubar" = 2,
	"menu" = 4,
	"toolbar" = 8,
	"placeholder" = 16,
	"popup" = 32,
	"menuitem" = 64,
	"toolitem" = 128,
	"separator" = 256,
	"accelerator" = 512)
storage.mode(GtkUIManagerItemType) <- 'numeric'
class(GtkUIManagerItemType) <- 'flags' 

GtkWidgetFlags<-c("toplevel" = 16,
	"no-window" = 32,
	"realized" = 64,
	"mapped" = 128,
	"visible" = 256,
	"sensitive" = 512,
	"parent-sensitive" = 1024,
	"can-focus" = 2048,
	"has-focus" = 4096,
	"can-default" = 8192,
	"has-default" = 16384,
	"has-grab" = 32768,
	"rc-style" = 16384,
	"composite-child" = 131072,
	"no-reparent" = 262144,
	"app-paintable" = 524288,
	"receives-default" = 1048576,
	"double-buffered" = 2097152,
	"no-show-all" = 4194304)
storage.mode(GtkWidgetFlags) <- 'numeric'
class(GtkWidgetFlags) <- 'flags' 

GtkRecentFilterFlags<-c("uri" = 1,
	"display-name" = 2,
	"mime-type" = 4,
	"application" = 8,
	"group" = 16,
	"age" = 32)
storage.mode(GtkRecentFilterFlags) <- 'numeric'
class(GtkRecentFilterFlags) <- 'flags' 

GtkToolPaletteDragTargets<-c("items" = 1,
	"groups" = 2)
storage.mode(GtkToolPaletteDragTargets) <- 'numeric'
class(GtkToolPaletteDragTargets) <- 'flags' 

