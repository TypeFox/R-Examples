gtkAboutDialog <- gtkAboutDialogNew

gtkAccelGroup <- gtkAccelGroupNew

gtkAccelLabel <- gtkAccelLabelNew

gtkAction <- gtkActionNew

gtkActionGroup <- gtkActionGroupNew

gtkAdjustment <- gtkAdjustmentNew

gtkAlignment <- gtkAlignmentNew

gtkArrow <- gtkArrowNew

gtkAspectFrame <- gtkAspectFrameNew

gtkButton <-
function(label, stock.id, show = TRUE)
{
  if (!missing(stock.id)) {
    gtkButtonNewFromStock(stock.id, show)
  }
  else {
    if (!missing(label)) {
      gtkButtonNewWithLabel(label, show)
    }
    else {
      gtkButtonNew(show)
    }
  }
}

gtkCalendar <- gtkCalendarNew

gtkCellRendererCombo <- gtkCellRendererComboNew

gtkCellRendererPixbuf <- gtkCellRendererPixbufNew

gtkCellRendererProgress <- gtkCellRendererProgressNew

gtkCellRendererText <- gtkCellRendererTextNew

gtkCellRendererToggle <- gtkCellRendererToggleNew

gtkCellView <- gtkCellViewNew

gtkCheckButton <-
function(label, show = TRUE)
{
  if (!missing(label)) {
    gtkCheckButtonNewWithLabel(label, show)
  }
  else {
    gtkCheckButtonNew(show)
  }
}

gtkCheckMenuItem <-
function(label, show = TRUE)
{
  if (!missing(label)) {
    gtkCheckMenuItemNewWithLabel(label, show)
  }
  else {
    gtkCheckMenuItemNew(show)
  }
}

gtkClipboard <- gtkClipboardGetForDisplay

gtkCList <-
function(columns = 1, titles, show = TRUE)
{
  if (!missing(titles)) {
    gtkCListNewWithTitles(columns, titles, show)
  }
  else {
    gtkCListNew(columns, show)
  }
}

gtkColorButton <-
function(color, show = TRUE)
{
  if (!missing(color)) {
    gtkColorButtonNewWithColor(color, show)
  }
  else {
    gtkColorButtonNew(show)
  }
}

gtkColorSelection <- gtkColorSelectionNew

gtkColorSelectionDialog <- gtkColorSelectionDialogNew

gtkCombo <- gtkComboNew

gtkComboBox <-
function(model, show = TRUE)
{
  if (!missing(model)) {
    gtkComboBoxNewWithModel(model, show)
  }
  else {
    gtkComboBoxNew(show)
  }
}

gtkComboBoxEntry <-
function(model, text.column, show = TRUE)
{
  if (!missing(model)) {
    gtkComboBoxEntryNewWithModel(model, text.column, show)
  }
  else {
    gtkComboBoxEntryNew(show)
  }
}

gtkCTree <-
function(columns = 1, tree.column = 0, titles, show = TRUE)
{
  if (!missing(titles)) {
    gtkCTreeNewWithTitles(columns, tree.column, titles, show)
  }
  else {
    gtkCTreeNew(columns, tree.column, show)
  }
}

gtkCurve <- gtkCurveNew

gtkDialog <-
function(title = NULL, parent = NULL, flags = 0, ..., show = TRUE)
{
  if (!missing(title)) {
    gtkDialogNewWithButtons(title, parent, flags, ..., show = show)
  }
  else {
    gtkDialogNew(show)
  }
}

gtkDrawingArea <- gtkDrawingAreaNew

gtkEntry <-
function(max = 0, buffer, show = TRUE)
{
  if (!missing(max)) {
    gtkEntryNewWithMaxLength(max, show)
  }
  else {
    if (!missing(buffer)) {
      gtkEntryNewWithBuffer(buffer, show)
    }
    else {
      gtkEntryNew(show)
    }
  }
}

gtkEntryCompletion <- gtkEntryCompletionNew

gtkEventBox <- gtkEventBoxNew

gtkExpander <- gtkExpanderNew

gtkFileChooserButton <-
function(title, action, backend, show = TRUE)
{
  if (!missing(backend)) {
    gtkFileChooserButtonNewWithBackend(title, action, backend, show)
  }
  else {
    gtkFileChooserButtonNew(title, action, show)
  }
}

gtkFileChooserDialog <-
function(title = NULL, parent = NULL, action, ..., backend, show = TRUE)
{
  if (!missing(backend)) {
    gtkFileChooserDialogNewWithBackend(title, parent, action, backend, ..., show = show)
  }
  else {
    gtkFileChooserDialogNew(title, parent, action, ..., show = show)
  }
}

gtkFileChooserWidget <-
function(action, backend, show = TRUE)
{
  if (!missing(backend)) {
    gtkFileChooserWidgetNewWithBackend(action, backend, show)
  }
  else {
    gtkFileChooserWidgetNew(action, show)
  }
}

gtkFileFilter <- gtkFileFilterNew

gtkFileSelection <- gtkFileSelectionNew

gtkFixed <- gtkFixedNew

gtkFontButton <- gtkFontButtonNew

gtkFontSelection <- gtkFontSelectionNew

gtkFontSelectionDialog <- gtkFontSelectionDialogNew

gtkFrame <- gtkFrameNew

gtkGammaCurve <- gtkGammaCurveNew

gtkHandleBox <- gtkHandleBoxNew

gtkHBox <- gtkHBoxNew

gtkHButtonBox <- gtkHButtonBoxNew

gtkHPaned <- gtkHPanedNew

gtkHRuler <- gtkHRulerNew

gtkHScale <-
function(adjustment = NULL, min, max, step, show = TRUE)
{
  if (!missing(adjustment)) {
    gtkHScaleNew(adjustment, show)
  }
  else {
    gtkHScaleNewWithRange(min, max, step, show)
  }
}

gtkHScrollbar <- gtkHScrollbarNew

gtkHSeparator <- gtkHSeparatorNew

gtkIconFactory <- gtkIconFactoryNew

gtkIconTheme <- gtkIconThemeNew

gtkIconView <-
function(model = NULL, show = TRUE)
{
  if (!missing(model)) {
    gtkIconViewNewWithModel(model, show)
  }
  else {
    gtkIconViewNew(show)
  }
}

gtkImage <-
function(size, mask = NULL, pixmap = NULL, image = NULL, filename, pixbuf = NULL, stock.id, icon.set, animation, icon, show = TRUE)
{
  if (!missing(pixmap)) {
    gtkImageNewFromPixmap(pixmap, mask, show)
  }
  else {
    if (!missing(mask)) {
      gtkImageNewFromImage(image, mask, show)
    }
    else {
      if (!missing(filename)) {
        gtkImageNewFromFile(filename, show)
      }
      else {
        if (!missing(pixbuf)) {
          gtkImageNewFromPixbuf(pixbuf, show)
        }
        else {
          if (!missing(stock.id)) {
            gtkImageNewFromStock(stock.id, size, show)
          }
          else {
            if (!missing(size)) {
              if (!missing(icon.set)) {
                gtkImageNewFromIconSet(icon.set, size, show)
              }
              else {
                gtkImageNewFromGicon(icon, size, show)
              }
            }
            else {
              if (!missing(animation)) {
                gtkImageNewFromAnimation(animation, show)
              }
              else {
                gtkImageNew(show)
              }
            }
          }
        }
      }
    }
  }
}

gtkImageMenuItem <-
function(label, stock.id, accel.group, show = TRUE)
{
  if (!missing(stock.id)) {
    gtkImageMenuItemNewFromStock(stock.id, accel.group, show)
  }
  else {
    if (!missing(label)) {
      gtkImageMenuItemNewWithLabel(label, show)
    }
    else {
      gtkImageMenuItemNew(show)
    }
  }
}

gtkIMContextSimple <- gtkIMContextSimpleNew

gtkIMMulticontext <- gtkIMMulticontextNew

gtkInputDialog <- gtkInputDialogNew

gtkInvisible <-
function(screen, show = TRUE)
{
  if (!missing(screen)) {
    gtkInvisibleNewForScreen(screen, show)
  }
  else {
    gtkInvisibleNew(show)
  }
}

gtkItemFactory <- gtkItemFactoryNew

gtkLabel <-
function(str = NULL, show = TRUE)
{
  gtkLabelNew(str, show)
}

gtkLayout <- gtkLayoutNew

gtkList <- gtkListNew

gtkListItem <-
function(label, show = TRUE)
{
  if (!missing(label)) {
    gtkListItemNewWithLabel(label, show)
  }
  else {
    gtkListItemNew(show)
  }
}

gtkListStore <-
function(..., value)
{
  if (!missing(...)) {
    gtkListStoreNew(...)
  }
  else {
    gtkListStoreNewv(value)
  }
}

gtkMenu <- gtkMenuNew

gtkMenuBar <- gtkMenuBarNew

gtkMenuItem <-
function(label, show = TRUE)
{
  if (!missing(label)) {
    gtkMenuItemNewWithLabel(label, show)
  }
  else {
    gtkMenuItemNew(show)
  }
}

gtkMenuToolButton <- gtkMenuToolButtonNew

gtkMessageDialog <-
function(parent, flags, type, buttons, ..., show = TRUE)
{
  gtkMessageDialogNew(parent, flags, type, buttons, ..., show = show)
}

gtkNotebook <- gtkNotebookNew

gtkOptionMenu <- gtkOptionMenuNew

gtkPixmap <- gtkPixmapNew

gtkPlug <- gtkPlugNew

gtkPreview <- gtkPreviewNew

gtkProgressBar <-
function(adjustment = NULL, show = TRUE)
{
  if (!missing(adjustment)) {
    gtkProgressBarNewWithAdjustment(adjustment, show)
  }
  else {
    gtkProgressBarNew(show)
  }
}

gtkRadioAction <- gtkRadioActionNew

gtkRadioButton <-
function(group = NULL, label, show = TRUE)
{
  if (!missing(label)) {
    gtkRadioButtonNewWithLabel(group, label, show)
  }
  else {
    gtkRadioButtonNew(group, show)
  }
}

gtkRadioMenuItem <-
function(group = NULL, label, show = TRUE)
{
  if (!missing(label)) {
    gtkRadioMenuItemNewWithLabel(group, label, show)
  }
  else {
    gtkRadioMenuItemNew(group, show)
  }
}

gtkRadioToolButton <-
function(group = NULL, stock.id, show = TRUE)
{
  if (!missing(stock.id)) {
    gtkRadioToolButtonNewFromStock(group, stock.id, show)
  }
  else {
    gtkRadioToolButtonNew(group, show)
  }
}

gtkScrolledWindow <- gtkScrolledWindowNew

gtkSeparatorMenuItem <- gtkSeparatorMenuItemNew

gtkSeparatorToolItem <- gtkSeparatorToolItemNew

gtkSizeGroup <- gtkSizeGroupNew

gtkSocket <- gtkSocketNew

gtkSpinButton <-
function(adjustment = NULL, climb.rate = NULL, digits = NULL, min, max, step, show = TRUE)
{
  if (!missing(adjustment)) {
    gtkSpinButtonNew(adjustment, climb.rate, digits, show)
  }
  else {
    gtkSpinButtonNewWithRange(min, max, step, show)
  }
}

gtkStatusbar <- gtkStatusbarNew

gtkStyle <- gtkStyleNew

gtkTable <- gtkTableNew

gtkTearoffMenuItem <- gtkTearoffMenuItemNew

gtkTextBuffer <- gtkTextBufferNew

gtkTextChildAnchor <- gtkTextChildAnchorNew

gtkTextMark <- gtkTextMarkNew

gtkTextTag <- gtkTextTagNew

gtkTextTagTable <- gtkTextTagTableNew

gtkTextView <-
function(buffer = NULL, show = TRUE)
{
  if (!missing(buffer)) {
    gtkTextViewNewWithBuffer(buffer, show)
  }
  else {
    gtkTextViewNew(show)
  }
}

gtkTipsQuery <- gtkTipsQueryNew

gtkToggleAction <- gtkToggleActionNew

gtkToggleButton <-
function(label, show = TRUE)
{
  if (!missing(label)) {
    gtkToggleButtonNewWithLabel(label, show)
  }
  else {
    gtkToggleButtonNew(show)
  }
}

gtkToggleToolButton <- gtkToggleToolButtonNew

gtkToolbar <- gtkToolbarNew

gtkToolButton <-
function(icon.widget = NULL, label = NULL, stock.id, show = TRUE)
{
  if (!missing(icon.widget)) {
    gtkToolButtonNew(icon.widget, label, show)
  }
  else {
    gtkToolButtonNewFromStock(stock.id, show)
  }
}

gtkToolItem <- gtkToolItemNew

gtkTooltips <- gtkTooltipsNew

gtkTreeModelFilter <- gtkTreeModelFilterNew

gtkTreeModelSort <- gtkTreeModelSortNewWithModel

gtkTreeStore <-
function(..., types)
{
  if (!missing(...)) {
    gtkTreeStoreNew(...)
  }
  else {
    gtkTreeStoreNewv(types)
  }
}

gtkTreeView <-
function(model = NULL, show = TRUE)
{
  if (!missing(model)) {
    gtkTreeViewNewWithModel(model, show)
  }
  else {
    gtkTreeViewNew(show)
  }
}

gtkTreeViewColumn <-
function(title, cell, ...)
{
  if (!missing(title)) {
    gtkTreeViewColumnNewWithAttributes(title, cell, ...)
  }
  else {
    gtkTreeViewColumnNew()
  }
}

gtkUIManager <- gtkUIManagerNew

gtkVBox <- gtkVBoxNew

gtkVButtonBox <- gtkVButtonBoxNew

gtkViewport <- gtkViewportNew

gtkVPaned <- gtkVPanedNew

gtkVRuler <- gtkVRulerNew

gtkVScale <-
function(adjustment = NULL, min, max, step, show = TRUE)
{
  if (!missing(adjustment)) {
    gtkVScaleNew(adjustment, show)
  }
  else {
    gtkVScaleNewWithRange(min, max, step, show)
  }
}

gtkVScrollbar <- gtkVScrollbarNew

gtkVSeparator <- gtkVSeparatorNew

gtkWidget <- gtkWidgetNew

gtkWindow <- gtkWindowNew

gtkWindowGroup <- gtkWindowGroupNew

gtkCellRendererAccel <- gtkCellRendererAccelNew

gtkCellRendererSpin <- gtkCellRendererSpinNew

gtkPageSetup <- gtkPageSetupNew

gtkPrintOperation <- gtkPrintOperationNew

gtkPrintSettings <- gtkPrintSettingsNew

gtkRecentFilter <- gtkRecentFilterNew

gtkRecentManager <- gtkRecentManagerNew

gtkStatusIcon <-
function(icon)
{
  if (!missing(icon)) {
    gtkStatusIconNewFromGicon(icon)
  }
  else {
    gtkStatusIconNew()
  }
}

gtkRecentChooserMenu <- gtkRecentChooserMenuNewForManager

gtkLinkButton <- gtkLinkButtonNewWithLabel

gtkRecentChooserWidget <- gtkRecentChooserWidgetNewForManager

gtkRecentChooserDialog <- gtkRecentChooserDialogNew

gtkAssistant <- gtkAssistantNew

gtkBuilder <- gtkBuilderNew

gtkRecentAction <- gtkRecentActionNew

gtkScaleButton <- gtkScaleButtonNew

gtkVolumeButton <- gtkVolumeButtonNew

gtkMountOperation <- gtkMountOperationNew

gtkEntryBuffer <- gtkEntryBufferNew

gtkInfoBar <-
function(first.button.text, ..., show = TRUE)
{
  if (!missing(first.button.text)) {
    gtkInfoBarNewWithButtons(first.button.text, ..., show = show)
  }
  else {
    gtkInfoBarNew(show)
  }
}

gtkToolItemGroup <- gtkToolItemGroupNew

gtkToolPalette <- gtkToolPaletteNew

gtkCellRendererSpinner <- gtkCellRendererSpinnerNew

gtkOffscreenWindow <- gtkOffscreenWindowNew

gtkSpinner <- gtkSpinnerNew

