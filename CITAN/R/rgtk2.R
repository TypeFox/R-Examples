## This file is part of the CITAN package for R
##
## Copyright 2011-2015 Marek Gagolewski
##
##
## CITAN is free software: you can redistribute it and/or modify
## it under the terms of the GNU Lesser General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## CITAN is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.See the
## GNU Lesser General Public License for more details.
##
## You should have received a copy of the GNU Lesser General Public License
## along with CITAN. If not, see <http://www.gnu.org/licenses/>.



# /internal/
.gtk2.selectDocuments <- function(conn, idDocuments, title="List of documents", remove=FALSE, parent=NULL)
{
   stopifnot(is.numeric(idDocuments))
   n <- length(idDocuments);
   if (remove) stopifnot(n > 1);

   if (n>2) {
      window <- .gtk2.progressBar(0, n, info="Preparing data...");
   } else window <- NULL;

   info <- list();
   for (i in 1:n)
   {
      info[[i]] <- lbsGetInfoDocuments(conn, idDocuments[i])[[1]];
      if (!is.null(window)) .gtk2.progressBar(i, n, window=window);
   }


   if (remove)
   {
      dialog <- gtkDialogNewWithButtons(NULL, parent, 0,
         "Remove selected", 1,
         "Do nothing", 0,
         "gtk-cancel", GtkResponseType["reject"], show=FALSE);
   } else
   {
         dialog <- gtkDialogNewWithButtons(NULL, parent, 0,
         "gtk-ok", GtkResponseType["ok"], show=FALSE);
   }
   dialog$setDefaultSize(900, 400);
   dialog$setTitle(title);

   data <- data.frame(
      "Del"       =FALSE,
      "Id"        =sapply(info, function(x) x$IdDocument),
      "Title"     =sapply(info, function(x) x$Title),
      "BibEntry"  =sapply(info, function(x) x$BibEntry),
      "Year"      =sapply(info, function(x) x$Year),
      "Type"      =sapply(info, function(x) x$Type),
      "Cit."      =sapply(info, function(x) x$Citations),
      "Authors"   =sapply(info, function(x) {
         paste(sapply(x$Authors,
            function(y) paste(y$Name, y$IdAuthor, sep="/")
         ),
         collapse=", ")
      }),
      "UniqueId"  =sapply(info, function(x) x$UniqueId)
   );

   if (!remove) data <- data[,-1];

   model <- rGtkDataFrame(data);
   tree_view <- gtkTreeView(model);

   sapply(1:ncol(model), function(j) {
      if (j == 1 && remove)
      {
         renderer <- gtkCellRendererToggle();
         renderer[["radio"]] <- FALSE;
         column <- gtkTreeViewColumn(colnames(model)[j], renderer, active = j-1)
         tree_view$appendColumn(column)
         gSignalConnect(renderer, "toggled", function(widget, path)
         {
            model[as.integer(path)+1][1]$Del <<- !model[as.integer(path)+1][1]$Del;
         })
      } else {
         renderer <- gtkCellRendererText();
         renderer[["wrap-width"]] <- 250;
         renderer[["wrap-mode"]] <- GtkWrapMode["word"];
         column <- gtkTreeViewColumn(colnames(model)[j], renderer, text = j-1)
         tree_view$appendColumn(column)
      }
   })
   if (is.null(gtkCheckVersion(2, 10, 0))) tree_view$setGridLines("both");



   swin <- gtkScrolledWindow()
   swin$add(tree_view);
   dialog[["vbox"]]$add(swin);

   doneVal <- dialog$run();

   dialog$destroy();

   data <- as.data.frame(model);


   if (doneVal==1) {
      return(data$Id[data$Del]);
   } else if (doneVal<0) {
      return(NULL);
   } else return(integer(0));
}


# /internal/
.gtk2.progressBar <- function(i, M, each=max(1,floor(M/100)), info=NULL, window=NULL)
{
   strTimeLeft <- "Estimated time left: %g secs.";

   if (!is.null(info))
   {
      stopifnot(is.null(window));

      window <- gtkWindowNew(NULL, FALSE);
      window$setDefaultSize(250, 60);
      window$setTitle("Operation progress");
      window$setModal(TRUE);

      box <- gtkVBox(FALSE);

      lblTitle <- gtkLabelNew(info);
      barProgress <- gtkProgressBarNew();
      gtkProgressBarSetText(barProgress, "0%");
      gtkProgressBarSetFraction(barProgress, 0.0);
      lblTimeLeft <- gtkLabelNew(sprintf(strTimeLeft, NA));

      box$add(lblTitle);
      box$add(barProgress);
      box$add(lblTimeLeft);

      gObjectSetData(window, "TimeStart", Sys.time());
      gObjectSetData(window, "lblTitle", lblTitle);
      gObjectSetData(window, "barProgress", barProgress);
      gObjectSetData(window, "lblTimeLeft", lblTimeLeft);

      window$add(box);
      window$setResizable(FALSE);
      window$show();

      return(window);
   } else
   {
      stopifnot(!is.null(window));
   }

   if (i==M)
   {
      window$destroy();
      return(NULL);
   } else if (i %% each == 0L)
   {
      TimeStart   <- gObjectGetData(window, "TimeStart");
      lblTitle    <- gObjectGetData(window, "lblTitle");
      barProgress <- gObjectGetData(window, "barProgress");
      lblTimeLeft <- gObjectGetData(window, "lblTimeLeft");

      gtkWidgetQueueDraw(lblTitle);
      gtkProgressBarSetFraction(barProgress, i/M);
      gtkProgressBarSetText(barProgress, sprintf("%g%%", floor(i/M*100)));
      gtkLabelSetText(lblTimeLeft, sprintf(strTimeLeft,
         round(as.numeric(Sys.time()-TimeStart, units="secs")*(M-i)/i, 0)
      ));

      warn <- getOption("warn");
      options("warn"=-1);
      gtkWidgetDraw(window, c(0,0,0,0)); # must update!
      options("warn"=warn);
   }

   return(window);
}
