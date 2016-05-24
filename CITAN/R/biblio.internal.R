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
# stops if connection is invalid or dead
.lbsCheckConnection <- function(conn)
{
   if (!class(conn) == "SQLiteConnection")
      stop("incorrect 'conn' given");

   dbGetInfo(conn); # check if conn is active
}


# /internal/
.lbs_SourceTypesFull  <- c("Book Series", "Conference Proceedings", "Journal");

# /internal/
.lbs_SourceTypesShort <- c("'bs'",        "'cp'",                   "'jo'");


# /internal/
.lbs_DocumentTypesFull  <- c("Article", "Article in Press", "Book", "Conference Paper", "Editorial", "Erratum", "Letter", "Note", "Report", "Review", "Short Survey");

# /internal/
.lbs_DocumentTypesShort <- c("'ar'",    "'ip'",             "'bk'", "'cp'",             "'ed'",      "'er'",    "'le'",   "'no'", "'rp'",   "'re'",   "'sh'");

# /internal/
.lbs_DocumentType_ShortToFull <- function(type)
{
   was.factor <- is.factor(type);
   type <- as.factor(type);
   lev <- sprintf("'%s'", levels(type));
   for (i in 1:length(.lbs_DocumentTypesFull))
   {
      lev[lev==(.lbs_DocumentTypesShort[i])] <- .lbs_DocumentTypesFull[i];
   }
   levels(type) <- lev;
   if (!was.factor) type <- as.character(type);
   return(type);
}


# /internal/
.lbs_PrepareRestriction_DocumentTypes <- function(conn, documentTypes)
{
   if (is.null(documentTypes)) return(NULL);

   if (!is.character(documentTypes))
      stop("incorrect 'documentTypes' given");

   documentTypesShort <- character(length(documentTypes));
   for (i in 1:length(documentTypes))
      documentTypesShort[i] <- sqlSwitchOrNULL(documentTypes[i],
         .lbs_DocumentTypesFull,
         .lbs_DocumentTypesShort
      );

   incorrect <- which(documentTypesShort == "NULL");

   if (length(incorrect)>0)
   {
      warning(sprintf("incorrect document types: %s. Ignoring.",
         paste(documentTypes[incorrect], collapse=", ")));
      documentTypesShort <- documentTypesShort[-incorrect];
   }

   if (length(documentTypesShort) == 0) stop("all given document types were incorrect.");

   return(documentTypesShort);
}



# /internal/
.lbs_PrepareRestriction_SurveyDescription <- function(conn, surveyDescription)
{
   if (is.null(surveyDescription)) return(NULL);

   if (!is.character(surveyDescription) || length(surveyDescription)!=1)
      stop("incorrect 'surveyDescription' given");

   surveyDescription <- sqlEscapeTrim(surveyDescription);

   res <- dbGetQuery(conn, sprintf("SELECT * FROM Biblio_Surveys WHERE Description='%s';",
      surveyDescription));
   if (nrow(res) == 0) stop("Survey not found.");

   return(surveyDescription);
}
