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


#' @description
#' Creates an empty Local Bibliometric Storage.
#'
#' @details
#' The function may be executed only if the database contains no tables
#' named \code{Biblio_*} and no views named \code{ViewBiblio_*}.
#'
#' The following SQL code is executed.
#' \preformatted{
#' CREATE TABLE Biblio_Categories (\cr
#'      -- Source classification codes (e.g. ASJC)\cr
#'    IdCategory         INTEGER PRIMARY KEY ASC,\cr
#'    IdCategoryParent   INTEGER NOT NULL,\cr
#'    Description        VARCHAR(63) NOT NULL,\cr
#'    FOREIGN KEY(IdCategoryParent) REFERENCES Biblio_Categories(IdCategory)\cr
#' );
#' }
#'
#'
#' \preformatted{
#' CREATE TABLE Biblio_Sources (
#'    IdSource      INTEGER PRIMARY KEY AUTOINCREMENT,
#'    AlternativeId VARCHAR(31)  UNIQUE NOT NULL,
#'    Title         VARCHAR(255) NOT NULL,
#'    IsActive      BOOLEAN,
#'    IsOpenAccess  BOOLEAN,
#'    Type          CHAR(2) CHECK (Type IN ('bs', 'cp', 'jo')),
#'        -- Book Series / Conference Proceedings / Journal
#'        -- or NULL in all other cases
#'    Impact1        REAL, -- value of an impact factor
#'    Impact2        REAL, -- value of an impact factor
#'    Impact3        REAL, -- value of an impact factor
#'    Impact4        REAL, -- value of an impact factor
#'    Impact5        REAL, -- value of an impact factor
#'    Impact6        REAL, -- value of an impact factor
#' );
#' }
#'
#' \preformatted{
#' CREATE TABLE Biblio_SourcesCategories (
#'      -- links Sources and Categories
#'    IdSource         INTEGER NOT NULL,
#'    IdCategory       INTEGER NOT NULL,
#'    PRIMARY KEY(IdSource, IdCategory),
#'    FOREIGN KEY(IdSource)     REFERENCES Biblio_Sources(IdSource),
#'    FOREIGN KEY(IdCategory)   REFERENCES Biblio_Categories(IdCategory)
#' );
#' }
#'
#'
#' \preformatted{
#' CREATE TABLE Biblio_Documents (
#'    IdDocument     INTEGER PRIMARY KEY AUTOINCREMENT,
#'    IdSource       INTEGER,
#'    AlternativeId  VARCHAR(31) UNIQUE NOT NULL,
#'    Title          VARCHAR(255) NOT NULL,
#'    BibEntry       TEXT,
#'        -- (e.g. Source Title,Year,Volume,Issue,Article Number,PageStart,PageEnd)
#'    Year           INTEGER,
#'    Pages          INTEGER,
#'    Citations      INTEGER NOT NULL,
#'    Type           CHAR(2) CHECK (Type IN ('ar', 'ip', 'bk',
#'        'cp', 'ed', 'er', 'le', 'no', 'rp', 're', 'sh')),
#'        -- Article-ar / Article in Press-ip / Book-bk /
#'        -- Conference Paper-cp / Editorial-ed / Erratum-er /
#'        -- Letter-le/ Note-no / Report-rp / Review-re / Short Survey-sh
#'        -- or NULL in all other cases
#'    FOREIGN KEY(IdSource)   REFERENCES Biblio_Sources(IdSource),
#'    FOREIGN KEY(IdLanguage) REFERENCES Biblio_Languages(IdLanguage)
#' );
#' }
#'
#'
#'  \preformatted{
#'  CREATE TABLE Biblio_Citations (
#'     IdDocumentParent     INTEGER NOT NULL,  # cited document
#'     IdDocumentChild      INTEGER NOT NULL,  # reference
#'     PRIMARY KEY(IdDocumentParent, IdDocumentChild),
#'     FOREIGN KEY(IdDocumentParent) REFERENCES Biblio_Documents(IdDocument),
#'     FOREIGN KEY(IdDocumentChild)  REFERENCES Biblio_Documents(IdDocument)
#' );
#' }
#'
#'
#' \preformatted{
#' CREATE TABLE Biblio_Surveys (
#'      -- each call to lbsImportDocuments() puts a new record here,
#'      -- they may be grouped into so-called 'Surveys' using 'Description' field
#'    IdSurvey       INTEGER PRIMARY KEY AUTOINCREMENT,
#'    Description    VARCHAR(63) NOT NULL,   -- survey group name
#'    FileName       VARCHAR(63),            -- original file name
#'    Timestamp      DATETIME                -- date of file import
#' );
#' }
#'
#'
#'
#' \preformatted{
#' CREATE TABLE Biblio_DocumentsSurveys (
#'    -- note that the one Document may often be found in many Surveys
#'    IdDocument     INTEGER NOT NULL,
#'    IdSurvey       INTEGER NOT NULL,
#'    PRIMARY KEY(IdDocument, IdSurvey),
#'    FOREIGN KEY(IdSurvey)   REFERENCES Biblio_Surveys(IdSurvey),
#'    FOREIGN KEY(IdDocument) REFERENCES Biblio_Documents(IdDocument)
#' );
#' }
#'
#' \preformatted{
#' CREATE TABLE Biblio_Authors (
#'    IdAuthor        INTEGER PRIMARY KEY AUTOINCREMENT,
#'    Name            VARCHAR(63) NOT NULL,
#'    AuthorGroup     VARCHAR(31), # used to merge authors with non-unique representations
#' );
#' }
#'
#' \preformatted{
#' CREATE TABLE Biblio_AuthorsDocuments (
#'      -- links Authors and Documents
#'    IdAuthor        INTEGER NOT NULL,
#'    IdDocument      INTEGER NOT NULL,
#'    PRIMARY KEY(IdAuthor, IdDocument),
#'    FOREIGN KEY(IdAuthor)   REFERENCES Biblio_Authors(IdAuthor),
#'    FOREIGN KEY(IdDocument) REFERENCES Biblio_Documents(IdDocument)
#' );
#' }
#'
#' In addition, the following views are created.
#' \preformatted{
#' CREATE VIEW ViewBiblio_DocumentsSurveys AS
#'    SELECT
#'       Biblio_DocumentsSurveys.IdDocument AS IdDocument,
#'       Biblio_DocumentsSurveys.IdSurvey AS IdSurvey,
#'       Biblio_Surveys.Description AS Description,
#'       Biblio_Surveys.Filename AS Filename,
#'       Biblio_Surveys.Timestamp AS Timestamp
#'    FROM Biblio_DocumentsSurveys
#'    JOIN Biblio_Surveys
#'       ON Biblio_DocumentsSurveys.IdSurvey=Biblio_Surveys.IdSurvey;
#' }
#'
#' \preformatted{
#' CREATE VIEW ViewBiblio_DocumentsCategories AS
#' SELECT
#'       IdDocument AS IdDocument,
#'       DocSrcCat.IdCategory AS IdCategory,
#'       DocSrcCat.Description AS Description,
#'       DocSrcCat.IdCategoryParent AS IdCategoryParent,
#'       Biblio_Categories.Description AS DescriptionParent
#'    FROM
#'    (
#'       SELECT
#'          Biblio_Documents.IdDocument AS IdDocument,
#'          Biblio_SourcesCategories.IdCategory AS IdCategory,
#'          Biblio_Categories.Description AS Description,
#'          Biblio_Categories.IdCategoryParent AS IdCategoryParent
#'       FROM Biblio_Documents
#'       JOIN Biblio_SourcesCategories
#'          ON Biblio_Documents.IdSource=Biblio_SourcesCategories.IdSource
#'       JOIN Biblio_Categories
#'          ON Biblio_SourcesCategories.IdCategory=Biblio_Categories.IdCategory
#'    ) AS DocSrcCat
#'    JOIN Biblio_Categories
#'          ON DocSrcCat.IdCategoryParent=Biblio_Categories.IdCategory;
#' }
#'
#' @title Create a Local Bibliometric Storage
#' @param conn a connection object, see \code{\link{lbsConnect}}.
#' @param verbose logical; \code{TRUE} to be more verbose.
#' @examples
#' \dontrun{
#' conn <- lbsConnect("Bibliometrics.db");
#' ## ...
#' lbsCreate(conn);
#' Scopus_ImportSources(conn);
#' ## ...
#' lbsDisconnect(conn);}
#'
#' @return \code{TRUE} on success.
#'
#' @export
#'
#' @seealso
#'  \code{\link{lbsConnect}},
#'  \code{\link{lbsClear}},
#'  \code{\link{Scopus_ImportSources}},
#'  \code{\link{lbsTidy}}
lbsCreate <- function(conn, verbose=TRUE)
{
   .lbsCheckConnection(conn); # will stop on invalid/dead connection


   ## --------- auxiliary function -------------------------------------------

   #' /internal/
   .lbsCreateTable <- function(conn, tablename, query, verbose)
   {
      if (verbose) cat(sprintf("Creating table '%s'... ", tablename));

      dbExecQuery(conn, query, FALSE);

      if (verbose) cat("Done.\n");
   }


   ## --------- auxiliary function -------------------------------------------

   #' /internal/
   .lbsCreateView <- function(conn, viewname, query, verbose)
   {
      if (verbose) cat(sprintf("Creating view '%s'... ", viewname));

      dbExecQuery(conn, query, FALSE);

      if (verbose) cat("Done.\n");
   }


   ## --------- auxiliary function -------------------------------------------

   #' /internal/
   .lbsCreateIndex <- function(conn, indexname, query, verbose)
   {
      if (verbose) cat(sprintf("Creating index for '%s'... ", indexname));

      dbExecQuery(conn, query, FALSE);

      if (verbose) cat("Done.\n");
   }



   ## ---- check whether LBS is empty ------------------------------------


   tablesviews <- dbListTables(conn);
   if (any(substr(tablesviews, 1, 7) == "Biblio_"))
      stop("Your Local Bibliometric Storage is not empty.");
   if (any(substr(tablesviews, 1, 11) == "ViewBiblio_"))
      stop("Your Local Bibliometric Storage is not empty.");


   ## -------------------------------------------------------------------

   query <- "CREATE TABLE Biblio_Categories (
      IdCategory        INTEGER PRIMARY KEY ASC,
      IdCategoryParent  INTEGER NOT NULL,
      Description       VARCHAR(63) NOT NULL,
      FOREIGN KEY(IdCategoryParent) REFERENCES Biblio_Categories(IdCategory)
   );"

   .lbsCreateTable(conn, "Biblio_Categories", query, verbose);




   query <- "CREATE TABLE Biblio_Sources (
      IdSource      INTEGER PRIMARY KEY NOT NULL,
      AlternativeId VARCHAR(31)  UNIQUE NOT NULL,
      Title         VARCHAR(255) NOT NULL,
      IsActive      BOOLEAN,
      IsOpenAccess  BOOLEAN,
      Type          CHAR(2) CHECK (Type IN ('bs', 'cp', 'jo')),
      Impact1       REAL,
      Impact2       REAL,
      Impact3       REAL,
      Impact4       REAL,
      Impact5       REAL,
      Impact6       REAL
   );"

   .lbsCreateTable(conn, "Biblio_Sources", query, verbose);



   query <- "CREATE INDEX IF NOT EXISTS Biblio_Sources_Title ON Biblio_Sources (Title ASC);";

   .lbsCreateIndex(conn, "Biblio_Sources", query, verbose);



   query <- "CREATE TABLE Biblio_SourcesCategories (
      IdSource          INTEGER NOT NULL,
      IdCategory        INTEGER NOT NULL,
      PRIMARY KEY(IdSource, IdCategory),
      FOREIGN KEY(IdSource)     REFERENCES Biblio_Sources(IdSource),
      FOREIGN KEY(IdCategory)   REFERENCES Biblio_Categories(IdCategory)
   );"

   .lbsCreateTable(conn, "Biblio_SourcesCategories", query, verbose);



   query <- "CREATE TABLE Biblio_Documents (
      IdDocument     INTEGER PRIMARY KEY AUTOINCREMENT,
      IdSource       INTEGER,
      AlternativeId  VARCHAR(31)  UNIQUE NOT NULL,
      Title          VARCHAR(255),
      BibEntry       TEXT,
      Year           INTEGER,
      Pages          INTEGER,
      Citations      INTEGER NOT NULL,
      Type           CHAR(2) CHECK (Type IN ('ar', 'ip', 'bk', 'cp', 'ed', 'er', 'le', 'no', 'rp', 're', 'sh')),
      FOREIGN KEY(IdSource)   REFERENCES Biblio_Sources(IdSource)
   );"

   .lbsCreateTable(conn, "Biblio_Documents", query, verbose);




   query <- "CREATE TABLE Biblio_Citations (
      IdDocumentParent     INTEGER NOT NULL,
      IdDocumentChild      INTEGER NOT NULL,
      PRIMARY KEY(IdDocumentParent, IdDocumentChild),
      FOREIGN KEY(IdDocumentParent) REFERENCES Biblio_Documents(IdDocument),
      FOREIGN KEY(IdDocumentChild)  REFERENCES Biblio_Documents(IdDocument)
   );"

   .lbsCreateTable(conn, "Biblio_Citations", query, verbose);



   query <- "CREATE TABLE Biblio_Surveys (
      IdSurvey     INTEGER PRIMARY KEY AUTOINCREMENT,
      Description  VARCHAR(63) NOT NULL,
      FileName     VARCHAR(63),
      Timestamp    DATETIME
   );"

   .lbsCreateTable(conn, "Biblio_Surveys", query, verbose);

   query <- "CREATE TABLE Biblio_DocumentsSurveys (
      IdDocument     INTEGER NOT NULL,
      IdSurvey       INTEGER NOT NULL,
      PRIMARY KEY(IdDocument, IdSurvey),
      FOREIGN KEY(IdSurvey)   REFERENCES Biblio_Surveys(IdSurvey),
      FOREIGN KEY(IdDocument) REFERENCES Biblio_Documents(IdDocument)
   );"

   .lbsCreateTable(conn, "Biblio_DocumentsSurveys", query, verbose);




   query <- "CREATE TABLE Biblio_Authors (
      IdAuthor     INTEGER PRIMARY KEY AUTOINCREMENT,
      Name         VARCHAR(63) NOT NULL,
      AuthorGroup  VARCHAR(31)
   );"

   .lbsCreateTable(conn, "Biblio_Authors", query, verbose);






   query <- "CREATE TABLE Biblio_AuthorsDocuments (
      IdAuthor     INTEGER NOT NULL,
      IdDocument   INTEGER NOT NULL,
      PRIMARY KEY(IdAuthor, IdDocument),
      FOREIGN KEY(IdAuthor)   REFERENCES Biblio_Authors(IdAuthor),
      FOREIGN KEY(IdDocument) REFERENCES Biblio_Documents(IdDocument)
   );"

   .lbsCreateTable(conn, "Biblio_AuthorsDocuments", query, verbose);



   query <- "CREATE VIEW ViewBiblio_DocumentsSurveys AS
      SELECT
         Biblio_DocumentsSurveys.IdDocument AS IdDocument,
         Biblio_DocumentsSurveys.IdSurvey AS IdSurvey,
         Biblio_Surveys.Description AS Description,
         Biblio_Surveys.Filename AS Filename,
         Biblio_Surveys.Timestamp AS Timestamp
      FROM Biblio_DocumentsSurveys
      JOIN Biblio_Surveys ON Biblio_DocumentsSurveys.IdSurvey=Biblio_Surveys.IdSurvey;";

   .lbsCreateView(conn, "ViewBiblio_DocumentsSurveys", query, verbose);


   query <- "CREATE VIEW ViewBiblio_DocumentsCategories AS
      SELECT
         IdDocument AS IdDocument,
         DocSrcCat.IdCategory AS IdCategory,
         DocSrcCat.Description AS Description,
         DocSrcCat.IdCategoryParent AS IdCategoryParent,
         Biblio_Categories.Description AS DescriptionParent
      FROM
      (
         SELECT
            Biblio_Documents.IdDocument AS IdDocument,
            Biblio_SourcesCategories.IdCategory AS IdCategory,
            Biblio_Categories.Description AS Description,
            Biblio_Categories.IdCategoryParent AS IdCategoryParent
         FROM Biblio_Documents
         JOIN Biblio_SourcesCategories ON Biblio_Documents.IdSource=Biblio_SourcesCategories.IdSource
         JOIN Biblio_Categories ON Biblio_SourcesCategories.IdCategory=Biblio_Categories.IdCategory
      ) AS DocSrcCat
      JOIN Biblio_Categories ON DocSrcCat.IdCategoryParent=Biblio_Categories.IdCategory;";

   .lbsCreateView(conn, "ViewBiblio_DocumentsCategories", query, verbose);

   # -------------------------------------------------------------------

   if (verbose) cat("Your Local Bibliometric Storage has been created.
   Perhaps now you may wish to use Scopus_ImportSources(...) to import source information.\n");

   return(TRUE);
}
