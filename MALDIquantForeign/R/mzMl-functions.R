## Copyright 2013-2015 Sebastian Gibb
## <mail@sebastiangibb.de>
##
## This file is part of MALDIquantForeign for R and related languages.
##
## MALDIquantForeign is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## MALDIquantForeign is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with MALDIquantForeign. If not, see <http://www.gnu.org/licenses/>

.writeMzMlDocument <- function(x, file, id, encoding="utf-8", imsArgs=list()) {
  ## stop if file isn't writeable
  if (file.exists(file) && file.access(file, 2) != 0) {
    stop("No permissions to write into ", sQuote(file), "!")
  }

  isIms <- length(imsArgs)

  ## file handle
  f <- file(file, open="wt", encoding=encoding)

  .writeXmlHeader(file=f, encoding=encoding)

  .writeXmlTag("mzML", attrs=c(xmlns="http://psi.hupo.org/ms/mzml",
    "xmlns:xsi"="http://www.w3.org/2001/XMLSchema-instance",
    "xsi:schemaLocation"=paste("http://psi.hupo.org/ms/mzml",
      "http://psidev.info/files/ms/mzML/xsd/mzML1.1.0.xsd"),
    id=.sanitize(ifelse(missing(id) || is.null(id),
                        deparse(substitute(x)), id)),
    version="1.1.0"), close=FALSE, file=f)

  .writeMzMlCvList(file=f, isIms=isIms)
  .writeMzMlFileDescription(x, file=f, isIms=isIms, imsArgs=imsArgs)
  .writeMzMlSoftwareList(x, file=f)
  if (isIms) {
    .writeImzMlReferenceableParamGroups(x, file=f)
    .writeImzMlScanSettings(x, file=f)
  }
  .writeMzMlInstrumentConfigurationList(x, file=f)
  .writeMzMlDataProcessingList(x, file=f)
  .writeMzMlRun(x, file=f, isIms=isIms)

  .writeCloseXmlTag("mzML", file=f)

  invisible(close(f))
}

.writeMzMlCvList <- function(file, isIms=FALSE) {
  items <- list(
    ms=list(id="MS",
      fullName="Proteomics Standards Initiative Mass Spectrometry Ontology",
      version="3.44.0",
      URI="http://psidev.cvs.sourceforge.net/*checkout*/psidev/psi/psi-ms/mzML/controlledVocabulary/psi-ms.obo"),
    uo=list(id="UO",
      fullName="Unit Ontology",
      version="12:10:2012",
      URI="http://obo.cvs.sourceforge.net/*checkout*/obo/obo/ontology/phenotype/unit.obo"))

  if (isIms) {
    items$imzml <- list(id="IMS",
                        fullName="Imaging MS Ontology",
                        version="0.9.1",
                        URI="http://www.maldi-msi.org/download/imzml/imagingMS.obo")
  }

  .writeXmlTag("cvList", attrs=c(count=2), intend=1, close=FALSE, file=file)
  for (i in seq(along=items)) {
    .writeXmlTag("cv", attrs=items[[i]], intend=2, file=file)
  }
  .writeCloseXmlTag("cvList", intend=1, file=file)
}

.writeMzMlFileDescription <- function(x, file, isIms=FALSE, imsArgs) {
  .writeXmlTag("fileDescription", intend=1, close=FALSE, file=file)
    .writeXmlTag("fileContent", intend=2, close=FALSE, file=file)
      .writeXmlTag("cvParam", intend=3,
                   attrs=c(cvRef="MS", accession="MS:1000579",
                           name="MS1 spectrum"), file=file)
      if (isIms) {
        .writeXmlTag("cvParam", intend=3,
                     attrs=c(cvRef="IMS", accession="IMS:1000080",
                             name="universally unique identifier",
                             value=paste0("{", imsArgs$uuid, "}")),
                     file=file)
        .writeXmlTag("cvParam", intend=3,
                     attrs=c(cvRef="IMS", accession="IMS:1000091",
                             name="ibd SHA-1", value=imsArgs$sha1),
                     file=file)

        if (imsArgs$processed) {
          .writeXmlTag("cvParam", intend=3,
                       attrs=c(cvRef="IMS", accession="IMS:1000031",
                               name="processed"), file=file)
        } else {
          .writeXmlTag("cvParam", intend=3,
                       attrs=c(cvRef="IMS", accession="IMS:1000030",
                               name="continuous"), file=file)
        }
      }

      .writeXmlTag("userParam", intend=3,
                   attrs=c(name="MALDIquantForeign",
                           value="MALDIquant object(s) exported to mzML"),
                   file=file)
    .writeCloseXmlTag("fileContent", intend=2, file=file)
    .writeMzMlSourceFileList(x, file=file)
  .writeCloseXmlTag("fileDescription", intend=1, file=file)
}

.writeMzMlSourceFileList <- function(x, file) {
  files <- unique(unlist(lapply(x, function(s)metaData(s)$file)))

  if (length(files)) {
    dname <- dirname(files)
    bname <- basename(files)
    ext <- tolower(.fileExtension(bname))

    .writeXmlTag("sourceFileList", attrs=c(count=length(files)), intend=2,
                 close=FALSE, file=file)
    for (i in seq(along=files)) {
      .writeXmlTag("sourceFile", intend=3,
                   attrs=c(id=paste0("SF", i), location=dname[i],
                           name=bname[i]), close=FALSE, file=file)

      if (ext[i] == "fid") {
        .writeXmlTag("cvParam", intend=4,
                     attrs=c(cvRef="MS", accession="MS:1000825",
                             name="Bruker FID file"), file=file)
        .writeXmlTag("cvParam", intend=4,
                     attrs=c(cvRef="MS", accession="MS:1000773",
                             name="Bruker FID nativeID format"), file=file)
      } else if (ext[i] == "mzxml") {
        .writeXmlTag("cvParam", intend=4,
                     attrs=c(cvRef="MS", accession="MS:1000566",
                             name="ISB mzXML file"), file=file)
      } else if (ext[i] == "mzml") {
        .writeXmlTag("cvParam", intend=4,
                     attrs=c(cvRef="MS", accession="MS:1000584",
                             name="mzML file"), file=file)
      }

      if (file.exists(files[i])) {
        .writeXmlTag("cvParam", intend=4,
                     attrs=c(cvRef="MS", accession="MS:1000569", name="SHA-1",
                             value=digest::digest(files[i], algo="sha1",
                                                  file=TRUE)), file=file)
      }
      .writeCloseXmlTag("sourceFile", intend=3, file=file)
    }
    .writeCloseXmlTag("sourceFileList", intend=2, file=file)
  }
}

.writeMzMlSoftwareList <- function(x, file) {
  .writeXmlTag("softwareList", attrs=c(count=1), intend=1, close=FALSE,
               file=file)
    .writeXmlTag("software", intend=2, attrs=c(id="MALDIquantForeign",
                  version=as.character(packageVersion("MALDIquantForeign"))),
                 file=file)
  .writeCloseXmlTag("softwareList", intend=1, file=file)
}

.writeMzMlInstrumentConfigurationList <- function(x, file) {
  .writeXmlTag("instrumentConfigurationList", attrs=c(count=1), intend=1,
               close=FALSE, file=file)
    .writeXmlTag("instrumentConfiguration", attrs=c(id="IC0"), intend=2,
                 file=file)
  .writeCloseXmlTag("instrumentConfigurationList", intend=1, file=file)
}

.writeMzMlDataProcessingList <- function(x, file) {
  .writeXmlTag("dataProcessingList", attrs=c(count=1), intend=1, close=FALSE,
               file=file)
    .writeXmlTag("dataProcessing", attrs=c(id="export"), intend=2, close=FALSE,
                 file=file)
      .writeXmlTag("processingMethod", intend=3,
                   attrs=c(order=1, softwareRef="MALDIquantForeign"),
                   close=FALSE, file=file)
        .writeXmlTag("userParam", intend=4,
                     attrs=c(name="MALDIquant object(s) exported to mzML",
                             value=""), file=file)
      .writeCloseXmlTag("processingMethod", intend=3, file=file)
    .writeCloseXmlTag("dataProcessing", intend=2, file=file)
  .writeCloseXmlTag("dataProcessingList", intend=1, file=file)
}

.writeMzMlRun <- function(x, file, isIms=FALSE) {
  .writeXmlTag("run", attrs=c(id="run0",
                              defaultInstrumentConfigurationRef="IC0"),
                intend=1, close=FALSE, file=file)
  .writeMzMlSpectrumList(x, file=file, isIms=isIms)
  .writeCloseXmlTag("run", intend=1, file=file)
}

.writeMzMlSpectrumList <- function(x, file, isIms=FALSE) {
  .writeXmlTag("spectrumList",
               attrs=c(count=length(x), defaultDataProcessingRef="export"),
               intend=2, close=FALSE, file=file)

  for (i in seq(along=x)) {
    id <- ifelse(is.null(metaData(x[[i]])$id), paste0("scan=", i-1L),
                         metaData(x[[i]])$id)

    .writeXmlTag("spectrum", intend=3,
                 attrs=c(index=i-1, id=id, defaultArrayLength=length(x[[i]]),
                         spotID=metaData(x[[i]])$fullName), close=FALSE,
                 file=file)

      msLevel <- ifelse(is.null(metaData(x[[i]])$msLevel), 1,
                        metaData(x[[i]])$msLevel)

      .writeXmlTag("cvParam", intend=4, attrs=c(cvRef="MS",
                     accession="MS:1000511", name="ms level", value=msLevel),
                   file=file)
      .writeXmlTag("cvParam", intend=4, attrs=c(cvRef="MS",
                     accession="MS:1000294", name="mass spectrum"), file=file)
      .writeXmlTag("cvParam", intend=4, attrs=c(cvRef="MS",
                     accession="MS:1000528", name="lowest observed m/z",
                     value=min(mass(x[[i]])), unitCvRef="MS",
                     unitAccession="MS:1000040", unitName="m/z"), file=file)
      .writeXmlTag("cvParam", intend=4, attrs=c(cvRef="MS",
                     accession="MS:1000527", name="highest observed m/z",
                     value=max(mass(x[[i]])), unitCvRef="MS",
                     unitAccession="MS:1000040", unitName="m/z"), file=file)
      .writeXmlTag("cvParam", intend=4, attrs=c(cvRef="MS",
                     accession="MS:1000285", name="total ion current",
                     value=totalIonCurrent(x[[i]])), file=file)

    if (MALDIquant::isMassSpectrum(x[[i]])) {
      .writeXmlTag("cvParam", intend=4, attrs=c(cvRef="MS",
                     accession="MS:1000128", name="profile spectrum"),
                   file=file)
    } else {
      .writeXmlTag("cvParam", intend=4, attrs=c(cvRef="MS",
                     accession="MS:1000127", name="centroid spectrum"),
                   file=file)
    }

    if (isIms) {
      .writeImzMlScanList(x[[i]], file=file)
      .writeImzMlBinaryDataArrayList(x[[i]], file=file)
    } else {
      .writeMzMlBinaryDataArrayList(x[[i]], file=file)
    }

    .writeCloseXmlTag("spectrum", intend=3, file=file)
  }
  .writeCloseXmlTag("spectrumList", intend=2, file=file)
}

.writeMzMlBinaryDataArrayList <- function(x, file) {
  count <- ifelse(MALDIquant::isMassSpectrum(x), 2, 3)
  .writeXmlTag("binaryDataArrayList", attrs=c(count=count), intend=4,
               close=FALSE, file=file)
  .writeMzMlBinaryData(mass(x), file=file, c(cvRef="MS", accession="MS:1000514",
                        name="m/z array", unitCvRef="MS",
                        unitAccession="MS:1000040", unitName="m/z"))

  .writeMzMlBinaryData(intensity(x), file=file, c(cvRef="MS",
                        accession="MS:1000515", name="intensity array",
                        unitCvRef="MS", unitAccession="MS:1000131",
                        unitName="number of counts"))

  if (MALDIquant::isMassPeaks(x)) {
    .writeMzMlBinaryData(snr(x), file=file, c(cvRef="MS",
                          accession="MS:1000517", name="signal to noise array"))
  }
  .writeCloseXmlTag("binaryDataArrayList", intend=4, file=file)
}

.writeMzMlBinaryData <- function(x, file, additionalAttrs) {
  binaryData <- .base64encode(x, size=8, endian="little",
                              compressionType="gzip")

  .writeXmlTag("binaryDataArray", attrs=c(encodedLength=nchar(binaryData)),
                intend=5, close=FALSE, file=file)
    .writeXmlTag("cvParam", intend=6, file=file,
                 attrs=c(cvRef="MS", accession="MS:1000574",
                         name="zlib compression"))
    .writeXmlTag("cvParam", intend=6, file=file,
                 attrs=c(cvRef="MS", accession="MS:1000523",
                         name="64-bit float"))
    if (!missing(additionalAttrs)) {
      .writeXmlTag("cvParam", attrs=additionalAttrs, intend=6, file=file)
    }
    .writeXmlTag("binary", text=binaryData, intend=6, file=file)
  .writeCloseXmlTag("binaryDataArray", intend=5, file=file)
}

.writeImzMlReferenceableParamGroups <- function(x, file) {
  .writeXmlTag("referenceableParamGroupList", attrs=c(count=2), intend=1,
               close=FALSE, file=file)

    .writeXmlTag("referenceableParamGroup", attrs=c(id="mzArray"), intend=2,
                 close=FALSE, file=file)

      .writeXmlTag("cvParam", intend=3, file=file,
                   attrs=c(cvRef="MS", accession="MS:1000514",
                           name="m/z array", unitCvRef="MS",
                           unitAccession="MS:1000040", unitName="m/z"))

      ref <- c("MS", "MS", "IMS")
      accession <- paste(ref, c(1000576, 1000523, 1000101), sep=":")
      name <- c("no compression", "64-bit float", "external data")
      value <- c("", "", "true")

      for (i in seq(along=accession)) {
        .writeXmlTag("cvParam", intend=3, file=file,
                     attrs=c(cvRef=ref[i], accession=accession[i], name=name[i],
                             value=value[i]))
      }

    .writeCloseXmlTag("referenceableParamGroup", intend=2, file=file)

    .writeXmlTag("referenceableParamGroup", attrs=c(id="intensityArray"),
                intend=2, close=FALSE, file=file)

      .writeXmlTag("cvParam", intend=3, file=file,
                   attrs=c(cvRef="MS", accession="MS:1000515",
                           name="intensity array", unitCvRef="MS",
                           unitAccession="MS:1000131"))

      for (i in seq(along=accession)) {
        .writeXmlTag("cvParam", intend=3, file=file,
                     attrs=c(cvRef=ref[i], accession=accession[i], name=name[i],
                             value=value[i]))
      }
    .writeCloseXmlTag("referenceableParamGroup", intend=2, file=file)
  .writeCloseXmlTag("referenceableParamGroupList", intend=1, file=file)
}

.writeImzMlScanSettings <- function(x, file) {
  .writeXmlTag("scanSettingsList", attrs=c(count=1), intend=1, close=FALSE,
               file=file)
    .writeXmlTag("scanSettings", attrs=c(id="scansetting1"), intend=2,
                 close=FALSE, file=file)

    accession <- paste("IMS", 1000042:1000043, sep=":")
    name <- paste("max count of pixel", c("x", "y"))
    value <- unname(metaData(x[[1L]])$imaging$size)

    for (i in seq(along=accession)) {
      .writeXmlTag("cvParam", intend=3, file=file,
                   attrs=c(cvRef="IMS", accession=accession[i], name=name[i],
                           value=value[i]))
    }

    accession <- paste("IMS", 1000044:1000047, sep=":")
    name <- c(paste("max dimension", c("x", "y")),
              paste("pixel size", c("x", "y")))
    value <- unname(c(metaData(x[[1L]])$imaging$dim,
                      metaData(x[[1L]])$imaging$pixelSize))

    for (i in seq(along=accession)) {
      .writeXmlTag("cvParam", intend=3, file=file,
                   attrs=c(cvRef="IMS", accession=accession[i], name=name[i],
                           value=value[i],
                           unitCvRef="UO", unitAccession="UO:0000017",
                           unitName="micrometer"))
    }
    .writeCloseXmlTag("scanSettings", intend=2, file=file)
  .writeCloseXmlTag("scanSettingsList", intend=1, file=file)
}

.writeImzMlScanList <- function(x, file) {
  .writeXmlTag("scanList", attrs=c(count=1), intend=4, close=FALSE, file=file)
    .writeXmlTag("scan", intend=5, close=FALSE, file=file)
    accession <- paste("IMS", 1000050:1000051, sep=":")
    name <- paste("position", c("x", "y"))
    value <- unname(metaData(x)$imaging$pos)

    for (i in seq(along=accession)) {
      .writeXmlTag("cvParam", intend=6, file=file,
                   attrs=c(cvRef="IMS", accession=accession[i], name=name[i],
                           value=value[i]))
    }
    .writeCloseXmlTag("scan", intend=5, file=file)
  .writeCloseXmlTag("scanList", intend=4, file=file)
}

.writeImzMlBinaryDataArrayList <- function(x, file) {
  .writeXmlTag("binaryDataArrayList", attrs=c(count=2), intend=4, close=FALSE,
               file=file)
  .writeImzMlBinaryData(metaData(x)$imaging$offsets["mass",], file=file,
                        ref="mzArray")
  .writeImzMlBinaryData(metaData(x)$imaging$offsets["intensity",], file=file,
                        ref="intensityArray")
  .writeCloseXmlTag("binaryDataArrayList", intend=4, file=file)
}

.writeImzMlBinaryData <- function(x, file, ref) {
  .writeXmlTag("binaryDataArray", attrs=c(encodedLength=0), intend=5,
               close=FALSE, file=file)

    .writeXmlTag("referenceableParamGroupRef", attrs=c(ref=ref),
                 intend=6, file=file)

    accession <- paste("IMS", 1000102:1000104, sep=":")
    name <- paste("external", c("offset", "array length", "encoded length"))
    value <- unname(x)

    for (i in seq(along=accession)) {
      .writeXmlTag("cvParam", intend=6, file=file,
                   attrs=c(cvRef="IMS", accession=accession[i], name=name[i],
                           value=value[i]))
    }
    .writeXmlTag("binary", intend=6, file=file)
  .writeCloseXmlTag("binaryDataArray", intend=5, file=file)
}
