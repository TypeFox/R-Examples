# @author "HB"
setConstructorS3("ASCRMAv2", function(...) {
  extend(AromaPipeline(..., .class="AffymetrixCelSet"), "ASCRMAv2");
})


setMethodS3("assertAnnotationData", "ASCRMAv2", function(this, ...) {
  ds <- getInputDataSet(this);

  # Assert that the CDF file exists
  cdf <- getCdf(ds);

  # Assert that an UGP annotation data file exists
  gi <- getGenomeInformation(cdf);

  # Assert that an UFL annotation data file exists
  si <- getSnpInformation(cdf);

  # Assert than an ACS (probe-sequence) annotation files
  acs <- AromaCellSequenceFile$byChipType(getChipType(cdf, fullname=FALSE));
});


setMethodS3("getSteps", "ASCRMAv2", function(this, ...) {
  list(
    "acc" = function(csR) {
      # Allelic cross-talk calibration tests
      acc <- AllelicCrosstalkCalibration(csR, model="CRMAv2");
      print(acc);
      csC <- process(acc, verbose=log);
      print(csC);
      csC;
    },

    "bpn" = function(csC) {
      # Base-position normalization
      bpn <- BasePositionNormalization(csC, target="zero");
      print(bpn);
      csN <- process(bpn, verbose=log);
      print(csN);
      csN;
    },

    "avg" = function(csN) {
      # Allele-specific probe summarization
      plm <- AvgCnPlm(csN, mergeStrands=TRUE, combineAlleles=FALSE);
      print(plm);
      if (length(findUnitsTodo(plm)) > 0) {
        # Fit CN probes quickly (~5-10s/array + some overhead)
        units <- fitCnProbes(plm, verbose=log);
        str(units);
        units <- fit(plm, verbose=log);
        str(units);
      }
      ces <- getChipEffectSet(plm);
      print(ces);
      ces;
    },

    "fln" = function(ces) {
      # Fragment-length normalization test
      fln <- FragmentLengthNormalization(ces, target="zero");
      print(fln);
      cesN <- process(fln, verbose=log);
      print(cesN);
      cesN;
    }
  );
})



##############################################################################
# HISTORY:
# 2009-12-21
# o Added class ASCRMAv2.
# o Created.
##############################################################################
