# @author "HB"
setConstructorS3("AromaPipeline", function(dataSet=NULL, ..., .class=NULL) {
  # Argument 'dataSet':
  if (!is.null(dataSet)) {
    if (!is.null(.class)) {
      dataSet <- Arguments$getInstanceOf(dataSet, .class);
    }
  }

  extend(Object(), "AromaPipeline",
    .dataSet = dataSet,
    .result = NULL
  );
})


setMethodS3("getSteps", "AromaPipeline", abstract=TRUE);

setMethodS3("nbrOfSteps", "AromaPipeline", function(this, ...) {
  steps <- getSteps(this);
  length(steps);
})


setMethodS3("process", "AromaPipeline", function(this, ..., verbose=FALSE) {
  assertAnnotationData(this);

  ds <- getInputDataSet(this);
  steps <- getSteps(this);
  for (kk in seq_along(steps)) {
    step <- steps[[kk]];
    ds <- step(ds, verbose=verbose);
  }

  this$.result <- ds;

  ds;
})


##############################################################################
# HISTORY:
# 2009-12-21
# o Added class AromaPipeline.
# o Created.
##############################################################################
