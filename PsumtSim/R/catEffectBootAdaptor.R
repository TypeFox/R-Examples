catEffectBootAdaptor<-function (df, index, testFnc = sumSqCat, useResp = TRUE, ...) {
    if (useResp) respVal <- df$resp
    else respVal <- df$bkg
    testFnc(respVal[index], df$cat, ...)
}
