opts= tdmOptsDefaultsSet();
tdm = list(umode="RSUB");
data(iris);
dataObj=tdmSplitTestData(opts,tdm,dset=iris);
