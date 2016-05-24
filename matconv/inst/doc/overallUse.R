## ----setup, include=FALSE------------------------------------------------
library(matconv)
library(knitr)

knitr::opts_chunk$set(fig.pos='center', echo=TRUE, comment='>')


matCode <- c("function [out] = csvReadPretty(csvPath, hd)",
 "\tfid = fopen(csvPath,'r');",
 "    parse = textscan(fid, '%s','delimiter','\\n');",
 "    parse = parse{1};",
 "    fclose(fid);",
 "    if(isempty(parse) || isempty(parse{1}))",
 "        out = parse;",
 "    end",
 "    ",
 "    for (row = (hd+1):length(parse))",
 "        line = [parse{row} ','];",
 "        commas = regexp(line,',');",
 "        col = 2;",
 "        bef = commas(1);",
 "        if bef~=1",
 "            out{row-hd,1} = line(1:bef-1);",
 "        end",
 "",
 "        for ca = commas(2:end)",
 "            %if the commas are next to each other just add col",
 "",
 "            if (bef+1 ~= ca)",
 "               out{row,col} = line(bef+1:ca-1);",
 "",
 "            end",
 "            col=col+1;",
 "            bef = ca;",
 "        end",
 "",
 "    end",
 "",
 "",
 "end")

## ----basic---------------------------------------------------------------

out <- mat2r(matCode)
names(out)

## ----functionMaps--------------------------------------------------------
hMaps <- makeFuncMaps(
	pathDict = system.file("extdata", "HiebelerDict.txt", package = "matconv"))

## ----dataConverters------------------------------------------------------
source(system.file("extdata", "defDataConv.R", package = "matconv"))

length(dataConvs)

## ----finish--------------------------------------------------------------
out <- mat2r(matCode, funcConverters = hMaps, dataConverters = dataConvs)


