# Copyright 2009-2014 Meik Michalke <meik.michalke@hhu.de>
#
# This file is part of the R package klausuR.
#
# klausuR is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# klausuR is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with klausuR.  If not, see <http://www.gnu.org/licenses/>.


# internal package description
pckg.dscrptn <- data.frame(
    Package="klausuR",
    Type="Package",
    Title="Multiple Choice Test Evaluation",
    Author="m.eik michalke",
    AuthorsR="c(person(given=\"Meik\", family=\"Michalke\", email=\"meik.michalke@hhu.de\", role=c(\"aut\", \"cre\")))",
    Maintainer="m.eik michalke <meik.michalke@hhu.de>",
    Depends="R (>= 2.9.0), xtable, psychometric, polycor, methods, graphics, tools",
    Enhances="rkward",
    Description="A set of functions designed to quickly generate results of
            a multiple choice test. Generates detailed global results, lists
            for anonymous feedback and personalised result feedback (in LaTeX
            and/or PDF format), as well as item statistics like Cronbach's alpha or
            disciminatory power. klausuR also includes a plugin for the R GUI and
            IDE RKWard, providing dialogs for basic features. To use them, install
            RKWard from http://rkward.sf.net (plugins are detected automatically).",
    License="GPL (>= 3)",
    Encoding="UTF-8",
    LazyLoad="yes",
    URL="http://r-forge.r-project.org/projects/klausur",
    stringsAsFactors=FALSE)
