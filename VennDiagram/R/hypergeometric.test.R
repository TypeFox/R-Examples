# The VennDiagram package is copyright (c) 2012 Ontario Institute for Cancer Research (OICR)
# This package and its accompanying libraries is free software; you can redistribute it and/or modify it under the terms of the GPL
# (either version 1, or at your option, any later version) or the Artistic License 2.0.  Refer to LICENSE for the full license text.
# OICR makes no representations whatsoever as to the SOFTWARE contained herein.  It is experimental in nature and is provided WITHOUT
# WARRANTY OF MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE OR ANY OTHER WARRANTY, EXPRESS OR IMPLIED. OICR MAKES NO REPRESENTATION
# OR WARRANTY THAT THE USE OF THIS SOFTWARE WILL NOT INFRINGE ANY PATENT OR OTHER PROPRIETARY RIGHT.
# By downloading this SOFTWARE, your Institution hereby indemnifies OICR against any loss, claim, damage or liability, of whatsoever kind or
# nature, which may arise from your Institution's respective use, handling or storage of the SOFTWARE.
# If publications result from research using this SOFTWARE, we ask that the Ontario Institute for Cancer Research be acknowledged and/or
# credit be given to OICR scientists, as scientifically appropriate.

#This function performs the hypergeometric test on the two categories. Taken from package BoutrosLab.statistics.general
calculate.overlap.and.pvalue = function(list1, list2, total.size, lower.tail = TRUE, adjust = FALSE) {

        # calculate actual overlap
        actual.overlap <- length(intersect(list1, list2));

        # calculate expected overlap
        # need to cast to avoid integer overflow when length(list1) * length(list2) is extremely large
        expected.overlap <- as.numeric(length(list1)) * length(list2) / total.size;

        adjust.value <- 0;

        # adjust actual.overlap to reflect P[X >= x]
        if (adjust & !lower.tail) {
                adjust.value <- 1;
                warning('Calculating P[X >= x]');
                }

        # calculate significance of the overlap
        overlap.pvalue <- phyper(
                q = actual.overlap - adjust.value,
                m = length(list1),
                n = total.size - length(list1),
                k = length(list2),
                lower.tail = lower.tail
                );

        # return values
        return( c(actual.overlap, expected.overlap, overlap.pvalue) );

        }

