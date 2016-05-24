##
##  s e g m e n t . R  Segment Functions
##


segm_intersect <- function(s1, s2) {
    stopifnot(is.numeric(s1), nrow(s1) == 2, ncol(s1) == 2,
              is.numeric(s2), nrow(s2) == 2, ncol(s2) == 2)

    bb <- function(s) {  # bounding box coordinates
        matrix(c(min(s[,1]), max(s[,1]), min(s[,2]), max(s[,2])), 2, 2)
    }
    # compute bounding boxes
    bb1 <- bb(s1)
    bb2 <- bb(s2)

    # bounding boxes do not intersect
    if (! all(rbind(bb1[2,], bb2[2,]) >= rbind(bb2[1,], bb1[1,])) )
        return(FALSE)

    # bounding boxes are intersecting
    p1 <- s1[1,]; p2 <- s1[2,]
    p3 <- s2[1,]; p4 <- s2[2,]
    sgn1 <- sign(cross(c(p3-p1, 0), c(p2-p1, 0))[3]) *
            sign(cross(c(p4-p1, 0), c(p2-p1, 0))[3])
    sgn2 <- sign(cross(c(p1-p3, 0), c(p4-p3, 0))[3]) *
            sign(cross(c(p2-p3, 0), c(p4-p3, 0))[3])

    if (sgn1 <= 0 && sgn2 <= 0) TRUE else FALSE
}


segm_distance <- function(p1, p2, p3, p4 = c()) {
    stopifnot(is.numeric(p1),  is.numeric(p2),  is.numeric(p3),
              length(p1) == 2, length(p2) == 2, length(p3) == 2)

    edist <- function(p1, p2) sqrt((p1[1]-p2[1])^2 + (p1[2]-p2[2])^2)

    if (is.null(p4)) {
        if (edist(p1, p2) == 0) return(list(d = edist(p1, p3), p = p1))

        p21 <- p2 - p1
        # det(A) = 0 only if p1 = p2
        A <- matrix(c(-p21[2], -p21[1],
                       p21[1], -p21[2]), 2, 2, byrow=TRUE)
        b <- as.matrix(p1 - p3)
        a <- solve(A, b)[2]     # crossing point at p1 + a*(p2-p1)

        if (a >= 0 && a <= 1) p <- p1 + a*(p2-p1)
        else if (a < 0)       p <- p1
        else                  p <- p2

        return(list(d = edist(p, p3), p = p, q = p3))

    } else if (is.numeric(p4) && length(p4) == 2) {
        if (segm_intersect(rbind(p1, p2), rbind(p3, p4)) &&
            cross(c(p2-p1, 0), c(p4-p3, 0))[3] != 0) {
            A <- cbind(p2 - p1, p4 - p3)
            b <- (p3 - p1)
            a <- solve(A, b)

            return(list(d = 0, p = p1+a[1]*(p2-p1), q = p3-a[2]*(p4-p3)))

        } else {
            P <- list(p3, p4, p1, p2)
            S <- list(segm_distance(p1, p2, p3), segm_distance(p1, p2, p4),
                      segm_distance(p3, p4, p1), segm_distance(p3, p4, p2))
            i <- which.min(lapply(S, function(s) s$d))

            return(list(d = S[[i]]$d, p = S[[i]]$p, q = P[[i]]))
        }
        
    } else
        stop("Argument 'p4' must be NULL or a vector of length 2.")
}
