"comval" <-
function() {
        ipsi <- integer(1)
        c <- single(1)
        h1 <- single(1)
        h2 <- single(1)
        h3 <- single(1)
        xk <- single(1)
        d <- single(1)
        bta <- single(1)
        bt0 <- single(1)
        iucv <- integer(1)
        a2 <- single(1)
        b2 <- single(1)
        chk <- single(1)
        ckw <- single(1)
        bb <- single(1)
        bt <- single(1)
        cw <- single(1)
        em <- single(1)
        cr <- single(1)
        vk <- single(1)
        np <- integer(1)
        enu <- single(1)
        v7 <- single(1)
        iwww <- integer(1)
        f.res <- .Fortran("comval",
                ipsi = as.integer(ipsi),
                c = as.single(c),
                h1 = as.single(h1),
                h2 = as.single(h2),
                h3 = as.single(h3),
                xk = as.single(xk),
                d = as.single(d),
                bta = as.single(bta),
                bt0 = as.single(bt0),
                iucv = as.integer(iucv),
                a2 = as.single(a2),
                b2 = as.single(b2),
                chk = as.single(chk),
                ckw = as.single(ckw),
                bb = as.single(bb),
                bt = as.single(bt),
                cw = as.single(cw),
                em = as.single(em),
                cr = as.single(cr),
                vk = as.single(vk),
                np = as.integer(np),
                enu = as.single(enu),
                v7 = as.single(v7),
                iwww = as.integer(iwww))
        list(ipsi = f.res$ipsi, c = f.res$c, h1 = f.res$h1, h2 = f.res$h2,
                h3 = f.res$h3, xk = f.res$xk, d = f.res$d, bta = f.res$bta,
                bt0 = f.res$bt0, iucv = f.res$iucv, a2 = f.res$a2, b2 = f.res$
                b2, chk = f.res$chk, ckw = f.res$ckw, bb = f.res$bb, bt = f.res
                $bt, cw = f.res$cw, em = f.res$em, cr = f.res$cr, vk = f.res$vk,
                np = f.res$np, iwww = f.res$iwww)
}

