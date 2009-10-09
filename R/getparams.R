getparams <-
function (d) 
{
    p <- list()
    p[["ns"]] <- length(levels(d$batch))
    p[["nc"]] <- range(d$cn)[2]
    p[["nind"]] <- dim(d)[1]/p[["nc"]]
    maxLike <- tapply(d$pr, d$subject, max)
    p[["lnL"]] <- sum(maxLike + log(tapply(exp(d$pr - maxLike[d$subject]), 
        FUN = sum, d$subject)))
    p[["alpha"]] <- matrix(0, nrow = p$nc, ncol = p$ns)
    p[["mean"]] <- matrix(0, nrow = p$nc, ncol = p$ns)
    p[["var"]] <- matrix(0, nrow = p$nc, ncol = p$ns)
    p[["nu"]] <- matrix(0, nrow = p$nc, ncol = p$ns)
    p[["pdc"]] <- matrix(0, nrow = p$nc, ncol = p$ns)
    lev <- levels(d$batch)
    for (j in 1:p$ns) {
        for (i in 1:p$nc) {
            p$mean[i, j] <- .Call("get_first_match", nrow(d), 
                as.numeric(d$cn), as.numeric(d$batch), as.numeric(i), 
                as.numeric(j), d$mean)
            p$alpha[i, j] <- .Call("get_first_match", nrow(d), 
                as.numeric(d$cn), as.numeric(d$batch), as.numeric(i), 
                as.numeric(j), d$alpha)
            p$var[i, j] <- .Call("get_first_match", nrow(d), 
                as.numeric(d$cn), as.numeric(d$batch), as.numeric(i), 
                as.numeric(j), d$var)
            p$nu[i, j] <- .Call("get_first_match", nrow(d), as.numeric(d$cn), 
                as.numeric(d$batch), as.numeric(i), as.numeric(j), 
                d$nu)
            p$pdc[i, j] <- .Call("get_first_match", nrow(d), 
                as.numeric(d$cn), as.numeric(d$batch), as.numeric(i), 
                as.numeric(j), d$pdc)
        }
    }
    return(p)
}
