apply.ldf <-
function (full.signal, posterior) 
{
    posterior <- posterior[dimnames(full.signal)[[1]], ]
    if (sum(dimnames(full.signal[[1]]) != dimnames(posterior)[[1]]) > 
        0) 
        stop("Names are not matching")
    can <- cancor(x = full.signal, y = posterior)
    ldf <- as.numeric(full.signal %*% can$xcoef[, 1])
    return(ldf/sd(ldf))
}
