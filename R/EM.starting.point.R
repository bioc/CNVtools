EM.starting.point <-
function (d, trait = "binary") 
{
    ncohort <- length(levels(d$batch))
    ncomp <- range(d$cn)[2]
    d$alpha <- c(1/ncomp)
    d$pr <- c(-10)
    d$posterior <- c(0)
    if (trait == "binary") {
        Ncontrol <- length(d[d$cn == 1 & d$d == "0", ]$subject)
        Ncase <- length(d[d$cn == 1 & d$d == "1", ]$subject)
        if (Ncase == 0) {
            d$pdc <- 0
        }
        else {
            d$pdc <- Ncase/(Ncontrol + Ncase)
        }
    }
    else {
        d$pdc <- mean(d$trait)
        d$rs <- sd(d$trait)
    }
    lev <- levels(d$batch)
    for (my.batch in lev) {
        ends <- range(d$signal)
        r <- diff(ends)
        pf <- 0.2
        if (ncomp > 1) {
            interval <- r/(ncomp - 1)
        }
        else {
            interval <- r
        }
        mean <- seq(from = ends[1], by = interval, length.out = ncomp)
        mean <- mean + c(runif(1:ncomp, -1 * pf * r, +1 * pf * 
            r))
        vars <- pmax(runif(1:ncomp, 0, diff(range(d$signal)/(6 * 
            ncomp))), 0.1)
        nus <- runif(ncomp, 5, 15)
        mean <- sort(mean)
        d$mean <- ifelse(d$batch == my.batch, mean[d$cn], d$mean)
        d$var <- ifelse(d$batch == my.batch, vars[d$cn], d$var)
        d$nu <- ifelse(d$batch == my.batch, nus[d$cn], d$nu)
    }
    return(d)
}
