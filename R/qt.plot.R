qt.plot <-
function (DataFrame.list, main = "", hist.or.dens = "histogram") 
{
    posterior.H0 <- DataFrame.list$posterior.H0
    LR <- -2 * (DataFrame.list$model.H0$lnL - DataFrame.list$model.H1$lnL)
    cohort.list <- unique(as.character(posterior.H0$batch))
    ncohort <- length(cohort.list)
    ncomp <- max(posterior.H0$cn)
    if (ncohort != 1) {
        stop("QT plotting supported for only one batch\n")
    }
    par(mfrow = c(1, 3))
    onecohort <- subset(posterior.H0, posterior.H0$batch == cohort.list[1])
    fit <- lm(posterior.H0$trait ~ posterior.H0$signal)
    p <- summary(fit)$coefficients[2, 4]
    plot(posterior.H0$signal, posterior.H0$trait, xlab = "Signal", 
        ylab = "Trait", main = paste("Signal vs Trait : p = ", 
            signif(p, 3)))
    abline(fit, col = "blue")
    fit <- lm(posterior.H0$trait ~ posterior.H0$cn)
    p <- summary(fit)$coefficients[2, 4]
    plot(posterior.H0$cn, posterior.H0$trait, xlab = "MAP (H0)", 
        ylab = "Trait", main = paste("MAP vs Trait : p = ", signif(p, 
            3)))
    abline(fit, col = "blue")
    loc.main <- paste("Full model : LR = ", signif(LR, 3))
    if (hist.or.dens == "density") {
        dens <- density(onecohort$signal)
        plot(dens, col = "red", xlim = range(onecohort$signal), 
            xlab = "Signal", ylab = "", main = loc.main)
        my.max <- max(dens$y)
    }
    if (hist.or.dens == "histogram") {
        my.hist <- hist(onecohort$signal, col = "red", xlim = range(onecohort$signal), 
            xlab = "Signal", ylab = "", main = loc.main, breaks = 30)
        my.max <- max(my.hist$counts)
    }
    onecohort <- onecohort[order(onecohort$signal), ]
    col <- 1
    for (i in c(1:ncomp)) {
        lines(onecohort$signal, my.max * onecohort[, paste("P", 
            i, sep = "")], type = "l", col = col)
        col <- col + 2
    }
    if (main != "") 
        title(main = main, outer = TRUE, line = -1)
}
