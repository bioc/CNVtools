cnv.plot <-
function (posterior, hist.or.dens = "histogram", batch = NULL, 
    freq = NULL, ...) 
{
    if (!is.null(batch)) 
        posterior <- posterior[posterior$batch %in% batch, ]
    posterior <- posterior[order(posterior$signal), ]
    if (hist.or.dens == "density") {
        dens <- density(posterior$signal)
        plot(dens, ...)
        my.max <- max(dens$y)
    }
    if (hist.or.dens == "histogram") {
        my.hist <- hist(posterior$signal, freq = freq, ...)
        my.max <- max(my.hist$counts)
        if (!is.null(freq)) {
            if (freq == FALSE) 
                my.max <- max(my.hist$density)
        }
    }
    col <- 1
    ncomp <- max(posterior$cn)
    for (i in c(1:ncomp)) {
        lines(posterior$signal, my.max * posterior[, paste("P", 
            i, sep = "")], type = "l", col = col)
        col <- col + 2
    }
}
