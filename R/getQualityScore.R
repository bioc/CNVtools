getQualityScore <-
function (posterior) 
{
    Qvec <- c()
    for (coh in levels(posterior$batch)) {
        posterior2 <- subset(posterior, posterior$batch == coh)
        sds <- tapply(posterior2$signal, FUN = sd, INDEX = posterior2$cn)
        means <- tapply(posterior2$signal, FUN = mean, INDEX = posterior2$cn)
        freq <- table(posterior2$cn)/length(posterior2$cn)
        sds <- subset(sds, table(posterior2$cn) > 4)
        means <- subset(means, table(posterior2$cn) > 4)
        freq <- subset(freq, table(posterior2$cn) > 4)
        l <- length(means)
        if (l == 1) {
            return(NA)
        }
        else {
            dmeans <- abs(means[1:(l - 1)] - means[2:l])
            av.sds <- (freq[1:(l - 1)] * sds[1:(l - 1)] + freq[2:l] * 
                sds[2:l])/(freq[1:(l - 1)] + freq[2:l])
            weights <- freq[1:(l - 1)] * freq[2:l]
            Q <- sum(weights * dmeans/av.sds)/sum(weights)
            Qvec <- append(Qvec, values = Q)
        }
    }
    return(min(Qvec))
}
