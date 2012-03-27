test.posterior <-
function (frame, ncomp, samples.by.disease = NULL) 
{
    bad.fit <- FALSE
    comp <- c(1:ncomp)
    col <- paste("P", comp, sep = "")
    frame$average <- apply(MARGIN = 1, frame[, col, drop = FALSE], 
        FUN = function(x) {
            sum(x * comp)
        })
    for (my.id in levels(frame$batch)) {
        small.frame <- subset(frame, frame$batch == my.id)
        if (is.null(samples.by.disease)) {
            small.frame <- small.frame[order(small.frame$signal), 
                ]
            di <- diff(small.frame$average)
            should.not.be <- subset(di, di < 0)
            if (sum(should.not.be) < -0.01) 
                bad.fit <- TRUE
        }
        else {
            for (t in c(1:2)) {
                smaller.frame <- subset(small.frame, small.frame$subject %in% 
                  samples.by.disease[[t]])
                smaller.frame <- smaller.frame[order(smaller.frame$signal), 
                  ]
                di <- diff(smaller.frame$average)
                should.not.be <- subset(di, di < 0)
                if (sum(should.not.be) < -0.01) 
                  bad.fit <- TRUE
            }
        }
    }
    return(bad.fit)
}
