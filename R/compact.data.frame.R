compact.data.frame <-
function (full.frame) 
{
    full.frame <- full.frame[order(full.frame$cn), ]
    full.frame.mod <- do.call(rbind.data.frame, split(x = full.frame$posterior, 
        f = full.frame$subject))
    names(full.frame.mod) <- paste("P", c(1:(dim(full.frame.mod)[2])), 
        sep = "")
    full.frame.mod$cn <- apply(full.frame.mod, FUN = which.max, 
        MAR = 1)
    full.frame.mod$subject <- row.names(full.frame.mod)
    full.frame <- subset(full.frame[, c("subject", "batch", "signal", 
        "trait")], full.frame$cn == 1)
    full.frame <- merge(full.frame, full.frame.mod)
    return(full.frame)
}
