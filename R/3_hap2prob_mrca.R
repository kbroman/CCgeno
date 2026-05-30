# convert haplotype files to a matrix of genotype calls
# this is for the MRCA files
# no M or Y chromosome information

library(here)

# read files
dir <- here("Hapfiles", "MRCA")
files <- list.files(dir, pattern=".hap$")
hap <- lapply(files, function(file) readLines(here(dir, file)))
names(hap) <- sapply(strsplit(files, "_"), "[", 1)

# haplotype file to data frame
hap2table <-
    function(h)
{
    strain <- strsplit(h[2], '["-]')[[1]][2]

    h <- grep("^chr[1-9X]", h, value=TRUE)

    result <- NULL
    for(i in seq_along(h)) {
        x <- strsplit(h[i], ",")[[1]]
        chr <- sub("chr", "", x[1])
        allele <- x[seq(3, length(x), by=3)]
        start <- as.numeric(x[seq(4, length(x), by=3)])
        end <- as.numeric(x[seq(5, length(x), by=3)])
        if((length(x) - 2) %% 3 != 0)
            warning("strain ", strain, " chr ", chr, " has unexpected number of values")

        if(length(start) > 1 && (any(diff(start)<0) || any(diff(end)<0) || any(start > end) || any(end[-length(end)] >= start[-1]))) {
            warning("out of order in strain ", strain, " chr ", chr)
        }

        if(length(start) > 1) {
            for(i in 1:(length(start)-1)) {
                if(end[i] >= start[i+1]) end[i] <- start[i+1]-1
            }
        }

        while(length(start) > 1 && (any(diff(start)<0) || any(diff(end)<0) || any(start > end) || any(end[-length(end)] >= start[-1]))) {
            omit <- NULL
            for(i in 1:(length(start)-1)) {
                if(end[i] >= start[i+1]) end[i] <- start[i+1]-1
                if(start[i] > end[i]) omit <- c(omit, i)
            }
            if(length(omit) > 0) {
                allele <- allele[-omit]
                start <- start[-omit]
                end <- end[-omit]
            }
        }

        (matpat <- ifelse((i %% 2)==1, "mat", "pat"))

        this <- data.frame(strain=rep(strain, length(allele)),
                           chr=rep(chr, length(allele)),
                           matpat=rep(matpat, length(allele)),
                           allele=allele,
                           start=as.numeric(start),
                           end=as.numeric(end))
        if(is.null(result)) result <- this
        else result <- rbind(result, this)
    }

    result
}


haptab <- lapply(hap, hap2table)
haptab <- do.call("rbind", haptab)
