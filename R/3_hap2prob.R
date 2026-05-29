# convert haplotype files to a set of genotype probabilities
#
# 1. read in the haplotypes (ignore chr Y and M)
#    turn into a data frame: chromosome, mat/pat, allele, start, end
# 2. for each chromosome:
#    a. find the unique positions
#    b. put markers at the endpoints and half-way between
#    c. determine allele on each chromosome at each marker
#    d. what to do with the missing pieces? linear interpolation
# 3. genoprob to alleleprobs

library(here)

# read files
dir <- here("Hapfiles")
files <- list.files(dir, pattern=".hap$")
hap <- lapply(files, function(file) readLines(here(dir, file)))
names(hap) <- sapply(strsplit(files, "-"), "[", 1)

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


# pull out strain, long name, Y and Mt genotypes
hap2df <-
    function(h)
{

    short_strain <- strsplit(h[2], '["-]')[[1]][2]
    long_strain <- strsplit(h[2], '["]')[[1]][2]

    yline <- grep("^chrY", h, value=TRUE)
    mline <- grep("^chrM", h, value=TRUE)

    yallele <- strsplit(yline, ",")[[1]][3]
    mallele <- strsplit(mline, ",")[[1]][3]

    data.frame(strain=short_strain,
               long_strain=long_strain,
               Y=yallele, M=mallele)
}

tab <- lapply(hap, hap2df)
tab <- do.call("rbind", tab)
write.table(tab, here("strains_info.csv"), quote=FALSE, sep=",", row.names=FALSE)

haptab <- lapply(hap, hap2table)
haptab <- do.call("rbind", haptab)

# find unique positions on each chromosome
# put a marker between each pair of positions
hap_upos <-
    function(ht)
{
    allchr <- c(1:19,"X")
    results <- setNames(vector("list", length(allchr)), allchr)

    for(chr in allchr) {
        this <- ht[ht$chr==chr,c("start","end")]
        results[[chr]] <- sort(unique(c(this[,1], this[,2], round((this[,1]+this[,2])/2))))
        names(results[[chr]]) <- paste0(chr, "@", results[[chr]])
    }
    results
}

upos <- hap_upos(haptab)

# get genotype at each map position
hap2geno <-
    function(ht, upos)
{
    strains <- unique(ht$strain)
    pos <- unlist(lapply(upos, names))

    result <- matrix("", nrow=length(strains), ncol=length(pos))
    dimnames(result) <- list(strains, pos)

    for(i in 1:nrow(ht)) {
        if(i==round(i,-1)) cat(i, " of ", nrow(ht), "\n")
        up <- upos[[ ht$chr[i] ]]
        col <- names(up[up >= ht$start[i] & up <= ht$end[i]])
        for(j in col) result[ht$strain[i], j] <- paste0(result[ht$strain[i], j], ht$allele[i])
    }

    result
}

# get genotype calls
# g <- hap2geno(haptab, upos)
