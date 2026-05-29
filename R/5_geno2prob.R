library(here)
library(broman) # for %wnin%

g <- readRDS(here("Hapfiles", "called_geno.rds"))

strains <- rownames(g)

# create map again
mname <- colnames(g)
mname_spl <- strsplit(mname, "@")
chr <- sapply(mname_spl, "[", 1)
pos <- as.numeric(sapply(mname_spl, "[", 2))

map <- split(pos, factor(chr, c(1:19,"X")))
names(map) <- c(1:19, "X")
for(i in seq_along(map)) names(map[[i]]) <- paste0(names(map)[i], "@", map[[i]])

# array to hold genoprobs
pr <- setNames(vector("list", 20), c(1:19,"X"))
attr(pr, "crosstype") <- "risib"
attr(pr, "is_x_chr") <- setNames( c(rep(FALSE, 19), TRUE), c(1:19,"X"))
attr(pr, "alleles") <- LETTERS[1:8]
attr(pr, "alleleprobs") <- TRUE
attr(pr, "class") <- c("calc_genoprob", "list")

for(i in 1:19) {
    thisg <- g[,names(map[[i]])]

    pr[[i]] <- array(dim=c(length(strains), 8, length(map[[i]])))
    dimnames(pr[[i]]) <- list(strains, LETTERS[1:8], names(map[[i]]))

    for(a in LETTERS[1:8]) {
        nota <- LETTERS[1:8] %wnin% a
        for(b in LETTERS[1:8]) {
            ab <- paste0(a, b)
            notb <- LETTERS[1:8] %wnin% b
            notab <- LETTERS[1:8] %wnin% c(a,b)

            if(a==b) {
                for(k in 1:nrow(thisg)) {
                    pr[[i]][k,a,thisg[k,]==ab] <- 1
                    pr[[i]][k,nota,thisg[k,]==ab] <- 0
                }
            } else {
                for(k in 1:nrow(thisg)) {
                    pr[[i]][k,a,thisg[k,]==ab] <- 0.5
                    pr[[i]][k,b,thisg[k,]==ab] <- 0.5
                    pr[[i]][k,notab,thisg[k,]==ab] <- 0
                }
            }

        } # loop over 2nd allele
    } # loop over 1st allele
} # loop over chromosome


i <- 20
thisg <- g[,names(map[[i]])]

pr[[i]] <- array(dim=c(length(strains), 8, length(map[[i]])))
dimnames(pr[[i]]) <- list(strains, LETTERS[1:8], names(map[[i]]))

for(a in LETTERS[1:8]) {
    nota <- LETTERS[1:8] %wnin% a
    for(k in 1:nrow(thisg)) {
        pr[[i]][k,a,thisg[k,]==a] <- 1
        pr[[i]][k,nota,thisg[k,]==a] <- 0
    }
} # loop over allele


# fill in missing values by linear interpolation
for(chr in 1:20) {
    pos <- map[[chr]]
    for(ind in 1:nrow(pr[[chr]])) {
        if(any(is.na(pr[[chr]][ind,1,]))) {
            wh <- which(is.na(pr[[chr]][ind,1,]))
            whnot <- which(!is.na(pr[[chr]][ind,1,]))
            left <- sapply(wh, function(x) max(whnot[whnot < x]))
            right <- sapply(wh, function(x) min(whnot[whnot > x]))

            wts <- sapply(seq_along(wh), function(i)
                (pos[wh[i]] - pos[left[i]])/(pos[right[i]]-pos[left[i]]) )


            for(k in seq(along=wh)) {
                pr[[chr]][ind,,wh[k]] <- pr[[chr]][ind,,left[k]]*(1-wts[k]) + pr[[chr]][ind,,right[k]]*wts[k]
            }

        }
    }
}
map <- lapply(map, function(a) a/1e6)

save(pr, map, file="cc_genoprob.RData")
