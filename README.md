## Collaborative Cross genotype information

Genotype information on the [Collaborative Cross](https://doi.org/10.1534/g3.111.001891) mouse strains
is available in several forms, but not immediately useful for analysis
with [R/qtl2](https://kbroman.org/qtl2).

GigaMUGA data on many strains are provided in standard R/qtl2 input at
the [qtl2data repository on github](https://github.com/rqtl/qtl2data/tree/main/CC), but not all strains are
represented.

Here, we converted the `.hap` files at
<https://csbio.unc.edu/CCstatus/CCGenomes/> to a genotype probability
array suitable for use with R/qtl2.

**Note**: There were a bunch of cases where
tha haplotype files were internally inconsistent, with sets of
overlapping haplotypes. If the end of one segment overlapped the start
of the next, we shifted it back so that they didn't overlap.


### Usage

The genotype probabilities as well as the map of locations is in
[`cc_genoprob.RData`](cc_genoprob.RData).

Load them into R by first downloading the file and then using `load()`.

```r
download.file("https://github.com/kbroman/CCgeno/raw/refs/heads/main/cc_genoprob.RData",
              "cc_genoprob.RData")
load("cc_genoprob.RData")
```

The objects are named `pr` (a list of 3-dimensional arrays, as
produced by
[`calc_genoprob()`](https://cran.r-project.org/web/packages/qtl2/refman/qtl2.html#calc_genoprob))
and `map` (a list of vectors of marker positions).


The [`strains_info.csv`](strains_info.csv) file contains the long
version of the strain names, as well as the allele on the Y chromosome
and the mitochondria.


### License

Code licensed under the [MIT license](License.md).
