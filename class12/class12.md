# Class 12 Lab Session
Malibu Slattery (A1844012)

- [Proportion of G/G in a population](#proportion-of-gg-in-a-population)
- [The variant that is associated with childhood asthma is more frequent
  in GBR population than in the MXL population….let’s explore this
  further…](#the-variant-that-is-associated-with-childhood-asthma-is-more-frequent-in-gbr-population-than-in-the-mxl-populationlets-explore-this-further)

## Proportion of G/G in a population

``` r
mxl<-read.csv("rs8067378.csv")
head(mxl)
```

      Sample..Male.Female.Unknown. Genotype..forward.strand. Population.s. Father
    1                  NA19648 (F)                       A|A ALL, AMR, MXL      -
    2                  NA19649 (M)                       G|G ALL, AMR, MXL      -
    3                  NA19651 (F)                       A|A ALL, AMR, MXL      -
    4                  NA19652 (M)                       G|G ALL, AMR, MXL      -
    5                  NA19654 (F)                       G|G ALL, AMR, MXL      -
    6                  NA19655 (M)                       A|G ALL, AMR, MXL      -
      Mother
    1      -
    2      -
    3      -
    4      -
    5      -
    6      -

``` r
table(mxl$Genotype..forward.strand.)
```


    A|A A|G G|A G|G 
     22  21  12   9 

``` r
table(mxl$Genotype..forward.strand.)/nrow(mxl)*100
```


        A|A     A|G     G|A     G|G 
    34.3750 32.8125 18.7500 14.0625 

``` r
gbr<-read.csv("gbr_rs8067378.csv")
head(gbr)
```

      Sample..Male.Female.Unknown. Genotype..forward.strand. Population.s. Father
    1                  HG00096 (M)                       A|A ALL, EUR, GBR      -
    2                  HG00097 (F)                       G|A ALL, EUR, GBR      -
    3                  HG00099 (F)                       G|G ALL, EUR, GBR      -
    4                  HG00100 (F)                       A|A ALL, EUR, GBR      -
    5                  HG00101 (M)                       A|A ALL, EUR, GBR      -
    6                  HG00102 (F)                       A|A ALL, EUR, GBR      -
      Mother
    1      -
    2      -
    3      -
    4      -
    5      -
    6      -

Find proportion of G\|G in GBR population

``` r
table(gbr$Genotype..forward.strand.)
```


    A|A A|G G|A G|G 
     23  17  24  27 

``` r
round(table(gbr$Genotype..forward.strand.)/nrow(gbr)*100,2)
```


      A|A   A|G   G|A   G|G 
    25.27 18.68 26.37 29.67 

## The variant that is associated with childhood asthma is more frequent in GBR population than in the MXL population….let’s explore this further…

``` r
expr<-read.table("rs8067378_ENSG00000172057.6.txt")
head(expr)
```

       sample geno      exp
    1 HG00367  A/G 28.96038
    2 NA20768  A/G 20.24449
    3 HG00361  A/A 31.32628
    4 HG00135  A/A 34.11169
    5 NA18870  G/G 18.25141
    6 NA11993  A/A 32.89721

``` r
nrow(expr)
```

    [1] 462

``` r
table(expr$geno)
```


    A/A A/G G/G 
    108 233 121 

``` r
library(ggplot2)
```

> Q 13 sample size: A/A (108, median = 31.25), A/G (233, 25.06) G/G
> (121, 20.07).

``` r
#I didn't really know how to use the median function that was being talked about, so I used the filter function we've been working with in datacamp.
library(dplyr)
```


    Attaching package: 'dplyr'

    The following objects are masked from 'package:stats':

        filter, lag

    The following objects are masked from 'package:base':

        intersect, setdiff, setequal, union

``` r
AA<-expr %>% filter(geno=="A/A")
head(AA)
```

        sample geno      exp
    3  HG00361  A/A 31.32628
    4  HG00135  A/A 34.11169
    6  NA11993  A/A 32.89721
    8  NA18498  A/A 47.64556
    13 NA20585  A/A 30.71355
    15 HG00235  A/A 25.44983

``` r
round(median(AA$exp),2)
```

    [1] 31.25

``` r
AG<-expr %>% filter(geno=="A/G")
head(AG)
```

        sample geno      exp
    1  HG00367  A/G 28.96038
    2  NA20768  A/G 20.24449
    7  HG00256  A/G 31.48736
    10 HG00115  A/G 33.85374
    11 NA20806  A/G 16.29854
    12 HG00278  A/G 19.73450

``` r
round(median(AG$exp),2)
```

    [1] 25.06

``` r
GG<-expr %>% filter(geno=="G/G")
head(GG)
```

        sample geno      exp
    5  NA18870  G/G 18.25141
    9  HG00327  G/G 17.67473
    17 NA12546  G/G 18.55622
    20 NA18488  G/G 23.10383
    23 NA19214  G/G 30.94554
    28 HG00112  G/G 21.14387

``` r
round(median(GG$exp),2)
```

    [1] 20.07

> Q.14 Generate a boxplot with a box per genotype, what could you infer
> from the relative expression value between A/A and G/G displayed in
> this plot? Does the SNP effect the expression of ORMDL3?

``` r
ggplot(expr, aes(x=geno, y=  exp, fill = geno)) + geom_boxplot(notch=TRUE) 
```

![](class12_files/figure-commonmark/unnamed-chunk-11-1.png)

> ans: A/A is more associated with the expression of ORMDL3 and G/G is
> less associated with the expression of ORMDL3. Thus, you are more
> likely to express ORMDL3 if you have the former, rather than the
> latter.
