# bimm143_lab5datavisR
Malibu Slattery (A18488012)

## Background

There are lots of ways to make figures in R. These include so-called
“base-R” graphics. (e.g. `plot()`) and tones of add-on packages like
**ggplot2**.

For example, here we make the same plot with both:

``` r
head(cars)
```

      speed dist
    1     4    2
    2     4   10
    3     7    4
    4     7   22
    5     8   16
    6     9   10

``` r
plot(cars)
```

![](bimm143_lab5datavisR_files/figure-commonmark/unnamed-chunk-1-1.png)

``` r
library(ggplot2)
p<- ggplot(cars) + aes(speed, dist) + geom_smooth(col="green", method = "lm", se = FALSE) + xlab("Speed") + geom_point() + ylab("Distance") + labs(title ="Speed vs Distance Plot", subtitle = "by Malibu Slattery") + theme_minimal() + theme(plot.title = element_text(hjust = 0.5)) + theme(plot.subtitle =  element_text(hjust = 0.5))+ theme(plot.title = element_text(face="bold"))
p
```

    `geom_smooth()` using formula = 'y ~ x'

![](bimm143_lab5datavisR_files/figure-commonmark/unnamed-chunk-2-1.png)

## Gene expression plot

``` r
library(ggplot2)
url <- "https://bioboot.github.io/bimm143_S20/class-material/up_down_expression.txt"
genes <- read.delim(url)
head(genes)
```

            Gene Condition1 Condition2      State
    1      A4GNT -3.6808610 -3.4401355 unchanging
    2       AAAS  4.5479580  4.3864126 unchanging
    3      AASDH  3.7190695  3.4787276 unchanging
    4       AATF  5.0784720  5.0151916 unchanging
    5       AATK  0.4711421  0.5598642 unchanging
    6 AB015752.4 -3.6808610 -3.5921390 unchanging

``` r
nrow(genes)
```

    [1] 5196

``` r
ncol(genes)
```

    [1] 4

``` r
table(genes$State)
```


          down unchanging         up 
            72       4997        127 

``` r
round( table(genes$State)/nrow(genes) * 100, 2 )
```


          down unchanging         up 
          1.39      96.17       2.44 

``` r
ggplot(genes) +
  aes(Condition1, Condition2) +
  geom_point() + xlab("Condition 1") + ylab("Condition 2")
```

![](bimm143_lab5datavisR_files/figure-commonmark/unnamed-chunk-4-1.png)

``` r
#genes$State
```

Condition 2: Let’s color by `State`, so we can see the “up” and “down”
significant genes compared to all the “unchanging” genes..

``` r
ggplot(genes) +
  aes(Condition1, Condition2, col = State) +
  geom_point() + xlab("Condition 1") + ylab("Condition 2")
```

![](bimm143_lab5datavisR_files/figure-commonmark/unnamed-chunk-5-1.png)

``` r
ggplot(genes) +
  aes(Condition1, Condition2, col = State) +
  geom_point() + xlab("Condition 1") + ylab("Condition 2") + scale_color_manual(values = c("darkblue", "grey", "darkgreen")) + labs(title="Gene Expression Changes upon GLP-1 Drug")+ theme(plot.title = element_text(hjust = 0.5))
```

![](bimm143_lab5datavisR_files/figure-commonmark/unnamed-chunk-6-1.png)

Let’s look at the famous **gapminder**!

``` r
url <- "https://raw.githubusercontent.com/jennybc/gapminder/master/inst/extdata/gapminder.tsv"

gapminder <- read.delim(url)
head(gapminder, 3)
```

          country continent year lifeExp      pop gdpPercap
    1 Afghanistan      Asia 1952  28.801  8425333  779.4453
    2 Afghanistan      Asia 1957  30.332  9240934  820.8530
    3 Afghanistan      Asia 1962  31.997 10267083  853.1007

``` r
ggplot(gapminder) + aes(gdpPercap, lifeExp, col = continent) + geom_point(alpha = 0.3) + xlab("GDP per Capita") + ylab("Life Expectancy") + guides(color = guide_legend(title = "Continent"))
```

![](bimm143_lab5datavisR_files/figure-commonmark/unnamed-chunk-7-1.png)

Let’s “facet” (i.e. make a separate plot) by continent rather than the
big hot mess above!

``` r
url <- "https://raw.githubusercontent.com/jennybc/gapminder/master/inst/extdata/gapminder.tsv"

gapminder <- read.delim(url)
head(gapminder, 3)
```

          country continent year lifeExp      pop gdpPercap
    1 Afghanistan      Asia 1952  28.801  8425333  779.4453
    2 Afghanistan      Asia 1957  30.332  9240934  820.8530
    3 Afghanistan      Asia 1962  31.997 10267083  853.1007

``` r
ggplot(gapminder) + aes(gdpPercap, lifeExp, col = continent) + geom_point(alpha = 0.3) + xlab("GDP per Capita") + ylab("Life Expectancy") + facet_wrap(~continent) + guides(color = guide_legend(title = "Continent"))
```

![](bimm143_lab5datavisR_files/figure-commonmark/unnamed-chunk-8-1.png)
