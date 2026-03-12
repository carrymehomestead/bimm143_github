# bimm143_lab5datavisR
Malibu Slattery (A18488012)

- [Background](#background)
- [Gene expression plot](#gene-expression-plot)
- [Custom Plots](#custom-plots)

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

![](bimm143_lab5datavisR.1_files/figure-commonmark/unnamed-chunk-1-1.png)

``` r
library(ggplot2)
p<- ggplot(cars) + aes(speed, dist) + geom_smooth(col="green", method = "lm", se = FALSE) + xlab("Speed") + geom_point() + ylab("Distance") + labs(title ="Speed vs Distance Plot", subtitle = "by Malibu Slattery") + theme_minimal() + theme(plot.title = element_text(hjust = 0.5)) + theme(plot.subtitle =  element_text(hjust = 0.5))+ theme(plot.title = element_text(face="bold"))
p
```

    `geom_smooth()` using formula = 'y ~ x'

![](bimm143_lab5datavisR.1_files/figure-commonmark/unnamed-chunk-2-1.png)

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

![](bimm143_lab5datavisR.1_files/figure-commonmark/unnamed-chunk-4-1.png)

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

![](bimm143_lab5datavisR.1_files/figure-commonmark/unnamed-chunk-5-1.png)

``` r
ggplot(genes) +
  aes(Condition1, Condition2, col = State) +
  geom_point() + xlab("Condition 1") + ylab("Condition 2") + scale_color_manual(values = c("darkblue", "grey", "darkgreen")) + labs(title="Gene Expression Changes upon GLP-1 Drug")+ theme(plot.title = element_text(hjust = 0.5))
```

![](bimm143_lab5datavisR.1_files/figure-commonmark/unnamed-chunk-6-1.png)

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

![](bimm143_lab5datavisR.1_files/figure-commonmark/unnamed-chunk-7-1.png)

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

![](bimm143_lab5datavisR.1_files/figure-commonmark/unnamed-chunk-8-1.png)

## Custom Plots

How big is this gapminder dataset?

``` r
nrow(gapminder)
```

    [1] 1704

``` r
ncol(gapminder)
```

    [1] 6

``` r
library(dplyr)
```


    Attaching package: 'dplyr'

    The following objects are masked from 'package:stats':

        filter, lag

    The following objects are masked from 'package:base':

        intersect, setdiff, setequal, union

``` r
gapminder_2007<-filter(gapminder, year == 2007)
head(gapminder_2007)
```

          country continent year lifeExp      pop  gdpPercap
    1 Afghanistan      Asia 2007  43.828 31889923   974.5803
    2     Albania    Europe 2007  76.423  3600523  5937.0295
    3     Algeria    Africa 2007  72.301 33333216  6223.3675
    4      Angola    Africa 2007  42.731 12420476  4797.2313
    5   Argentina  Americas 2007  75.320 40301927 12779.3796
    6   Australia   Oceania 2007  81.235 20434176 34435.3674

``` r
filter(gapminder, year==2007, country == "Ireland")
```

      country continent year lifeExp     pop gdpPercap
    1 Ireland    Europe 2007  78.885 4109086     40676

``` r
filter(gapminder, year==1977, country == "Ireland")
```

      country continent year lifeExp     pop gdpPercap
    1 Ireland    Europe 1977   72.03 3271900  11150.98

``` r
filter(gapminder, year==1977, country == "United States")
```

            country continent year lifeExp       pop gdpPercap
    1 United States  Americas 1977   73.38 220239000  24072.63

``` r
filter(gapminder, year==1977, country == "United States")$lifeExp
```

    [1] 73.38

Make a plot comparing 1977 and 2007 for all countries!

``` r
year_1977_2007<- filter(gapminder, year %in% c(1977, 2007))
head(year_1977_2007)
```

          country continent year lifeExp      pop gdpPercap
    1 Afghanistan      Asia 1977  38.438 14880372  786.1134
    2 Afghanistan      Asia 2007  43.828 31889923  974.5803
    3     Albania    Europe 1977  68.930  2509048 3533.0039
    4     Albania    Europe 2007  76.423  3600523 5937.0295
    5     Algeria    Africa 1977  58.014 17152804 4910.4168
    6     Algeria    Africa 2007  72.301 33333216 6223.3675

``` r
ggplot(year_1977_2007) + aes(gdpPercap, lifeExp, col = continent) + geom_point() + facet_wrap(~year) + guides(color = guide_legend(title = "Continent")) + xlab("GDP per Capita") + ylab("Life Expectancy")
```

![](bimm143_lab5datavisR.1_files/figure-commonmark/unnamed-chunk-11-1.png)
