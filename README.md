
<!-- README.md is generated from README.Rmd. Please edit that file -->

# grantham <img src='man/figures/logo.svg' align="right" height="139" />

<!-- badges: start -->

[![CRAN
status](https://www.r-pkg.org/badges/version/grantham)](https://CRAN.R-project.org/package=grantham)
<!-- badges: end -->

The goal of `{grantham}` is to provide a minimal set of routines to
calculate the Grantham distance<sup>[1](#1)</sup>.

The Grantham distance attempts to provide a proxy for the evolutionary
distance between two amino acids based on three key side chain chemical
properties: composition, polarity and molecular volume. In turn,
evolutionary distance is used as a proxy for the impact of missense
substitutions. The higher the distance, the more deleterious the
substitution is expected to be.

## Installation

Install `{grantham}` from CRAN:

``` r
install.packages("grantham")
```

You can install the development version of `{grantham}` like so:

``` r
# install.packages("remotes")
remotes::install_github("maialab/grantham")
```

## Usage

Grantham distance between two amino acids:

``` r
library(grantham)

grantham_distance(x = 'Ser', y = 'Phe')
#> # A tibble: 1 × 3
#>   x     y         d
#>   <chr> <chr> <dbl>
#> 1 Ser   Phe     155
```

The function `grantham_distance()` is vectorised with amino acids being
matched element-wise to form pairs for comparison:

``` r
grantham_distance(x = c('Ser', 'Arg'), y = c('Phe', 'Leu'))
#> # A tibble: 2 × 3
#>   x     y         d
#>   <chr> <chr> <dbl>
#> 1 Ser   Phe     155
#> 2 Arg   Leu     102
```

The two vectors of amino acids must have compatible sizes in the sense
of [vec_recycle()](https://vctrs.r-lib.org/reference/vec_recycle.html)
for element recycling to be possible, i.e., either the two vectors have
the same length, or one of them is of length one, and it is recycled up
to the length of the other.

``` r
# `'Ser'` is recycled to match the length of the second vector, i.e. 3.
grantham_distance(x = 'Ser', y = c('Phe', 'Leu', 'Arg'))
#> # A tibble: 3 × 3
#>   x     y         d
#>   <chr> <chr> <dbl>
#> 1 Ser   Phe     155
#> 2 Ser   Leu     145
#> 3 Ser   Arg     110
```

Use the function `amino_acid_pairs()` to generate all 20 x 20 amino acid
pairs:

``` r
aa_pairs <- amino_acid_pairs()
aa_pairs
#> # A tibble: 400 × 2
#>    x     y    
#>    <chr> <chr>
#>  1 Ser   Ser  
#>  2 Ser   Arg  
#>  3 Ser   Leu  
#>  4 Ser   Pro  
#>  5 Ser   Thr  
#>  6 Ser   Ala  
#>  7 Ser   Val  
#>  8 Ser   Gly  
#>  9 Ser   Ile  
#> 10 Ser   Phe  
#> # … with 390 more rows
#> # ℹ Use `print(n = ...)` to see more rows
```

And now calculate all Grantham distances for all pairs `aa_pairs`:

``` r
grantham_distance(x = aa_pairs$x, y = aa_pairs$y)
#> # A tibble: 400 × 3
#>    x     y         d
#>    <chr> <chr> <dbl>
#>  1 Ser   Ser       0
#>  2 Ser   Arg     110
#>  3 Ser   Leu     145
#>  4 Ser   Pro      74
#>  5 Ser   Thr      58
#>  6 Ser   Ala      99
#>  7 Ser   Val     124
#>  8 Ser   Gly      56
#>  9 Ser   Ile     142
#> 10 Ser   Phe     155
#> # … with 390 more rows
#> # ℹ Use `print(n = ...)` to see more rows
```

Because distances are symmetric, and for pairs formed by the same amino
acid are trivially zero, you might want to exclude these pairs:

``` r
# `keep_self = FALSE`: excludes pairs such as ("Ser", "Ser")
# `keep_reverses = FALSE`: excludes reversed pairs, e.g. ("Arg", "Ser") will be
# removed because ("Ser", "Arg") already exists.
aa_pairs <- amino_acid_pairs(keep_self = FALSE, keep_reverses = FALSE)

# These amino acid pairs are the 190 pairs shown in Table 2 of Grantham's
# original publication.
aa_pairs
#> # A tibble: 190 × 2
#>    x     y    
#>    <chr> <chr>
#>  1 Ser   Arg  
#>  2 Ser   Leu  
#>  3 Ser   Pro  
#>  4 Ser   Thr  
#>  5 Ser   Ala  
#>  6 Ser   Val  
#>  7 Ser   Gly  
#>  8 Ser   Ile  
#>  9 Ser   Phe  
#> 10 Ser   Tyr  
#> # … with 180 more rows
#> # ℹ Use `print(n = ...)` to see more rows

# Grantham distance for the 190 unique amino acid pairs
grantham_distance(x = aa_pairs$x, y = aa_pairs$y)
#> # A tibble: 190 × 3
#>    x     y         d
#>    <chr> <chr> <dbl>
#>  1 Ser   Arg     110
#>  2 Ser   Leu     145
#>  3 Ser   Pro      74
#>  4 Ser   Thr      58
#>  5 Ser   Ala      99
#>  6 Ser   Val     124
#>  7 Ser   Gly      56
#>  8 Ser   Ile     142
#>  9 Ser   Phe     155
#> 10 Ser   Tyr     144
#> # … with 180 more rows
#> # ℹ Use `print(n = ...)` to see more rows
```

The Grantham distance
![d\_{i,j}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;d_%7Bi%2Cj%7D "d_{i,j}")
for two amino acids
![i](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;i "i")
and
![j](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;j "j")
is:

![d\_{i,j} = \rho (\alpha (c_i-c_j)^2+\beta (p_i-p_j)^2+ \gamma (v_i-v_j)^2)^{1/2}\\ .](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;d_%7Bi%2Cj%7D%20%3D%20%5Crho%20%28%5Calpha%20%28c_i-c_j%29%5E2%2B%5Cbeta%20%28p_i-p_j%29%5E2%2B%20%5Cgamma%20%28v_i-v_j%29%5E2%29%5E%7B1%2F2%7D%5C%20. "d_{i,j} = \rho (\alpha (c_i-c_j)^2+\beta (p_i-p_j)^2+ \gamma (v_i-v_j)^2)^{1/2}\ .")

The distance is based on three chemical properties of amino acid side
chains:

-   composition
    (![c](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;c "c"))
-   polarity
    (![p](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;p "p"))
-   molecular volume
    (![v](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;v "v"))

We provide a data set with these properties:

``` r
amino_acids_properties
#> # A tibble: 20 × 4
#>    amino_acid     c     p     v
#>    <chr>      <dbl> <dbl> <dbl>
#>  1 Ser         1.42   9.2  32  
#>  2 Arg         0.65  10.5 124  
#>  3 Leu         0      4.9 111  
#>  4 Pro         0.39   8    32.5
#>  5 Thr         0.71   8.6  61  
#>  6 Ala         0      8.1  31  
#>  7 Val         0      5.9  84  
#>  8 Gly         0.74   9     3  
#>  9 Ile         0      5.2 111  
#> 10 Phe         0      5.2 132  
#> 11 Tyr         0.2    6.2 136  
#> 12 Cys         2.75   5.5  55  
#> 13 His         0.58  10.4  96  
#> 14 Gln         0.89  10.5  85  
#> 15 Asn         1.33  11.6  56  
#> 16 Lys         0.33  11.3 119  
#> 17 Asp         1.38  13    54  
#> 18 Glu         0.92  12.3  83  
#> 19 Met         0      5.7 105  
#> 20 Trp         0.13   5.4 170
```

If you want to calculate the Grantham distance from these property
values you may use the function `grantham_equation()`.

## Related software

Other sources we’ve found in the R ecosystem that also provide code for
calculation of the Grantham distance:

-   A GitHub Gist by Daniel E Cook provides the function
    `calculate_grantham()`, see
    [Fetch_Grantham.R](https://gist.github.com/danielecook/501f03650bca6a3db31ff3af2d413d2a).
-   The `{midasHLA}` package includes the unexported function
    `distGrantham()` in
    [utils.R](https://github.com/Genentech/midasHLA/blob/ec29296f9bfd7c4fae9e2040592b618e5f2a99a1/R/utils.R).
-   The `{HLAdivR}` package exports a data set with the Grantham
    distances in the format of a matrix, see
    [data.R](https://github.com/rbentham/HLAdivR/blob/master/R/data.R).
-   The Bioconductor package `{MSA2dist}` by Kristian K. Ullrich
    provides the function
    [`aastring2dist()`](https://www.bioconductor.org/packages/devel/bioc/vignettes/MSA2dist/inst/doc/MSA2dist.html#granthams-distance).

## Code of Conduct

Please note that the `{grantham}` package is released with a
[Contributor Code of
Conduct](https://contributor-covenant.org/version/2/0/CODE_OF_CONDUCT.html).
By contributing to this project, you agree to abide by its terms.

## References

<a id="1">1.</a> Grantham, R. *Amino acid difference formula to help
explain protein evolution*. Science 185, 862–864 (1974). doi:
[10.1126/science.185.4154.862](https://doi.org/10.1126/science.185.4154.862).
