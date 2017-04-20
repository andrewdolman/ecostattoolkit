
ecostattoolkit
==============

An R package to implement and test the ECOSTAT toolkit for estimating nutrient boundaries.

Implements and tests the ECOSTAT toolkit for estimating nutrient boundaries. R functions have been created that replicate the original scripts by Geoff Phillips.

Installation
------------

ecostattoolkit can be installed directly from GitHub using the `install_github` function in package `devtools`

``` r
# install.packages("devtools")
devtools::install_github("andrewdolman/ecostattoolkit")
```

Usage
-----

### Required packages

``` r
library(ecostattoolkit)
library(dplyr)
library(tidyr)
```

### Read in and prepare data

``` r
filename <- system.file("extdata", "DataTemplate_Example1.csv", package = "ecostattoolkit")
dat <- read.csv(filename)

status_boundaries <- c("HG" = 0.8, "GM" = 0.6, "MP" = 0.4)

# min and max nutrient value used for model
nut.min <- 1    
nut.max <- 138  

dat$status <- cut(dat$EQR,
                  breaks = rev(c(Inf, status_boundaries, -Inf)),
                  labels = rev(c("High", "Good", "Moderate", "Poor")))

dat$x.u <- log10(dat$P)
dat$y.u <- dat$EQR

knitr::kable(head(dat))
```

|  Record|  Unique\_ID|    EQR|     P|  Exclude\_P|     N|  Exclude\_N|  BioClass| Group | status   |       x.u|    y.u|
|-------:|-----------:|------:|-----:|-----------:|-----:|-----------:|---------:|:------|:---------|---------:|------:|
|       1|           1|  0.781|  55.3|           0|  1.39|           0|         4| BE    | Good     |  1.742725|  0.781|
|       2|           5|  0.771|  76.0|           0|  4.33|           0|         4| BE    | Good     |  1.880814|  0.771|
|       3|          50|  0.505|  57.2|           0|  1.22|           0|         3| DE    | Moderate |  1.757396|  0.505|
|       4|          59|  1.015|  44.4|           0|  1.09|           0|         5| DE    | High     |  1.647383|  1.015|
|       5|          60|  0.661|  28.6|           0|  0.64|           0|         4| DE    | Good     |  1.456366|  0.661|
|       6|          63|  0.438|  93.3|           0|  1.19|           0|         3| DE    | Moderate |  1.969882|  0.438|

### Estimate boundaries with all methods

``` r
dat_sub <- filter(dat, Exclude_P == 0,
                  P <= nut.max,
                  P >= nut.min)

boundaries <- boundaries_all(dat)


boundaries <- dplyr::mutate_if(boundaries, is.numeric, dplyr::funs(round(10^., 0))) %>% 
  select(-MP) 

knitr::kable(boundaries)
```

| Method     |   HG|   GM|
|:-----------|----:|----:|
| OLS1\_YonX |   27|   61|
| OLS2\_XonY |   31|   46|
| RMA        |   29|   53|
| MM         |   29|   46|
| AdjQ       |   24|   49|
| Median     |   31|   51|
| q75        |   27|   53|

### Compare with results from scripts in *TKit\_20170208\_CORRECTED*

Run script

``` r
source(system.file("extdata",
                   "TKit_20170208_CORRECTED/TKit_fit_lin_mod1.R",
                   package = "ecostattoolkit"))
```

![](README_files/figure-markdown_github/unnamed-chunk-4-1.png)![](README_files/figure-markdown_github/unnamed-chunk-4-2.png)![](README_files/figure-markdown_github/unnamed-chunk-4-3.png)

``` r
# Collect results
script.res <- data.frame(
  Method = c(rep("OLS1_YonX", 2), rep("OLS2_XonY", 2),
             rep("RMA", 2)),
  
  Boundary = c("GM", "HG", "GM", "HG",
               "GM", "HG"),
  
  Value = c(GM.mod1, HG.mod1, GM.mod2, HG.mod2,
            GM.mod4, HG.mod4),
  stringsAsFactors = TRUE)

script.res %>% 
  tidyr::spread(Boundary, Value) %>% 
  dplyr::select(Method, HG, GM) %>% 
  arrange(Method) %>% 
  knitr::kable(.)
```

| Method     |   HG|   GM|
|:-----------|----:|----:|
| OLS1\_YonX |   27|   53|
| OLS2\_XonY |   30|   43|
| RMA        |   28|   49|
