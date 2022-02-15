
<!-- README.md is generated from README.Rmd. Please edit that file -->

# QuESTr

<!-- badges: start -->
<!-- badges: end -->

The goal of QuESTr is to provide the functions to construct a prediction
model of environments using transcriptomes linked with the environments
based on QuEST algorithm. This package can quest for candidate genes for
the model construction even in non-model organinsms’ transcriptomes
without any genetic information.

## Installation

You can install the development version of QuESTr from
[GitHub](https://github.com/takakoizumi/QuESTr) with:  
\#install.packages(“devtools”)  
devtools::install\_github(“takakoizumi/QuESTr”)

## Example

install.packages(“QuESTr”)  
library(QuESTr)  
\# basic example code  
data(Pinus)  
x.train &lt;- Pinus\[\[1\]\]  
x.test &lt;- Pinus\[\[2\]\]  
y &lt;- Pinus\[\[3\]\]  
cor(y, quest(x.train, y, newx = x.test, n.gene = 100))
