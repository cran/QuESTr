##' Sort and truncate genes according to the strength of gene-environment interaction
##'
##' @title Sort and truncate genes according to the strength of gene-environment interaction
##' @param x A data matrix (raw: samples, col: genes).
##' @param y A vector of an environment in which the samples were collected.
##' @param method A string to specify the method of regression for calculating R-squared values.
##' "linear" (default), "quadratic" or "cubic" regression model can be specified.
##' @param n.gene The number of genes to be included in QuEST model (default: ncol(x)).
##' @param trunc a threshold to be truncated (default: 1).
##' @importFrom stats lm
##' @return A data matrix (raw: samples, col: sorted genes)
##' @examples
##' data(Pinus)
##' train <- q.clean(Pinus$train)
##' target <- Pinus$target
##' cor(target, train[, 1])
##'
##' train <- q.sort(train, target, trunc = 0.5)
##' cor(target, train[, 1])
##' @author Takahiko Koizumi
##' @export
q.sort <- function(x, y, method = "linear", n.gene = ncol(x), trunc = 1){
  degree <- switch(method,
                   "linear" = 1,
                   "quadratic" = 2,
                   "cubic" = 3,
                   stop("Select the <method> linear, quadratic, or cubic")
  )

  if(n.gene < ncol(x) & trunc < 1){
    stop("Don't specify <n.gene> and <trunc> at a time")
  }
  if(n.gene < 0){
    stop("<n.gene> should not be a negative value")
  }else if(n.gene > ncol(x)){
    stop(paste("<n.gene> must not exceed", ncol(x), sep = " "))
  }
  if(trunc < 0 | trunc > 1){
    stop("<trunc> should be within the range of 0-1")
  }

  ## calculate R2 values
  result <- rep(NA, ncol(x))
  for(i in 1:ncol(x)){
    result[i] <- summary(lm(y ~ poly(x[, i], degree = degree, raw = TRUE)))$r.squared
  }
  ## sort genes in descending order of R2
  x <- x[, order(result, decreasing = TRUE)]

  ## extract genes with higher R2
  if(n.gene <= ncol(x)){
    x <- x[, 1:n.gene]
  }else if(trunc < 1){
    x <- x[, 1:length(result[result >= trunc])]
  }else{
    x <- x
  }
}
