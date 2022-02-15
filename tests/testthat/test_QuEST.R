test_that("sort", {

  y <- runif(30)
  x <- matrix(runif(30 * 100), 30, 100)

  expect_equal(nrow(q.sort(x, y)), 30)
  expect_equal(ncol(q.sort(x, y)), 100)
})
