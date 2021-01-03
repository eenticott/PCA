test_that("multiplication works", {
  expect_equal(2 * 2, 4)
})

test_that("centre works", {
  expect_equal(centre(matrix(seq(1:9), nrow = 3)), centre(3 * matrix(seq(1:9), nrow = 3), T))
})
