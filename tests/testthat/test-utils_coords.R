test_coords <- matrix(c(1:9, 1:9), ncol = 2)

test_that("coordinates can be specified as a data frame", {
  expect_equal(coord_matrix(data.frame(x = 1:9, y = 1:9)), test_coords)
})

test_that("coordinates can be specified as a matrix", {
  expect_equal(coord_matrix(test_coords), test_coords)
})

test_that("coordinates can be specified as a list", {
  expect_equal(coord_matrix(list(1:9, 1:9)), test_coords)
})

test_that("non-numeric coordinates are an error", {
  expect_error(coord_matrix(list(letters, LETTERS)), "is.numeric")
})

test_that("coordinates of different lengths are an error", {
  expect_error(
    coord_matrix(list(1:2, 1:3)),
    "arguments imply differing number of rows"
  )
})

test_that("coordinates with NAs are an error", {
  expect_error(coord_matrix(list(c(1, NA, 3), c(1, 2, NA))), "all.*is.na")
})
