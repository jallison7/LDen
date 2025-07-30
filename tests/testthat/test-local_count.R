source(file = "../testdata.R")

# copied from the file "LDEN1.out"
LDen1_out <- data.frame(
  x = c(1, 1, 2, 3, 8, 9, 9, 10, 6, 6, 5, 7, 12, 13, 14, 15),
  y = c(1, 2, 1, 3, 2, 1, 3, 2, 10, 11, 12, 12, 8, 7, 9, 7),
  type = c(0, 0, 0, 0, 2, 0, 2, 0, 10, 11, 12, 13, 4, 2, 2, 4),
  #var = 1,
  count_0 = c(3, 2, 2, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0),
  count_2 = rep(c(0, 1, 0, 1, 0, 1, 0), c(4, 1, 1, 1, 6, 2, 1)),
  count_4 = rep(rep(0:1, 2), c(12, 1, 2, 1)),
  count_10 = rep(c(0, 1, 0), c(8, 2, 6)),
  count_11 = rep(c(0, 1, 0), c(8, 2, 6)),
  count_12 = rep(c(0, 1, 0), c(10, 1, 5)),
  count_13 = rep(c(0, 1, 0), c(11, 1, 4)),
  count_total = rep(c(3, 2, 1, 2, 1), c(1, 2, 5, 2, 6))
)

test_that("local_count() matches TFQA'S LDEN with radius 1", {
  lcount <- local_count(test_locs, radius = 1)
  expect_equal(lcount, LDen1_out$count_total, ignore_attr = TRUE)
})


# copied from the file "LDEN2.out"
LDen2_out <- data.frame(
  type = c(0, 0, 0, 0, 2, 0, 2, 0, 10, 11, 12, 13, 4, 2, 2, 4),
  #var = 1,
  x = c(1, 1, 2, 3, 8, 9, 9, 10, 6, 6, 5, 7, 12, 13, 14, 15),
  y = c(1, 2, 1, 3, 2, 1, 3, 2, 10, 11, 12, 12, 8, 7, 9, 7),
  count_0 = rep(c(3, 1, 2, 0), c(3, 1, 4, 8)),
  count_2 = rep(c(0, 2, 0, 1), each = 4),
  count_4 = rep(c(0, 1, 2, 0, 1), rep(c(12, 1), c(1, 4))),
  count_10 = rep(c(0, 1, 0), c(8, 2, 6)),
  count_11 = rep(c(0, 1, 0), c(8, 4, 4)),
  count_12 = rep(c(0, 1, 0), c(9, 3, 4)),
  count_13 = rep(c(0, 1, 0), c(9, 3, 4)),
  count_total = c(3, 3, 3, 1, 4, 4, 4, 4, 2, 4, 3, 3, 2, 3, 1, 2)
)

test_that("local_count() matches TFQA'S LDEN with radius 2", {
  lcount <- local_count(test_locs, radius = 2)
  expect_equal(lcount, LDen2_out$count_total, ignore_attr = TRUE)
})


# copied from the file "LDEN5.out"
LDen5_out <- data.frame(
  type = c(0, 0, 0, 0, 2, 0, 2, 0, 10, 11, 12, 13, 4, 2, 2, 4),
  #var = 1,
  x = c(1, 1, 2, 3, 8, 9, 9, 10, 6, 6, 5, 7, 12, 13, 14, 15),
  y = c(1, 2, 1, 3, 2, 1, 3, 2, 10, 11, 12, 12, 8, 7, 9, 7),
  count_0 = rep(c(4, 2, 0), c(4, 4, 8)),
  count_2 = rep(rep(c(0, 2), 2), each = 4),
  count_4 = rep(c(0, 2), c(12, 4)),
  count_10 = rep(c(0, 1, 0), c(8, 4, 4)),
  count_11 = rep(c(0, 1, 0), c(8, 4, 4)),
  count_12 = rep(c(0, 1, 0), c(8, 4, 4)),
  count_13 = rep(c(0, 1, 0), c(8, 4, 4)),
  count_total = 4
)

test_that("local_count() matches TFQA'S LDEN with radius 5", {
  lcount <- local_count(test_locs, radius = 5)
  expect_equal(lcount, LDen5_out$count_total, ignore_attr = TRUE)
})
