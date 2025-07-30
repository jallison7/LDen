source(file = "../testdata.R")

test_that("global_density() calculates density correctly", {
  expect_equal(
    global_density(test_locs, area = 10),
    nrow(test_locs) / 10
  )
})

test_that("global_density() infers area from bounding box if unspecified", {
  expect_equal(
    global_density(test_locs),
    nrow(test_locs) / area_bbox(test_locs)
  )
})

test_that("area_bbox() calculates area of bounding box correctly", {
  expect_equal(area_bbox(data.frame(c(0, 1), c(0, 1))), 1)
})
