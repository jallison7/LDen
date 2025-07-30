source(file = "../testdata.R")

test_that("lda() with LDA.TXT data matches TFQA output with radius 1", {
  ldens <- lda(test_locs, type = test_locs$type, radius = 1, area = 154)
  expect_equal(test_lda1, round(ldens, 2)) # TFQA rounds, we don't
})

test_that("lda() with LDA.TXT data matches TFQA output with radius 2", {
  ldens <- lda(test_locs, type = test_locs$type, radius = 2, area = 154)
  expect_equal(test_lda2, round(ldens, 2)) # TFQA rounds, we don't
})

test_that("lda() with LDA.TXT data matches TFQA output with radius 5", {
  ldens <- lda(test_locs, type = test_locs$type, radius = 5, area = 154)
  expect_equal(test_lda5, round(ldens, 2)) # TFQA rounds, we don't
})
