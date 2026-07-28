# Test data from TFQA's HOA.TXT
hoa_test <- read.table(test_path("test_data", "HOA.TXT"),
                       col.names = c("x", "y", "type"))

test_that("association() matches TFQA's HOA", {
  hoa_test_type2 <- hoa_test[hoa_test$type == 2,]
  hoa_test_type4 <- hoa_test[hoa_test$type == 4,]
  expect_equal(
    round(association(hoa_test_type2, hoa_test_type4), 2), # TFQA rounds
    0.83
  )
})

test_that("dispersion() matches TFQA's HOA", {
  hoa_test_type2 <- hoa_test[hoa_test$type == 2,]
  hoa_test_type4 <- hoa_test[hoa_test$type == 4,]
  expect_equal(
    round(dispersion(hoa_test_type2, hoa_test_type4), 2), # TFQA rounds
    1.76
  )
})
