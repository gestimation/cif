test_that("Surv() produced the same outputs as Surv() of survfit", {
  testthat::skip_if_not_installed("survival")
  testdata <- createTestData(100, 1, first_zero=TRUE, last_zero=TRUE, subset_present=FALSE, logical_strata=FALSE, na_strata=FALSE)
  tested <- Surv(testdata$t, testdata$d)
  expected <- survival::Surv(testdata$t, testdata$d)
  expect_equal(expected, tested)
})

test_that("Surv() produced the same outputs as Surv() of survfit with NA", {
  testthat::skip_if_not_installed("survival")
  testdata <- createTestData(100, 1, first_zero=TRUE, last_zero=TRUE, subset_present=FALSE, logical_strata=FALSE, na_strata=FALSE)
  testdata$t[1] <- NA
  testdata$d[1] <- NA
  tested <- Surv(testdata$t, testdata$d)
  expected <- survival::Surv(testdata$t, testdata$d)
  expect_equal(expected, tested)
})

test_that("Surv() produced the same outputs as Surv() of survfit with a factor", {
  testthat::skip_if_not_installed("survival")
  testdata <- createTestData(100, 1, first_zero=TRUE, last_zero=TRUE, subset_present=FALSE, logical_strata=FALSE, na_strata=FALSE)
  f <- as.factor(testdata$d)
  tested <- Surv(testdata$t, f)
  expected <- survival::Surv(testdata$t, f)
  expect_equal(expected[,2], tested[,2])
})
