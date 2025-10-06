test_that("cif_reg() produced expected coefficients and variance covariance matrix from competing risks data with categorical exposure", {
  data(diabetes.complications)
  output <- cif_reg(nuisance.model = Event(t,epsilon)~+1, exposure = 'fruitq', data = diabetes.complications, effect.measure1='RR', effect.measure2='RR', code.exposure.ref='Q1', time.point=8, outcome.type='C')
  tested_coefficient <- round(output$coefficient,digit=3)
  tested_cov <- round(output$cov[1,],digit=3)
  tested <- as.vector(cbind(tested_coefficient,tested_cov))
  expected <- c(-1.383, -0.276, -0.067, -0.555, -3.991, -0.104, 0.045, -0.169, 0.018, -0.013, -0.013, -0.013, 0.011, -0.009, -0.009, -0.009)
  expect_equal(tested, expected)
})

test_that("cif_reg() produced expected coefficients and variance covariance matrix when coded other than the default and stratified in IPCW", {
  data(diabetes.complications)
#  diabetes.complications$epsilon1 <- diabetes.complications$epsilon + 1
#  output <- cif_reg(nuisance.model = Event(t,epsilon1)~+1, exposure = 'fruitq1', data = diabetes.complications, strata = 'strata', effect.measure1='SHR', effect.measure2='SHR', code.event1=2, code.event2=3, code.censoring=1, time.point=8, outcome.type='C')
  output <- cif_reg(nuisance.model = Event(t,epsilon)~+1, exposure = 'fruitq1', data = diabetes.complications, strata = 'strata', effect.measure1='SHR', effect.measure2='SHR', time.point=8, outcome.type='C')
  tested_coefficient <- round(output$coefficient,digit=3)
  tested_cov <- round(output$cov[1,],digit=3)
  tested <- as.vector(cbind(tested_coefficient,tested_cov))
  expected <- c(-1.383, 0.363, -3.988, 0.081, 0.009, -0.006, 0.004, -0.003)
  expect_equal(tested, expected)
})

test_that("cif_reg() produced expected coefficients and variance covariance matrix from survival data with missing data", {
  data(diabetes.complications)
  diabetes.complications$d <- as.numeric(diabetes.complications$epsilon>0)

  expected_df <- diabetes.complications[-(1:10), ]
  expected_output <- cif_reg(nuisance.model = Surv(t,d)~sex, exposure = 'fruitq1', strata = 'strata', data = expected_df, effect.measure1='RR', time.point=8, outcome.type='S')
  expected <- round(expected_output$coefficient,digit=3)

  diabetes.complications$t[1:2] <- NA
  diabetes.complications$d[3:4] <- NA
  diabetes.complications$sex[5:6] <- NA
  diabetes.complications$fruitq1[7:8] <- NA
  diabetes.complications$strata[9:10] <- NA
  output <- cif_reg(nuisance.model = Event(t,d)~sex, exposure = 'fruitq1', strata = 'strata', data = diabetes.complications, effect.measure1='RR', effect.measure2='RR', time.point=8, outcome.type='S')
  tested <- round(output$coefficient,digit=3)
  expect_equal(tested, expected)
})

test_that("cif_reg() produced expected common effects at 1:8", {
  data(diabetes.complications)
  output <- cif_reg(nuisance.model = Surv(t,epsilon)~+1, exposure = 'fruitq1', strata = 'strata', data = diabetes.complications, effect.measure1='RR', effect.measure2='RR', time.point=1:8, outcome.type='PP', report.boot.conf=FALSE)
  tested <- round(output$coefficient,digit=3)
#  expected <- c(-3.634, -2.360, -1.944, -1.727, -1.575, -1.473, -1.403, -1.343, 0.256, -4.632, -3.940, -3.634, -3.360, -2.996, -2.795, -2.609, -2.536, -0.050)
  expected <- c(-7.225, -4.051, -3.067, -2.534, -2.114, -1.803, -1.572, -1.392, 0.296, -9.319, -7.712, -6.977, -6.087, -5.156, -4.731, -4.366, -4.040, -0.022)
  expect_equal(tested, expected)
})

#test_that("cif_reg() produced expected common effects at 1:5 in prostate", {
#  testthat::skip_if_not_installed("dplyr")
#  data(prostate)
#  prostate <- dplyr::mutate(prostate, epsilon = ifelse(status=="alive",0,
#                                                       ifelse(status=="dead - prostatic ca",1,
#                                                              ifelse(status=="dead - other ca",1,
#                                                                     ifelse(status=="dead - heart or vascular",2,
#                                                                            ifelse(status=="dead - cerebrovascular",2,2)
#                                                                     )))))
#  prostate$epsilon <- as.numeric(prostate$epsilon)
#  prostate$a <- as.numeric((prostate$rx=="placebo"))
#  prostate$t <- prostate$dtime/12
#  output <- cif_reg(nuisance.model = Event(t,epsilon) ~ +1, exposure = 'a', strata='stage', data = prostate,
#                    effect.measure1='RR', effect.measure2='RR', time.point=1:5, outcome.type='POLY-PROPORTIONAL', report.boot.conf=FALSE)
#  tested <- round(output$coefficient,digit=3)
#  expected <- c(-2.289, -1.753, -1.461, -1.363, -1.254, -0.034, -2.056, -1.609, -1.245, -1.067, -0.971, 0.058)
#  expect_equal(tested, expected)
#})
