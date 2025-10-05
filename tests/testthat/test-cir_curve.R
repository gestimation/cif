test_that("cif_curve() produced the same outputs as survfit() in survival in survival data", {
  library(survival)
  testdata <- createTestData(20, 1, first_zero=FALSE, last_zero=TRUE, subset_present=FALSE, logical_strata=TRUE, na_strata=FALSE)
  e <- survfit(Surv(t, d)~strata, testdata, weight=w, conf.type = "none")
  t <- cif_curve(Surv(t, d) ~ strata, data = testdata, weight="w", report.ggsurvfit = FALSE, conf.type = "none", outcome.type = "SURVIVAL")
  #  expected <- as.numeric(c(e$time, round(e$surv,digit=5), e$n, e$n.risk, e$n.event, e$n.censor, round(e$std.err,digit=5), e$strata))
  #  tested <- as.numeric(c(t$time, round(t$surv,digit=5), t$n, t$n.risk, t$n.event, t$n.censor, round(t$std.err,digit=5), t$strata))
  expected <- as.numeric(c(e$time, round(e$surv,digit=5), e$n, e$n.risk, e$n.event, e$n.censor, e$lower, e$strata))
  tested <- as.numeric(c(t$time, round(t$surv,digit=5), t$n, t$n.risk, t$n.event, t$n.censor, t$lower, t$strata))
  expect_equal(expected, tested)
})

test_that("cif_curve() produced the same estimates as cif() in mets in competing risks data", {
  library(mets)
  data(diabetes.complications)
  cif_fit <- cif(Event(t,epsilon) ~ +1, data=diabetes.complications, cause=1)
  surv_fit <- cif_curve(Surv(t, epsilon) ~ +1, data = diabetes.complications, report.ggsurvfit = FALSE, conf.type = "none", outcome.type = "C", error="delta")
  index<-c(1,3,5,6,8,10,12,13,15,16,17,18)
  expected <- round(1-surv_fit$surv[index],digit=5)
  tested <- round(cif_fit$mu[1:12],digit=5)
  expect_equal(expected, tested)
})


