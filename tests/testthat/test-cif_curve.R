test_that("cif_curve() produced the same outputs as survfit() in survival in survival data", {
  testdata <- createTestData(20, 1, first_zero=FALSE, last_zero=TRUE, subset_present=FALSE, logical_strata=TRUE, na_strata=FALSE)
  t <- cif_curve(Surv(t, d) ~ strata, data = testdata, weight="w", report.ggsurvfit = FALSE, conf.type = "none", outcome.type = "SURVIVAL")
  tested <- as.numeric(c(t$time, round(t$surv,digit=5), t$n, t$n.risk, t$n.event, t$n.censor, t$lower, t$strata))
  #  e <- survival::survfit(Surv(t, d)~strata, testdata, weight=w, conf.type = "none")
  #  expected <- as.numeric(c(e$time, round(e$surv,digit=5), e$n, e$n.risk, e$n.event, e$n.censor, round(e$std.err,digit=5), e$strata))
  #  tested <- as.numeric(c(t$time, round(t$surv,digit=5), t$n, t$n.risk, t$n.event, t$n.censor, round(t$std.err,digit=5), t$strata))
  #  expected <- as.numeric(c(e$time, round(e$surv,digit=5), e$n, e$n.risk, e$n.event, e$n.censor, e$lower, e$strata))
  expected <- c(1.0,3.0,5.0,7.0,9.0,2.0,4.0,6.0,8.0,10.0,0.8,0.6,0.4,0.2,0.0,0.8,0.6,0.4,0.2,0.2,10.0,10.0,10.0,8.0,6.0,
                4.0,2.0,10.0,8.0,6.0,4.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,2.0,5.0,5.0)
  expect_equal(expected, tested)
})

test_that("cif_curve() produced the same estimates as cif() in mets in competing risks data", {
  data(diabetes.complications)
  surv_fit <- cif_curve(Surv(t, epsilon) ~ +1, data = diabetes.complications, report.ggsurvfit = FALSE, conf.type = "none", outcome.type = "C", error="delta")
  index<-c(1,3,5,6,8,10,12,13,15,16,17,18)
  tested <- round(1-surv_fit$surv[index],digit=5)
  #  library(mets)
  #  cif_fit <- mets::cif(Event(t,epsilon) ~ +1, data=diabetes.complications, cause=1)
  #  expected <- round(cif_fit$mu[1:12],digit=5)
  expected <- c(0.00102,0.00204,0.00307,0.00409,0.00511,0.00613,0.00716,0.00818,0.00920,0.01022,0.01125,0.01227)
  expect_equal(expected, tested)
})

