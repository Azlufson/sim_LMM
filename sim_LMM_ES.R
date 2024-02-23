##Stärke des Effekts

library(future.apply)
library(parallel)
library(tidyverse)
library(lme4)
library(lmerTest)
library(pbkrtest)
library(afex)

sim_data.n_int <- function(n.subj, n.obs, b0, beta_obs, beta_cond, sd.int_subj, sd_eps) {
  subj <- rep(1:n.subj, each = n.obs)
  obs <- rep(rep(c(0,1), each = n.obs/2), n.subj)
  cond <- rep(c(0,1), n.subj*n.obs/2)
  subj_int <- rep(rnorm(n.subj, 0, sd.int_subj), each = n.obs)
  y <- b0 + beta_obs * obs + beta_cond * cond + subj_int + rnorm(n.obs*n.subj, 0, sd_eps)
  return(data.frame(subj, obs, cond, y))
}

##LRT:
#alt (keine modellspezifikation möglich):
# test_lrtstat.fixed <- function(n.subj = 6, n.obs = 10, beta_obs = 0, REML = TRUE) {
#   data <- sim_data.n_int(n.subj = n.subj, n.obs = n.obs, b0 = 10, beta_obs = beta_obs, beta_cond = 5, sd.int_subj = 5, sd_eps = 1)
#   full <- lmer(y ~ obs + cond + (1|subj), data = data, REML = REML)
#   null <- lmer(y ~ cond + (1|subj), data = data, REML = REML)
#   return(pchisq(as.numeric(2 * (logLik(full) - logLik(null))), 1, lower = FALSE))
# }

test_lrtstat <- function(data, m.full, m.null, REML = TRUE) {
  full <- lmer(m.full, data = data, REML = REML)
  null <- lmer(m.null, data = data, REML = REML)
  return(pchisq(as.numeric(2 * (logLik(full) - logLik(null))), 1, lower = FALSE))
}

##t-as-z
test_TasZ.fixed <- function(data, m.full, REML = TRUE) {
  return(summary(lmer(model, data = data, REML = REML))$coefficients[2,4])
}

##KR, SW
#ddf ... Art der Approximation
test_approx.fixed <- function(data, model, REML = TRUE, ddf = "Satterthwaite") {
  return(anova(lmer(model, data = data, REML = REML), ddf = ddf)$`Pr(>F)`[1])
}

##Funktion zur Ausgabe des p-Wertes des fixed effects via parametric bootstrap
#mixed auf afex (nutzt pbmodcomp)
#nsim.pb bestimmt anzahl an bootstrap-simulationen von pbmodcomp
#cl erlaubt multicore nutzung (via package parallel)
test_PB.fixed <- function(mode, data, nsim.pb = 1000, cl = NULL) {
  return(suppressMessages(mixed(model, data = data, method = "PB", progress = FALSE, cl = cl, args_test = list(nsim = nsim.pb, cl = cl))$anova_table$`Pr(>PB)`[1]))
}
#suppressMessages: "mixed" will throw a message if numerical variables are not centered on 0, as main effects (of other variables then the numeric one) can be hard to interpret if numerical variables appear in interactions. See Dalal & Zickar (2012).
#kleine Tests haben ergeben, dass es am effizientesten ist, fürs fitten und fürs bootstrap multicore zu nutzen


#full and null model (LRT):
m.full <- y ~ obs + cond + (1|subj)
m.null <- y ~ cond + (1|subj)

#model
model <- y ~ obs + cond + (1|subj)

#Parameter für Simulationen
beta_obs <- 0 #auf diesen fixed effect wird jeweils getestet
n.subj <- c(4, 6, 10, 16)
n.obs <- c(4, 6, 10, 16)
grid <- expand.grid(n.subj, n.obs)
colnames(grid) <- c("n.subj", "n.obs")

#future_apply
plan("multisession", workers = detectCores())

#Parameter für parametric bootstrap
nsim.mixed <- 100 #niedriger, weil pro iteration auch noch gebootstrapped wird (mit nsim.pb)
nsim.pb <- 500 